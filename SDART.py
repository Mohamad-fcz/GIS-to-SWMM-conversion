# ==============================================================================
# SDART.py
# SDART: Stepwise GIS → SWMM preparation toolkit for ArcGIS Pro / ArcPy
#
# Version: 2.0.0 (2025-09-20)
# Author: Mohamad Ahmadzadeh
# Copyright (c) 2025 Mohamad Ahmadzadeh
# License: GNU General Public License v3.0 (GPL-3.0)
#
# Project intro
# -------------
# SDART automates common GIS-to-SWMM preprocessing: verifying dataset inventory,
# making safe backups, preparing line topology (unsplit → optional snap → split),
# aligning source layers by XY shift, robust spatial join of points and lines
# (multi-radius and multi-method sweeps with geometry-aware scoring), writing
# audit-ready Excel reports, filling missing names, and assigning From/To nodes.
#
# Quick prerequisites
# -------------------
# • ArcGIS Pro with ArcPy available in the current Python environment
# • GeoPandas (used for reading shapefiles into DataFrames for analysis)
# • pandas
#
# How to cite in your repo
# ------------------------
# SDART — Stepwise GIS → SWMM preparation toolkit, v2.0.0 (2025-09-20).
# Copyright (c) 2025 Mohamad Ahmadzadeh. Licensed under GPL-3.0.
#
# GPL-3.0 summary
# ---------------
# You CAN: use, study, run, and modify this code, and redistribute original or
# modified versions, including commercially, as long as you comply with GPL-3.0.
#
# You MUST: keep the copyright and license notices, make complete corresponding
# source code available when you distribute binaries or modified versions, and
# license the entire derivative work under GPL-3.0.
#
# You CANNOT: relicense the code (or any derivative work you distribute) under a
# non-GPL license, impose further restrictions beyond GPL-3.0, or remove notices.
#
# Warranty
# --------
# This software is provided “as is” without any warranty. Use at your own risk.
# ==============================================================================

from __future__ import annotations

# Standard library imports
import logging
import os
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# Third-party analytics
import pandas as pd

# ArcPy import with graceful fallback so callers see a clear message if ArcPy
# is not available. All ArcPy-using functions call _require_arcpy() before use.
try:
    import arcpy
except Exception as e:  # pragma: no cover
    arcpy = None
    _ARCPY_IMPORT_ERROR = e

# ------------------------------------------------------------------------------
# Naming conventions and canonical dataset keys
# ------------------------------------------------------------------------------
# Layer keys are the stable identifiers used across inventory, topology, and
# spatial join routines. They intentionally match SWMM concepts.

# Layer keys for points, lines, polygons
POINT_KEYS = ("junctions", "storages", "dividers", "outfalls")
LINE_KEYS  = ("conduits", "orifices", "outlets", "pumps", "weirs")
POLY_KEYS  = ("subcatchments",)  # polygons if present

# Canonical feature class and shapefile names that SDART expects and resolves
UPDATED_FC = {
    "points": {k: f"Updated_{k}" for k in POINT_KEYS},
    "lines":  {k: f"Updated_{k}" for k in LINE_KEYS},
    "polys":  {k: f"Updated_{k}" for k in POLY_KEYS},
}
SOURCE_SHP = {
    "points": {k: f"Source_{k}.shp" for k in POINT_KEYS},
    "lines":  {k: f"Source_{k}.shp" for k in LINE_KEYS},
    "polys":  {k: f"Source_{k}.shp" for k in POLY_KEYS},
}

# ------------------------------------------------------------------------------
# Configuration object
# ------------------------------------------------------------------------------
@dataclass
class SDARTConfig:
    """
    SDART configuration
    -------------------
    source_shapefiles: folder with Source_*.shp
    updated_shapefiles: FileGDB containing Updated_* feature classes
    scratch_gdb: FileGDB for temp outputs (created if missing)
    """
    source_shapefiles: Path
    updated_shapefiles: Path
    scratch_gdb: Path

# ------------------------------------------------------------------------------
# Utility helpers — minimal and focused, reused across the class
# ------------------------------------------------------------------------------
def _require_arcpy() -> None:
    # Centralized check so errors are clear and actionable for users
    if arcpy is None:
        raise RuntimeError(
            "ArcPy is not available. Activate ArcGIS Pro's Python environment.\n"
            f"Original import error: {_ARCPY_IMPORT_ERROR}"
        )

def _ensure_gdb(path: Path) -> None:
    # Create a file geodatabase if it does not exist yet
    _require_arcpy()
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        arcpy.CreateFileGDB_management(path.parent.as_posix(), path.name)

def _set_workspace(gdb: Path) -> None:
    # Make ArcPy tools write into the given workspace and allow overwrites
    _require_arcpy()
    arcpy.env.workspace = gdb.as_posix()
    arcpy.env.overwriteOutput = True

def _exists_fc(dataset: str | Path) -> bool:
    # ArcGIS-safe existence check for geodatabase feature classes or shapefiles
    _require_arcpy()
    return bool(arcpy.Exists(str(dataset)))

def _safe_delete(dataset: str | Path) -> None:
    # Delete a dataset only if it exists to avoid noisy exceptions
    _require_arcpy()
    if _exists_fc(dataset):
        arcpy.management.Delete(str(dataset))

def _unique_name(base: str, parent: Path) -> str:
    """Create a unique dataset/shapefile name inside a GDB or folder."""
    _require_arcpy()
    try:
        return arcpy.management.CreateUniqueName(base, parent.as_posix()).getOutput(0)
    except Exception:
        # Fallback name uses a timestamp to guarantee uniqueness
        ts = time.strftime("%Y%m%d_%H%M%S")
        stem = f"{Path(base).stem}_{ts}"
        if parent.suffix.lower() == ".gdb":
            return (parent / stem).as_posix()
        else:
            return (parent / f"{stem}{Path(base).suffix}").as_posix()

# ------------------------------------------------------------------------------
# Main class — all core workflows live here
# ------------------------------------------------------------------------------
class SDART:
    def __init__(self, cfg: SDARTConfig, logger: logging.Logger | None = None):
        # Initialize environment and a simple INFO-level logger
        _require_arcpy()
        self.cfg = cfg
        _ensure_gdb(self.cfg.scratch_gdb)
        if logger is None:
            logging.basicConfig(level=logging.INFO, format="INFO: %(message)s")
            logger = logging.getLogger("SDART")
        self.log = logger

        # Internal caches populated by inventory_table() for quick path lookups
        self.updated_found: dict[str, str] = {}  # feature class name → path
        self.source_found: dict[str, str]  = {}  # shapefile name → path

        # Legacy mapping of required fields for line layers preserved from prior
        # versions for compatibility. Some routines compute transfer fields
        # dynamically, this stays as a reference without altering behavior.
        self._LINE_REQ_FIELDS = {
            "conduits": [
                "Name", "FromNode", "ToNode", "Length", "Roughness", "InOffset", "OutOffset",
                "InitFlow", "MaxFlow", "Shape_1", "Geom1", "Geom2", "Geom3", "Geom4", "Barrels",
                "Culvert", "Shp_Trnsct", "Kentry", "Kexit", "Kavg", "FlapGate", "Seepage",
                "Shape_Leng", "Angle", "N_Curve", "Vertices"
            ],
            "weirs": [
                "Name", "FromNode", "ToNode", "Type", "CrestHeigh", "Qcoeff", "FlapGate",
                "EndContrac", "EndCoeff", "Surcharge", "RoadWidth", "RoadSurf", "CoeffCurve",
                "Height", "Length", "SideSlope",
                "Shape_Leng", "Angle", "N_Curve", "Vertices"
            ],
            "pumps": [
                "Name", "FromNode", "ToNode", "PumpCurve", "Status", "Startup", "Shutoff",
                "Shape_Leng", "Angle", "N_Curve", "Vertices"
            ],
            "outlets": [
                "Name", "FromNode", "ToNode", "InOffset", "RateCurve", "Qcoeff", "Qexpon",
                "FlapGate", "CurveName",
                "Shape_Leng", "Angle", "N_Curve", "Vertices"
            ],
            "orifices": [
                "Name", "FromNode", "ToNode", "Type", "InOffset", "Qcoeff", "FlapGate",
                "CloseTime", "Shape_1", "Height", "Width",
                "Shape_Leng", "Angle", "N_Curve", "Vertices"
            ],
        }

    # ------------------------------------------------------------------
    # Inventory — scan Updated_* in GDB and Source_* in folder
    # ------------------------------------------------------------------
    def inventory_table(self) -> pd.DataFrame:
        """
        Build a concise table of expected datasets, their presence, and resolved
        paths. Also populates self.updated_found and self.source_found caches.
        """
        _set_workspace(self.cfg.updated_shapefiles)

        rows: list[list[str]] = []

        def push(group: str, key: str):
            # Resolve expected Updated_* feature class in the target GDB
            fc_name = UPDATED_FC[group][key]
            fc_path = self.cfg.updated_shapefiles / fc_name
            fc_found = _exists_fc(fc_path)
            if fc_found:
                self.updated_found[fc_name] = fc_path.as_posix()

            # Resolve expected Source_* shapefile in the source folder
            shp_name = SOURCE_SHP[group][key]
            shp_path = self.cfg.source_shapefiles / shp_name
            shp_found = shp_path.exists()
            if shp_found:
                self.source_found[shp_name] = shp_path.as_posix()

            rows.append([
                key,
                fc_name, "Found" if fc_found else "Missing",
                shp_name, "Found" if shp_found else "Missing",
                fc_path.as_posix() if fc_found else "",
                shp_path.as_posix() if shp_found else "",
            ])

        # Enumerate all layer keys by type
        for k in POINT_KEYS: push("points", k)
        for k in LINE_KEYS:  push("lines",  k)
        for k in POLY_KEYS:  push("polys",  k)

        df = pd.DataFrame(
            rows,
            columns=[
                "Layer Key",
                "Updated Expected FC", "Updated Status",
                "Source Expected Shapefile", "Source Status",
                "Resolved Updated Path", "Resolved Source Path",
            ],
        )
        return df

    # ------------------------------------------------------------------
    # Backups — copy Updated GDB and Source folder with lock skipping
    # ------------------------------------------------------------------
    def backup(self, out_root: Path | None = None, overwrite: bool = False) -> tuple[Path, Path]:
        """
        Copy Updated_Shapefiles (GDB) and Source_Shapefiles (folder) into backups.
        Skips *.lock files to avoid permission errors.
        Returns (gis_backup_path, source_backup_path).
        """
        if out_root is None:
            out_root = self.cfg.updated_shapefiles.parent / "_sdart_backups"
        out_root.mkdir(parents=True, exist_ok=True)

        def _copytree_skip_locks(src: Path, dst: Path):
            if dst.exists():
                if overwrite:
                    shutil.rmtree(dst, ignore_errors=True)
                else:
                    ts = time.strftime("%Y%m%d_%H%M%S")
                    dst = dst.parent / f"{dst.name}_{ts}"
            def _ignore(dirpath, names):
                return [n for n in names if n.lower().endswith(".lock")]
            shutil.copytree(src, dst, ignore=_ignore)
            return dst

        gis_dst   = out_root / "GIS_Backup"
        source_dst= out_root / "Source_Backup"

        used_gis    = _copytree_skip_locks(self.cfg.updated_shapefiles, gis_dst)
        used_source = _copytree_skip_locks(self.cfg.source_shapefiles, source_dst)
        self.log.info(f"Backups created → {used_gis} , {used_source}")
        return used_gis, used_source

    # ------------------------------------------------------------------
    # Topology preparation — unsplit, optional snap, chain split at points
    # ------------------------------------------------------------------
    def prepare_topology(self, *, snap: bool = True, snap_radius_ft: float = 5.0) -> None:
        """
        For each present Updated_* line layer:
          1) Unsplit (merge) into scratch,
          2) If snap=True, Snap each Updated_* point layer to that line (radius in FEET),
          3) Split the merged lines at each point layer sequentially,
          4) Write final output back to Updated GDB overwriting the same line feature class.
        Notes:
        - The Snap & Split radius is passed as feet (ArcGIS accepts linear units).
        - We split once per point layer, chaining outputs so ALL points are honored.
        """
        _set_workspace(self.cfg.updated_shapefiles)

        # Refresh caches if needed so downstream steps resolve datasets
        if not self.updated_found:
            _ = self.inventory_table()

        # Collect point and line feature classes that actually exist
        updated_points = [self.updated_found.get(f"Updated_{k}") for k in POINT_KEYS]
        updated_points = [p for p in updated_points if p]

        updated_lines = [self.updated_found.get(f"Updated_{k}") for k in LINE_KEYS]
        updated_lines = [l for l in updated_lines if l]

        if not updated_lines:
            self.log.warning("No Updated_* line layers found; skipping topology preparation.")
            return

        # Prefer to split with junctions first, then storages, dividers, outfalls
        order = []
        for k in ("junctions", "storages", "dividers", "outfalls"):
            p = self.updated_found.get(f"Updated_{k}")
            if p: order.append(p)
        # Append any other point layers if present
        order += [p for p in updated_points if p not in order]

        for line_fc in updated_lines:
            line_name = Path(line_fc).name
            self.log.info(f"[{line_name}] Unsplit → {self.cfg.scratch_gdb}/<{line_name}_unsplit>")

            tmp_unsplit = (self.cfg.scratch_gdb / f"{line_name}_unsplit").as_posix()
            _safe_delete(tmp_unsplit)
            arcpy.management.UnsplitLine(line_fc, tmp_unsplit)

            # Optional point snapping prior to splitting to catch near-misses
            if snap and order:
                for pt_fc in order:
                    self.log.info(f"[{line_name}] Snap points → unsplit @ {snap_radius_ft} ft")
                    try:
                        arcpy.Snap_edit(pt_fc, [[tmp_unsplit, "EDGE", f"{snap_radius_ft} Feet"]])
                    except Exception:
                        self.log.warning(f"[{line_name}] Snap failed for {Path(pt_fc).name}: {arcpy.GetMessages(2)}")

            # Sequential split at each point layer; chain outputs so all point types are honored
            prev = tmp_unsplit
            step = 0
            for pt_fc in order:
                step += 1
                out_step = (self.cfg.scratch_gdb / f"{line_name}_split_{step}").as_posix()
                _safe_delete(out_step)
                self.log.info(f"[{line_name}] Split at {Path(pt_fc).name} → step {step}")
                arcpy.management.SplitLineAtPoint(prev, pt_fc, out_step, f"{snap_radius_ft} Feet")
                prev = out_step

            # Final output overwrites the original Updated_* feature class in place
            final_fc = line_fc  # overwrite in place
            _safe_delete(final_fc)
            arcpy.management.CopyFeatures(prev, final_fc)
            self.log.info(f"[{line_name}] Done → {final_fc}")

    # ------------------------------------------------------------------
    # XY shift utilities — compute dx,dy from 2 features, shift all Source_*.shp
    # ------------------------------------------------------------------
    # These helpers perform geometry-preserving translations on points, lines,
    # and polygons so a misaligned source package can be aligned to Updated_*.

    def _translate_points_fc(self, shp_path: Path, dx: float, dy: float) -> None:
        """Translate every point by dx, dy (dataset units)."""
        _require_arcpy()
        with arcpy.da.UpdateCursor(shp_path.as_posix(), ["SHAPE@XY"]) as cur:
            for (xy,) in cur:  # xy is a (x, y) tuple
                cur.updateRow([(xy[0] + dx, xy[1] + dy)])

    def _iter_all_source_shapefiles(self):
        """
        Yield Path objects for every shapefile in Source_Shapefiles whose name
        begins with 'Source_' this covers points, lines, polygons uniformly.
        """
        src = self.cfg.source_shapefiles
        for p in src.glob("Source_*.shp"):
            if p.is_file():
                yield p

    def _translate_polyline_fc(self, shp_path: Path, dx: float, dy: float) -> None:
        """Translate every polyline vertex by dx, dy."""
        _require_arcpy()
        sr = arcpy.Describe(shp_path.as_posix()).spatialReference
        with arcpy.da.UpdateCursor(shp_path.as_posix(), ["SHAPE@"]) as cur:
            for (geom,) in cur:
                parts = []
                for part in geom:  # part is an array of points
                    if part is None:
                        continue
                    new_part = [arcpy.Point(p.X + dx, p.Y + dy, p.Z, p.M) for p in part if p]
                    parts.append(arcpy.Array(new_part))
                new_geom = arcpy.Polyline(arcpy.Array(parts), sr)
                cur.updateRow([new_geom])

    def _translate_polygon_fc(self, shp_path: Path, dx: float, dy: float) -> None:
        """Translate every polygon including multipart and holes by dx, dy."""
        _require_arcpy()
        import json
        sr = arcpy.Describe(shp_path.as_posix()).spatialReference
        with arcpy.da.UpdateCursor(shp_path.as_posix(), ["SHAPE@"]) as cur:
            for (geom,) in cur:
                # Use Esri JSON to preserve rings and multipart structure
                data = json.loads(geom.JSON)
                for ring in data.get("rings", []):
                    for i in range(len(ring)):
                        ring[i][0] += dx
                        ring[i][1] += dy
                new_geom = arcpy.AsShape(data, True)  # True => esriJSON
                cur.updateRow([new_geom])

    def _xy_from_feature_label(self, dataset: str, label: str) -> tuple[float, float]:
        """Return centroid/point XY for a label like 'Feature 12' using OID."""
        _require_arcpy()
        try:
            oid = int(str(label).split()[-1])
        except Exception:
            raise ValueError(f"Bad feature label: {label!r}. Expected like 'Feature 12'.")
        with arcpy.da.SearchCursor(dataset, ["OID@", "SHAPE@TRUECENTROID"]) as cur:
            for row_oid, xy in cur:
                if row_oid == oid:
                    return float(xy[0]), float(xy[1])
        raise ValueError(f"OID {oid} not found in {dataset}.")

    def shift_sources_by_two_oids(
            self,
            *,
            moving_feature_label: str,  # e.g. "Feature 0" in Source_junctions.shp
            anchor_feature_label: str,  # e.g. "Feature 1" in Updated_junctions
            moving_points: str = "Source_junctions.shp",
            anchor_points_fc: str = "Updated_junctions",
    ) -> tuple[float, float]:
        """
        Compute dx,dy from (moving_points, moving_feature_label) to
        (anchor_points_fc, anchor_feature_label), then shift every Source_*.shp
        by that same (dx,dy). Works across geometry types.
        """
        _require_arcpy()

        # Resolve paths for moving and anchor point datasets
        moving_pts_path = (self.cfg.source_shapefiles / moving_points).as_posix()
        anchor_fc_path = (self.cfg.updated_shapefiles / anchor_points_fc).as_posix()

        # Compute translation vector from centroids so it is robust to type
        mx, my = self._xy_from_feature_label(moving_pts_path, moving_feature_label)
        ax, ay = self._xy_from_feature_label(anchor_fc_path, anchor_feature_label)
        dx, dy = (ax - mx), (ay - my)
        self.log.info(f"Computed shift dx={dx:.3f}, dy={dy:.3f}")

        # Apply translation to every Source_*.shp by geometry type
        for shp in self._iter_all_source_shapefiles():
            desc = arcpy.Describe(shp.as_posix())
            gtype = desc.shapeType.lower()  # 'point', 'polyline', 'polygon'
            self.log.info(f"Shifting {shp.name} ({gtype}) by dx={dx:.3f}, dy={dy:.3f}")
            if gtype == "point":
                self._translate_points_fc(shp, dx, dy)
            elif gtype == "polyline":
                self._translate_polyline_fc(shp, dx, dy)
            elif gtype == "polygon":
                self._translate_polygon_fc(shp, dx, dy)
            else:
                self.log.warning(f"Unknown shape type for {shp.name}: {gtype} — skipped.")
        self.log.info("All Source_* layers shifted.")
        return dx, dy

    # ------------------------------------------------------------------
    # Spatial Join suite — helpers and two main runners (points, lines)
    # ------------------------------------------------------------------
    # Notes:
    # • Paths are normalized for ArcPy calls.
    # • Geometry attributes are computed once per sweep to enable scoring.
    # • Shapefile exports are used to get deterministic names for GeoPandas.

    # ---- tiny utilities ----
    def _abs(self, p: str | Path) -> str:
        """Absolute POSIX path for ArcPy tools."""
        return Path(p).resolve().as_posix()

    def _set_ws_for(self, path_str: str) -> None:
        """Set arcpy.env.workspace appropriate to the input path."""
        if path_str.lower().endswith(".shp"):
            arcpy.env.workspace = Path(path_str).parent.as_posix()
        else:
            _set_workspace(self.cfg.updated_shapefiles)

    def _calc_geom_attrs(self, fc_path: str | Path, for_lines: bool = True) -> None:
        """Compute geometry attributes in place for a feature class or shapefile."""
        _require_arcpy()
        fc_abs = self._abs(fc_path)
        self._set_ws_for(fc_abs)
        if for_lines:
            arcpy.management.CalculateGeometryAttributes(
                fc_abs,
                [["Shape_Leng", "LENGTH"],
                 ["Angle", "LINE_BEARING"],
                 ["N_Curve", "CURVE_COUNT"],
                 ["Vertices", "POINT_COUNT"]]
            )

    def _to_shp_tmp(self, src_fc: str | Path, shp_basename: str) -> Path:
        """
        Export a feature class to a shapefile in CWD with the exact base name.
        Uses FeatureClassToFeatureClass so naming is deterministic.
        """
        _require_arcpy()
        src_abs = self._abs(src_fc)
        out_folder = Path(os.getcwd())
        out_name = Path(shp_basename).stem  # sanitize
        out_shp = out_folder / f"{out_name}.shp"
        # Clean any previous leftovers with the same name
        try:
            if out_shp.exists():
                arcpy.management.Delete(out_shp.as_posix())
        except Exception:
            pass
        self._set_ws_for(src_abs)
        # If already the exact shapefile, just return it
        if src_abs.lower().endswith(".shp") and Path(src_abs).name == out_shp.name and Path(src_abs).parent == out_folder:
            return Path(src_abs)
        # Export with deterministic name
        arcpy.conversion.FeatureClassToFeatureClass(src_abs, out_folder.as_posix(), out_name)
        return out_shp

    def _read_gdf(self, shp_path: Path):
        # GeoPandas is used for flexible DataFrame operations across sweeps
        import geopandas as gpd
        return gpd.read_file(shp_path)

    def _pick_final_candidate(self, cand_rows, prefer_order=("HAVE_THEIR_CENTER_IN", "LARGEST_OVERLAP", "CLOSEST")):
        """Choose first candidate with Pass_All == 1 by method preference."""
        if not cand_rows:
            return None
        passing = [r for r in cand_rows if r.get("Pass_All", 0) == 1]
        if not passing:
            return None
        for m in prefer_order:
            for r in passing:
                if r.get("Method") == m:
                    return r
        return passing[0]

    # Mapping retained from your original logic. Note that the line matcher
    # can also infer fields dynamically; this preserves prior expectations.
    _LINE_REQ_FIELDS = {
        "conduits": ['Name', 'FromNode', 'ToNode', 'Length', 'Roughness', 'InOffset', 'OutOffset',
                     'InitFlow', 'MaxFlow', 'Shape_1', 'Geom1', 'Geom2', 'Geom3', 'Geom4',
                     'Barrels', 'Culvert', 'Shp_Trnsct', 'Kentry', 'Kexit', 'Kavg', 'FlapGate', 'Seepage'],
        "orifices": ['Name', 'FromNode', 'ToNode', 'Type', 'InOffset', 'Qcoeff', 'FlapGate', 'CloseTime',
                     'Shape_1', 'Height', 'Width'],
        "outlets":  ['Name', 'FromNode', 'ToNode', 'InOffset', 'RateCurve', 'Qcoeff', 'Qexpon',
                     'FlapGate', 'CurveName'],
        "pumps":    ['Name', 'FromNode', 'ToNode', 'PumpCurve', 'Status', 'Startup', 'Shutoff'],
        "weirs":    ['Name', 'FromNode', 'ToNode', 'Type', 'CrestHeigh', 'Qcoeff', 'FlapGate',
                     'EndContrac', 'EndCoeff', 'Surcharge', 'RoadWidth', 'RoadSurf',
                     'CoeffCurve', 'Height', 'Length', 'SideSlope'],
    }

    def _ensure_fields(self, target_fc: str, src_fc: str, fields_to_copy: list[str]) -> list[str]:
        """
        Ensure target has fields to receive attributes copied from src_fc.
        Returns the list of fields that exist or were created in the target.
        """
        _require_arcpy()
        target_abs = self._abs(target_fc)
        src_abs = self._abs(src_fc)
        self._set_ws_for(target_abs)
        tgt_fields = {f.name.lower(): f for f in arcpy.ListFields(target_abs)}
        src_fields = {f.name: f for f in arcpy.ListFields(src_abs)}

        ensured = []
        for name in fields_to_copy:
            lname = name.lower()
            ensured.append(name)
            if lname in tgt_fields:
                continue
            # Create missing fields by inferring suitable ArcGIS types
            sf = src_fields.get(name)
            if sf is None:
                arcpy.management.AddField(target_abs, name, "TEXT", field_length=255)
                continue
            if sf.type in ("String", "Guid"):
                arcpy.management.AddField(target_abs, name, "TEXT", field_length=min(sf.length or 255, 255))
            elif sf.type in ("SmallInteger", "Integer"):
                arcpy.management.AddField(target_abs, name, "LONG")
            elif sf.type in ("Single", "Double"):
                arcpy.management.AddField(target_abs, name, "DOUBLE")
            elif sf.type in ("Date",):
                arcpy.management.AddField(target_abs, name, "DATE")
            else:
                arcpy.management.AddField(target_abs, name, "TEXT", field_length=255)
        return ensured

    # -------------------------------- POINTS --------------------------------
    def spatial_join_points_report(
        self,
        *,
        layer_key: str,
        radii_ft=(1, 2, 3, 4),
        xlsx_book: dict | None = None,
        apply_update: bool = True
    ) -> pd.DataFrame:
        """
        Grid-search radii for a point layer (Updated_<layer_key> vs Source_<layer_key>.shp).
        Writes summary to xlsx_book and, if apply_update, overwrites Updated_* with the chosen join.
        Choice heuristic: maximize unique 'Name' coverage and prefer smaller radius on ties.
        """
        _require_arcpy()
        _ensure_gdb(self.cfg.scratch_gdb)

        if not self.updated_found:
            _ = self.inventory_table()

        target = self.updated_found.get(f"Updated_{layer_key}")
        join   = (self.cfg.source_shapefiles / f"Source_{layer_key}.shp")

        if not target or not join.exists():
            self.log.warning(f"Skip points layer {layer_key} because inputs are missing")
            return pd.DataFrame()

        target_abs = self._abs(target)
        join_abs   = self._abs(join)

        try:
            tgt_count = int(arcpy.management.GetCount(target_abs).getOutput(0))
        except Exception:
            tgt_count = 0

        rows = []
        best_radius = None
        best_unique = -1

        for r in radii_ft:
            out_fc = (self.cfg.scratch_gdb / f"tmp_{layer_key}_spa_{int(r)}").as_posix()
            _safe_delete(out_fc)
            try:
                self._set_ws_for(target_abs)
                arcpy.analysis.SpatialJoin(
                    target_abs, join_abs, out_fc,
                    "JOIN_ONE_TO_ONE", "KEEP_ALL", "#",
                    "CLOSEST", f"{int(r)} Feet", "#"
                )
                if not arcpy.Exists(out_fc):
                    raise RuntimeError(f"SpatialJoin produced no output for radius {r} ft. {arcpy.GetMessages(2)}")
            except Exception as e:
                self.log.warning(f"Points SpatialJoin failed for {layer_key} at {r} ft: {e} | {arcpy.GetMessages(2)}")
                continue

            # Unique non-empty 'Name' values indicate successful matches
            uniq_names = set()
            if arcpy.ListFields(out_fc, "Name"):
                with arcpy.da.SearchCursor(out_fc, ["Name"]) as cur:
                    for (nm,) in cur:
                        if nm not in (None, "", " "):
                            uniq_names.add(nm)
            uniq = len(uniq_names)
            cov = round((uniq / tgt_count) * 100, 1) if tgt_count else 0.0

            rows.append({
                "Layer": layer_key,
                "GIS_Target": Path(target_abs).name,
                "SWMM_Join": Path(join_abs).name,
                "Test_Radius_ft": int(r),
                "Target_Count": tgt_count,
                "Unique_Names": uniq,
                "Coverage_pct": cov
            })

            # Keep the radius with the highest unique coverage, then smaller radius
            if uniq > best_unique or (uniq == best_unique and (best_radius is None or r < best_radius)):
                best_unique = uniq
                best_radius = r

        df = pd.DataFrame(rows)
        if not df.empty and best_radius is not None:
            df["Chosen"] = (df["Test_Radius_ft"] == int(best_radius)).astype(int)

            if apply_update:
                # Apply chosen result in place on the Updated_* target
                final_fc = (self.cfg.scratch_gdb / f"tmp_{layer_key}_final").as_posix()
                _safe_delete(final_fc)
                self._set_ws_for(target_abs)
                arcpy.analysis.SpatialJoin(
                    target_abs, join_abs, final_fc,
                    "JOIN_ONE_TO_ONE", "KEEP_ALL", "#",
                    "CLOSEST", f"{int(best_radius)} Feet", "#"
                )
                _safe_delete(target_abs)
                arcpy.management.CopyFeatures(final_fc, target_abs)
                self.log.info(f"[{layer_key}] Points SpatialJoin applied in-place at {best_radius} ft → {target_abs}")

        if xlsx_book is not None:
            xlsx_book[f"Points_Summary_{layer_key}"] = df

        return df

    # -------------------------------- LINES --------------------------------
    def _list_transfer_fields(self, src_fc: str | Path) -> list[str]:
        """Return attribute fields to transfer from Source_* (filters system/geom)."""
        _require_arcpy()
        src_abs = self._abs(src_fc)
        self._set_ws_for(src_abs)
        exclude_exact = {"objectid", "fid", "join_count", "target_fid"}
        exclude_prefix = {"oid_"}
        exclude_geom = {"shape", "shape_length", "shape_area"}  # keep 'Shape_1'
        fields = []
        for f in arcpy.ListFields(src_abs):
            n = f.name.lower()
            if f.type in ("Geometry", "OID", "Blob", "Raster"):  continue
            if n in exclude_exact:                               continue
            if any(n.startswith(p) for p in exclude_prefix):     continue
            if n in exclude_geom:                                continue
            fields.append(f.name)
        return fields

    def _detect_swmm_name_col(self, df) -> str | None:
        """
        Identify which 'Name' column in the Spatial Join output came from the
        JOIN side. Prefer suffixed variants over plain 'Name' to avoid mixing
        target-side names with SWMM names.
        """
        if df is None or df.empty:
            return None
        suffixed = [c for c in df.columns if c.lower().startswith("name_")]
        if suffixed:
            best = max(suffixed, key=lambda c: df[c].notna().sum())
            if df[best].notna().sum() > 0:
                return best
        plain = [c for c in df.columns if c.lower() == "name"]
        if plain:
            best = max(plain, key=lambda c: df[c].notna().sum())
            return best if df[best].notna().sum() > 0 else None
        return None

    def _angle_diff_allow_reverse(self, a, b) -> float:
        """Minimal angular difference in degrees, treating 180° reversal as 0°."""
        try:
            a = float(a) % 360.0
            b = float(b) % 360.0
        except Exception:
            return float("nan")
        d = abs(a - b)
        if d > 180.0:
            d = 360.0 - d
        return min(d, abs(d - 180.0))

    def _ensure_track_id(self, target_fc: str | Path, field_name: str = "SDART_TID") -> str:
        """Ensure a stable integer key exists on the target FC to align sweeps."""
        _require_arcpy()
        t_abs = self._abs(target_fc)
        self._set_ws_for(t_abs)
        present = {f.name for f in arcpy.ListFields(t_abs)}
        if field_name not in present:
            arcpy.management.AddField(t_abs, field_name, "LONG")
            oid = arcpy.Describe(t_abs).OIDFieldName
            arcpy.management.CalculateField(t_abs, field_name, f"!{oid}!", "PYTHON3")
        return field_name

    def spatial_join_lines_matrix(
            self,
            *,
            layer_key: str,  # 'conduits' 'orifices' 'outlets' 'pumps' 'weirs'
            radius_ft: int = 1,  # kept for API compatibility, sweep controls radii
            tolerances: dict | None = None,
            xlsx_book: dict | None = None,
            apply_update: bool = True
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Multi-method, multi-radius line matcher.
        Methods: CENTER/HAVE_THEIR_CENTER_IN, LARGEST_OVERLAP, CLOSEST.
        Radii searched: 1,5,10,20,30,50,100,150 ft plus unlimited for CLOSEST.
        Winner per target is chosen among geometry-pass rows, else best-scoring
        candidate. Enforces one-to-one by SWMM_Name, then applies attributes to
        Updated_* in place from the exact winning sweep.
        Returns three DataFrames: method summary, candidates, finals.
        """
        _require_arcpy()
        import math
        import pandas as pd

        # Tolerances are forgiving by default, with length absolute and relative
        tol = {"len": 10.0, "angle": 12.0, "curves": 10.0, "verts": 9999.0, "len_ratio": 0.30}
        if tolerances:
            tol.update(tolerances)

        if not self.updated_found:
            _ = self.inventory_table()

        target = self.updated_found.get(f"Updated_{layer_key}")
        join = (self.cfg.source_shapefiles / f"Source_{layer_key}.shp")
        if not target or not join.exists():
            self.log.warning(f"Skip line layer {layer_key} because inputs are missing")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

        target_abs = self._abs(target)
        join_abs = self._abs(join)

        # Stable key on target for cross-sweep correlation
        tid_field = self._ensure_track_id(target_abs, "SDART_TID")

        # Compute SWMM geometry once for scoring consistency
        self._calc_geom_attrs(join_abs, for_lines=True)
        swmm_shp = self._to_shp_tmp(join_abs, f"__swmm_{layer_key}")
        gdf_swmm = self._read_gdf(swmm_shp)
        keep_swmm = ["Name", "Shape_Leng", "Angle", "N_Curve", "Vertices"]
        for c in keep_swmm:
            if c not in gdf_swmm.columns:
                gdf_swmm[c] = None
        gdf_swmm = gdf_swmm[keep_swmm].copy()
        gdf_swmm["SWMM_OID"] = gdf_swmm.index.astype(int)

        methods = ["HAVE_THEIR_CENTER_IN", "LARGEST_OVERLAP", "CLOSEST"]
        labels = ["CENTER", "LARGEST", "CLOSEST"]
        radii_ft_list = [1, 5, 10, 20, 30, 50, 100, 150]  # searched radii

        method_counts = []
        candidates = []
        out_fcs = {}  # (label, radius_key) -> FC ; radius_key=-1 => unlimited

        # Sweep across methods × radii and record outputs and metrics
        for m, label in zip(methods, labels):
            radii_this = list(radii_ft_list)
            if m == "CLOSEST":
                radii_this.append(None)  # unlimited search

            for r in radii_this:
                out_fc = (
                            self.cfg.scratch_gdb / f"tmp_{layer_key}_{label.lower()}_{('all' if r is None else int(r))}").as_posix()
                _safe_delete(out_fc)
                try:
                    self._set_ws_for(target_abs)
                    search = "#" if r is None else f"{int(r)} Feet"
                    arcpy.analysis.SpatialJoin(
                        target_abs, join_abs, out_fc,
                        "JOIN_ONE_TO_ONE", "KEEP_ALL", "#",
                        m, search, "#"
                    )
                    if not arcpy.Exists(out_fc):
                        raise RuntimeError(
                            f"SpatialJoin produced no output for {label} @ {search}. {arcpy.GetMessages(2)}")
                    self._calc_geom_attrs(out_fc, for_lines=True)
                except Exception as e:
                    self.log.warning(
                        f"Line SpatialJoin failed for {layer_key} {label}@{r}: {e} | {arcpy.GetMessages(2)}")
                    continue

                radius_key = -1 if r is None else int(r)  # ranking; -1 => unlimited (worst)
                radius_for_report = 0 if r is None else int(r)  # 0 only for reporting
                out_fcs[(label, radius_key)] = out_fc

                shp = self._to_shp_tmp(out_fc, f"__{layer_key}_{label.lower()}_{('all' if r is None else int(r))}")
                gdf = self._read_gdf(shp)

                swmm_name_col = self._detect_swmm_name_col(gdf)

                # Per-sweep summary for Excel reporting
                uniq = int(gdf[swmm_name_col].nunique()) if swmm_name_col and swmm_name_col in gdf.columns else 0
                method_counts.append({
                    "SWMM_Layer": f"Updated_{layer_key}",
                    "GIS_Target": Path(target_abs).name,
                    "Radius_ft": radius_for_report,
                    label: uniq
                })

                # Left side: GIS geometry and stable key
                left_cols = [tid_field, "Shape_Leng", "Angle", "N_Curve", "Vertices"]
                for c in left_cols:
                    if c not in gdf.columns:
                        gdf[c] = None
                left = gdf[left_cols].copy().rename(columns={
                    tid_field: "GIS_OID",
                    "Shape_Leng": "GIS_Shape_Leng",
                    "Angle": "GIS_Angle",
                    "N_Curve": "GIS_N_Curve",
                    "Vertices": "GIS_Vertices"
                })
                left["SWMM_Name"] = gdf[swmm_name_col].astype(
                    object) if swmm_name_col and swmm_name_col in gdf.columns else None

                # Right side: SWMM geometry for comparison
                right = gdf_swmm.rename(columns={
                    "Shape_Leng": "SWMM_Shape_Leng",
                    "Angle": "SWMM_Angle",
                    "N_Curve": "SWMM_N_Curve",
                    "Vertices": "SWMM_Vertices"
                })
                merged = left.merge(
                    right[["Name", "SWMM_Shape_Leng", "SWMM_Angle", "SWMM_N_Curve", "SWMM_Vertices", "SWMM_OID"]],
                    how="left", left_on="SWMM_Name", right_on="Name"
                )

                # Compute deltas and pass flags; treat both-NaN as pass for that metric
                merged["Len_Diff"] = (merged["GIS_Shape_Leng"] - merged["SWMM_Shape_Leng"]).abs()
                merged["Angle_Diff"] = merged.apply(
                    lambda z: self._angle_diff_allow_reverse(z["GIS_Angle"], z["SWMM_Angle"]), axis=1)
                merged["Curves_Diff"] = (merged["GIS_N_Curve"] - merged["SWMM_N_Curve"]).abs()
                merged["Verts_Diff"] = (merged["GIS_Vertices"] - merged["SWMM_Vertices"]).abs()

                both_nan = merged["GIS_N_Curve"].isna() & merged["SWMM_N_Curve"].isna()
                merged.loc[both_nan, "Curves_Diff"] = 0.0
                both_nan = merged["GIS_Vertices"].isna() & merged["SWMM_Vertices"].isna()
                merged.loc[both_nan, "Verts_Diff"] = 0.0
                both_nan = merged["GIS_Angle"].isna() & merged["SWMM_Angle"].isna()
                merged.loc[both_nan, "Angle_Diff"] = 0.0

                # Length pass supports absolute and relative thresholds
                len_abs = merged["Len_Diff"]
                denom = (merged[["GIS_Shape_Leng", "SWMM_Shape_Leng"]].abs()).max(axis=1).replace(0, 1e-9)
                len_rel = (len_abs / denom).fillna(1.0)
                merged["Pass_Len"] = ((len_abs <= tol["len"]) | (len_rel <= tol["len_ratio"])).astype(int)
                merged["Pass_Angle"] = (merged["Angle_Diff"] <= tol["angle"]).fillna(False).astype(int)
                merged["Pass_Curves"] = (merged["Curves_Diff"] <= tol["curves"]).fillna(False).astype(int)
                merged["Pass_Verts"] = (merged["Verts_Diff"] <= tol["verts"]).fillna(False).astype(int)
                merged["Pass_All"] = (
                            merged[["Pass_Len", "Pass_Angle", "Pass_Curves", "Pass_Verts"]].sum(axis=1) == 4).astype(
                    int)

                # Sweep metadata for later selection and tie-breaking
                merged["Method"] = m
                merged["Method_lbl"] = label
                merged["Radius_key"] = radius_key  # int; -1 for unlimited
                merged["Radius_rep"] = radius_for_report  # for reporting only

                # Collect candidates for scoring
                for _, row in merged.iterrows():
                    candidates.append({
                        "GIS_OID": int(row["GIS_OID"]) if pd.notna(row.get("GIS_OID")) else None,
                        "SWMM_OID": int(row["SWMM_OID"]) if pd.notna(row.get("SWMM_OID")) else None,
                        "SWMM_Name": row.get("SWMM_Name"),
                        "Method": m,
                        "Method_lbl": label,
                        "Radius_key": radius_key,
                        "Radius_rep": int(row.get("Radius_rep", 0)) if pd.notna(row.get("Radius_rep")) else 0,
                        "Len_Diff": float(row["Len_Diff"]) if pd.notna(row.get("Len_Diff")) else math.inf,
                        "Angle_Diff": float(row["Angle_Diff"]) if pd.notna(row.get("Angle_Diff")) else math.inf,
                        "Curves_Diff": float(row["Curves_Diff"]) if pd.notna(row.get("Curves_Diff")) else math.inf,
                        "Verts_Diff": float(row["Verts_Diff"]) if pd.notna(row.get("Verts_Diff")) else math.inf,
                        "Pass_All": int(row["Pass_All"]) if pd.notna(row.get("Pass_All")) else 0
                    })

                # Clean shapefile view for this sweep
                try:
                    arcpy.management.Delete(shp.as_posix())
                except Exception:
                    pass

        # Build summary and candidate frames for reporting and selection
        df_counts = pd.DataFrame(method_counts)
        if not df_counts.empty:
            df_counts = df_counts.groupby(["SWMM_Layer", "GIS_Target", "Radius_ft"], as_index=False).sum(
                numeric_only=True)

        df_cand = pd.DataFrame(candidates)
        if df_cand.empty:
            df_final = pd.DataFrame()
            if xlsx_book is not None:
                xlsx_book[f"Lines_Method_Summary_{layer_key}"] = df_counts
                xlsx_book[f"Lines_Candidates_{layer_key}"] = df_cand
                xlsx_book[f"Lines_Final_{layer_key}"] = df_final
            try:
                arcpy.management.Delete(swmm_shp.as_posix())
            except Exception:
                pass
            return df_counts, df_cand, df_final

        # ---------- pick winners per target (fallback to best if none pass) ----------
        finals = []
        best_pick_for_target = {}  # GIS_OID -> (Method_lbl, radius_key)
        method_rank = {"HAVE_THEIR_CENTER_IN": 0, "LARGEST_OVERLAP": 1, "CLOSEST": 2}
        eps = 1e-9

        # Frequency helps tie-breaking when the same pair is seen across sweeps
        freq_all = (df_cand[df_cand["SWMM_Name"].notna()]
                    .groupby(["GIS_OID", "SWMM_Name"]).size()
                    .rename("Seen_Count").reset_index())

        for gid, sub in df_cand.groupby("GIS_OID", dropna=True):
            sub = sub.copy()
            sub["Score"] = (
                    (sub["Len_Diff"] / max(tol["len"], eps)) ** 2 +
                    (sub["Angle_Diff"] / max(tol["angle"], eps)) ** 2 +
                    (sub["Curves_Diff"] / max(tol["curves"], eps)) ** 2 +
                    (sub["Verts_Diff"] / max(tol["verts"], eps)) ** 2
            )
            sub["_mrank"] = sub["Method"].map(method_rank).fillna(99)
            sub["_rpref"] = sub["Radius_key"].apply(lambda k: 1_000_000 if int(k) == -1 else int(k))
            sub = sub.merge(freq_all, on=["GIS_OID", "SWMM_Name"], how="left")
            sub["Seen_Count"] = sub["Seen_Count"].fillna(0)

            pool = sub[(sub["Pass_All"] == 1) & sub["SWMM_Name"].notna()]
            if pool.empty:
                pool = sub[sub["SWMM_Name"].notna()]
            if pool.empty:
                continue

            pick = pool.sort_values(by=["Score", "_rpref", "_mrank", "Seen_Count"],
                                    ascending=[True, True, True, False]).iloc[0].to_dict()
            finals.append({
                "GIS_OID": int(gid),
                "SWMM_Name": pick.get("SWMM_Name"),
                "Chosen_Method": pick.get("Method"),
                "Radius_ft": int(pick.get("Radius_rep", 0)),
                "Score": float(pick.get("Score", 0.0))
            })
            best_pick_for_target[int(gid)] = (pick.get("Method_lbl"), int(pick.get("Radius_key", -1)))

        # Mark winners in the candidate table for transparency
        if finals:
            keyset = {(int(r["GIS_OID"]), r["SWMM_Name"], best_pick_for_target[int(r["GIS_OID"])][0],
                       best_pick_for_target[int(r["GIS_OID"])][1]) for r in finals}

            def _is_final(x):
                if pd.isna(x.get("GIS_OID")) or pd.isna(x.get("SWMM_Name")):
                    return 0
                tup = (int(x["GIS_OID"]), x["SWMM_Name"], x["Method_lbl"], int(x["Radius_key"]))
                return 1 if tup in keyset else 0

            df_cand["Final_Choice"] = df_cand.apply(_is_final, axis=1)
        else:
            df_cand["Final_Choice"] = 0

        df_final = pd.DataFrame(finals)

        # ---------- enforce one-to-one by SWMM_Name (deduplicate winners) ----------
        losers_for_clear = set()
        if not df_final.empty and "SWMM_Name" in df_final.columns:
            def _rad_pref(r):
                rf = int(r.get("Radius_ft", 0))
                return 1_000_000 if rf == 0 else rf  # 0 = unlimited (worst)
            for swmm_name, grp in df_final.groupby("SWMM_Name"):
                if pd.isna(swmm_name):
                    continue
                if len(grp) <= 1:
                    continue
                grp = grp.copy()
                grp["_rp"] = grp.apply(_rad_pref, axis=1)
                grp["_mr"] = grp["Chosen_Method"].map(method_rank).fillna(99)
                best_row = grp.sort_values(by=["Score", "_rp", "_mr"], ascending=[True, True, True]).iloc[0]
                keep_gid = int(best_row["GIS_OID"])
                for gid in grp["GIS_OID"]:
                    gid = int(gid)
                    if gid != keep_gid:
                        losers_for_clear.add(gid)
                        if gid in best_pick_for_target:
                            del best_pick_for_target[gid]
            if losers_for_clear:
                self.log.info(
                    f"[{layer_key}] De-duplicated {len(losers_for_clear)} features mapped to the same SWMM element; will clear them to NULL.")

        # ---- Apply attributes in place from the exact winning sweep output ----
        if apply_update and (best_pick_for_target or losers_for_clear) and out_fcs:
            fields_to_copy = self._list_transfer_fields(join_abs)
            ensured = self._ensure_fields(target_abs, join_abs, fields_to_copy)

            # Helper to map joined field names that ArcGIS may suffix back to base names
            def pick_join_field_name(present_names_set: set[str], base: str) -> str | None:
                base_l = base.lower()
                present_l = {n.lower(): n for n in present_names_set}
                if base_l in present_l:
                    return present_l[base_l]
                for suf in ("_1", "_2", "_r", "_j", "_join", "_join1"):
                    cand = base_l + suf
                    if cand in present_l:
                        return present_l[cand]
                return None

            # Build lookup tables indexed by SDART_TID for the chosen sweep per winner
            lookups = {}
            for (label, rkey), fc in out_fcs.items():
                present = {f.name for f in arcpy.ListFields(fc)}
                if "SDART_TID" not in present:
                    continue  # lost the key; skip this sweep output
                chosen_out_fields, base_cols = [], []
                for base in ensured:
                    out_name = pick_join_field_name(present, base)
                    if out_name:
                        chosen_out_fields.append(out_name)
                        base_cols.append(base)
                use_fields = ["SDART_TID"] + chosen_out_fields
                data = {}
                if len(use_fields) > 1:
                    with arcpy.da.SearchCursor(fc, use_fields) as cur:
                        for row in cur:
                            key = row[0]
                            if key is None:
                                continue
                            d = {}
                            for i, base in enumerate(base_cols, start=1):
                                d[base] = row[i]
                            data[int(key)] = d
                lookups[(label, rkey)] = data

            update_fields = ["OBJECTID"] + ensured
            updated_rows = 0
            cleared_rows = 0
            self._set_ws_for(target_abs)
            with arcpy.da.UpdateCursor(target_abs, update_fields) as ucur:
                for row in ucur:
                    oid = int(row[0])

                    # Losers in dedup step → clear attributes to NULL
                    if oid in losers_for_clear:
                        changed = False
                        for i in range(1, len(row)):
                            if row[i] is not None:
                                row[i] = None
                                changed = True
                        if changed:
                            ucur.updateRow(row)
                            cleared_rows += 1
                        continue

                    # Winners → copy values from the selected sweep
                    choice = best_pick_for_target.get(oid)  # (Method_lbl, radius_key)
                    if not choice:
                        continue
                    srcvals = lookups.get(choice, {}).get(oid)
                    if not srcvals:
                        continue
                    changed = False
                    for i, fld in enumerate(ensured, start=1):
                        val = srcvals.get(fld)
                        if row[i] != val:
                            row[i] = val
                            changed = True
                    if changed:
                        ucur.updateRow(row)
                        updated_rows += 1

            self.log.info(
                f"[{layer_key}] Line attributes applied in-place from sweep winners → {target_abs} "
                f"(updated {updated_rows} features; cleared {cleared_rows} duplicates)"
            )

        # Cleanup temp SWMM shapefile
        try:
            arcpy.management.Delete(swmm_shp.as_posix())
        except Exception:
            pass

        if xlsx_book is not None:
            xlsx_book[f"Lines_Method_Summary_{layer_key}"] = df_counts
            xlsx_book[f"Lines_Candidates_{layer_key}"] = df_cand
            xlsx_book[f"Lines_Final_{layer_key}"] = df_final

        return df_counts, df_cand, df_final

    # -------------------------------- PACKAGE RUNNER --------------------------------
    def spatial_join_package_to_excel(
        self,
        *,
        points_layers=("junctions","storages","dividers","outfalls"),
        line_layers=("conduits","orifices","outlets","pumps","weirs"),
        points_radii=(1,2,3,4),
        line_radius=1,
        tolerances=None,
        xlsx_path: Path | None = None,
        apply_updates: bool = True
    ) -> Path:
        """
        Run the full package:
        • Points grid search and apply chosen join in place
        • Lines method×radius matrix, pick winners, enforce one-to-one, apply attributes
        • Write SDART_Sjoin_Result.xlsx with multiple sheets and a Run_Metadata tab
        Always writes metadata so the workbook is valid even if others are empty.
        """
        _require_arcpy()
        from datetime import datetime

        if xlsx_path is None:
            xlsx_path = Path(os.getcwd()) / "SDART_Sjoin_Result.xlsx"

        book: dict[str, pd.DataFrame] = {}
        run_notes = []

        # Points sweep per layer with coverage-based selection
        for k in points_layers:
            try:
                df = self.spatial_join_points_report(
                    layer_key=k, radii_ft=points_radii, xlsx_book=book, apply_update=apply_updates
                )
                if df.empty:
                    run_notes.append(f"Points layer {k} produced no rows")
            except Exception as e:
                msg = f"Points layer {k} failed in report: {e}"
                self.log.warning(msg)
                run_notes.append(msg)

        # Lines sweep per layer with scoring and deduplication
        for k in line_layers:
            try:
                cnts, cands, final = self.spatial_join_lines_matrix(
                    layer_key=k, radius_ft=line_radius, tolerances=tolerances, xlsx_book=book, apply_update=apply_updates
                )
                if cnts.empty and cands.empty and final.empty:
                    run_notes.append(f"Line layer {k} produced no rows")
            except Exception as e:
                msg = f"Line layer {k} failed in matrix: {e}"
                self.log.warning(msg)
                run_notes.append(msg)

        # Metadata sheet for reproducibility and audit
        meta = pd.DataFrame({
            "When": [datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
            "Updated_GDB": [self.cfg.updated_shapefiles.as_posix()],
            "Source_Folder": [self.cfg.source_shapefiles.as_posix()],
            "Scratch_GDB": [self.cfg.scratch_gdb.as_posix()],
            "Points_Layers": [", ".join(points_layers) if points_layers else "<none>"],
            "Line_Layers": [", ".join(line_layers) if line_layers else "<none>"],
            "Notes": ["; ".join(run_notes) if run_notes else "OK"]
        })
        book["Run_Metadata"] = meta

        with pd.ExcelWriter(xlsx_path.as_posix()) as xw:
            for sheet, df in book.items():
                df.to_excel(xw, index=False, sheet_name=sheet[:31])

        self.log.info(f"Spatial Join workbook written → {xlsx_path}")
        return xlsx_path

    # ------------------------------------------------------------------
    # Name hygiene — ensure fields exist and fill NULL/blank Names
    # ------------------------------------------------------------------
    def _ensure_text_field(self, fc_path: str | Path, field_name: str, length: int = 64) -> None:
        _require_arcpy()
        fc_abs = self._abs(fc_path)
        self._set_ws_for(fc_abs)
        fields = {f.name: f for f in arcpy.ListFields(fc_abs)}
        if field_name not in fields:
            arcpy.management.AddField(fc_abs, field_name, "TEXT", field_length=min(length, 255))
        elif getattr(fields[field_name], "type", "").lower() != "string":
            self.log.warning(f"{Path(fc_abs).name}: field '{field_name}' exists but is not TEXT.")

    def _ensure_long_field(self, fc_path: str | Path, field_name: str) -> None:
        _require_arcpy()
        fc_abs = self._abs(fc_path)
        self._set_ws_for(fc_abs)
        if not arcpy.ListFields(fc_abs, field_name):
            arcpy.management.AddField(fc_abs, field_name, "LONG")

    # Fill placeholders for missing Names on points and optionally conduits
    def fill_null_names_unknown(
        self,
        *,
        include_points: bool = True,
        include_conduits: bool = True,
        use_objectid: bool = True,
        points_layers: tuple[str, ...] | None = None,   # optional: limit to specific point layers
        line_layers: tuple[str, ...] | None = None,     # optional: limit to specific line layers
    ) -> None:
        """
        Fills NULL/blank Name for Updated_*:
          points: junctions, storages, dividers, outfalls
          lines:  conduits only
        Names become 'u_<type>_<n>' with <n>=OBJECTID if use_objectid=True, else sequential.
        """
        _require_arcpy()
        if not self.updated_found:
            _ = self.inventory_table()

        default_points = ("junctions","storages","dividers","outfalls")
        default_lines  = ("conduits",)

        p_layers = points_layers if (include_points and points_layers) else (default_points if include_points else ())
        l_layers = line_layers  if (include_conduits and line_layers)  else (default_lines  if include_conduits else ())

        def _prefix(k: str) -> str:
            return {
                "junctions": "u_junction",
                "storages":  "u_storage",
                "dividers":  "u_divider",
                "outfalls":  "u_outfall",
                "conduits":  "u_conduit",
            }[k]

        def _fill_for(layer_key: str) -> None:
            fc = self.updated_found.get(f"Updated_{layer_key}")
            if not fc:
                return
            fc_abs = self._abs(fc)
            self._set_ws_for(fc_abs)
            self._ensure_text_field(fc_abs, "Name", 64)

            oidf = arcpy.Describe(fc_abs).OIDFieldName
            pref = _prefix(layer_key)

            # Sequential mode continues after the highest existing placeholder
            next_idx = 1
            if not use_objectid:
                existing = []
                with arcpy.da.SearchCursor(fc_abs, ["Name"]) as cur:
                    for (nm,) in cur:
                        if isinstance(nm, str) and nm.startswith(pref + "_"):
                            try:
                                existing.append(int(nm.split("_")[-1]))
                            except Exception:
                                pass
                if existing:
                    next_idx = max(existing) + 1

            filled = 0
            with arcpy.da.UpdateCursor(fc_abs, [oidf, "Name"]) as ucur:
                for oid, nm in ucur:
                    if nm in (None, "", " "):
                        new_nm = f"{pref}_{int(oid)}" if use_objectid else f"{pref}_{next_idx}"
                        if not use_objectid:
                            next_idx += 1
                        ucur.updateRow([oid, new_nm])
                        filled += 1
            self.log.info(f"[{layer_key}] Filled {filled} NULL Names → {fc_abs}")

        for k in p_layers:
            _fill_for(k)
        for k in l_layers:
            _fill_for(k)

    # ------------------------------------------------------------------
    # From/To assignment — nearest nodes within a search radius
    # ------------------------------------------------------------------
    def assign_from_to_nodes(
        self,
        *,
        line_layers=("conduits","orifices","outlets","pumps","weirs"),
        node_layers=("junctions","storages","dividers","outfalls"),
        search_radius_ft: float = 2.0
    ) -> None:
        """
        For each Updated_<line_layer>:
        • Create start and end points
        • Find nearest node among merged Updated_<node_layers> within search_radius_ft
        • Assign FromNode and ToNode text fields accordingly
        Leaves NULL if no node is found within the radius.
        """
        _require_arcpy()
        _ensure_gdb(self.cfg.scratch_gdb)

        # Build merged node candidate set; ensure Names exist for every node FC
        if not self.updated_found:
            _ = self.inventory_table()

        inputs = []
        for nk in node_layers:
            fc = self.updated_found.get(f"Updated_{nk}")
            if not fc:
                continue
            # Fill placeholders only for this node layer if it has blanks
            self.fill_null_names_unknown(include_points=True, include_conduits=False, points_layers=(nk,))
            inputs.append(self._abs(fc))

        if not inputs:
            self.log.warning("No node layers found for From/To assignment.")
            return

        merged_nodes = (self.cfg.scratch_gdb / "SDART_Merged_Nodes").as_posix()
        _safe_delete(merged_nodes)
        arcpy.management.Merge(inputs, merged_nodes)
        self._ensure_text_field(merged_nodes, "Name", 64)

        # Build OID -> Name lookup for the merged nodes
        node_oid = arcpy.Describe(merged_nodes).OIDFieldName
        node_name_lookup = {}
        with arcpy.da.SearchCursor(merged_nodes, [node_oid, "Name"]) as cur:
            for oid, nm in cur:
                node_name_lookup[int(oid)] = nm

        radius = f"{float(search_radius_ft)} Feet"

        def _nearest_names(pts_fc: str, label: str) -> dict[int, str | None]:
            """Run Near(pts_fc -> merged_nodes) and return ORIG_FID -> Name."""
            pts_copy = (self.cfg.scratch_gdb / f"{Path(pts_fc).name}_near_{label}").as_posix()
            _safe_delete(pts_copy)
            arcpy.management.CopyFeatures(pts_fc, pts_copy)

            arcpy.analysis.Near(pts_copy, merged_nodes, radius, "NO_LOCATION", "NO_ANGLE", "PLANAR")

            best: dict[int, tuple[float, int | None]] = {}
            with arcpy.da.SearchCursor(pts_copy, ["ORIG_FID", "NEAR_FID", "NEAR_DIST"]) as cur:
                for ofid, nfid, ndist in cur:
                    if ofid is None:
                        continue
                    if nfid in (None, -1):
                        if ofid not in best:
                            best[ofid] = (float("inf"), None)
                        continue
                    dist = float(ndist) if ndist is not None else float("inf")
                    prev = best.get(ofid, (float("inf"), None))
                    if dist < prev[0]:
                        best[ofid] = (dist, int(nfid))

            out = {}
            for ofid, (_, nfid) in best.items():
                out[int(ofid)] = node_name_lookup.get(int(nfid)) if nfid is not None else None

            _safe_delete(pts_copy)
            return out

        # Process each requested line layer independently
        for lk in line_layers:
            lines = self.updated_found.get(f"Updated_{lk}")
            if not lines:
                continue

            lines_abs = self._abs(lines)
            self._set_ws_for(lines_abs)
            self._ensure_text_field(lines_abs, "FromNode", 64)
            self._ensure_text_field(lines_abs, "ToNode", 64)
            self._ensure_long_field(lines_abs, "SDART_TID")

            oidf = arcpy.Describe(lines_abs).OIDFieldName
            arcpy.management.CalculateField(lines_abs, "SDART_TID", f"!{oidf}!", "PYTHON3")

            start_fc = (self.cfg.scratch_gdb / f"SDART_{lk}_start").as_posix()
            end_fc   = (self.cfg.scratch_gdb / f"SDART_{lk}_end").as_posix()
            _safe_delete(start_fc); _safe_delete(end_fc)
            arcpy.management.FeatureVerticesToPoints(lines_abs, start_fc, "START")
            arcpy.management.FeatureVerticesToPoints(lines_abs, end_fc,   "END")

            from_map = _nearest_names(start_fc, "from")
            to_map   = _nearest_names(end_fc, "to")

            updated = 0
            with arcpy.da.UpdateCursor(lines_abs, [oidf, "FromNode", "ToNode"]) as ucur:
                for oid, fnode, tnode in ucur:
                    new_f = from_map.get(int(oid), fnode)
                    new_t = to_map.get(int(oid), tnode)
                    if new_f != fnode or new_t != tnode:
                        ucur.updateRow([oid, new_f, new_t])
                        updated += 1

            _safe_delete(start_fc); _safe_delete(end_fc)
            self.log.info(f"[{lk}] From/To assigned within {radius} → {lines_abs} (updated {updated} features)")

# ------------------------------------------------------------------------------
# End of SDART
# ------------------------------------------------------------------------------
