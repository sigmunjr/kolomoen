#!/usr/bin/env python3
"""Lowest terrain passage in a GeoTIFF, with a token-free CesiumJS viewer.

Scientific objective: minimize the largest raster elevation on a connected
path. Among paths meeting that exact ceiling, minimize planimetric grid length.
An optional, explicitly heuristic, valley-biased representative is also saved.

Usage:
  python lowest_passage.py 6706_3_10m_z32.tif --out results
  python lowest_passage.py dem.tif --start 60.82 11.63 --end 60.70 11.27
  python lowest_passage.py dem.tif --via LAT LON --search-buffer-km 5
  python lowest_passage.py --self-test

Dependencies: numpy scipy numba rasterio pyproj matplotlib pillow
All scientific calculations use the original raster, not the web-view grid.
Heights retain the input vertical datum. No vertical datum is guessed.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import heapq
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from numba import njit
from pyproj import CRS, Transformer
from rasterio.features import shapes
from rasterio.transform import from_bounds, rowcol, xy
from rasterio.warp import reproject, Resampling, transform_bounds, transform_geom
from rasterio.windows import Window, from_bounds as window_from_bounds
from scipy import ndimage as ndi
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_START = (60.820232979269655, 11.635771799665232)
DEFAULT_END = (60.69954161460697, 11.265232590579124)
FOUR = ndi.generate_binary_structure(2, 1)
EIGHT = ndi.generate_binary_structure(2, 2)
VERSION = "1.0.0"


def log(message: str) -> None:
    print(message, flush=True)


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def write_standalone_viewer(path: Path, html: str, output_dir: Path,
                            data: dict, plateaus: dict) -> None:
    """Embed every local viewer asset so the HTML can be opened via file://."""
    terrain = data["terrain"]
    image_files = ["hillshade.png", *(level["file"] for level in terrain["levels"])]
    payload = dict(
        data=data,
        plateaus=plateaus,
        terrain_base64=base64.b64encode((output_dir / terrain["raw_file"]).read_bytes()).decode("ascii"),
        images={name: "data:image/png;base64," + base64.b64encode((output_dir / name).read_bytes()).decode("ascii")
                for name in image_files},
    )
    embedded = json.dumps(payload, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    bootstrap = "<script>window.PASSAGE_EMBEDDED=" + embedded.replace("<", "\\u003c") + ";</script>\n"
    marker = '<script>\n"use strict";'
    if marker not in html:
        raise ValueError("Cannot find the viewer script insertion point")
    path.write_text(html.replace(marker, bootstrap + marker, 1), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def labels_at(z: np.ndarray, level: float, diagonal_only: bool = False):
    """4-connectivity is also the topology of 8-neighbours without corner cutting."""
    return ndi.label(np.isfinite(z) & (z <= np.float32(level)),
                     structure=EIGHT if diagonal_only else FOUR)[0]


def anchors_connected(labels: np.ndarray, anchors: list[tuple[int, int]]) -> bool:
    first = labels[anchors[0]]
    return bool(first != 0 and all(labels[p] == first for p in anchors[1:]))


def minimum_ceiling(z: np.ndarray, anchors: list[tuple[int, int]], diagonal_only=False) -> float:
    """Exact float32 ceiling, found by monotone connectivity bisection.

    Termination is between adjacent representable float32 values, not an
    arbitrary decimal tolerance. The selected ceiling is an actual input value.
    """
    finite = np.isfinite(z)
    if not finite.any() or not all(finite[p] for p in anchors):
        raise ValueError("An endpoint/via point falls on NoData; it is not silently moved.")
    lo = np.float32(max(float(z[p]) for p in anchors))
    if anchors_connected(labels_at(z, lo, diagonal_only), anchors):
        return float(lo)
    hi = np.float32(z[finite].max())
    if not anchors_connected(labels_at(z, hi, diagonal_only), anchors):
        raise ValueError("No path exists inside the valid raster, even at its maximum elevation.")
    for _ in range(256):
        mid = np.float32((float(lo) + float(hi)) * 0.5)
        if mid <= lo or mid >= hi:
            break
        if anchors_connected(labels_at(z, mid, diagonal_only), anchors):
            hi = mid
        else:
            lo = mid
    else:
        raise RuntimeError("Float32 ceiling search did not converge.")
    # Independent proof certificate at the two sides of the selected level.
    below = np.nextafter(hi, np.float32(-np.inf))
    assert anchors_connected(labels_at(z, hi, diagonal_only), anchors)
    assert not anchors_connected(labels_at(z, below, diagonal_only), anchors)
    assert np.any(z == hi)
    return float(hi)


@njit(cache=True)
def _astar(z, level, start, goal, dx, dy, valley_weight=0.0, z_floor=0.0):
    """Shortest admissible 8-neighbour path; both side cells must be admissible.

    A* heuristic is the exact obstacle-free rectangular-grid metric and remains
    admissible with nonnegative valley penalties. No smoothing of the solution.
    """
    nr, nc = z.shape
    n = nr * nc
    distance = np.full(n, np.inf, np.float64)
    parent = np.full(n, -1, np.int32)
    closed = np.zeros(n, np.uint8)
    si = start[0] * nc + start[1]
    ti = goal[0] * nc + goal[1]
    distance[si] = 0.0
    queue = [(0.0, 0.0, si)]
    diagonal = math.hypot(dx, dy)
    denominator = max(level - z_floor, 1e-6)
    expanded = 0
    while queue:
        _, g, i = heapq.heappop(queue)
        if closed[i] or g > distance[i]:
            continue
        closed[i] = 1
        expanded += 1
        if i == ti:
            break
        r, c = i // nc, i % nc
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                if dr == 0 and dc == 0:
                    continue
                rr, cc = r + dr, c + dc
                if rr < 0 or cc < 0 or rr >= nr or cc >= nc:
                    continue
                zz = z[rr, cc]
                if not np.isfinite(zz) or zz > level:
                    continue
                if dr != 0 and dc != 0:
                    if (not np.isfinite(z[r, cc]) or not np.isfinite(z[rr, c])
                            or z[r, cc] > level or z[rr, c] > level):
                        continue
                j = rr * nc + cc
                if closed[j]:
                    continue
                step = math.hypot(dc * dx, dr * dy)
                if valley_weight > 0:
                    relative_z = max(0., (max(z[r, c], zz) - z_floor) / denominator)
                    step *= 1.0 + valley_weight * relative_z ** 4
                candidate = g + step
                if candidate < distance[j]:
                    distance[j] = candidate
                    parent[j] = i
                    hr, hc = abs(rr - goal[0]), abs(cc - goal[1])
                    shared = min(hr, hc)
                    heuristic = shared * diagonal + (hr-shared)*dy + (hc-shared)*dx
                    heapq.heappush(queue, (candidate + heuristic, candidate, j))
    if not np.isfinite(distance[ti]):
        raise ValueError("A* found no admissible path.")
    reverse_path = [ti]
    i = ti
    while i != si:
        i = parent[i]
        if i < 0:
            raise RuntimeError("Broken predecessor chain.")
        reverse_path.append(i)
    return np.array(reverse_path[::-1], np.int64), distance[ti], expanded


def solve_path(z, level, anchors, dx, dy, weight=0.0):
    pieces, total_cost, expanded = [], 0.0, 0
    floor = float(z[np.isfinite(z)].min())
    for i, (s, t) in enumerate(zip(anchors, anchors[1:])):
        p, cost, count = _astar(z, level, s, t, dx, dy, weight, floor)
        pieces.append(p if i == 0 else p[1:])
        total_cost += cost
        expanded += count
    path = np.concatenate(pieces)
    rr, cc = np.divmod(path, z.shape[1])
    assert float(z[rr, cc].max()) <= level
    # Explicit no-corner-cutting certificate for every diagonal edge.
    for r, c, r2, c2 in zip(rr[:-1], cc[:-1], rr[1:], cc[1:]):
        assert abs(r2-r) <= 1 and abs(c2-c) <= 1
        if r2 != r and c2 != c:
            assert np.isfinite(z[r, c2]) and np.isfinite(z[r2, c])
            assert z[r, c2] <= level and z[r2, c] <= level
    return rr, cc, float(total_cost), int(expanded)


class Grid:
    def __init__(self, source: Path, points: list[tuple[float, float]], buffer_km, max_cells):
        self.source = source
        with rasterio.open(source) as ds:
            if ds.count != 1:
                raise ValueError("Expected one elevation band.")
            if not ds.crs or not ds.crs.is_projected:
                raise ValueError("Input must have a projected CRS in metres; reproject first.")
            crs = CRS.from_user_input(ds.crs)
            if any(abs(a.unit_conversion_factor - 1.) > 1e-9 for a in crs.axis_info[:2]):
                raise ValueError("Projected horizontal units must be metres.")
            if ds.transform.b != 0 or ds.transform.d != 0 or ds.transform.a <= 0 or ds.transform.e >= 0:
                raise ValueError("Expected a north-up, unrotated raster.")
            if ds.units[0] and ds.units[0].lower() not in ("m", "meter", "metre", "meters", "metres"):
                raise ValueError(f"Unsupported vertical unit: {ds.units[0]}")
            self.crs = ds.crs
            self.forward = Transformer.from_crs(4326, ds.crs, always_xy=True)
            self.inverse = Transformer.from_crs(ds.crs, 4326, always_xy=True)
            points_xy = [self.forward.transform(lon, lat) for lat, lon in points]
            w = Window(0, 0, ds.width, ds.height)
            if buffer_km is not None:
                margin = 1000. * buffer_km
                xs, ys = np.array(points_xy).T
                raw = window_from_bounds(xs.min()-margin, ys.min()-margin,
                                         xs.max()+margin, ys.max()+margin, ds.transform)
                c0 = max(0, math.floor(raw.col_off)); r0 = max(0, math.floor(raw.row_off))
                c1 = min(ds.width, math.ceil(raw.col_off+raw.width))
                r1 = min(ds.height, math.ceil(raw.row_off+raw.height))
                if c1 <= c0 or r1 <= r0:
                    raise ValueError("Study window does not overlap the input raster.")
                w = Window(c0, r0, c1-c0, r1-r0)
            cells = int(w.width * w.height)
            if cells > max_cells:
                raise ValueError(f"{cells:,} cells exceed --max-cells={max_cells:,}. "
                                 "Crop with --search-buffer-km or use a larger memory limit. "
                                 "The program will not silently coarsen the analysis.")
            arr = ds.read(1, window=w, masked=True).astype(np.float32)
            self.z = np.asarray(arr.filled(np.nan), dtype=np.float32)
            self.z = self.z * np.float32(ds.scales[0]) + np.float32(ds.offsets[0])
            self.z[~np.isfinite(self.z)] = np.nan
            self.transform = ds.window_transform(w)
            self.shape = self.z.shape
            self.dx, self.dy = abs(self.transform.a), abs(self.transform.e)
            self.bounds = rasterio.transform.array_bounds(*self.shape, self.transform)
            self.window = w
            self.metadata = dict(
                source_filename=source.name,
                sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                source_shape=[ds.height, ds.width], analysis_shape=list(self.shape),
                source_crs_wkt=ds.crs.to_wkt(), source_crs_authority=ds.crs.to_authority(),
                source_crs_declared_id=crs.to_json_dict().get("id"),
                resolved_horizontal_authority=crs.to_authority(),
                source_crs_name=crs.name,
                crs_note="The embedded CRS identifier and axis declaration may differ; retain the WKT and use its actual X/Y projection with always_xy=True.",
                raster_axis_order="X=easting, Y=northing; transforms always_xy=True",
                transform=list(self.transform)[:6], pixel_size_m=[self.dx, self.dy],
                source_bounds=list(ds.bounds), analysis_bounds=list(self.bounds),
                source_area_or_point=ds.tags().get("AREA_OR_POINT"),
                source_nodata=(float(ds.nodata) if ds.nodata is not None and np.isfinite(ds.nodata) else None), band_unit=ds.units[0],
                source_scale=ds.scales[0], source_offset=ds.offsets[0],
                source_tags=ds.tags(), band_tags=ds.tags(1),
                vertical_datum="NOT VERIFIED: no vertical CRS is encoded in the provided file",
                vertical_values="Native raster heights, metres; no geoid/vertical transform applied",
                valid_cells=int(np.isfinite(self.z).sum()),
                search_buffer_km=buffer_km,
                analysis_window=[int(w.col_off), int(w.row_off), int(w.width), int(w.height)],
            )
        self.anchors = []
        self.point_records = []
        for i, ((lat, lon), (x0, y0)) in enumerate(zip(points, points_xy)):
            r, c = self.rc(x0, y0)
            r, c = int(r), int(c)
            if not (0 <= r < self.shape[0] and 0 <= c < self.shape[1]):
                raise ValueError(f"Point {i} lies outside the analysis window.")
            if not np.isfinite(self.z[r, c]):
                raise ValueError(f"Point {i} lies on NoData; choose another explicit point.")
            self.anchors.append((r, c))
            xx, yy = self.xy(r, c)
            llon, llat = self.inverse.transform(xx, yy)
            self.point_records.append(dict(
                name="Heradsbygd / start" if i == 0 else ("Kolomoen / end" if i == len(points)-1 else f"Via {i}"),
                requested_latitude=lat, requested_longitude=lon,
                requested_easting_m=x0, requested_northing_m=y0,
                cell_center_latitude=float(llat), cell_center_longitude=float(llon),
                cell_center_easting_m=float(xx), cell_center_northing_m=float(yy),
                snap_distance_m=float(math.hypot(xx-x0, yy-y0)),
                row=r, col=c, height_native_m=float(self.z[r, c]),
            ))
        self.points_xy = points_xy

    def rc(self, x0, y0):
        # Rasterio/GDAL already accounts for PixelIsPoint in the affine transform.
        col = np.floor((np.asarray(x0) - self.transform.c) / self.transform.a).astype(np.int64)
        row = np.floor((np.asarray(y0) - self.transform.f) / self.transform.e).astype(np.int64)
        return row, col

    def xy(self, r, c):
        return (self.transform.c + (np.asarray(c)+.5)*self.transform.a,
                self.transform.f + (np.asarray(r)+.5)*self.transform.e)

    def table(self, x0, y0, distances):
        x0, y0, distances = np.asarray(x0), np.asarray(y0), np.asarray(distances)
        rr, cc = self.rc(x0, y0)
        lon, lat = self.inverse.transform(x0, y0)
        rows = []
        for i in range(len(x0)):
            inside = 0 <= rr[i] < self.shape[0] and 0 <= cc[i] < self.shape[1]
            zz = self.z[rr[i], cc[i]] if inside else np.nan
            rows.append(dict(index=i, distance_m=float(distances[i]),
                             easting_m=float(x0[i]), northing_m=float(y0[i]),
                             longitude=float(lon[i]), latitude=float(lat[i]),
                             height_native_m=float(zz) if np.isfinite(zz) else None,
                             raster_row=int(rr[i]+self.window.row_off),
                             raster_col=int(cc[i]+self.window.col_off)))
        return rows

    def save_mask(self, path, mask):
        with rasterio.open(path, "w", driver="GTiff", height=self.shape[0], width=self.shape[1],
                           count=1, dtype="uint8", crs=self.crs, transform=self.transform,
                           compress="deflate", tiled=True) as ds:
            ds.write(mask.astype(np.uint8), 1)
            ds.update_tags(description="1=connected to start at ceiling; 0=not in component; 255=NoData")


def route_table(grid, rr, cc):
    x0, y0 = grid.xy(rr, cc)
    distances = np.r_[0., np.cumsum(np.hypot(np.diff(x0), np.diff(y0)))]
    return grid.table(x0, y0, distances)


def sampled_route(grid, rows, spacing):
    distances = np.array([p["distance_m"] for p in rows])
    x0 = np.array([p["easting_m"] for p in rows]); y0 = np.array([p["northing_m"] for p in rows])
    ts = np.r_[np.arange(0, distances[-1], spacing), distances[-1]]
    return grid.table(np.interp(ts, distances, x0), np.interp(ts, distances, y0), ts)


def straight_profile(grid, count=None, spacing=None):
    a, b = np.array(grid.points_xy[0]), np.array(grid.points_xy[-1])
    length = float(np.linalg.norm(b-a))
    ts = np.linspace(0., length, count) if count else np.r_[np.arange(0., length, spacing), length]
    us = ts / max(length, 1e-15)
    points = a[None, :] + us[:, None]*(b-a)[None, :]
    return grid.table(points[:, 0], points[:, 1], ts)


def straight_cell_profile(grid):
    """Visit every cell crossed by the projected straight line (positive-length intersections).

    This is not a 50 m sampling approximation: each row/column grid crossing
    contributes a break point. Exact corner contacts use the adjacent segments.
    """
    a, b = np.array(grid.points_xy[0]), np.array(grid.points_xy[-1])
    ta = ~grid.transform
    ca, ra = ta * tuple(a); cb, rb = ta * tuple(b)
    breaks = [0., 1.]
    for v0, v1 in ((ca, cb), (ra, rb)):
        if abs(v1-v0) > 1e-12:
            boundaries = np.arange(math.floor(min(v0,v1))+1, math.ceil(max(v0,v1)))
            breaks.extend(((boundaries-v0)/(v1-v0)).tolist())
    breaks = np.unique(np.clip(breaks, 0, 1))
    us = np.r_[0., .5*(breaks[:-1]+breaks[1:]), 1.]
    p = a[None, :] + us[:, None]*(b-a)[None, :]
    return grid.table(p[:,0], p[:,1], us*np.linalg.norm(b-a))


def boundary_status(labels, anchor, finite):
    label = int(labels[anchor])
    if label == 0:
        return dict(present=False, touches_edge=False, touches_nodata=False)
    at_edge = bool(np.any(labels[0,:] == label) or np.any(labels[-1,:] == label)
                   or np.any(labels[:,0] == label) or np.any(labels[:,-1] == label))
    at_nodata = False
    if not finite.all():
        fringe = ndi.binary_dilation(~finite, structure=FOUR)
        at_nodata = bool(np.any((labels == label) & fringe))
    return dict(present=True, touches_edge=at_edge, touches_nodata=at_nodata)


def find_plateaus(grid, level, below_labels, path_rc):
    """Critical-height plateaus, with a distinction between real bridges and contour touches.

    A 'direct bridge' touches both the start and goal subcritical components.
    There can be more complex chains through intermediate components; route-
    intersecting candidates are retained, but not called proven unique saddles.
    """
    critical, count = ndi.label(grid.z == np.float32(level), structure=FOUR)
    start_id, end_id = int(below_labels[grid.anchors[0]]), int(below_labels[grid.anchors[-1]])
    route_ids = set()
    for rr, cc in path_rc:
        route_ids.update(int(x) for x in np.unique(critical[rr,cc]) if x)
    table, features = [], []
    for k, sl0 in enumerate(ndi.find_objects(critical), start=1):
        if sl0 is None:
            continue
        sl = (slice(max(sl0[0].start-1,0), min(sl0[0].stop+1,grid.shape[0])),
              slice(max(sl0[1].start-1,0), min(sl0[1].stop+1,grid.shape[1])))
        mask = critical[sl] == k
        ring = ndi.binary_dilation(mask, structure=FOUR) & ~mask
        adjacent = np.unique(below_labels[sl][ring])
        bridge = bool(start_id != end_id and start_id != 0 and end_id != 0
                      and start_id in adjacent and end_id in adjacent)
        if not bridge and k not in route_ids:
            continue
        rr, cc = np.nonzero(mask); rr += sl[0].start; cc += sl[1].start
        xs, ys = grid.xy(rr, cc)
        centre = int(np.argmin((xs-xs.mean())**2+(ys-ys.mean())**2))
        lon, lat = grid.inverse.transform(xs[centre], ys[centre])
        necessary = None
        if bridge:
            allowed = np.isfinite(grid.z) & (grid.z <= np.float32(level))
            allowed[rr, cc] = False
            without_plateau = ndi.label(allowed, structure=FOUR)[0]
            necessary = not anchors_connected(without_plateau, grid.anchors)
            del allowed, without_plateau
        record = dict(plateau_id=k, direct_bridge=bridge, necessary_at_ceiling=necessary, route_intersects=k in route_ids,
                      cells=len(rr), raster_support_area_m2=float(len(rr)*grid.dx*grid.dy),
                      height_native_m=level, longitude=float(lon), latitude=float(lat),
                      easting_m=float(xs[centre]), northing_m=float(ys[centre]),
                      extent_east_west_m=float(xs.max()-xs.min()+grid.dx),
                      extent_north_south_m=float(ys.max()-ys.min()+grid.dy),
                      adjacent_subcritical_components=";".join(str(int(v)) for v in adjacent if v))
        table.append(record)
        tf = grid.transform * rasterio.Affine.translation(sl[1].start, sl[0].start)
        for geom, value in shapes(mask.astype(np.uint8), mask=mask, transform=tf):
            features.append(dict(type="Feature", geometry=transform_geom(grid.crs, "EPSG:4326", geom),
                                 properties=record.copy()))
    table.sort(key=lambda p: (not p["direct_bridge"], not p["route_intersects"], p["plateau_id"]))
    return table, features


def cross_section(grid, saddle, radius=1000., spacing=10.):
    # Fixed north-south section: reproducible, not an inferred optimal saddle axis.
    ts = np.arange(-radius, radius+.5*spacing, spacing)
    rows = grid.table(np.full_like(ts, saddle["easting_m"]), saddle["northing_m"]+ts, ts)
    for row in rows:
        row["north_offset_m"] = row.pop("distance_m")
    return rows


def line_feature(rows, name, **properties):
    """RFC 7946: keep geometry 2D; native, unconverted heights live in properties."""
    return dict(type="Feature", geometry=dict(type="LineString", coordinates=[
        [p["longitude"], p["latitude"]] for p in rows]), properties=dict(
        name=name, native_elevations_m=[p["height_native_m"] for p in rows],
        distances_m=[p["distance_m"] for p in rows],
        vertical_datum="input raster datum, unverified", **properties))


def route_stats(rows):
    z = np.array([p["height_native_m"] for p in rows], dtype=float)
    return dict(length_m=rows[-1]["distance_m"], vertices=len(rows),
                min_height_native_m=float(z.min()), max_height_native_m=float(z.max()))


def create_plots(grid, routes, direct, level, saddles, out):
    primary = routes[0]
    # One separate figure per chart, using Matplotlib's default palette.
    allx = [p["easting_m"] for p in primary]; ally = [p["northing_m"] for p in primary]
    xmin, xmax = min(allx)-2200., max(allx)+2200.
    ymin, ymax = min(ally)-2200., max(ally)+2200.
    r0, c0 = grid.rc(xmin, ymax); r1, c1 = grid.rc(xmax, ymin)
    r0=max(0,int(r0)); r1=min(grid.shape[0],int(r1)+1)
    c0=max(0,int(c0)); c1=min(grid.shape[1],int(c1)+1)
    step=max(1,math.ceil(max(r1-r0,c1-c0)/1600))
    zview=grid.z[r0:r1:step,c0:c1:step]
    left=grid.transform.c+c0*grid.dx; right=grid.transform.c+c1*grid.dx
    top=grid.transform.f-r0*grid.dy; bottom=grid.transform.f-r1*grid.dy
    fig, ax=plt.subplots(figsize=(13,9),layout="constrained")
    im=ax.imshow(zview,extent=[left/1000,right/1000,bottom/1000,top/1000],origin="upper")
    fig.colorbar(im,ax=ax,shrink=.8,label="Høyde i originalrasteret (m)")
    ax.plot(np.array(allx)/1000,np.array(ally)/1000,linewidth=2.1,label="Laveste terskel, korteste rute")
    ax.plot([direct[0]["easting_m"]/1000,direct[-1]["easting_m"]/1000],
            [direct[0]["northing_m"]/1000,direct[-1]["northing_m"]/1000],
            linestyle="--",linewidth=1.8,label="Rett referanselinje")
    if len(routes)>1:
        ax.plot([p["easting_m"]/1000 for p in routes[1]], [p["northing_m"]/1000 for p in routes[1]],
                linewidth=1.1,label="Dalvektet alternativ, samme terskel")
    for p,label in [(primary[0],"Heradsbygd"),(primary[-1],"Kolomoen")]:
        ax.scatter(p["easting_m"]/1000,p["northing_m"]/1000,s=35)
        ax.annotate(label,(p["easting_m"]/1000,p["northing_m"]/1000),xytext=(8,10),textcoords="offset points",fontsize=11,bbox=dict(boxstyle="round,pad=0.22",facecolor=ax.get_facecolor(),alpha=.8,linewidth=0))
    if saddles:
        s=saddles[0]
        ax.scatter(s["easting_m"]/1000,s["northing_m"]/1000,s=90,marker="x")
        ax.annotate(f"Terskelområde: {level:.1f} m",(s["easting_m"]/1000,s["northing_m"]/1000),
                    xytext=(18,35),textcoords="offset points",arrowprops=dict(arrowstyle="->"),fontsize=11,bbox=dict(boxstyle="round,pad=0.22",facecolor=ax.get_facecolor(),alpha=.8,linewidth=0))
    ax.set(xlabel="Øst (km, rasterets UTM-sone)",ylabel="Nord (km, rasterets UTM-sone)",
           title="Heradsbygd–Kolomoen | laveste sammenhengende terrengpassasje")
    ax.legend(loc="lower right",fontsize=9)
    fig.savefig(out/"overview.png",dpi=170); plt.close(fig)

    fig,ax=plt.subplots(figsize=(13,4.8),layout="constrained")
    ax.plot([p["distance_m"]/1000 for p in primary],[p["height_native_m"] for p in primary],
            label="Laveste passasje (alle ruteceller)",linewidth=1.2)
    ax.plot([p["distance_m"]/1000 for p in direct],[p["height_native_m"] for p in direct],
            label="Rett linje (alle kryssede celler)",linewidth=1.1)
    ax.axhline(level,linestyle="--",label=f"Minste sammenhengende terskel: {level:.1f} m")
    ax.set(xlabel="Avstand langs den enkelte linjen (km)",ylabel="Høyde i originalrasteret (m)",
           title="Terrengprofiler | avstandsaksen følger hver sin rute")
    ax.grid(alpha=.2); ax.legend(); fig.savefig(out/"profiles.png",dpi=170);plt.close(fig)

    if saddles:
        section=cross_section(grid,saddles[0],spacing=min(grid.dx,grid.dy))
        write_csv(out/"saddle_north_south_2km.csv",section)
        fig,ax=plt.subplots(figsize=(11,4.4),layout="constrained")
        ax.plot([p["north_offset_m"] for p in section],[p["height_native_m"] for p in section])
        ax.axhline(level,linestyle="--",label=f"Terskelnivå {level:.1f} m")
        ax.set(xlabel="Avstand fra representativ terskelcelle (m; positiv = nord)",
               ylabel="Høyde i originalrasteret (m)",
               title="Fast nord–sør-profil gjennom terskelområdet (ikke utledet sadelakse)")
        ax.grid(alpha=.2);ax.legend();fig.savefig(out/"saddle_section.png",dpi=170);plt.close(fig)


def export_web_grid(grid, routes, levels, out, max_side=2048):
    """Reproject a DISPLAY-ONLY crop. Scientific connectivity precedes reprojection."""
    points=[p for route in routes for p in route]
    xs=[p["easting_m"] for p in points];ys=[p["northing_m"] for p in points]
    l,b,r,t=grid.bounds
    box=(max(l,min(xs)-4500),max(b,min(ys)-4500),min(r,max(xs)+4500),min(t,max(ys)+4500))
    west,south,east,north=transform_bounds(grid.crs,4326,*box,densify_pts=41)
    mean_lat=(south+north)/2
    width_m=(east-west)*111320.*math.cos(math.radians(mean_lat))
    height_m=(north-south)*111320.
    resolution=max(grid.dx,grid.dy,max(width_m,height_m)/max_side)
    width=max(64,math.ceil(width_m/resolution));height=max(64,math.ceil(height_m/resolution))
    transform=from_bounds(west,south,east,north,width,height)
    dem=np.full((height,width),np.nan,np.float32)
    reproject(grid.z,dem,src_transform=grid.transform,src_crs=grid.crs,
              dst_transform=transform,dst_crs="EPSG:4326",src_nodata=np.nan,dst_nodata=np.nan,
              resampling=Resampling.bilinear,num_threads=2)
    dem.astype("<f4").tofile(out/"terrain.f32")
    finite=np.isfinite(dem)
    # Hillshade texture is an illustrative grayscale image, not another DEM source.
    filled=np.where(finite,dem,np.nanmedian(dem))
    gy,gx=np.gradient(filled,resolution,resolution)
    normal=np.dstack((-gx,gy,np.ones_like(gx)))
    normal/=np.maximum(np.linalg.norm(normal,axis=2,keepdims=True),1e-12)
    sun=np.array([-.5,.5,math.sqrt(.5)])
    light=np.clip(np.sum(normal*sun,axis=2),0.,1.)
    shade=np.uint8(np.clip(70+170*light,0,255))
    rgba=np.dstack((shade,shade,shade,np.where(finite,255,0).astype(np.uint8)))
    Image.fromarray(rgba).save(out/"hillshade.png")
    web_levels=[]
    for i,item in enumerate(levels):
        lab=labels_at(grid.z,item["height_m"])
        sid=int(lab[grid.anchors[0]]);tid=int(lab[grid.anchors[-1]])
        classification=np.zeros(grid.shape,np.uint8)
        if tid:classification[lab==tid]=2
        if sid:classification[lab==sid]=1
        dst=np.zeros((height,width),np.uint8)
        reproject(classification,dst,src_transform=grid.transform,src_crs=grid.crs,
                  dst_transform=transform,dst_crs="EPSG:4326",resampling=Resampling.nearest,num_threads=2)
        # These are classified connected sublevel sets, NOT hydraulic flood footprints.
        colors=np.array([[0,0,0,0],[39,167,223,135],[235,165,70,135]],np.uint8)
        filename=f"connected_{i}.png"
        Image.fromarray(colors[dst]).save(out/filename)
        web_levels.append(dict(**item,file=filename,
                               connected=anchors_connected(lab,grid.anchors)))
    return dict(width=width,height=height,west=west,south=south,east=east,north=north,
                raw_file="terrain.f32",raw_type="little-endian float32, row-major, north to south",
                approximate_display_spacing_m=resolution,
                pixel_registration="outer bounds; samples at pixel centres (subtract 0.5 on lookup)",
                vertical="native DEM heights, NOT converted to ellipsoidal height",levels=web_levels)


def self_test():
    def h(z,s,t):return minimum_ceiling(np.array(z,dtype=np.float32),[s,t])
    assert h([[1,9,1],[1,3,1],[1,4,1]],(0,0),(0,2))==3.
    assert h([[1,9],[9,1]],(0,0),(1,1))==9.  # disallow diagonal corner leak
    assert minimum_ceiling(np.array([[1,9],[9,1]],np.float32),[(0,0),(1,1)],True)==1.
    assert h([[-4,-2,-3]],(0,0),(0,2))==-2.
    assert h([[1,2,7]],(0,0),(0,2))==7.      # endpoint-dominated case
    z=np.array([[1,np.nan,1],[1,2,1]],np.float32)
    assert minimum_ceiling(z,[(0,0),(0,2)])==2.
    try:h([[1,np.nan,1]],(0,0),(0,2))
    except ValueError:pass
    else:raise AssertionError("NoData barrier was crossed")
    # Brute-force height enumeration validates 120 random small cases, including ties and negatives.
    rng=np.random.default_rng(72)
    for _ in range(120):
        z=rng.integers(-5,20,size=(7,9)).astype(np.float32)
        anchors=[(0,0),(6,8)]
        expected=next(float(v) for v in np.unique(z) if anchors_connected(labels_at(z,v),anchors))
        level=minimum_ceiling(z,anchors)
        assert level==expected
        rr,cc,cost,expanded=solve_path(z,level,anchors,10.,10.)
        assert z[rr,cc].max()<=level
    # Shortest path on a clear 8-neighbour rectangular grid, exact analytic length.
    z=np.zeros((5,7),np.float32)
    rr,cc,cost,_=solve_path(z,0.,[(0,0),(4,6)],10.,20.)
    assert abs(cost-(4*math.hypot(10,20)+2*10))<1e-8
    # Via points are solved at the GLOBAL minimax ceiling, not at separate leg ceilings.
    anchors=[(0,0),(4,3),(0,6)]
    rr,cc,_,_=solve_path(z,0.,anchors,10.,20.)
    assert (4,3) in set(zip(rr,cc))
    log("PASS: 120 random minimax comparisons + corner, NoData, negative-height, endpoint, shortest-path and via tests.")


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dem",nargs="?",type=Path,default=Path("6706_3_10m_z32.tif"))
    parser.add_argument("--out",type=Path,default=Path("results"))
    parser.add_argument("--start",nargs=2,type=float,metavar=("LAT","LON"),default=DEFAULT_START)
    parser.add_argument("--end",nargs=2,type=float,metavar=("LAT","LON"),default=DEFAULT_END)
    parser.add_argument("--via",nargs=2,type=float,action="append",default=[],metavar=("LAT","LON"))
    parser.add_argument("--search-buffer-km",type=float,default=None,help="Crop around anchors; default: entire supplied tile")
    parser.add_argument("--max-cells",type=int,default=50_000_000,help="Guard against excessive RAM; 25M cells need roughly 1-2 GB")
    parser.add_argument("--profile-step-m",type=float,default=50.)
    parser.add_argument("--valley-weight",type=float,default=2.,help="0 disables the heuristic alternative; does not change the ceiling")
    parser.add_argument("--web-max-side",type=int,default=2048,help="Display-only grid maximum side")
    parser.add_argument("--skip-web",action="store_true")
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:
        self_test();return
    if args.profile_step_m<=0 or args.valley_weight<0 or args.web_max_side<64:
        parser.error("Profile spacing must be positive, valley weight nonnegative, web size >=64.")
    if args.search_buffer_km is not None and args.search_buffer_km<=0:
        parser.error("Search buffer must be positive.")
    points=[tuple(args.start),*[tuple(p) for p in args.via],tuple(args.end)]
    if any(not (-90<=lat<=90 and -180<=lon<=180) for lat,lon in points):
        parser.error("Use --start/--end/--via LAT LON in decimal degrees.")
    if not args.dem.is_file():
        parser.error(f"DEM not found: {args.dem}")
    tic=time.time();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    log("Reading native raster (no analysis resampling) ...")
    grid=Grid(args.dem.resolve(),points,args.search_buffer_km,args.max_cells)
    log(f"Raster {grid.shape[1]} x {grid.shape[0]}, {grid.dx:g} x {grid.dy:g} m; resolved {grid.metadata['resolved_horizontal_authority']}")
    log("Finding the minimum connecting height ...")
    level=minimum_ceiling(grid.z,grid.anchors)
    loose8=minimum_ceiling(grid.z,grid.anchors,diagonal_only=True)
    below=np.nextafter(np.float32(level),np.float32(-np.inf))
    lower_values=grid.z[np.isfinite(grid.z)&(grid.z<np.float32(level))]
    preceding=float(lower_values.max()) if len(lower_values) else float(below)
    log(f"Minimum ceiling = {level:.6f} native metres; previous raster level = {preceding:.6f}")
    log("Choosing the shortest path below that ceiling ...")
    rr,cc,cost,expanded=solve_path(grid.z,level,grid.anchors,grid.dx,grid.dy)
    primary=route_table(grid,rr,cc);routes=[primary];path_rc=[(rr,cc)]
    write_csv(out/"route_all_cells.csv",primary)
    write_csv(out/f"route_{args.profile_step_m:g}m.csv",sampled_route(grid,primary,args.profile_step_m))
    alt_stats=None
    if args.valley_weight>0:
        log("Computing a valley-biased alternative, without changing the height ceiling ...")
        r2,c2,cost2,expanded2=solve_path(grid.z,level,grid.anchors,grid.dx,grid.dy,args.valley_weight)
        alternative=route_table(grid,r2,c2);routes.append(alternative);path_rc.append((r2,c2))
        alt_stats=route_stats(alternative);alt_stats.update(weight=args.valley_weight,weighted_cost=cost2)
        write_csv(out/"valley_route_all_cells.csv",alternative)
    log("Sampling the reference line and identifying the critical plateau ...")
    control=straight_profile(grid,count=27)
    reference50=straight_profile(grid,spacing=args.profile_step_m)
    direct=straight_cell_profile(grid)
    write_csv(out/"reference_27_NEW_points.csv",control)
    write_csv(out/f"reference_{args.profile_step_m:g}m.csv",reference50)
    write_csv(out/"reference_all_crossed_cells.csv",direct)
    write_csv(out/"input_points.csv",grid.point_records)
    below_labels=labels_at(grid.z,float(below))
    statuses=[boundary_status(below_labels,p,np.isfinite(grid.z)) for p in (grid.anchors[0],grid.anchors[-1])]
    saddles,plateau_features=find_plateaus(grid,level,below_labels,path_rc)
    plateau_geojson=dict(type="FeatureCollection",features=plateau_features)
    del below_labels
    write_csv(out/"critical_plateaus.csv",saddles)
    write_json(out/"critical_plateaus.geojson",plateau_geojson)
    lab=labels_at(grid.z,level)
    source_id=int(lab[grid.anchors[0]])
    mask=(lab==source_id).astype(np.uint8)
    mask[~np.isfinite(grid.z)]=255
    grid.save_mask(out/"connected_at_threshold.tif",mask)
    area=float(np.count_nonzero(mask==1)*grid.dx*grid.dy)
    del lab,mask
    levels=[dict(label=f"{level-1:.1f} m",height_m=level-1),
            dict(label=f"Under {level:.1f} m (siste rasterverdi {preceding:.1f})",height_m=preceding),
            dict(label=f"Terskel: {level:.1f} m",height_m=level),
            dict(label=f"Terskel + 1 m: {level+1:.1f}",height_m=level+1),
            dict(label=f"Terskel + 5 m: {level+5:.1f}",height_m=level+5)]
    checks=[]
    for item in levels:
        lab=labels_at(grid.z,item["height_m"])
        checks.append(dict(height_native_m=item["height_m"],connected=anchors_connected(lab,grid.anchors),
                           start_component=int(lab[grid.anchors[0]]),end_component=int(lab[grid.anchors[-1]])))
    write_csv(out/"connectivity_checks.csv",checks)
    del lab
    geo=[line_feature(primary,"Lowest ceiling / shortest grid route",ceiling_native_m=level),
         line_feature(direct,"Straight reference line")]
    if len(routes)>1:geo.append(line_feature(routes[1],"Valley-biased representative",ceiling_native_m=level,
                                            valley_weight=args.valley_weight))
    for p in grid.point_records:
        geo.append(dict(type="Feature",geometry=dict(type="Point",coordinates=[p["requested_longitude"],p["requested_latitude"]]),
                        properties=p))
    write_json(out/"routes.geojson",dict(type="FeatureCollection",features=geo))
    # H*>both endpoint elevations is sufficient to disprove an everywhere-downhill bed path in either direction.
    no_downhill_start_to_end=level>grid.point_records[0]["height_native_m"]
    no_downhill_end_to_start=level>grid.point_records[-1]["height_native_m"]
    result=dict(
        software_version=VERSION,created_utc=datetime.now(timezone.utc).isoformat(),
        objective="Minimize maximum native raster height; then minimize 8-neighbour planimetric path length, no corner cutting",
        input=grid.metadata,points=grid.point_records,
        minimum_ceiling_native_m=level,previous_distinct_raster_height_m=preceding,
        raster_4_connected_ceiling_m=level,permissive_8_connected_ceiling_m=loose8,
        path=route_stats(primary),valley_alternative=alt_stats,
        straight_reference=dict(length_m=direct[-1]["distance_m"],
            max_height_native_m=max(p["height_native_m"] for p in direct if p["height_native_m"] is not None),
            control_points=len(control),sampled_profile_points=len(reference50),profile_spacing_m=args.profile_step_m,
            crossed_cell_samples=len(direct),
            note="NEW straight UTM line through user coordinates, NOT the unavailable historical 27 points"),
        critical_plateaus=saddles,connectivity_checks=checks,
        critical_component_area_within_tile_m2=area,
        boundary_status_below_ceiling=dict(start=statuses[0],end=statuses[1]),
        boundary_warning=("Both subcritical endpoint components touch the analysis boundary or NoData. "
                          "A lower route outside the supplied/selected domain is not ruled out." if
            all(s["touches_edge"] or s["touches_nodata"] for s in statuses) else
            "At least one endpoint is enclosed below the ceiling, subject to the raster topology and input validity."),
        downhill_bed_test=dict(
            start_to_end="ruled_out_by_height_ceiling" if no_downhill_start_to_end else "not_tested",
            end_to_start="ruled_out_by_height_ceiling" if no_downhill_end_to_start else "not_tested"),
        caveats=[
            "This is a topographic sublevel-connectivity test, not a hydraulic model or proof of a former river.",
            "The threshold permits zero-depth contact; positive flow requires additional head and a finite cross-section.",
            "Water direction, ice dams, erosion history, deposits, and time-varying topography are not reconstructed.",
            "No sink filling, stream burning, smoothing, or manual lowering of roads/embankments is applied.",
            "A centreline is not unique and its bends depend on raster topology and secondary cost.",
            "10 m horizontal spacing and decimal raster values do not imply decimetre vertical accuracy.",
            "No bathymetry is inferred from apparent lake surfaces in the terrain raster.",
            "GeoJSON geometry is 2D; native, non-ellipsoidal elevations are stored as explicit properties.",
            "Viewer terrain is display-only and uses native heights as local display heights; absolute globe altitude is not geodetically corrected.",
        ],
        python=platform.python_version(),libraries=dict(numpy=np.__version__,rasterio=rasterio.__version__),
        astar_expanded_nodes=expanded,
    )
    log("Rendering overview and terrain profiles ...")
    create_plots(grid,routes,direct,level,saddles,out)
    if not args.skip_web:
        log("Preparing the CesiumJS display-only height grid and connectivity overlays ...")
        webgrid=export_web_grid(grid,routes,levels,out,args.web_max_side)
        display=dict(summary=result,terrain=webgrid,routes=[
            dict(name="Laveste terskel / korteste rute",points=[[p["longitude"],p["latitude"],p["height_native_m"],p["distance_m"]] for p in primary]),
            dict(name="Rett referanselinje",points=[[p["longitude"],p["latitude"],p["height_native_m"],p["distance_m"]] for p in direct])])
        if len(routes)>1:
            display["routes"].append(dict(name="Dalvektet alternativ",points=[[p["longitude"],p["latitude"],p["height_native_m"],p["distance_m"]] for p in routes[1]]))
        write_json(out/"viewer_data.json",display)
        (out/"index.html").write_text(VIEWER_HTML,encoding="utf-8")
        write_standalone_viewer(out/"Heradsbygd-Kolomoen.html",VIEWER_HTML,out,display,plateau_geojson)
        result["viewer"]=dict(approximate_display_spacing_m=webgrid["approximate_display_spacing_m"],
                               cesium_version="1.145",ion_token_required=False,
                               note="Standalone HTML embeds local results; CDN internet access is still required for CesiumJS")
    result["elapsed_seconds"]=time.time()-tic
    write_json(out/"summary.json",result)
    (out/"RESULTAT.md").write_text(make_report(result),encoding="utf-8")
    log(f"Done. Ceiling {level:.1f} m; route {primary[-1]['distance_m']/1000:.3f} km; {time.time()-tic:.1f} s.")
    log(f"Results: {out}")
    log(f"Viewer: python -m http.server 8000 --bind 127.0.0.1 --directory '{out}'")


def make_report(s):
    critical=[x for x in s["critical_plateaus"] if x["direct_bridge"]]
    plateau_text="Ingen direkte bro mellom de to subkritiske bassengene ble identifisert; se kandidatlisten."
    if critical:
        p=critical[0]
        plateau_text=(f"Den kritiske sammenkoblingen skjer i et rasterplatå med {p['cells']} celler på "
                      f"{p['height_native_m']:.1f} m. Representativ celle: {p['latitude']:.8f} N, "
                      f"{p['longitude']:.8f} E; øst {p['easting_m']:.1f}, nord {p['northing_m']:.1f}. "
                      f"Rasterstøtten er {p['extent_east_west_m']:.0f} × {p['extent_north_south_m']:.0f} m "
                      "i utstrekning. Dette er ikke en målt kanalbredde eller en endelig geologisk sadel.")
    return f'''# Heradsbygd–Kolomoen: laveste sammenhengende terrengpassasje

## Beregnet resultat

**Minste forbindende rasterhøyde: {s['minimum_ceiling_native_m']:.1f} m.**
Korteste rute under dette taket: **{s['path']['length_m']/1000:.3f} km**.
Ved forrige lavere rasterhøyde ({s['previous_distinct_raster_height_m']:.1f} m) er endepunktene ikke forbundet.
Dette er beregnet fra alle relevante celler i originalrasteret, ikke fra et anslag eller en 50-metersprofil.

| Størrelse | Verdi |
|---|---:|
| Startcellens høyde | {s['points'][0]['height_native_m']:.1f} m |
| Sluttcellens høyde | {s['points'][-1]['height_native_m']:.1f} m |
| Laveste forbindende terskel | {s['minimum_ceiling_native_m']:.1f} m |
| Korteste rute ved denne terskelen | {s['path']['length_m']/1000:.3f} km |
| Rett referanselinje, lengde | {s['straight_reference']['length_m']/1000:.3f} km |
| Høyeste kryssede rastercelle langs rett linje | {s['straight_reference']['max_height_native_m']:.1f} m |

{plateau_text}

## Kilde og referansesystem

Fil: `{s['input']['source_filename']}`. SHA256: `{s['input']['sha256']}`.
Originalraster: {s['input']['source_shape'][1]} × {s['input']['source_shape'][0]} celler.
Cellestørrelse: {s['input']['pixel_size_m'][0]:g} × {s['input']['pixel_size_m'][1]:g} m.
CRS-navn: {s['input']['source_crs_name']}. Innebygd CRS-ID: {s['input']['source_crs_declared_id']}.
Pyproj identifiserer den faktiske projeksjonen/akseordenen som {s['input']['resolved_horizontal_authority']}.
Rasterio returnerer ingen entydig autoritetskode for den opprinnelige WKT-en. Dette avviket
bevares i metadata i stedet for å skjules. Rasterets X/Y behandles som øst/nord;
koordinattransformasjoner bruker `always_xy=True`. `AREA_OR_POINT` er `{s['input']['source_area_or_point']}`.
Rasterio/GDAL håndterer dette allerede: ingen ekstra halv-celle-forskyvning legges til.

**Vertikalt datum er ikke bekreftet av filens metadata.** Tallene er originalrasterets høyder i meter.
De må ikke automatisk omtales som kvalitetssikrede NN2000- eller ellipsoidehøyder.
Desimalen i en rasterverdi er ikke en påstand om tilsvarende målenøyaktighet.

## Definisjon og algoritme

For en sammenhengende rute P minimeres `max(z for cell in P)`.
Vi øker et høydenivå og finner det første nivået som kobler sammen endepunktene
(og eventuelle eksplisitte `--via`-punkter). Deretter søker A* den korteste ruten
innenfor dette nivået. Denne totrinnsmetoden unngår å blande høyde og lengde i ett
vilkårlig hovedmål. Alle sluttverdier valideres med sammenheng like under og ved terskelen.

Sammenheng beregnes med felles celleside. Ruten kan også gå diagonalt når **begge**
sidestøttende celler ligger under høydegrensen. Det hindrer forbindelser som bare
berører hverandre i et hjørne. Som følsomhetskontroll ga også ukritisk 8-naboskap
{s['permissive_8_connected_ceiling_m']:.1f} m. Dette er et resultat på rastergrafen,
ikke en eksakt kontinuerlig optimal rute over et ukjent terreng mellom målepunktene.

Det dalvektede alternativet bruker samme terskel og en sekundær kostnad
`ds * (1 + weight * ((max(z1,z2) - raster_min)/(ceiling-raster_min))**4)`.
Dette er en illustrativ preferanse for lavere terreng, ikke en hydraulisk strømmodell.
Ruter er ikke glattet eller forenklet før analyse/eksport.

## Avgrensning

{s['boundary_warning']}

Søkedomenet er {s['input']['analysis_shape'][1]} × {s['input']['analysis_shape'][0]} celler.
Ingen lavere passasje **innenfor dette domenet** er funnet eller tillatt av rastertesten.
Det er ikke satt et krav om passasje gjennom Rokosjøen, siden ingen koordinat eller
avgrenset passasjesone der ble oppgitt. `--via LAT LON` kan legge inn et slikt eksplisitt krav.

## Hva dette sier – og ikke sier – om vann

Terskelen er høyere enn begge endepunktene. Dermed finnes ingen hele veien
ikke-stigende terrengbunn mellom akkurat disse punktene i denne rastermodellen.
Det utelukker ikke vanntransport over en terskel ved høyere vannstand eller oppdemming.
Minimax-testen er symmetrisk og avgjør ikke hvilken retning vannet ville beveget seg.

Nivålagene i kartet er sammenhengende **terrengceller under et nivå**, ikke dokumenterte
innsjøer eller beregnede flomsoner. Ved selve minimumsterskelen kan forbindelsen ha
null vanndybde; vannføring krever et reelt tverrsnitt, hydraulisk fall og tilstrekkelig vann.
Tidligere isdammer, innløp, andre utløp, senere erosjon og avsetninger er ikke modellert.
Vannoverflater i DTM er ikke innsjøbunn. Vegfyllinger og små naturlige terskler kan ha
avgjørende betydning og må kontrolleres i bedre data og eventuelt i felt.

## Referansepunkter og profiler

`reference_27_NEW_points.csv` inneholder **27 nye punkter** jevnt fordelt langs en
rett UTM-linje mellom koordinatene i brukerens melding. Disse er ikke de historiske
27 punktene som PDF-en sier må kvalitetssikres. PDF-ens omtrentlige endepunkter
avviker også fra brukerens nye endepunkter. Se `input_points.csv` for faktiske innparametre.

`reference_50m.csv` og `route_50m.csv` gir profiler med valgt intervall (standard 50 m).
Høydene er hentet fra cellen ved prøvepunktet; de er ikke nye feltmålinger.
`reference_all_crossed_cells.csv` og `route_all_cells.csv` bevarer tettere rasterkontroll.
Høyeste verdi på referanselinjen er kontrollert mot kryssede celler, ikke bare 50 m-prøver.
Rutens lengde måles mellom cellesentre. Endepunktene kobles til sin inneholdende celle;
små forskyvninger er eksplisitt oppgitt som `snap_distance_m` i `input_points.csv`.

`critical_plateaus.csv` skiller det forbindende terskelområdet fra andre steder der
ruten bare berører samme høydekote. `saddle_north_south_2km.csv` er en fast nord–sør-
profil gjennom en representativ terskelcelle, ikke en ferdig tolkning av landskapsformen.

## CesiumJS

Kjør `python -m http.server 8000 --bind 127.0.0.1 --directory results` fra prosjektmappen,
og åpne `http://localhost:8000`.
Se `README.md` i prosjektmappen for oppsett.

Viseren bruker egne terrengdata fra TIFF-en, ikke Cesium World Terrain.
Analyse og profiler bruker originaloppløsningen; nettleserterreng er reprojisert og
forenklet kun for visning. Terreng, ruter og terskelmarkører har konsistent lokal
visningshøyde, men originalhøydene er **ikke** transformert til ellipsoidehøyder.
En korrekt geodetisk overlagring krever bekreftet vertikalt datum og riktig HREF/geoid-
transformasjon; en vilkårlig konstant forskyvning er ikke lagt inn som om den var korrekt.

Kartet kan vise en valgt, felles høydeforstørrelse; avlest analysehøyde er uforandret.
Utenfor det eksporterte utsnittet er det ingen terrengdata (visningsellipsoide ved 0).
GeoJSON-geometri er 2D for ikke å feilmerke rå rasterhøyder som WGS84-ellipsoidehøyder.

## Reproduserbarhet

Programversjon {s['software_version']}; Python {s['python']}.
`summary.json` inneholder metadata, filhash, valideringsresultater, koordinater og forbehold.
Kjør `python lowest_passage.py --self-test` for syntetiske enhetstester.
'''

# The self-contained viewer is appended below. It needs only the generated files
# and the explicitly version-pinned CesiumJS CDN (or a local copy of Cesium assets).
VIEWER_HTML = r'''<!doctype html>
<html lang="nb">
<head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Heradsbygd–Kolomoen | terrengpassasje</title>
<link rel="stylesheet" href="vendor/Cesium/Widgets/widgets.css">
<script src="vendor/Cesium/Cesium.js"></script>
<style>
:root{--bg:#111b27;--panel:#172431;--line:#314453;--text:#e7eef2;--muted:#a8bcc8;--route:#29d4c5;--reference:#ffb65c;--alternate:#ce9dff}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}
#side{position:absolute;top:0;left:0;bottom:0;width:345px;padding:25px 22px 30px;overflow:auto;background:var(--panel);border-right:1px solid var(--line);z-index:2}
.eyebrow{font-size:10px;letter-spacing:2.1px;text-transform:uppercase;color:var(--route);font-weight:700}
h1{font-size:28px;letter-spacing:-.8px;line-height:1.14;margin:10px 0 8px}h2{font-size:13px;margin:22px 0 9px;letter-spacing:.2px}.intro{color:var(--muted);font-size:12px;margin:0 0 20px}
.metrics{display:grid;grid-template-columns:1fr 1fr;gap:9px}.metric{padding:12px;background:#10202b;border:1px solid var(--line);border-radius:6px}.metric strong{display:block;font-size:23px;font-weight:650;letter-spacing:-.6px}.metric span{display:block;font-size:10px;color:var(--muted);margin-top:2px}
select,button{background:#213444;color:var(--text);border:1px solid #445c6c;border-radius:5px;padding:8px;font:inherit}select{width:100%;font-size:12px}button{cursor:pointer;font-size:12px}button:hover{background:#304b5e}.buttons{display:flex;gap:7px;margin:12px 0}label.toggle{display:flex;align-items:center;gap:8px;margin:9px 0;font-size:12px}input[type=checkbox]{accent-color:var(--route)}.swatch{width:17px;height:3px;display:inline-block}.route{background:var(--route)}.reference{background:var(--reference)}.alternate{background:var(--alternate)}.map-controls{padding:10px;background:#10202b;border:1px solid var(--line);border-radius:6px;margin-top:12px}.map-controls select{margin:5px 0 9px}.range-label{display:flex;justify-content:space-between;font-size:11px;color:var(--muted)}input[type=range]{width:100%;margin:5px 0 0;accent-color:var(--route)}
#connectivity{font-size:12px;padding:9px 10px;border-left:3px solid var(--route);background:#10202b;margin:9px 0}.small{font-size:11px;color:var(--muted)}.warn{margin-top:18px;padding:11px 12px;border:1px solid #65553a;background:#302a20;color:#e9d5b1;font-size:11px;border-radius:5px}.warn strong{display:block;margin-bottom:4px}.links{display:flex;flex-wrap:wrap;gap:7px 13px;margin-top:12px}a{color:#a7dbe8;text-decoration:none;font-size:11px}a:hover{text-decoration:underline}#picked{white-space:pre-line;font-size:12px;min-height:58px;color:var(--muted)}
#globe{position:absolute;left:345px;right:0;top:0;bottom:204px}#profilePanel{position:absolute;left:345px;right:0;bottom:0;height:204px;background:#101d29;border-top:1px solid var(--line);padding:8px 22px 8px}.profilehead{font-size:11px;display:flex;align-items:center;gap:12px;margin-bottom:3px;color:var(--muted)}.profile-picker{display:flex;align-items:center;gap:7px;white-space:nowrap}.profile-picker select{width:auto;min-width:230px;padding:4px 7px}.profilehead #profileReadout{margin-left:auto;text-align:right}#profile{width:100%;height:164px;display:block}#status{position:absolute;left:365px;top:18px;max-width:560px;z-index:5;padding:10px 14px;background:#10202bea;border:1px solid var(--line);border-radius:5px;font-size:12px}#status:empty{display:none}#legend{position:absolute;left:365px;bottom:222px;background:#10202be0;padding:7px 10px;font-size:10px;border-radius:4px;pointer-events:none;color:#bdd1dc}.cesium-viewer-bottom{margin-left:0!important}
@media(max-width:850px){#side{width:270px;padding:17px 15px}#globe,#profilePanel{left:270px}#status,#legend{left:283px}h1{font-size:23px}.metric strong{font-size:19px}}@media(max-width:560px){#side{width:210px;padding:13px 10px}#globe,#profilePanel{left:210px}#status,#legend{left:220px}.metrics{grid-template-columns:1fr}h1{font-size:20px}#profilePanel{padding:10px 5px}.profilehead span:last-child{display:none}}
</style>
</head>
<body>
<aside id="side">
<div class="eyebrow">Rasteranalyse · lokalt terreng</div>
<h1>Heradsbygd<br>til Kolomoen</h1>
<p class="intro">Laveste sammenhengende terrengterskel.<br>Ikke en beregnet eller historisk dokumentert elv.</p>
<div class="metrics"><div class="metric"><strong id="ceiling">—</strong><span>MINSTE TERSKEL · NATIVE METER</span></div><div class="metric"><strong id="length">—</strong><span>KORTESTE RUTE VED TERSKELEN</span></div><div class="metric"><strong id="directMax">—</strong><span>MAKSIMUM LANGS RETT LINJE</span></div><div class="metric"><strong id="resolution">—</strong><span>ORIGINALRASTERETS CELLESTØRRELSE</span></div></div>
<div class="buttons"><button id="overview">Hele ruten</button><button id="saddle">Til terskelen</button><button id="top">Ovenfra</button></div>
<h2>Terreng under valgt nivå</h2>
<select id="level" aria-label="Høydenivå"></select><div id="connectivity">Laster beregning …</div>
<label class="toggle"><input type="checkbox" id="showConnected" checked> Vis sammenhengende lavområder</label>
<p class="small">Blått: forbindelse fra start. Oransje: fra slutt, når områdene er adskilt. Dette er terrengforbindelse, ikke en flomsimulering.</p>
<h2>Ruter og kartlag</h2>
<label class="toggle"><input type="checkbox" id="showRoute" checked><i class="swatch route"></i> Laveste terskel / korteste rute</label>
<label class="toggle"><input type="checkbox" id="showReference" checked><i class="swatch reference"></i> Rett referanselinje</label>
<label class="toggle"><input type="checkbox" id="showAlternate"><i class="swatch alternate"></i> Dal nær (Indigo Valley Line)</label>
<label class="toggle"><input type="checkbox" id="showPlateau" checked> Kritisk terskelområde</label>
<div class="map-controls">
<label class="toggle"><input type="checkbox" id="showBasemap" checked> Vis bakgrunnskart (nett)</label>
<label for="basemap" class="small">Velg bakgrunnskart</label>
<select id="basemap" aria-label="Bakgrunnskart"><option value="arcgis">ArcGIS World Imagery</option><option value="osm">OpenStreetMap</option></select>
<label for="basemapOpacity" class="range-label"><span>Gjennomsiktighet</span><span id="basemapOpacityValue">65 %</span></label>
<input type="range" id="basemapOpacity" min="0" max="100" value="65" aria-label="Bakgrunnskartets gjennomsiktighet">
</div>
<h2>Høydeforstørrelse – kun visning</h2>
<select id="exaggeration" aria-label="Høydeforstørrelse"><option value="1">1× – riktig vertikal skala</option><option value="3" selected>3× – tydeligere terreng</option><option value="6">6× – kraftig forstørrelse</option></select>
<p class="small" id="displayInfo"></p>
<h2>Punktinformasjon</h2><div id="picked">Klikk i terrenget eller beveg pekeren langs profilen nederst.</div>
<div class="warn"><strong>Avgrensning og høyder</strong>Resultatet gjelder rasterutsnittet. Rå rasterhøyder er ikke transformert til ellipsoidehøyder; absolutt høydeplassering på jordkloden er derfor ikke geodetisk korrigert. En lav terrengpassasje er ikke i seg selv bevis for tidligere vannføring.</div>
<div class="links"><a href="RESULTAT.md" target="_blank">Resultat / metode</a><a href="summary.json" target="_blank">Metadata</a><a href="routes.geojson" download>GeoJSON</a><a href="route_all_cells.csv" download>Alle rutepunkter</a><a href="reference_27_NEW_points.csv" download>27 nye kontrollpunkter</a><a href="profiles.png" target="_blank">Profilfigur</a><a href="saddle_section.png" target="_blank">Tverrprofil</a></div>
</aside>
<div id="globe"></div><div id="status">Laster Cesium og lokale analyseresultater …</div><div id="legend">Eget DTM · ingen Cesium World Terrain · høydeforstørrelse 3×</div>
<section id="profilePanel"><div class="profilehead"><label class="profile-picker" for="profileRoute"><span>Høydeprofil</span><select id="profileRoute" aria-label="Velg høydekurve"></select></label><span id="profileReadout">Pek på profilen for å finne et sted på kartet</span></div><canvas id="profile"></canvas></section>
<script>
"use strict";
(async()=>{
 const $=id=>document.getElementById(id), status=$('status');
 const fail=e=>{console.error(e);status.textContent='Kunne ikke starte visningen: '+(e.message||e)+'.';};
 try{
  if(!window.Cesium)throw new Error('CesiumJS ble ikke lastet fra CDN');
  const embedded=window.PASSAGE_EMBEDDED;
  const json=async url=>{const r=await fetch(url);if(!r.ok)throw new Error(url+': HTTP '+r.status);return r.json();};
  const data=embedded?embedded.data:await json('viewer_data.json'), summary=data.summary, meta=data.terrain;
  let bytes;
  if(embedded){
   const binary=atob(embedded.terrain_base64),decoded=new Uint8Array(binary.length);
   for(let i=0;i<binary.length;i++)decoded[i]=binary.charCodeAt(i);bytes=decoded.buffer;
  }else{
   const rawResponse=await fetch(meta.raw_file);if(!rawResponse.ok)throw new Error('Mangler terrengfil');bytes=await rawResponse.arrayBuffer();
  }
  if(bytes.byteLength!==meta.width*meta.height*4)throw new Error('Terrengfilen har feil lengde');
  let heights=new Float32Array(bytes);
  if(new Uint8Array(new Uint16Array([1]).buffer)[0]!==1){const dv=new DataView(bytes);heights=Float32Array.from({length:meta.width*meta.height},(_,i)=>dv.getFloat32(4*i,true));}
  const radians=Cesium.Math.toRadians, degrees=Cesium.Math.toDegrees;
  const rectangle=Cesium.Rectangle.fromDegrees(meta.west,meta.south,meta.east,meta.north);
  function sample(lon,lat){
   if(lon<meta.west||lon>meta.east||lat<meta.south||lat>meta.north)return undefined;
   const col=Math.min(meta.width-1,Math.max(0,(lon-meta.west)/(meta.east-meta.west)*meta.width-.5));
   const row=Math.min(meta.height-1,Math.max(0,(meta.north-lat)/(meta.north-meta.south)*meta.height-.5));
   const c=Math.floor(col),r=Math.floor(row),c1=Math.min(c+1,meta.width-1),r1=Math.min(r+1,meta.height-1);
   const fc=col-c,fr=row-r;
   const vals=[heights[r*meta.width+c],heights[r*meta.width+c1],heights[r1*meta.width+c],heights[r1*meta.width+c1]];
   const weights=[(1-fc)*(1-fr),fc*(1-fr),(1-fc)*fr,fc*fr];
   let sum=0,w=0;for(let i=0;i<4;i++){if(Number.isFinite(vals[i])){sum+=weights[i]*vals[i];w+=weights[i];}}
   return w>0?sum/w:undefined;
  }
  let exaggeration=3;
  const scheme=new Cesium.GeographicTilingScheme();
  function makeTerrain(){
   const factor=exaggeration, tileSide=65;
   return new Cesium.CustomHeightmapTerrainProvider({width:tileSide,height:tileSide,tilingScheme:scheme,
    credit:'Terreng: opplastet '+summary.input.source_filename+' · rå rasterhøyder, lokal visningsreferanse',
    callback:(x,y,level)=>{
     const rect=scheme.tileXYToRectangle(x,y,level), h=new Float32Array(tileSide*tileSide);
     if(!Cesium.Rectangle.intersection(rect,rectangle))return h;
     for(let r=0;r<tileSide;r++){
      const lat=degrees(rect.north-r/(tileSide-1)*(rect.north-rect.south));
      for(let c=0;c<tileSide;c++){
       const lon=degrees(rect.west+c/(tileSide-1)*(rect.east-rect.west));
       const z=sample(lon,lat);h[r*tileSide+c]=Number.isFinite(z)?z*factor:0;
      }
     }return h;
    }});
  }
  Cesium.Ion.defaultAccessToken='';
  const viewer=new Cesium.Viewer('globe',{terrainProvider:makeTerrain(),baseLayer:false,
   baseLayerPicker:false,geocoder:false,animation:false,timeline:false,sceneModePicker:false,
   navigationHelpButton:false,homeButton:false,fullscreenButton:true,infoBox:false,
   selectionIndicator:false,skyBox:false,requestRenderMode:true,maximumRenderTimeChange:Infinity});
  // Expose inspection handles for integration/testing; no API tokens or private data.
  window.passageViewer=viewer;window.passageData=data;window.passageSample=sample;
  viewer.scene.globe.enableLighting=false;
  viewer.scene.globe.depthTestAgainstTerrain=false;
  viewer.scene.globe.maximumScreenSpaceError=2;
  viewer.scene.screenSpaceCameraController.minimumZoomDistance=60;
  viewer.scene.backgroundColor=Cesium.Color.fromCssColorString('#111b27');
  const embeddedAssetUrls={};
  function asset(name){
   if(!embedded)return name;if(embeddedAssetUrls[name])return embeddedAssetUrls[name];
   const encoded=embedded.images[name].split(',')[1],binary=atob(encoded),bytes=new Uint8Array(binary.length);
   for(let i=0;i<binary.length;i++)bytes[i]=binary.charCodeAt(i);
   return embeddedAssetUrls[name]=URL.createObjectURL(new Blob([bytes],{type:'image/png'}));
  }
  const hill=viewer.imageryLayers.addImageryProvider(await Cesium.SingleTileImageryProvider.fromUrl(asset('hillshade.png'),{rectangle}));
  const connectedLayers=[];
  for(let i=0;i<meta.levels.length;i++){
   const p=await Cesium.SingleTileImageryProvider.fromUrl(asset(meta.levels[i].file),{rectangle});
   const layer=viewer.imageryLayers.addImageryProvider(p);layer.alpha=.65;layer.show=i===2;connectedLayers.push(layer);
   const option=document.createElement('option');option.value=i;option.textContent=meta.levels[i].label;$('level').appendChild(option);
  }
  $('level').value='2';
  function updateLevel(){
   const selected=Number($('level').value), show=$('showConnected').checked;
   connectedLayers.forEach((layer,i)=>layer.show=show&&i===selected);
   const connected=meta.levels[selected].connected;
   $('connectivity').textContent=connected?'Forbindelse mellom alle oppgitte punkter ved dette nivået.':'Ingen sammenhengende forbindelse mellom alle punktene ved dette nivået.';
   $('connectivity').style.borderLeftColor=connected?'#29d4c5':'#ffb65c';viewer.scene.requestRender();
  }
  $('level').onchange=updateLevel;$('showConnected').onchange=updateLevel;updateLevel();
  $('ceiling').textContent=summary.minimum_ceiling_native_m.toFixed(1)+' m';
  $('length').textContent=(summary.path.length_m/1000).toFixed(1)+' km';
  $('directMax').textContent=summary.straight_reference.max_height_native_m.toFixed(1)+' m';
  $('resolution').textContent=summary.input.pixel_size_m[0]+' m';
  $('displayInfo').textContent='Analysen: '+summary.input.pixel_size_m[0]+' m originalraster. 3D-visning: ca. '+meta.approximate_display_spacing_m.toFixed(0)+' m reprojisert visningsraster. Profiltallene påvirkes ikke av forstørrelsen.';
  const colors=[Cesium.Color.fromCssColorString('#29d4c5'),Cesium.Color.fromCssColorString('#ffb65c'),Cesium.Color.fromCssColorString('#ce9dff')];
  const routeEntities=[],plateauEntities=[],markers=[];
  function positions(points,lift=8){return points.map(p=>Cesium.Cartesian3.fromDegrees(p[0],p[1],(p[2]||0)*exaggeration+lift));}
  for(let i=0;i<data.routes.length;i++){
   const route=data.routes[i];
   const ent=viewer.entities.add({name:route.name,show:i!==2,polyline:{positions:positions(route.points),
    width:i===0?4:2.5,arcType:Cesium.ArcType.NONE,material:colors[i],depthFailMaterial:colors[i]}});
   routeEntities.push(ent);
  }
  const plateauJSON=embedded?embedded.plateaus:await json('critical_plateaus.geojson');
  const critical=summary.critical_plateaus.find(p=>p.direct_bridge)||summary.critical_plateaus[0];
  for(const f of plateauJSON.features){
   if(!f.properties.direct_bridge)continue;
   const rings=f.geometry.type==='Polygon'?[f.geometry.coordinates]:f.geometry.coordinates;
   for(const polygon of rings){
    const hierarchy=new Cesium.PolygonHierarchy(Cesium.Cartesian3.fromDegreesArray(polygon[0].flat()));
    const ent=viewer.entities.add({name:'Kritisk terskelplatå',polygon:{hierarchy,height:summary.minimum_ceiling_native_m*exaggeration+4,
      material:Cesium.Color.fromCssColorString('#ffe675').withAlpha(.8),outline:true,outlineColor:Cesium.Color.WHITE}});
    ent._payload=f.properties;plateauEntities.push(ent);
   }
  }
  function addMarker(lon,lat,z,label,payload,color=Cesium.Color.WHITE){
   const ent=viewer.entities.add({name:label,position:Cesium.Cartesian3.fromDegrees(lon,lat,z*exaggeration+16),
    point:{pixelSize:9,color,outlineColor:Cesium.Color.BLACK,outlineWidth:2,disableDepthTestDistance:Infinity},
    label:{text:label,font:'14px sans-serif',fillColor:color,style:Cesium.LabelStyle.FILL_AND_OUTLINE,
     outlineWidth:4,outlineColor:Cesium.Color.BLACK,pixelOffset:new Cesium.Cartesian2(10,-18),
     horizontalOrigin:Cesium.HorizontalOrigin.LEFT,disableDepthTestDistance:Infinity}});
   ent._payload=payload;markers.push({ent,lon,lat,z});return ent;
  }
  summary.points.forEach((p,i)=>addMarker(p.cell_center_longitude,p.cell_center_latitude,p.height_native_m,
    i===0?'Heradsbygd':(i===summary.points.length-1?'Kolomoen':'Via '+i),p));
  let saddleMarker;
  if(critical)saddleMarker=addMarker(critical.longitude,critical.latitude,critical.height_native_m,
    'Terskel '+critical.height_native_m.toFixed(1)+' m',critical,Cesium.Color.fromCssColorString('#ffe675'));
  const hover=viewer.entities.add({show:false,position:Cesium.Cartesian3.ZERO,
    point:{pixelSize:10,color:Cesium.Color.WHITE,outlineColor:Cesium.Color.fromCssColorString('#29d4c5'),outlineWidth:3,disableDepthTestDistance:Infinity}});
  const toggleRoute=(id,index)=>{$(id).onchange=()=>{if(routeEntities[index])routeEntities[index].show=$(id).checked;viewer.scene.requestRender();};};
  toggleRoute('showRoute',0);toggleRoute('showReference',1);toggleRoute('showAlternate',2);
  if(!routeEntities[2]){$('showAlternate').disabled=true;}
  $('showPlateau').onchange=()=>{plateauEntities.forEach(e=>e.show=$('showPlateau').checked);if(saddleMarker)saddleMarker.show=$('showPlateau').checked;viewer.scene.requestRender();};
  $('exaggeration').onchange=()=>{
   exaggeration=Number($('exaggeration').value);viewer.terrainProvider=makeTerrain();
   routeEntities.forEach((e,i)=>e.polyline.positions=positions(data.routes[i].points));
   markers.forEach(m=>m.ent.position=Cesium.Cartesian3.fromDegrees(m.lon,m.lat,m.z*exaggeration+16));
   plateauEntities.forEach(e=>e.polygon.height=summary.minimum_ceiling_native_m*exaggeration+4);
   hover.show=false;$('legend').textContent='Eget DTM · ingen Cesium World Terrain · høydeforstørrelse '+exaggeration+'×';viewer.scene.requestRender();
  };
  const basemapLayers={};
  const basemapOptions={
   arcgis:{url:'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',maximumLevel:19,credit:'Esri, Maxar, Earthstar Geographics, and the GIS User Community',label:'ArcGIS World Imagery'},
   osm:{url:'https://tile.openstreetmap.org/{z}/{x}/{y}.png',maximumLevel:19,credit:'© OpenStreetMap contributors',label:'OpenStreetMap'}
  };
  function getBasemap(name){
   if(basemapLayers[name])return basemapLayers[name];
   const option=basemapOptions[name];
   const provider=new Cesium.UrlTemplateImageryProvider({url:option.url,maximumLevel:option.maximumLevel,credit:option.credit});
   provider.errorEvent.addEventListener(()=>{status.textContent=option.label+' kunne ikke lastes. Egne terrengdata er fortsatt tilgjengelige.';});
   const layer=viewer.imageryLayers.addImageryProvider(provider,0);layer.show=false;basemapLayers[name]=layer;return layer;
  }
  function updateBasemap(){
   const enabled=$('showBasemap').checked,name=$('basemap').value,opacity=Number($('basemapOpacity').value)/100;
   if(enabled)getBasemap(name);
   Object.entries(basemapLayers).forEach(([key,layer])=>{layer.show=enabled&&key===name;layer.alpha=opacity;});
   hill.alpha=enabled?1-opacity*.65:1;
   $('basemap').disabled=!enabled;$('basemapOpacity').disabled=!enabled;
   $('basemapOpacityValue').textContent=Math.round(opacity*100)+' %';viewer.scene.requestRender();
  }
  $('showBasemap').onchange=updateBasemap;$('basemap').onchange=updateBasemap;$('basemapOpacity').oninput=updateBasemap;updateBasemap();
  const overviewPoints=data.routes[0].points;
  const profileNames=['Laveste passasje','Rett referanselinje','Dal nær (Indigo Valley Line)'];
  data.routes.forEach((route,i)=>{const option=document.createElement('option');option.value=i;option.textContent=profileNames[i]||route.name;$('profileRoute').appendChild(option);});
  let profileIndex=0,p=overviewPoints;
  const lonMin=Math.min(...overviewPoints.map(v=>v[0])),lonMax=Math.max(...overviewPoints.map(v=>v[0]));
  const latMin=Math.min(...overviewPoints.map(v=>v[1])),latMax=Math.max(...overviewPoints.map(v=>v[1]));
  const centreLon=(lonMin+lonMax)/2,centreLat=(latMin+latMax)/2;
  const sceneSpan=Math.max((lonMax-lonMin)*111320*Math.cos(radians(centreLat)),(latMax-latMin)*111320);
  function look(lon,lat,range,pitch){
   const h=sample(lon,lat)||0;
   viewer.camera.lookAt(Cesium.Cartesian3.fromDegrees(lon,lat,h*exaggeration),new Cesium.HeadingPitchRange(0,radians(pitch),range));
   viewer.camera.lookAtTransform(Cesium.Matrix4.IDENTITY);viewer.scene.requestRender();
  }
  $('overview').onclick=()=>look(centreLon,centreLat,sceneSpan*1.65,-58);
  $('top').onclick=()=>look(centreLon,centreLat,sceneSpan*1.65,-90);
  $('saddle').onclick=()=>{if(critical)look(critical.longitude,critical.latitude,2200,-58);};
  $('overview').click();
  const clicker=new Cesium.ScreenSpaceEventHandler(viewer.scene.canvas);
  clicker.setInputAction(event=>{
   const picked=viewer.scene.pick(event.position),payload=picked?.id?._payload;
   if(payload){$('picked').textContent=JSON.stringify(payload,null,2).replace(/[{}\"]/g,'');return;}
   const ray=viewer.camera.getPickRay(event.position),hit=ray&&viewer.scene.globe.pick(ray,viewer.scene);
   if(!hit)return;const carto=Cesium.Cartographic.fromCartesian(hit),lon=degrees(carto.longitude),lat=degrees(carto.latitude),z=sample(lon,lat);
   $('picked').textContent=lat.toFixed(6)+' N, '+lon.toFixed(6)+' E\n'+(Number.isFinite(z)?'Visningsraster, interpolert: ca. '+z.toFixed(1)+' native m\nIkke en ny prøve av originalrasteret.':'Utenfor gyldig terrengutsnitt.');
  },Cesium.ScreenSpaceEventType.LEFT_CLICK);
  // Interactive profile is plotted from every original route cell; no display resampling of heights.
  const canvas=$('profile'),ctx=canvas.getContext('2d');let chartGeometry,hoverIndex=-1;
  function drawProfile(){
   const dpr=window.devicePixelRatio||1,w=canvas.clientWidth,h=canvas.clientHeight;
   canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);
   const pad={l:43,r:14,t:9,b:25},plotW=w-pad.l-pad.r,plotH=h-pad.t-pad.b;
   const zMin=Math.floor((Math.min(...p.map(v=>v[2]))-5)/10)*10,zMax=Math.ceil((Math.max(...p.map(v=>v[2]))+8)/10)*10;
   const xmax=p[p.length-1][3],xx=d=>pad.l+d/xmax*plotW,yy=z=>pad.t+(zMax-z)/(zMax-zMin)*plotH;
   chartGeometry={xx,yy,xmax,left:pad.l,width:plotW};
   ctx.clearRect(0,0,w,h);ctx.font='10px system-ui';ctx.textAlign='right';ctx.fillStyle='#a8bcc8';ctx.strokeStyle='#314453';ctx.lineWidth=1;
   for(let z=zMin;z<=zMax;z+=20){ctx.beginPath();ctx.moveTo(pad.l,yy(z));ctx.lineTo(w-pad.r,yy(z));ctx.stroke();ctx.fillText(z+' m',pad.l-6,yy(z)+3);}
   ctx.textAlign='center';for(let d=0;d<=xmax;d+=5000)ctx.fillText((d/1000)+' km',xx(d),h-6);
   ctx.setLineDash([5,4]);ctx.strokeStyle='#ffb65c';ctx.beginPath();ctx.moveTo(pad.l,yy(summary.minimum_ceiling_native_m));ctx.lineTo(w-pad.r,yy(summary.minimum_ceiling_native_m));ctx.stroke();ctx.setLineDash([]);
   ctx.strokeStyle=['#29d4c5','#ffb65c','#8f7cff'][profileIndex]||'#29d4c5';ctx.lineWidth=1.7;ctx.beginPath();p.forEach((v,i)=>i?ctx.lineTo(xx(v[3]),yy(v[2])):ctx.moveTo(xx(v[3]),yy(v[2])));ctx.stroke();
   if(hoverIndex>=0){const v=p[hoverIndex];ctx.fillStyle='#ffffff';ctx.beginPath();ctx.arc(xx(v[3]),yy(v[2]),4,0,2*Math.PI);ctx.fill();}
  }
  canvas.addEventListener('mousemove',event=>{
   if(!chartGeometry)return;const rect=canvas.getBoundingClientRect(),g=chartGeometry;
   const distance=Math.max(0,Math.min(g.xmax,(event.clientX-rect.left-g.left)/g.width*g.xmax));
   let lo=0,hi=p.length-1;while(lo<hi){const mid=(lo+hi)>>1;if(p[mid][3]<distance)lo=mid+1;else hi=mid;}
   hoverIndex=lo;const v=p[lo];$('profileReadout').textContent=(v[3]/1000).toFixed(2)+' km · '+v[2].toFixed(1)+' m · '+v[1].toFixed(5)+' N, '+v[0].toFixed(5)+' E';
   hover.position=Cesium.Cartesian3.fromDegrees(v[0],v[1],v[2]*exaggeration+22);hover.show=true;drawProfile();viewer.scene.requestRender();
  });
  canvas.addEventListener('mouseleave',()=>{hoverIndex=-1;hover.show=false;drawProfile();viewer.scene.requestRender();});
  $('profileRoute').onchange=()=>{profileIndex=Number($('profileRoute').value);p=data.routes[profileIndex].points;hoverIndex=-1;hover.show=false;$('profileReadout').textContent='Pek på profilen for å finne et sted på kartet';drawProfile();viewer.scene.requestRender();};
  new ResizeObserver(drawProfile).observe(canvas);drawProfile();
  status.textContent='';window.passageReady=true;
 }catch(error){fail(error);window.passageError=String(error);}
})();
</script>
</body></html>
'''

if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, rasterio.errors.RasterioError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        sys.exit(2)
