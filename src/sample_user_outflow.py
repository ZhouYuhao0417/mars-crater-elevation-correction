from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np


MARS_RADIUS_KM = 3396.19


def project(point: list[float], lon0: float, lat0: float) -> tuple[float, float]:
    lon, lat = point
    return (
        math.radians(lon - lon0) * MARS_RADIUS_KM * math.cos(math.radians(lat0)),
        math.radians(lat - lat0) * MARS_RADIUS_KM,
    )


def unproject(point: tuple[float, float], lon0: float, lat0: float) -> tuple[float, float]:
    x, y = point
    return (
        lon0 + math.degrees(x / (MARS_RADIUS_KM * math.cos(math.radians(lat0)))),
        lat0 + math.degrees(y / MARS_RADIUS_KM),
    )


def densify(coordinates: list[list[float]], spacing_km: float = 0.1):
    lon0 = sum(point[0] for point in coordinates) / len(coordinates)
    lat0 = sum(point[1] for point in coordinates) / len(coordinates)
    xy = [project(point, lon0, lat0) for point in coordinates]
    segment_lengths = [math.dist(a, b) for a, b in zip(xy, xy[1:])]
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    samples = np.arange(0.0, cumulative[-1], spacing_km)
    if samples[-1] < cumulative[-1]:
        samples = np.append(samples, cumulative[-1])
    sampled_xy = []
    for distance in samples:
        index = min(np.searchsorted(cumulative, distance, side="right") - 1, len(xy) - 2)
        fraction = (distance - cumulative[index]) / max(segment_lengths[index], 1e-12)
        sampled_xy.append(
            (
                xy[index][0] + fraction * (xy[index + 1][0] - xy[index][0]),
                xy[index][1] + fraction * (xy[index + 1][1] - xy[index][1]),
            )
        )
    lonlat = [unproject(point, lon0, lat0) for point in sampled_xy]
    return samples, np.asarray([p[0] for p in lonlat]), np.asarray([p[1] for p in lonlat])


def interp_regular(lons, lats, grid, query_lon, query_lat):
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    values = np.asarray(grid, dtype=float)
    if lats[0] > lats[-1]:
        lats = lats[::-1]
        values = values[::-1, :]
    result = np.full(query_lon.shape, np.nan, dtype=float)
    for k, (lon, lat) in enumerate(zip(query_lon, query_lat)):
        if lon < lons[0] or lon > lons[-1] or lat < lats[0] or lat > lats[-1]:
            continue
        j = min(max(np.searchsorted(lons, lon) - 1, 0), len(lons) - 2)
        i = min(max(np.searchsorted(lats, lat) - 1, 0), len(lats) - 2)
        x0, x1 = lons[j], lons[j + 1]
        y0, y1 = lats[i], lats[i + 1]
        block = values[i : i + 2, j : j + 2]
        if not np.isfinite(block).all():
            continue
        tx = (lon - x0) / (x1 - x0)
        ty = (lat - y0) / (y1 - y0)
        result[k] = (
            block[0, 0] * (1 - tx) * (1 - ty)
            + block[0, 1] * tx * (1 - ty)
            + block[1, 0] * (1 - tx) * ty
            + block[1, 1] * tx * ty
        )
    return result


def load_mola(csv_name: str):
    rows = []
    with Path(csv_name).open("r", encoding="utf-8-sig", newline="") as handle:
        for record in csv.DictReader(handle):
            rows.append(record)
    lons = sorted({float(row["longitude_e"]) for row in rows})
    lats = sorted({float(row["latitude_n"]) for row in rows})
    lon_index = {value: index for index, value in enumerate(lons)}
    lat_index = {value: index for index, value in enumerate(lats)}
    raw = np.full((len(lats), len(lons)), np.nan)
    corrected = np.full_like(raw, np.nan)
    for row in rows:
        i = lat_index[float(row["latitude_n"])]
        j = lon_index[float(row["longitude_e"])]
        raw[i, j] = float(row["raw_elevation_m"])
        corrected[i, j] = float(row["corrected_residual_m"])
    return np.asarray(lons), np.asarray(lats), raw, corrected


def main(outflow_geojson_name: str, intersection_json_name: str, mola_csv_name: str, ctx_npz_name: str, output_dir_name: str):
    output_dir = Path(output_dir_name)
    output_dir.mkdir(parents=True, exist_ok=True)
    outflow = json.loads(Path(outflow_geojson_name).read_text(encoding="utf-8"))
    coordinates = outflow["features"][0]["geometry"]["coordinates"]
    distance, lon, lat = densify(coordinates)
    intersection = json.loads(Path(intersection_json_name).read_text(encoding="utf-8"))
    rim_distance = intersection["rim_intersections"][0]["distance_along_outflow_km"]

    mola_lon, mola_lat, mola_raw_grid, mola_corrected_grid = load_mola(mola_csv_name)
    mola_raw = interp_regular(mola_lon, mola_lat, mola_raw_grid, lon, lat)
    mola_corrected = interp_regular(mola_lon, mola_lat, mola_corrected_grid, lon, lat)

    ctx = np.load(ctx_npz_name)
    ctx_raw = interp_regular(ctx["longitude"], ctx["latitude"], ctx["elevation"], lon, lat)
    ctx_corrected = interp_regular(ctx["longitude"], ctx["latitude"], ctx["corrected"], lon, lat)

    output_csv = output_dir / "JMARS_Profile2_outflow_axis_elevation_profile.csv"
    with output_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "distance_from_start_km",
                "distance_from_user_rim_km",
                "longitude_e_deg",
                "latitude_n_deg",
                "mola_raw_m",
                "mola_corrected_residual_m",
                "ctx_raw_m",
                "ctx_tilt_corrected_m",
            ]
        )
        for values in zip(distance, distance - rim_distance, lon, lat, mola_raw, mola_corrected, ctx_raw, ctx_corrected):
            writer.writerow(["nan" if not np.isfinite(value) else f"{value:.6f}" for value in values])

    valid_ctx = np.isfinite(ctx_raw)
    downstream_ctx = valid_ctx & (distance >= rim_distance)
    if downstream_ctx.any():
        first_ctx_distance = float(np.min(distance[downstream_ctx] - rim_distance))
        last_ctx_distance = float(np.max(distance[downstream_ctx] - rim_distance))
    else:
        first_ctx_distance = last_ctx_distance = math.nan
    summary = {
        "sample_spacing_km": 0.1,
        "sample_count": int(len(distance)),
        "user_rim_crossing_distance_from_start_km": rim_distance,
        "user_rim_crossing_lon_lat": [
            intersection["rim_intersections"][0]["longitude_e_deg"],
            intersection["rim_intersections"][0]["latitude_n_deg"],
        ],
        "ctx_valid_samples": int(valid_ctx.sum()),
        "ctx_valid_fraction": float(valid_ctx.mean()),
        "ctx_coverage_downstream_from_rim_km": [first_ctx_distance, last_ctx_distance],
        "limitation": "The available CTX DTM begins downstream of the user-mapped rim crossing; it cannot resolve the breach crest itself.",
    }
    (output_dir / "JMARS_Profile2_outflow_axis_elevation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:])
