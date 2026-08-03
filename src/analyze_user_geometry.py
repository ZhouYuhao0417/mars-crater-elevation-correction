from __future__ import annotations

import json
import math
from pathlib import Path
import sys


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


def cross(a: tuple[float, float], b: tuple[float, float]) -> float:
    return a[0] * b[1] - a[1] * b[0]


def intersection(
    p: tuple[float, float],
    p2: tuple[float, float],
    q: tuple[float, float],
    q2: tuple[float, float],
) -> tuple[float, float, tuple[float, float]] | None:
    r = (p2[0] - p[0], p2[1] - p[1])
    s = (q2[0] - q[0], q2[1] - q[1])
    denominator = cross(r, s)
    if abs(denominator) < 1e-12:
        return None
    qp = (q[0] - p[0], q[1] - p[1])
    t = cross(qp, s) / denominator
    u = cross(qp, r) / denominator
    if -1e-10 <= t <= 1 + 1e-10 and -1e-10 <= u <= 1 + 1e-10:
        return t, u, (p[0] + t * r[0], p[1] + t * r[1])
    return None


def point_in_polygon(point: tuple[float, float], polygon: list[tuple[float, float]]) -> bool:
    x, y = point
    inside = False
    for i in range(len(polygon) - 1):
        x1, y1 = polygon[i]
        x2, y2 = polygon[i + 1]
        if (y1 > y) != (y2 > y):
            x_cross = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < x_cross:
                inside = not inside
    return inside


def main(rim_name: str, outflow_name: str, output_name: str) -> None:
    rim_geojson = json.loads(Path(rim_name).read_text(encoding="utf-8"))
    outflow_geojson = json.loads(Path(outflow_name).read_text(encoding="utf-8"))
    rim_lonlat = rim_geojson["features"][0]["geometry"]["coordinates"]
    outflow_lonlat = outflow_geojson["features"][0]["geometry"]["coordinates"]
    lon0 = sum(point[0] for point in rim_lonlat) / len(rim_lonlat)
    lat0 = sum(point[1] for point in rim_lonlat) / len(rim_lonlat)
    rim = [project(point, lon0, lat0) for point in rim_lonlat]
    if rim[0] != rim[-1]:
        rim.append(rim[0])
    outflow = [project(point, lon0, lat0) for point in outflow_lonlat]

    cumulative = [0.0]
    for a, b in zip(outflow, outflow[1:]):
        cumulative.append(cumulative[-1] + math.dist(a, b))

    hits = []
    for out_index, (a, b) in enumerate(zip(outflow, outflow[1:])):
        for rim_index, (c, d) in enumerate(zip(rim, rim[1:])):
            hit = intersection(a, b, c, d)
            if hit is None:
                continue
            t, u, xy = hit
            lon, lat = unproject(xy, lon0, lat0)
            along = cumulative[out_index] + t * math.dist(a, b)
            hits.append(
                {
                    "longitude_e_deg": lon,
                    "latitude_n_deg": lat,
                    "distance_along_outflow_km": along,
                    "outflow_segment": [out_index + 1, out_index + 2],
                    "rim_segment": [rim_index + 1, rim_index + 2],
                }
            )

    hits.sort(key=lambda item: item["distance_along_outflow_km"])
    result = {
        "rim": str(rim_name),
        "outflow": str(outflow_name),
        "outflow_length_km": cumulative[-1],
        "start_inside_user_rim": point_in_polygon(outflow[0], rim),
        "end_inside_user_rim": point_in_polygon(outflow[-1], rim),
        "rim_intersection_count": len(hits),
        "rim_intersections": hits,
        "note": "Intersections use the user-traced rim closed by joining its last vertex to its first vertex.",
    }
    Path(output_name).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:])
