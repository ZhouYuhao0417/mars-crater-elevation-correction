from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re
import sys


MARS_RADIUS_KM = 3396.19
POINT_RE = re.compile(
    r"^Point(?P<index>\d+):\s*\[(?P<lat>-?\d+(?:\.\d+)?)°N,\s*(?P<lon>-?\d+(?:\.\d+)?)°E\]"
)


def distance_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    lon1, lat1 = map(math.radians, a)
    lon2, lat2 = map(math.radians, b)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    h = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * MARS_RADIUS_KM * math.asin(min(1.0, math.sqrt(h)))


def main(
    source_name: str,
    output_dir_name: str,
    output_stem: str = "JMARS_Profile1_crater_rim",
    feature_name: str = "crater_rim_user_trace",
    interpretation: str = "User-traced provisional crater rim; do not assume a circle or use as an outlet centerline.",
) -> None:
    source = Path(source_name)
    output_dir = Path(output_dir_name)
    output_dir.mkdir(parents=True, exist_ok=True)

    points: list[tuple[int, float, float]] = []
    for line in source.read_text(encoding="utf-8-sig").splitlines():
        match = POINT_RE.match(line.strip())
        if match:
            points.append(
                (
                    int(match.group("index")),
                    float(match.group("lon")),
                    float(match.group("lat")),
                )
            )
    if len(points) < 2:
        raise RuntimeError("No JMARS profile points were found")

    coordinates = [[lon, lat] for _, lon, lat in points]
    segment_lengths = [
        distance_km((points[i - 1][1], points[i - 1][2]), (points[i][1], points[i][2]))
        for i in range(1, len(points))
    ]
    cumulative = [0.0]
    for length in segment_lengths:
        cumulative.append(cumulative[-1] + length)

    geojson = {
        "type": "FeatureCollection",
        "name": output_stem,
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "name": feature_name,
                    "source": source.stem,
                    "status": "user-traced provisional geometry",
                    "longitude": "degrees east",
                    "latitude": "planetocentric north",
                    "vertex_count": len(points),
                    "length_km": round(cumulative[-1], 3),
                },
                "geometry": {"type": "LineString", "coordinates": coordinates},
            }
        ],
    }
    (output_dir / f"{output_stem}.geojson").write_text(
        json.dumps(geojson, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    with (output_dir / f"{output_stem}_vertices.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["vertex_index", "longitude_e_deg", "latitude_n_planetocentric_deg", "distance_along_trace_km"])
        for (index, lon, lat), along in zip(points, cumulative):
            writer.writerow([index, f"{lon:.6f}", f"{lat:.6f}", f"{along:.4f}"])

    summary = {
        "source": str(source),
        "vertex_count": len(points),
        "trace_length_km": cumulative[-1],
        "start_end_separation_km": distance_km(
            (points[0][1], points[0][2]), (points[-1][1], points[-1][2])
        ),
        "bbox_lon_lat": [
            min(p[1] for p in points),
            min(p[2] for p in points),
            max(p[1] for p in points),
            max(p[2] for p in points),
        ],
        "coordinate_convention": "east-positive longitude, planetocentric north latitude",
        "interpretation": interpretation,
    }
    (output_dir / f"{output_stem}_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:])
