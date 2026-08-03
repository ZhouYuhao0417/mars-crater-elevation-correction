# Mars crater elevation correction and outlet hypothesis

This repository contains the reproducible analysis package for an irregular-rim Martian crater elevation correction and suspected outflow-path assessment near 18.3°N, 75.8°E.

The circular-rim approximation used in early exploratory work has been discarded. All final geometry uses the crater rim and suspected outflow axis manually traced in JMARS.

## Main result

The user-mapped outflow axis starts inside the crater, ends outside it, and crosses the irregular rim once at 18.106545°N, 76.040457°E. Across three robust regional detrending models, the crossing is 29.7–95.5 m below the regional background surface and falls within the lowest 5.46–8.73% of the mapped rim elevations. The available CTX DTM begins 8.45 km downstream from the crossing, so it cannot resolve the breach crest or independently prove a dam-breach origin.

## Repository layout

- `data/geometry/` — JMARS-exported rim and outflow geometry.
- `data/derived/` — lightweight, analysis-ready MOLA subset and extracted CTX profile.
- `data/mola/` — partial NASA PDS MOLA MEGDR rows used in the original extraction.
- `src/` — reproducible analysis, figure and interactive-3D builders.
- `results/` — sensitivity tables, figures and machine-readable summary.
- `interactive/` — standalone interactive 3D terrain viewer.

## Reproduce

Create an environment and install dependencies:

```bash
python -m pip install -r requirements.txt
python src/paper_grade_analysis.py
python src/generate_paper_figures.py
python src/build_updated_3d.py
```

The analysis rebuilds `results/analysis_summary.json`, the sensitivity table and profile tables. Figure generation recreates the two publication figures. The 3D builder writes `interactive/crater-corrected-3d-paper.html`; `interactive/crater-corrected-3d.html` is a ready-to-open standalone viewer.

## Data provenance and limits

- Regional elevations: [NASA PDS MOLA MEGDR 128 pixel/degree](https://pds-geosciences.wustl.edu/missions/mgs/megdr.html).
- Local DTM: [USGS Astrogeology CTX DTMs](https://stac.astrogeology.usgs.gov/docs/data/mars/ctxdtms/), product `P05_002809_1975_XI_17N283W__P13_006000_1974_XI_17N283W`.
- Geometry: manual JMARS Profile traces; see `data/geometry/`.

The source CTX raster and processing cache are intentionally not redistributed because they are large derived downloads. The extracted profile, QA metrics and product identifier are included. Treat the 3D model as a regional MOLA-resolution interpretation: it is useful for visualizing rim asymmetry and the mapped outflow relationship, but not for measuring the unresolved outlet morphology.

## Interpretation boundary

These results support a low-rim outlet and outward-descending terrain configuration consistent with an outflow or breach hypothesis. They do **not** establish dam breach as the unique origin. A continuous high-resolution DTM over the mapped crossing and its immediate upstream/downstream reaches is required for breach morphology, paleolake level, storage or discharge estimates.

## License

Code is released under the MIT License. NASA PDS and USGS data retain their original provenance and terms.
