from pathlib import Path
import json
import math

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
VIS = ROOT / "interactive"
OUT = VIS / "crater-corrected-3d-paper.html"


def load(path):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    data.dtype.names = tuple(name.lstrip("\ufeff") for name in data.dtype.names)
    return data


def xy(lon, lat, lon0, lat0):
    radius_km = 3396.19
    east = (np.asarray(lon) - lon0) * math.pi / 180 * radius_km * math.cos(math.radians(lat0))
    north = (np.asarray(lat) - lat0) * math.pi / 180 * radius_km
    return east, north


def bilinear(lons, lats, grid, lon, lat):
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    ix = np.clip(np.searchsorted(lons, lon) - 1, 0, len(lons) - 2)
    iy = np.clip(np.searchsorted(lats, lat) - 1, 0, len(lats) - 2)
    tx = (lon - lons[ix]) / (lons[ix + 1] - lons[ix])
    ty = (lat - lats[iy]) / (lats[iy + 1] - lats[iy])
    return (
        (1 - tx) * (1 - ty) * grid[iy, ix]
        + tx * (1 - ty) * grid[iy, ix + 1]
        + (1 - tx) * ty * grid[iy + 1, ix]
        + tx * ty * grid[iy + 1, ix + 1]
    )


def arr(values, decimals=3):
    return [
        None if not np.isfinite(value) else round(float(value), decimals)
        for value in np.asarray(values).ravel()
    ]


mola = load(ROOT / "data" / "derived" / "MOLA_subset_points.csv")
rim = load(ROOT / "data" / "geometry" / "JMARS_Profile1_crater_rim_vertices.csv")
flow = load(ROOT / "results" / "tables" / "outflow_corrected_profile.csv")
summary = json.loads((ROOT / "results" / "analysis_summary.json").read_text(encoding="utf-8"))

lon0 = float(np.mean(rim["longitude_e_deg"]))
lat0 = float(np.mean(rim["latitude_n_planetocentric_deg"]))

lons = np.unique(mola["longitude_e"])
lats = np.unique(mola["latitude_n"])
lons.sort()
lats.sort()

raw = np.full((len(lats), len(lons)), np.nan)
ix = np.searchsorted(lons, mola["longitude_e"])
iy = np.searchsorted(lats, mola["latitude_n"])
raw[iy, ix] = mola["raw_elevation_m"]

gx, gy = np.meshgrid(lons, lats)
xx, yy = xy(gx, gy, lon0, lat0)
residual_models = []
for model in summary["models"]:
    plane = (
        model["intercept_m"]
        + model["east_gradient_m_per_km"] * xx
        + model["north_gradient_m_per_km"] * yy
    )
    residual_models.append(raw - plane)
corr = np.nanmedian(np.stack(residual_models), axis=0)

step = 2
xs = xx[0, ::step]
ys = yy[::step, 0]
raw_s = raw[::step, ::step]
corr_s = corr[::step, ::step]

rlon = np.r_[rim["longitude_e_deg"], rim["longitude_e_deg"][0]]
rlat = np.r_[
    rim["latitude_n_planetocentric_deg"],
    rim["latitude_n_planetocentric_deg"][0],
]
rx, ry = xy(rlon, rlat, lon0, lat0)
rraw = bilinear(lons, lats, raw, rlon, rlat)
rim_residuals = []
for model in summary["models"]:
    plane = (
        model["intercept_m"]
        + model["east_gradient_m_per_km"] * rx
        + model["north_gradient_m_per_km"] * ry
    )
    rim_residuals.append(rraw - plane)
rcorr = np.nanmedian(np.stack(rim_residuals), axis=0)

fx, fy = xy(flow["longitude_e_deg"], flow["latitude_n_deg"], lon0, lat0)
clon = summary["outlet_lon_e_deg"]
clat = summary["outlet_lat_n_deg"]
ox, oy = xy([clon], [clat], lon0, lat0)
oraw = float(bilinear(lons, lats, raw, [clon], [clat])[0])
outlet_residuals = []
for model in summary["models"]:
    plane = (
        model["intercept_m"]
        + model["east_gradient_m_per_km"] * ox[0]
        + model["north_gradient_m_per_km"] * oy[0]
    )
    outlet_residuals.append(oraw - plane)

data = {
    "x": arr(xs),
    "y": arr(ys),
    "raw": [arr(row, 1) for row in raw_s],
    "corrected": [arr(row, 1) for row in corr_s],
    "rimX": arr(rx),
    "rimY": arr(ry),
    "rimRaw": arr(rraw, 1),
    "rimCorrected": arr(rcorr, 1),
    "flowX": arr(fx),
    "flowY": arr(fy),
    "flowRaw": arr(flow["mola_raw_m"], 1),
    "flowCorrected": arr(flow["corrected_median_m"], 1),
    "outlet": {
        "x": round(float(ox[0]), 3),
        "y": round(float(oy[0]), 3),
        "raw": round(oraw, 1),
        "corrected": round(float(np.median(outlet_residuals)), 1),
        "lon": clon,
        "lat": clat,
    },
    "stats": {
        "residualMin": round(summary["outlet_residual_range_m"][0], 1),
        "residualMax": round(summary["outlet_residual_range_m"][1], 1),
        "percentileMin": round(summary["outlet_percentile_range"][0], 2),
        "percentileMax": round(summary["outlet_percentile_range"][1], 2),
    },
}

j = json.dumps(data, ensure_ascii=True, separators=(",", ":"))
fragment = f'''<div id="crater-paper-3d" style="width:100%;">
  <div class="viz-controls" aria-label="3D terrain controls">
    <button type="button" class="btn btn-primary" id="cp3-corrected" aria-pressed="true">Detrended</button>
    <button type="button" class="btn" id="cp3-raw" aria-pressed="false">Raw elevation</button>
    <label class="form-label" for="cp3-exaggeration">Vertical exaggeration <span id="cp3-exaggeration-value">6x</span>
      <input class="form-range" id="cp3-exaggeration" type="range" min="1" max="12" step="1" value="6">
    </label>
    <label class="form-label" for="cp3-speed">Orbit speed <span id="cp3-speed-value">6x</span>
      <input class="form-range" id="cp3-speed" type="range" min="1" max="8" step="1" value="6">
    </label>
    <button type="button" class="btn" id="cp3-orbit" aria-pressed="true"><i data-lucide="pause" aria-hidden="true"></i> Pause orbit</button>
    <button type="button" class="btn btn-ghost" id="cp3-reset">Reset view</button>
  </div>
  <div id="cp3-detail" class="text-small text-muted" aria-live="polite">Mapped irregular rim, suspected outflow axis, and projected outlet. Auto-orbit is about 10 seconds per revolution.</div>
  <div id="cp3-plot" role="img" aria-label="Three-dimensional detrended Mars crater terrain. The white line is the user-mapped irregular rim, the colored line is the suspected outflow axis, and the diamond marks the projected outlet."></div>
  <div class="sr-only">Across three background-surface models, the projected outlet falls within the lowest 5.46 to 8.73 percent of mapped rim elevations. Switch between raw elevation and detrended terrain, drag to rotate, zoom, or pause the orbit.</div>
</div>
<style>
  #crater-paper-3d #cp3-plot {{ width:100%; height:600px; min-height:420px; }}
  #crater-paper-3d .viz-controls {{ margin-bottom:var(--spacing-sm); }}
  #crater-paper-3d #cp3-detail {{ margin-bottom:var(--spacing-xs); }}
  @media (max-width:520px) {{ #crater-paper-3d #cp3-plot {{ height:460px; }} }}
</style>
<script src="/vendor/plotly.min.js"></script>
<script>
(() => {{
  const root=document.getElementById('crater-paper-3d');
  const plot=root.querySelector('#cp3-plot');
  const correctedButton=root.querySelector('#cp3-corrected');
  const rawButton=root.querySelector('#cp3-raw');
  const exaggeration=root.querySelector('#cp3-exaggeration');
  const exaggerationValue=root.querySelector('#cp3-exaggeration-value');
  const speed=root.querySelector('#cp3-speed');
  const speedValue=root.querySelector('#cp3-speed-value');
  const orbitButton=root.querySelector('#cp3-orbit');
  const resetButton=root.querySelector('#cp3-reset');
  const detail=root.querySelector('#cp3-detail');
  const terrain={j};
  const themeColor=name=>{{ const p=document.createElement('span'); p.style.color=`var(${{name}})`; p.style.display='none'; root.appendChild(p); const v=getComputedStyle(p).color; p.remove(); return v; }};
  const colorscale=[[0,themeColor('--viz-series-1')],[0.3,themeColor('--viz-series-2')],[0.52,themeColor('--muted')],[0.76,themeColor('--viz-series-3')],[1,themeColor('--viz-series-4')]];
  const foreground=themeColor('--foreground'), muted=themeColor('--muted-foreground'), border=themeColor('--border'), background=themeColor('--background'), flowColor=themeColor('--viz-series-2'), outletColor=themeColor('--viz-series-5');
  const toKm=m=>m.map(row=>row.map(v=>v===null?null:v/1000));
  const zSurface={{raw:toKm(terrain.raw),corrected:toKm(terrain.corrected)}};
  const zRim={{raw:terrain.rimRaw.map(v=>v/1000),corrected:terrain.rimCorrected.map(v=>v/1000)}};
  const zFlow={{raw:terrain.flowRaw.map(v=>v/1000),corrected:terrain.flowCorrected.map(v=>v/1000)}};
  const finiteRange=m=>{{const a=m.flat().filter(Number.isFinite).sort((a,b)=>a-b);return [a[Math.floor(a.length*.02)],a[Math.floor(a.length*.98)]];}};
  const ranges={{raw:finiteRange(terrain.raw),corrected:finiteRange(terrain.corrected)}};
  let mode='corrected';
  const defaultCamera={{eye:{{x:1.55,y:-1.62,z:1.05}},center:{{x:0,y:0,z:-.08}},up:{{x:0,y:0,z:1}}}};
  let camera=JSON.parse(JSON.stringify(defaultCamera));
  let angle=Math.atan2(camera.eye.y,camera.eye.x), radius=Math.hypot(camera.eye.x,camera.eye.y), height=camera.eye.z;
  const reducedMotion=window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  let paused=reducedMotion, interactionUntil=0, last=performance.now(), lastUpdate=0;
  const surface={{type:'surface',x:terrain.x,y:terrain.y,z:zSurface.corrected,colorscale,cmin:ranges.corrected[0]/1000,cmax:ranges.corrected[1]/1000,colorbar:{{title:{{text:'Detrended residual (km)',font:{{color:foreground}}}},thickness:14,len:.66,tickfont:{{color:muted}}}},hovertemplate:'East offset %{{x:.1f}} km<br>North offset %{{y:.1f}} km<br>Detrended residual %{{z:.3f}} km<extra></extra>',lighting:{{ambient:.7,diffuse:.78,specular:.14,roughness:.82,fresnel:.08}},lightposition:{{x:-80,y:-110,z:160}},showscale:true}};
  const rimTrace={{type:'scatter3d',mode:'lines',x:terrain.rimX,y:terrain.rimY,z:zRim.corrected,line:{{color:foreground,width:7}},hovertemplate:'Mapped irregular rim<br>East offset %{{x:.1f}} km<br>North offset %{{y:.1f}} km<extra></extra>',name:'Mapped rim',showlegend:true}};
  const flowTrace={{type:'scatter3d',mode:'lines',x:terrain.flowX,y:terrain.flowY,z:zFlow.corrected,line:{{color:flowColor,width:8}},hovertemplate:'Suspected outflow axis<br>East offset %{{x:.1f}} km<br>North offset %{{y:.1f}} km<extra></extra>',name:'Outflow trace',showlegend:true}};
  const outletTrace={{type:'scatter3d',mode:'markers+text',x:[terrain.outlet.x],y:[terrain.outlet.y],z:[terrain.outlet.corrected/1000],marker:{{color:outletColor,size:7,symbol:'diamond',line:{{color:foreground,width:2}}}},text:['Projected outlet'],textposition:'top center',textfont:{{color:foreground}},hovertemplate:`Projected outlet<br>${{terrain.outlet.lat.toFixed(6)}} deg N, ${{terrain.outlet.lon.toFixed(6)}} deg E<br>Median detrended residual ${{terrain.outlet.corrected.toFixed(1)}} m<extra></extra>`,name:'Projected outlet',showlegend:true}};
  const scene={{xaxis:{{title:'East offset (km)',color:muted,gridcolor:border,zerolinecolor:border,backgroundcolor:background}},yaxis:{{title:'North offset (km)',color:muted,gridcolor:border,zerolinecolor:border,backgroundcolor:background}},zaxis:{{title:'Detrended residual (km)',color:muted,gridcolor:border,zerolinecolor:border,backgroundcolor:background}},aspectmode:'manual',aspectratio:{{x:1.2,y:1,z:.25}},camera}};
  const layout={{margin:{{l:0,r:0,t:8,b:0}},paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',font:{{color:foreground}},scene,legend:{{orientation:'h',x:.5,xanchor:'center',y:.01,yanchor:'bottom',font:{{color:foreground}},bgcolor:'rgba(0,0,0,0)'}},uirevision:'paper-crater-3d'}};
  const config={{responsive:true,displaylogo:false,modeBarButtonsToRemove:['sendDataToCloud','lasso2d','select2d']}};
  Plotly.newPlot(plot,[surface,rimTrace,flowTrace,outletTrace],layout,config).then(()=>applyExaggeration());
  function applyExaggeration(){{const v=Number(exaggeration.value);exaggerationValue.textContent=`${{v}}x`;const spanX=terrain.x.at(-1)-terrain.x[0];const r=ranges[mode];const z=Math.max(.04,Math.min(1.2,(r[1]-r[0])/1000/spanX*v));Plotly.relayout(plot,{{'scene.aspectratio.z':z}});}}
  function setMode(next){{mode=next;const c=mode==='corrected';correctedButton.classList.toggle('btn-primary',c);rawButton.classList.toggle('btn-primary',!c);correctedButton.setAttribute('aria-pressed',String(c));rawButton.setAttribute('aria-pressed',String(!c));const label=c?'Detrended residual (km)':'Raw elevation (km)';const hover=c?'East offset %{{x:.1f}} km<br>North offset %{{y:.1f}} km<br>Detrended residual %{{z:.3f}} km<extra></extra>':'East offset %{{x:.1f}} km<br>North offset %{{y:.1f}} km<br>Raw elevation %{{z:.3f}} km<extra></extra>';surface.z=zSurface[mode];surface.cmin=ranges[mode][0]/1000;surface.cmax=ranges[mode][1]/1000;surface.hovertemplate=hover;surface.colorbar={{title:{{text:label,font:{{color:foreground}}}},thickness:14,len:.66,tickfont:{{color:muted}}}};rimTrace.z=zRim[mode];flowTrace.z=zFlow[mode];outletTrace.z=[terrain.outlet[mode]/1000];scene.zaxis.title=label;detail.textContent=c?`Projected outlet residual: ${{terrain.stats.residualMin}} to ${{terrain.stats.residualMax}} m; lowest ${{terrain.stats.percentileMin}}% to ${{terrain.stats.percentileMax}}% of the mapped rim.`:'Raw MOLA elevation is shown; the regional background slope has not been removed.';Plotly.react(plot,[surface,rimTrace,flowTrace,outletTrace],layout,config).then(()=>applyExaggeration());}}
  function updateOrbitButton(){{orbitButton.setAttribute('aria-pressed',String(!paused));orbitButton.innerHTML=paused?'<i data-lucide="play" aria-hidden="true"></i> Resume orbit':'<i data-lucide="pause" aria-hidden="true"></i> Pause orbit';if(window.lucide)lucide.createIcons({{attrs:{{width:16,height:16}}}});}}
  function animate(now){{const dt=Math.min(50,now-last);last=now;if(!paused&&now>interactionUntil&&now-lastUpdate>38){{angle+=dt*(.0001*Number(speed.value));camera.eye.x=radius*Math.cos(angle);camera.eye.y=radius*Math.sin(angle);camera.eye.z=height;Plotly.relayout(plot,{{'scene.camera':camera}});lastUpdate=now;}}requestAnimationFrame(animate);}}
  plot.on('plotly_relayout',e=>{{if(e['scene.camera']){{camera=e['scene.camera'];angle=Math.atan2(camera.eye.y,camera.eye.x);radius=Math.hypot(camera.eye.x,camera.eye.y);height=camera.eye.z;interactionUntil=performance.now()+2200;}}}});
  correctedButton.addEventListener('click',()=>setMode('corrected'));rawButton.addEventListener('click',()=>setMode('raw'));exaggeration.addEventListener('input',applyExaggeration);speed.addEventListener('input',()=>{{speedValue.textContent=`${{speed.value}}x`;}});orbitButton.addEventListener('click',()=>{{paused=!paused;updateOrbitButton();}});resetButton.addEventListener('click',()=>{{camera=JSON.parse(JSON.stringify(defaultCamera));angle=Math.atan2(camera.eye.y,camera.eye.x);radius=Math.hypot(camera.eye.x,camera.eye.y);height=camera.eye.z;Plotly.relayout(plot,{{'scene.camera':camera}});}});
  if(reducedMotion)detail.textContent='Auto-orbit is paused because the system prefers reduced motion; it can be resumed manually.';updateOrbitButton();requestAnimationFrame(animate);
}})();
</script>
'''

VIS.mkdir(parents=True, exist_ok=True)
OUT.write_text(fragment, encoding="utf-8")
print(OUT, OUT.stat().st_size)
