from pathlib import Path
import csv, json, math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

MARS_R_KM = 3396.19

def read_csv(path):
    d=np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    d.dtype.names=tuple(n.lstrip('\ufeff') for n in d.dtype.names)
    return d

def xy_km(lon, lat, lon0, lat0):
    return ((np.asarray(lon)-lon0)*math.pi/180*MARS_R_KM*math.cos(math.radians(lat0)),
            (np.asarray(lat)-lat0)*math.pi/180*MARS_R_KM)

def lonlat_from_xy(x, y, lon0, lat0):
    return (lon0 + np.asarray(x)/(MARS_R_KM*math.cos(math.radians(lat0)))*180/math.pi,
            lat0 + np.asarray(y)/MARS_R_KM*180/math.pi)

def point_segment_min_distance(px, py, sx, sy):
    p = np.column_stack([px, py])
    best = np.full(len(p), np.inf)
    for i in range(len(sx)-1):
        a = np.array([sx[i], sy[i]])
        b = np.array([sx[i+1], sy[i+1]])
        ab = b-a
        den = float(ab@ab)
        t = np.zeros(len(p)) if den == 0 else np.clip(((p-a)@ab)/den, 0, 1)
        d = np.sqrt(np.sum((p-(a+t[:,None]*ab))**2, axis=1))
        best = np.minimum(best, d)
    return best

def points_in_polygon(px, py, vx, vy):
    px=np.asarray(px); py=np.asarray(py); inside=np.zeros(px.shape,dtype=bool)
    j=len(vx)-1
    for i in range(len(vx)):
        cross=((vy[i]>py)!=(vy[j]>py)) & (px < (vx[j]-vx[i])*(py-vy[i])/(vy[j]-vy[i]+1e-15)+vx[i])
        inside ^= cross
        j=i
    return inside

def robust_plane(x, y, z, max_iter=30):
    A = np.column_stack([np.ones(len(x)), x, y])
    coef = np.linalg.lstsq(A, z, rcond=None)[0]
    for _ in range(max_iter):
        r = z-A@coef
        scale = 1.4826*np.median(np.abs(r-np.median(r))) + 1e-9
        u = np.abs(r)/(1.345*scale)
        w = np.ones_like(u)
        w[u>1] = 1/u[u>1]
        new = np.linalg.lstsq(A*w[:,None], z*w, rcond=None)[0]
        if np.max(np.abs(new-coef)) < 1e-7:
            coef = new; break
        coef = new
    rmse = float(np.sqrt(np.mean((z-A@coef)**2)))
    return coef, rmse

def densify_closed(lon, lat, spacing_km, lon0, lat0):
    x,y = xy_km(lon,lat,lon0,lat0)
    if np.hypot(x[-1]-x[0],y[-1]-y[0]) > 1e-6:
        x=np.r_[x,x[0]]; y=np.r_[y,y[0]]
    seg=np.hypot(np.diff(x),np.diff(y)); cum=np.r_[0,np.cumsum(seg)]
    s=np.arange(0,cum[-1],spacing_km)
    xo=np.interp(s,cum,x); yo=np.interp(s,cum,y)
    lo,la=lonlat_from_xy(xo,yo,lon0,lat0)
    return s,lo,la,xo,yo,float(cum[-1])

def make_grid(d):
    lons=np.unique(d['longitude_e']); lats=np.unique(d['latitude_n'])
    lons.sort(); lats.sort()
    grid=np.full((len(lats),len(lons)),np.nan)
    ix=np.searchsorted(lons,d['longitude_e']); iy=np.searchsorted(lats,d['latitude_n'])
    grid[iy,ix]=d['raw_elevation_m']
    return lons,lats,grid

def bilinear(lons,lats,grid,lon,lat):
    lon=np.asarray(lon); lat=np.asarray(lat)
    ix=np.clip(np.searchsorted(lons,lon)-1,0,len(lons)-2)
    iy=np.clip(np.searchsorted(lats,lat)-1,0,len(lats)-2)
    x0=lons[ix]; x1=lons[ix+1]; y0=lats[iy]; y1=lats[iy+1]
    tx=(lon-x0)/(x1-x0); ty=(lat-y0)/(y1-y0)
    z00=grid[iy,ix]; z10=grid[iy,ix+1]; z01=grid[iy+1,ix]; z11=grid[iy+1,ix+1]
    return (1-tx)*(1-ty)*z00+tx*(1-ty)*z10+(1-tx)*ty*z01+tx*ty*z11

def plane_eval(coef,x,y): return coef[0]+coef[1]*x+coef[2]*y

mola=read_csv(ROOT/'data'/'derived'/'MOLA_subset_points.csv')
rim=read_csv(ROOT/'data'/'geometry'/'JMARS_Profile1_crater_rim_vertices.csv')
flow=read_csv(ROOT/'data'/'geometry'/'JMARS_Profile2_outflow_axis_vertices.csv')
prof=read_csv(ROOT/'data'/'derived'/'JMARS_Profile2_outflow_axis_elevation_profile.csv')
inter=json.loads((ROOT/'data'/'derived'/'geometry_intersection.json').read_text(encoding='utf-8'))
cross=inter['rim_intersections'][0]
clon,clat=cross['longitude_e_deg'],cross['latitude_n_deg']
lon0=float(np.mean(rim['longitude_e_deg'])); lat0=float(np.mean(rim['latitude_n_planetocentric_deg']))

mx,my=xy_km(mola['longitude_e'],mola['latitude_n'],lon0,lat0)
rx,ry=xy_km(rim['longitude_e_deg'],rim['latitude_n_planetocentric_deg'],lon0,lat0)
rx=np.r_[rx,rx[0]]; ry=np.r_[ry,ry[0]]
inside=points_in_polygon(mx,my,rx,ry)
dist_rim=point_segment_min_distance(mx,my,rx,ry)
fx,fy=xy_km(flow['longitude_e_deg'],flow['latitude_n_planetocentric_deg'],lon0,lat0)
dist_flow=point_segment_min_distance(mx,my,fx,fy)
cx,cy=xy_km(clon,clat,lon0,lat0)
downstream=(mx>=cx-1.0)

variants=[
    ('A: 3–15 km annulus',3,15,None),
    ('B: 3–20 km; 5 km corridor excluded',3,20,5),
    ('C: 5–20 km; 8 km corridor excluded',5,20,8),
]
models=[]
for name,dmin,dmax,corr in variants:
    mask=(~inside)&(dist_rim>=dmin)&(dist_rim<=dmax)&np.isfinite(mola['raw_elevation_m'])
    if corr is not None:
        mask &= ~((dist_flow<corr)&downstream)
    coef,rmse=robust_plane(mx[mask],my[mask],mola['raw_elevation_m'][mask])
    models.append(dict(name=name,n=int(mask.sum()),coef=coef,rmse=rmse,mask=mask))

lons,lats,rawgrid=make_grid(mola)
s_rim,rlo,rla,rsx,rsy,perim=densify_closed(rim['longitude_e_deg'],rim['latitude_n_planetocentric_deg'],.25,lon0,lat0)
rim_raw=bilinear(lons,lats,rawgrid,rlo,rla)
cross_s=float(s_rim[np.argmin((rlo-clon)**2+(rla-clat)**2)])
aligned=(s_rim-cross_s)%perim
order=np.argsort(aligned)
aligned=aligned[order]; rlo=rlo[order]; rla=rla[order]; rsx=rsx[order]; rsy=rsy[order]; rim_raw=rim_raw[order]

outx,outy=xy_km(prof['longitude_e_deg'],prof['latitude_n_deg'],lon0,lat0)
resid_rim=[]; resid_map=[]; resid_prof=[]; stats=[]
for m in models:
    rr=rim_raw-plane_eval(m['coef'],rsx,rsy)
    rp=prof['mola_raw_m']-plane_eval(m['coef'],outx,outy)
    ro=float(bilinear(lons,lats,rawgrid,clon,clat)-plane_eval(m['coef'],cx,cy))
    pct=float(100*np.mean(rr<=ro))
    allres=mola['raw_elevation_m']-plane_eval(m['coef'],mx,my)
    resid_rim.append(rr); resid_prof.append(rp); resid_map.append(allres)
    stats.append(dict(model=m['name'],reference_points=m['n'],robust_fit_rmse_m=m['rmse'],
                      intercept_m=float(m['coef'][0]),east_gradient_m_per_km=float(m['coef'][1]),
                      north_gradient_m_per_km=float(m['coef'][2]),outlet_residual_m=ro,
                      rim_median_residual_m=float(np.nanmedian(rr)),rim_min_residual_m=float(np.nanmin(rr)),
                      rim_max_residual_m=float(np.nanmax(rr)),outlet_low_tail_percentile=pct))
resid_rim=np.vstack(resid_rim); resid_prof=np.vstack(resid_prof); resid_map=np.vstack(resid_map)
rim_med=np.nanmedian(resid_rim,axis=0); rim_lo=np.nanmin(resid_rim,axis=0); rim_hi=np.nanmax(resid_rim,axis=0)
prof_med=np.nanmedian(resid_prof,axis=0); prof_lo=np.nanmin(resid_prof,axis=0); prof_hi=np.nanmax(resid_prof,axis=0)
map_med=np.nanmedian(resid_map,axis=0)

ctx=np.isfinite(prof['ctx_raw_m'])
xd=prof['distance_from_user_rim_km'][ctx]; zd=prof['ctx_raw_m'][ctx]
coef_ctx=np.polyfit(xd,zd,1); pred=np.polyval(coef_ctx,xd)
se=float(np.sqrt(np.sum((zd-pred)**2)/(len(xd)-2)/np.sum((xd-xd.mean())**2)))

with (TAB/'detrending_sensitivity.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=stats[0].keys()); w.writeheader(); w.writerows(stats)
with (TAB/'rim_corrected_profile.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f); w.writerow(['distance_from_outlet_along_rim_km','longitude_e_deg','latitude_n_deg','mola_raw_m','corrected_min_m','corrected_median_m','corrected_max_m'])
    w.writerows(zip(aligned,rlo,rla,rim_raw,rim_lo,rim_med,rim_hi))
with (TAB/'outflow_corrected_profile.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f); w.writerow(['distance_from_rim_km','longitude_e_deg','latitude_n_deg','mola_raw_m','corrected_min_m','corrected_median_m','corrected_max_m','ctx_raw_m'])
    w.writerows(zip(prof['distance_from_user_rim_km'],prof['longitude_e_deg'],prof['latitude_n_deg'],prof['mola_raw_m'],prof_lo,prof_med,prof_hi,prof['ctx_raw_m']))

summary={
 'method':'Robust planar regional-trend removal using three external-reference sensitivity masks around the user-traced irregular rim.',
 'rim_perimeter_km':perim,'rim_sampling_km':0.25,'outlet_lon_e_deg':clon,'outlet_lat_n_deg':clat,
 'outflow_length_km':inter['outflow_length_km'],'rim_intersection_count':inter['rim_intersection_count'],
 'outlet_residual_range_m':[float(min(s['outlet_residual_m'] for s in stats)),float(max(s['outlet_residual_m'] for s in stats))],
 'outlet_percentile_range':[float(min(s['outlet_low_tail_percentile'] for s in stats)),float(max(s['outlet_low_tail_percentile'] for s in stats))],
 'ctx_coverage_from_rim_km':[float(xd.min()),float(xd.max())],
 'ctx_observed_slope_m_per_km':float(coef_ctx[0]),'ctx_slope_standard_error_m_per_km':se,
 'interpretation_guardrail':'Geometry and corrected topography may test consistency with a breach, but the present CTX DTM does not cover the rim crossing and therefore cannot resolve a breach crest or prove a dam-breach origin.',
 'models':stats}
(OUT/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

try:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
except ModuleNotFoundError:
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    print('Matplotlib unavailable; numerical outputs completed, figures skipped.')
    raise SystemExit(0)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.linewidth':0.7,'savefig.facecolor':'white'})
extent=[lons.min(),lons.max(),lats.min(),lats.max()]
grid_res=np.full_like(rawgrid,np.nan)
ix=np.searchsorted(lons,mola['longitude_e']); iy=np.searchsorted(lats,mola['latitude_n']); grid_res[iy,ix]=map_med
fig,axs=plt.subplots(1,2,figsize=(7.2,3.55),constrained_layout=True)
for ax,data,cmap,label in [(axs[0],rawgrid,'terrain','Elevation (m)'),(axs[1],grid_res,'RdBu_r','Detrended elevation (m)')]:
    if label.startswith('Detrended'):
        lim=np.nanpercentile(np.abs(data),98); im=ax.imshow(data,origin='lower',extent=extent,cmap=cmap,vmin=-lim,vmax=lim,aspect='equal')
    else: im=ax.imshow(data,origin='lower',extent=extent,cmap=cmap,aspect='equal')
    ax.plot(np.r_[rim['longitude_e_deg'],rim['longitude_e_deg'][0]],np.r_[rim['latitude_n_planetocentric_deg'],rim['latitude_n_planetocentric_deg'][0]],color='white',lw=1.4,label='User-traced rim')
    ax.plot(flow['longitude_e_deg'],flow['latitude_n_planetocentric_deg'],color='#00e5ff',lw=1.6,label='Mapped outflow axis')
    ax.scatter([clon],[clat],s=30,marker='*',c='#ffd54f',edgecolor='black',linewidth=.5,zorder=5,label='Rim crossing')
    ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.2f}')); ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.2f}'))
    cb=fig.colorbar(im,ax=ax,shrink=.82,pad=.02); cb.set_label(label)
axs[0].text(.02,.97,'a',transform=axs[0].transAxes,va='top',fontweight='bold',fontsize=11,color='white')
axs[1].text(.02,.97,'b',transform=axs[1].transAxes,va='top',fontweight='bold',fontsize=11)
axs[1].legend(loc='lower right',fontsize=6.5,frameon=True)
for ext in ['png','pdf','svg']:
    fig.savefig(FIG/f'Figure_1_corrected_topography.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

fig,axs=plt.subplots(2,1,figsize=(7.2,5.3),constrained_layout=True)
ax=axs[0]
ax.fill_between(aligned,rim_lo,rim_hi,color='#90a4ae',alpha=.35,label='Model sensitivity envelope')
ax.plot(aligned,rim_med,color='#263238',lw=1.1,label='Median corrected rim')
ax.axvline(0,color='#d32f2f',lw=1,ls='--'); ax.scatter([0],[np.interp(0,aligned,rim_med)],c='#d32f2f',s=22,zorder=4,label='Mapped crossing')
ax.axhline(0,color='0.6',lw=.6); ax.set_xlim(0,perim); ax.set_ylabel('Detrended rim elevation (m)'); ax.set_xlabel('Distance along rim from mapped crossing (km)'); ax.legend(ncol=3,fontsize=7,loc='upper right'); ax.text(.01,.94,'a',transform=ax.transAxes,fontweight='bold',fontsize=11)
ax=axs[1]
xprof=prof['distance_from_user_rim_km']
ax.fill_between(xprof,prof_lo,prof_hi,color='#90caf9',alpha=.35,label='MOLA correction sensitivity')
ax.plot(xprof,prof_med,color='#1565c0',lw=1.2,label='MOLA, median detrended')
ax.axvline(0,color='#d32f2f',lw=1,ls='--',label='Rim crossing')
ax.axvspan(0,float(xd.min()),color='#ffcc80',alpha=.25,label='No CTX coverage')
ax.set_xlabel('Distance from mapped rim crossing (km)'); ax.set_ylabel('MOLA detrended elevation (m)',color='#1565c0'); ax.tick_params(axis='y',labelcolor='#1565c0')
ax2=ax.twinx(); ax2.plot(xd,zd,color='#4a148c',lw=.8,alpha=.8,label='CTX raw elevation (covered reach)'); ax2.plot(xd,np.polyval(coef_ctx,xd),color='#4a148c',lw=1.4,ls='--',label=f'CTX linear trend: {coef_ctx[0]:.1f} m km⁻¹'); ax2.set_ylabel('CTX elevation (m)',color='#4a148c'); ax2.tick_params(axis='y',labelcolor='#4a148c')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,ncol=2,fontsize=6.6,loc='best'); ax.text(.01,.94,'b',transform=ax.transAxes,fontweight='bold',fontsize=11)
for ext in ['png','pdf','svg']:
    fig.savefig(FIG/f'Figure_2_profiles_and_coverage.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)
print(json.dumps(summary,ensure_ascii=False,indent=2))
