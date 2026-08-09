from pathlib import Path
import json, math
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from reportlab.lib.pagesizes import landscape

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; FIG=OUT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
REG='C:/Windows/Fonts/arial.ttf'; BOLD='C:/Windows/Fonts/arialbd.ttf'

def font(n,b=False): return ImageFont.truetype(BOLD if b else REG,n)
def load(path):
    d=np.genfromtxt(path,delimiter=',',names=True,dtype=float,encoding='utf-8')
    d.dtype.names=tuple(n.lstrip('\ufeff') for n in d.dtype.names); return d
def lerp(a,b,t): return tuple(int(a[i]+(b[i]-a[i])*t) for i in range(3))
def palette(vals,stops,vrange):
    z=np.asarray(vals,float); u=np.clip((z-vrange[0])/(vrange[1]-vrange[0]),0,1); rgb=np.zeros(z.shape+(3,),dtype=np.uint8)
    for k in range(len(stops)-1):
        p0,c0=stops[k]; p1,c1=stops[k+1]; m=(u>=p0)&(u<=p1); t=np.zeros_like(u); t[m]=(u[m]-p0)/(p1-p0)
        for j in range(3): rgb[...,j][m]=(c0[j]+(c1[j]-c0[j])*t[m]).astype(np.uint8)
    rgb[~np.isfinite(z)]=245; return rgb
TERR=[(0,(32,65,83)),(.2,(54,116,116)),(.42,(117,140,91)),(.62,(181,161,111)),(.82,(173,123,86)),(1,(245,242,231))]
DIV=[(0,(33,102,172)),(.25,(103,169,207)),(.5,(247,247,247)),(.75,(239,138,98)),(1,(178,24,43))]
INK=(33,43,54); MUTED=(91,103,112); GRID=(215,220,224); CYAN=(0,188,212); GOLD=(255,193,7); RED=(198,40,40); PURPLE=(74,20,140); BLUE=(21,101,192)

def save_pdf(img,path,w_in,h_in):
    c=canvas.Canvas(str(path),pagesize=(w_in*72,h_in*72)); c.drawImage(ImageReader(img),0,0,width=w_in*72,height=h_in*72); c.showPage(); c.save()
def text_center(d,xy,s,f,fill=INK):
    b=d.textbbox((0,0),s,font=f); d.text((xy[0]-(b[2]-b[0])/2,xy[1]-(b[3]-b[1])/2),s,font=f,fill=fill)
def rotated_text(img,xy,s,f,fill=INK):
    b=f.getbbox(s); layer=Image.new('RGBA',(b[2]-b[0]+20,b[3]-b[1]+20),(255,255,255,0)); ImageDraw.Draw(layer).text((10,10),s,font=f,fill=fill); layer=layer.rotate(90,expand=True); img.alpha_composite(layer,(int(xy[0]-layer.width/2),int(xy[1]-layer.height/2)))
def mapxy(lon,lat,box,extent):
    l,t,r,b=box; x=l+(np.asarray(lon)-extent[0])/(extent[1]-extent[0])*(r-l); y=b-(np.asarray(lat)-extent[2])/(extent[3]-extent[2])*(b-t); return x,y
def draw_axes(d,box,extent,fs=34):
    l,t,r,b=box; d.rectangle(box,outline=INK,width=3)
    for lon in np.linspace(extent[0],extent[1],5):
        x,_=mapxy(lon,extent[2],box,extent); d.line((x,b,x,b+12),fill=INK,width=2); text_center(d,(x,b+34),f'{lon:.2f}',font(fs),INK)
    for lat in np.linspace(extent[2],extent[3],5):
        _,y=mapxy(extent[0],lat,box,extent); d.line((l-12,y,l,y),fill=INK,width=2); s=f'{lat:.2f}'; bb=d.textbbox((0,0),s,font=font(fs)); d.text((l-20-(bb[2]-bb[0]),y-(bb[3]-bb[1])/2),s,font=font(fs),fill=INK)

mola=load(ROOT/'data'/'derived'/'MOLA_subset_points.csv'); rim=load(ROOT/'data'/'geometry'/'JMARS_Profile1_crater_rim_vertices.csv'); flow=load(ROOT/'data'/'geometry'/'JMARS_Profile2_outflow_axis_vertices.csv')
summary=json.loads((OUT/'analysis_summary.json').read_text(encoding='utf-8'))
models=summary['models']; lon0=float(np.mean(rim['longitude_e_deg'])); lat0=float(np.mean(rim['latitude_n_planetocentric_deg'])); R=3396.19
mx=(mola['longitude_e']-lon0)*math.pi/180*R*math.cos(math.radians(lat0)); my=(mola['latitude_n']-lat0)*math.pi/180*R
res=[]
for m in models: res.append(mola['raw_elevation_m']-(m['intercept_m']+m['east_gradient_m_per_km']*mx+m['north_gradient_m_per_km']*my))
res=np.nanmedian(np.vstack(res),axis=0)
lons=np.unique(mola['longitude_e']); lats=np.unique(mola['latitude_n']); lons.sort(); lats.sort(); raw=np.full((len(lats),len(lons)),np.nan); rg=np.full_like(raw,np.nan)
ix=np.searchsorted(lons,mola['longitude_e']); iy=np.searchsorted(lats,mola['latitude_n']); raw[iy,ix]=mola['raw_elevation_m']; rg[iy,ix]=res
extent=(float(lons.min()),float(lons.max()),float(lats.min()),float(lats.max())); cross=(summary['outlet_lon_e_deg'],summary['outlet_lat_n_deg'])

W,H=4320,2160; img=Image.new('RGBA',(W,H),'white'); d=ImageDraw.Draw(img); boxes=[(350,180,2040,1830),(2450,180,4140,1830)]
for n,(box,data,stops,label) in enumerate(zip(boxes,[raw,rg],[TERR,DIV],['Elevation (m)','Detrended elevation (m)'])):
    if n==0: vr=(float(np.nanpercentile(data,2)),float(np.nanpercentile(data,98)))
    else: lim=float(np.nanpercentile(np.abs(data),98)); vr=(-lim,lim)
    rgb=palette(data,stops,vr); heat=Image.fromarray(rgb[::-1,:,:],'RGB').resize((box[2]-box[0],box[3]-box[1]),Image.Resampling.BILINEAR); img.paste(heat,(box[0],box[1])); draw_axes(d,box,extent)
    rlon=np.r_[rim['longitude_e_deg'],rim['longitude_e_deg'][0]]; rlat=np.r_[rim['latitude_n_planetocentric_deg'],rim['latitude_n_planetocentric_deg'][0]]; x,y=mapxy(rlon,rlat,box,extent); d.line(list(zip(x,y)),fill='white',width=9,joint='curve')
    x,y=mapxy(flow['longitude_e_deg'],flow['latitude_n_planetocentric_deg'],box,extent); d.line(list(zip(x,y)),fill=CYAN,width=10,joint='curve')
    x,y=mapxy(*cross,box,extent); rad=18; d.ellipse((x-rad,y-rad,x+rad,y+rad),fill=GOLD,outline=INK,width=4)
    d.text((box[0]+20,box[1]+10),'ab'[n],font=font(62,True),fill='white' if n==0 else INK)
    # north arrow and 10 km scale
    nx,ny=box[2]-80,box[1]+120; d.polygon([(nx,ny-55),(nx-20,ny+15),(nx,ny),(nx+20,ny+15)],fill=INK); text_center(d,(nx,ny-90),'N',font(34,True))
    km_per_lon=R*math.pi/180*math.cos(math.radians(lat0)); px10=10/(extent[1]-extent[0])/km_per_lon*(box[2]-box[0]); sx=box[0]+60; sy=box[3]-70; d.line((sx,sy,sx+px10,sy),fill='white',width=10); d.line((sx,sy-12,sx,sy+12),fill='white',width=5); d.line((sx+px10,sy-12,sx+px10,sy+12),fill='white',width=5); text_center(d,(sx+px10/2,sy-32),'10 km',font(30,True),fill='white')
    # colorbar
    cb=(box[0],box[3]+95,box[2],box[3]+135); grad=np.linspace(vr[0],vr[1],box[2]-box[0])[None,:]; cbi=Image.fromarray(palette(grad,stops,vr),'RGB').resize((box[2]-box[0],40)); img.paste(cbi,(cb[0],cb[1])); d.rectangle(cb,outline=INK,width=2); d.text((cb[0],cb[3]+10),f'{vr[0]:.0f}',font=font(27),fill=INK); s=f'{vr[1]:.0f}'; bb=d.textbbox((0,0),s,font=font(27)); d.text((cb[2]-(bb[2]-bb[0]),cb[3]+10),s,font=font(27),fill=INK); text_center(d,((cb[0]+cb[2])/2,cb[3]+28),label,font(29))
    text_center(d,((box[0]+box[2])/2,box[3]+78),'Longitude (°E)',font(36)); rotated_text(img,(box[0]-170,(box[1]+box[3])/2),'Latitude (°N)',font(36))
# legend
ly=2050; x=2450
d.line((x,ly,x+90,ly),fill='white',width=10); d.text((x+110,ly-20),'User-traced rim',font=font(31),fill=INK); x+=520
d.line((x,ly,x+90,ly),fill=CYAN,width=10); d.text((x+110,ly-20),'Mapped outflow axis',font=font(31),fill=INK); x+=620
d.ellipse((x,ly-14,x+28,ly+14),fill=GOLD,outline=INK,width=3); d.text((x+48,ly-20),'Mapped rim crossing',font=font(31),fill=INK)
png=FIG/'Figure_1_corrected_topography.png'; img.convert('RGB').save(png,dpi=(600,600),quality=95); save_pdf(img.convert('RGB'),FIG/'Figure_1_corrected_topography.pdf',7.2,3.6)

rp=load(OUT/'tables'/'rim_corrected_profile.csv'); op=load(OUT/'tables'/'outflow_corrected_profile.csv')
W,H=4320,3180; img=Image.new('RGBA',(W,H),'white'); d=ImageDraw.Draw(img)
def plot_panel(box,x,ylo,ymed,yhi,xlabel,ylabel,panel):
    l,t,r,b=box; xmin,xmax=float(np.nanmin(x)),float(np.nanmax(x)); ymin=float(np.nanmin(ylo)); ymax=float(np.nanmax(yhi)); pad=.06*(ymax-ymin); ymin-=pad; ymax+=pad
    X=lambda q:l+(np.asarray(q)-xmin)/(xmax-xmin)*(r-l); Y=lambda q:b-(np.asarray(q)-ymin)/(ymax-ymin)*(b-t)
    for yy in np.linspace(ymin,ymax,5):
        py=Y(yy); d.line((l,py,r,py),fill=GRID,width=2); s=f'{yy:.0f}'; bb=d.textbbox((0,0),s,font=font(30)); d.text((l-18-(bb[2]-bb[0]),py-15),s,font=font(30),fill=INK)
    for xx in np.linspace(xmin,xmax,6):
        px=X(xx); d.line((px,t,px,b),fill=GRID,width=2); text_center(d,(px,b+38),f'{xx:.0f}',font(30))
    pts=list(zip(X(x),Y(ylo)))+list(zip(X(x[::-1]),Y(yhi[::-1]))); d.polygon(pts,fill=(144,164,174,85)); d.line(list(zip(X(x),Y(ymed))),fill=INK,width=6,joint='curve'); d.rectangle(box,outline=INK,width=3)
    d.text((l+18,t+10),panel,font=font(58,True),fill=INK); text_center(d,((l+r)/2,b+100),xlabel,font(36)); rotated_text(img,(l-180,(t+b)/2),ylabel,font(36)); return X,Y,(xmin,xmax,ymin,ymax)
box1=(400,160,4050,1370); X1,Y1,_=plot_panel(box1,rp['distance_from_outlet_along_rim_km'],rp['corrected_min_m'],rp['corrected_median_m'],rp['corrected_max_m'],'Distance along rim from mapped crossing (km)','Detrended rim elevation (m)','a')
d.line((X1(0),box1[1],X1(0),box1[3]),fill=RED,width=6); d.ellipse((X1(0)-12,Y1(rp['corrected_median_m'][0])-12,X1(0)+12,Y1(rp['corrected_median_m'][0])+12),fill=RED)
d.rectangle((2450,220,2500,260),fill=(144,164,174)); d.text((2520,215),'Model sensitivity envelope',font=font(30),fill=INK); d.line((3150,240,3250,240),fill=INK,width=6); d.text((3270,215),'Median corrected rim',font=font(30),fill=INK)

box2=(400,1750,4050,2960); x=op['distance_from_rim_km']; ylo=op['corrected_min_m']; ymed=op['corrected_median_m']; yhi=op['corrected_max_m']; X2,Y2,lims=plot_panel(box2,x,ylo,ymed,yhi,'Distance from mapped rim crossing (km)','MOLA detrended elevation (m)','b')
d.line((X2(0),box2[1],X2(0),box2[3]),fill=RED,width=6)
ctx=np.isfinite(op['ctx_raw_m']); xc=x[ctx]; zc=op['ctx_raw_m'][ctx]; zmin,zmax=float(zc.min()),float(zc.max()); zpad=.08*(zmax-zmin); zmin-=zpad; zmax+=zpad; Yc=lambda q:box2[3]-(np.asarray(q)-zmin)/(zmax-zmin)*(box2[3]-box2[1])
xstart=float(xc.min()); d.rectangle((X2(0),box2[1],X2(xstart),box2[3]),fill=(255,247,232));
# redraw the MOLA sensitivity band and median over the coverage-warning background
pts=list(zip(X2(x),Y2(ylo)))+list(zip(X2(x[::-1]),Y2(yhi[::-1]))); d.polygon(pts,fill=(144,164,174)); d.line(list(zip(X2(x),Y2(ymed))),fill=BLUE,width=7,joint='curve')
d.line((X2(0),box2[1],X2(0),box2[3]),fill=RED,width=6); d.line((X2(xstart),box2[1],X2(xstart),box2[3]),fill=(239,150,40),width=4)
d.line(list(zip(X2(xc),Yc(zc))),fill=PURPLE,width=5,joint='curve'); slope=summary['ctx_observed_slope_m_per_km']; intercept=float(np.nanmean(zc)-slope*np.nanmean(xc)); d.line((X2(xc.min()),Yc(slope*xc.min()+intercept),X2(xc.max()),Yc(slope*xc.max()+intercept)),fill=PURPLE,width=8)
for yy in np.linspace(zmin,zmax,5):
    py=Yc(yy); d.line((box2[2],py,box2[2]+12,py),fill=PURPLE,width=2); d.text((box2[2]+18,py-15),f'{yy:.0f}',font=font(30),fill=PURPLE)
rotated_text(img,(4230,(box2[1]+box2[3])/2),'CTX elevation (m)',font(36),PURPLE)
d.rectangle((470,1810,520,1850),fill=(255,204,128)); d.text((540,1805),f'No CTX coverage to {xstart:.2f} km',font=font(29),fill=INK); d.line((1260,1830,1360,1830),fill=BLUE,width=6); d.text((1380,1805),'MOLA median detrended',font=font(29),fill=INK); d.line((2110,1830,2210,1830),fill=PURPLE,width=7); d.text((2230,1805),f'CTX trend: {slope:.1f} ± {summary["ctx_slope_standard_error_m_per_km"]:.1f} m km⁻¹',font=font(29),fill=INK)
png=FIG/'Figure_2_profiles_and_coverage.png'; img.convert('RGB').save(png,dpi=(600,600),quality=95); save_pdf(img.convert('RGB'),FIG/'Figure_2_profiles_and_coverage.pdf',7.2,5.3)
print(png)
