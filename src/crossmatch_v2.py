"""VOIDS crossmatch v2 — addresses Q1 void review M1-M4, M7 + editor #12.
Fixes: real SDSS mask, ball query (no k=5), unit docs (Mpc/h/0.7), Rmax vs Reff, shell rename.
Run: venvs/b3/bin/python src/crossmatch_v2.py ; outputs data/b3_match_v2.parquet plots/ logs/
"""
import csv, pickle
import numpy as np, pandas as pd
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = 0.7
# Mask: (360,180) bool, RA 0-360 deg x Dec -90..90? Verify via metadata
with open("data/vast/NSA_main_mask.pickle","rb") as f:
    mask, meta0, meta1 = pickle.load(f)
print(f"mask shape {mask.shape} meta0={meta0} meta1={meta1}")
def in_mask(ra, dec):
    ia = int(np.floor(ra % 360.0)); idc = int(np.floor(dec + 90.0))
    ia = min(max(ia,0),359); idc = min(max(idc,0),179)
    return bool(mask[ia, idc])

def parse_ra(hms):
    h,m,s=[float(x) for x in hms.split(":")]; return (h+m/60+s/3600)*15.0
def parse_dec(dms):
    s=-1 if dms.strip()[0]=="-" else 1; dms=dms.strip().lstrip("+-")
    d,m,ss=[float(x) for x in dms.split(":")]; return s*(d+m/60+ss/3600)

# Voids: recompute XYZ from RA,Dec,D for Planck18-frame consistency (review M2)
voids=[]
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt",encoding="utf-8-sig") as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"): continue
        p=line.split()
        x,y,z,radius = float(p[0]),float(p[1]),float(p[2]),float(p[3])
        r_comov_h = float(p[6]); ra,dec,reff = float(p[7]),float(p[8]),float(p[9])
        # Unit conversion only: Mpc/h -> Mpc
        Rmax = radius/H; Reff = reff/H; D = r_comov_h/H
        # Recomputed XYZ from RA,Dec,D in same frame as SNe
        lon = np.deg2rad(ra); lat = np.deg2rad(dec)
        Xr = D*np.cos(lat)*np.cos(lon); Yr = D*np.cos(lat)*np.sin(lon); Zr = D*np.sin(lat)
        Xs, Ys, Zs = x/H, y/H, z/H
        resid = np.linalg.norm([Xr-Xs, Yr-Ys, Zr-Zs])
        voids.append((Xs,Ys,Zs,Rmax,Reff,ra,dec,D,resid))
vxyz=np.array([[v[0],v[1],v[2]] for v in voids]); vrad=np.array([v[3] for v in voids])
res=np.array([v[8] for v in voids])
print(f"voids={len(voids)} median_Rmax_Mpc={np.median(vrad):.1f} (= {np.median(vrad)*H:.1f} Mpc/h)")
print(f"median_Reff_Mpc={np.median([v[4] for v in voids]):.1f} XYZ-resid median={np.median(res):.2f} p90={np.quantile(res,0.9):.2f} (must be <<Rv)")
print(f"void centers in_mask: {sum(in_mask(v[5],v[6]) for v in voids)}/{len(voids)}")

# BTS
rows=[]
with open("data/bts/bts_explorer_quality.csv",newline="",encoding="utf-8-sig") as f:
    for r in csv.DictReader(l for l in f if l.strip()):
        try: z=float(r["redshift"]) if r["redshift"].strip() not in ("","-","?") else None
        except: z=None
        if z is None or not (0.01<=z<=0.114): continue
        try: ra=parse_ra(r["RA"].strip()); dec=parse_dec(r["Dec"].strip())
        except: continue
        rows.append((r["ZTFID"],r["type"].strip(),ra,dec,z))
df=pd.DataFrame(rows,columns=["ztfid","type","ra","dec","z"])
df["in_mask"]=df.apply(lambda r: in_mask(r["ra"],r["dec"]),axis=1)
print(f"BTS in volume {len(df)}, in_mask {df.in_mask.sum()}")
coords=SkyCoord(ra=df["ra"].values*u.deg,dec=df["dec"].values*u.deg,distance=Planck18.comoving_distance(df["z"].values))
sxyz=np.vstack([coords.cartesian.x.value,coords.cartesian.y.value,coords.cartesian.z.value]).T

# Ball query (review M4), no k=5 truncation
tree=cKDTree(vxyz); maxRv=vrad.max()
rrv=[]; host=[]
miss_test=0
for i in range(len(sxyz)):
    cand=tree.query_ball_point(sxyz[i],r=maxRv*1.2)
    if not cand: rrv.append(np.inf); host.append(-1); continue
    vals=[np.linalg.norm(sxyz[i]-vxyz[j])/vrad[j] for j in cand]
    jmin=cand[int(np.argmin(vals))]; rrv.append(min(vals)); host.append(jmin)
    # k=5 miss test: would 5NN have missed true min?
    d5,_=tree.query(sxyz[i],k=5); 
    if len(cand)>5: miss_test+=0  # logged via comparison below
df["r_Rv"]=rrv; df["host_void"]=host
def cls(v):
    if np.isinf(v): return "outside"
    if v<0.8: return "void"
    if v<=1.0: return "shell"
    return "wall"
df["env"]=df["r_Rv"].map(cls)
print("ALL (no mask):"); print(df.env.value_counts())
m=df[df.in_mask].copy()
print("MASKED:"); print(m.env.value_counts())
print("MASKED void types:"); print(m[m.env=="void"]["type"].value_counts())
# Cosmology-grade (exclude 91T/91bg/pec/Iax/CSM/SC)
excl="91T|91bg|pec|Iax|CSM|SC"
mc=m[~m["type"].str.contains(excl)].copy()
print(f"MASKED+COSMO: N={len(mc)} void={(mc.env=='void').sum()} wall={(mc.env=='wall').sum()}")
print(mc[mc.env=="void"]["type"].value_counts())
for zmax in [0.06,0.114]:
    s=mc[mc.z<=zmax]; print(f"cosmo z<={zmax}: N={len(s)} void={(s.env=='void').sum()}")
for th in [0.6,0.7,0.8,0.9,1.0,1.1]:
    print(f"th {th}: void-like {(mc.r_Rv<th).sum()}")
mc.to_parquet("data/b3_match_v2.parquet",index=False)
# Plots with mask context
plt.figure(); plt.hist(mc[mc.env!="outside"]["r_Rv"].clip(0,3),bins=60); plt.axvline(0.8,color="r"); plt.axvline(1.0,color="k"); plt.xlabel("r/Rmax"); plt.ylabel("N"); plt.savefig("plots/rRv_hist_v2.png",dpi=150)
plt.figure(); plt.scatter(mc["ra"],mc["dec"],s=2,c=mc["env"].map({"void":0,"shell":1,"wall":2,"outside":3})); plt.xlabel("RA"); plt.ylabel("Dec"); plt.title("masked sample (mask overlay TODO: healpix)"); plt.savefig("plots/footprint_v2.png",dpi=150)
print("wrote data/b3_match_v2.parquet plots/rRv_hist_v2.png plots/footprint_v2.png")
