"""VOIDS randoms + Veff — addresses Q1 M6/editor repro.
Uniform comoving randoms in z 0.01-0.114 x SDSS mask (NSA_main_mask.pickle),
same ball-query pipeline as crossmatch_v2. Outputs data/randoms_summary.txt + Veff table.
Run: venvs/b3/bin/python src/randoms_veff.py
"""
import pickle
import numpy as np
from astropy.cosmology import Planck18
import astropy.units as u
from scipy.spatial import cKDTree

H = 0.7
rng = np.random.default_rng(42)
with open("data/vast/NSA_main_mask.pickle","rb") as f:
    mask,_,_ = pickle.load(f)
def in_mask(ra, dec):
    return bool(mask[min(max(int(np.floor(ra%360)),0),359), min(max(int(np.floor(dec+90)),0),179)])

voids=[]
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt",encoding="utf-8-sig") as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"): continue
        p=line.split()
        voids.append((float(p[0])/H,float(p[1])/H,float(p[2])/H,float(p[3])/H))
vxyz=np.array([[v[0],v[1],v[2]] for v in voids]); vrad=np.array([v[3] for v in voids])
tree=cKDTree(vxyz); maxRv=vrad.max()

# Uniform comoving: sample z with p(z)~dV/dz, RA uniform, Dec uniform on sphere, keep in-mask
N=200000
zmin,zmax=0.01,0.114
# Inverse-CDF via Planck18 comoving volume
from astropy.cosmology import z_at_value
Vmin=Planck18.comoving_volume(zmin).value; Vmax=Planck18.comoving_volume(zmax).value
V=rng.uniform(Vmin,Vmax,N)
# approximate z from V via interpolation grid
zg=np.linspace(zmin,zmax,500); Vg=Planck18.comoving_volume(zg).value
z=np.interp(V,Vg,zg)
ra=rng.uniform(0,360,N); dec=np.rad2deg(np.arcsin(rng.uniform(-1,1,N)))
keep=np.array([in_mask(r,d) for r,d in zip(ra,dec)])
ra,dec,z=ra[keep],dec[keep],z[keep]
print(f"randoms in-mask {len(ra)}/{N}")
from astropy.coordinates import SkyCoord
c=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,distance=Planck18.comoving_distance(z))
sxyz=np.vstack([c.cartesian.x.value,c.cartesian.y.value,c.cartesian.z.value]).T
rr=[]
for i in range(len(sxyz)):
    cand=tree.query_ball_point(sxyz[i],r=maxRv*1.2)
    if not cand: rr.append(np.inf); continue
    rr.append(min(np.linalg.norm(sxyz[i]-vxyz[j])/vrad[j] for j in cand))
rr=np.array(rr)
f_void=(rr<0.8).mean(); f_shell=((rr>=0.8)&(rr<=1.0)).mean(); f_wall=(rr>1.0).mean()
print(f"f_rand void={f_void:.3f} shell={f_shell:.3f} wall={f_wall:.3f}")
# Veff: survey comoving volume in-mask (approx full-sky V * mask sky fraction * radial fraction)
# Sky fraction from mask grid
sky_frac=mask.mean()
Vshell_full=(Vmax-Vmin)*(4*np.pi/3)/(4*np.pi/3)  # Vmax-Vmin already full-sky comoving volume in Mpc^3
Veff_tot=(Vmax-Vmin)*sky_frac
print(f"sky_frac={sky_frac:.3f} Veff_tot={Veff_tot:.3e} Mpc^3")
print(f"Veff_void={Veff_tot*f_void:.3e} Veff_wall={Veff_tot*f_wall:.3e}")
# Observed BTS in-mask cosmo counts for rate ratios (from v2 log)
N_void_obs, N_wall_obs = 176, 1883
R_ratio=(N_void_obs/(Veff_tot*f_void))/(N_wall_obs/(Veff_tot*f_wall))
print(f"raw rate ratio void/wall={(N_void_obs/N_wall_obs):.4f} volume-corrected={R_ratio:.3f}")
print(f"Poisson N_void={N_void_obs} => rel err {1/np.sqrt(N_void_obs):.2f}")
with open("data/randoms_summary.txt","w") as o:
    o.write(f"N_rand_inmask={len(ra)}\nsky_frac={sky_frac:.4f}\n")
    o.write(f"f_void={f_void:.4f}\nf_shell={f_shell:.4f}\nf_wall={f_wall:.4f}\n")
    o.write(f"Veff_tot={Veff_tot:.6e}\nR_ratio_volcorr={R_ratio:.4f}\nseed=42\n")
print("wrote data/randoms_summary.txt")
