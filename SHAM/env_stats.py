import numpy as np
import matplotlib.pyplot as plt
from utils import *
import plotparams
#plotparams.default()
import sys

env = sys.argv[1]#'_lo'#'_hi'#'low'#'high'
which = sys.argv[2]#"infall"#"peak"#"infall"
sat_str = "_sats"#"_sats"

dm_inds_env = np.load("dm_inds"+env+sat_str+".npy")
fp_inds_env = np.load("fp_inds"+env+sat_str+".npy")

def get_frac(v_dm,v_fp):
    frac = (v_fp-v_dm)/v_dm
    #frac[np.abs(v_dm-0.) < 1.e-6] = -2
    return frac


peakstats_fp = np.load("peakstats"+env+sat_str+"_fp.npy")
peakstats_dm = np.load("peakstats"+env+sat_str+"_dm.npy")
print(peakstats_dm.shape)

peakstats_fp = peakstats_fp[fp_inds_env]
peakstats_dm = peakstats_dm[dm_inds_env]
    
zpeak_dm = peakstats_dm[:,0]
zinfall_dm = peakstats_dm[:,1]
mpeak_dm = peakstats_dm[:,2]*1.e10
minfall_dm = peakstats_dm[:,3]*1.e10
dmpeak_dm = peakstats_dm[:,4]*1.e10
dminfall_dm = peakstats_dm[:,5]*1.e10
print(np.sum(mpeak_dm > 0.))
print(len(mpeak_dm))

zpeak_fp = peakstats_fp[:,0]
zinfall_fp = peakstats_fp[:,1]
mpeak_fp = peakstats_fp[:,2]*1.e10
minfall_fp = peakstats_fp[:,3]*1.e10
dmpeak_fp = peakstats_fp[:,4]*1.e10
dminfall_fp = peakstats_fp[:,5]*1.e10
print(np.sum(mpeak_fp > 0.));

if which == "peak":
    frac_z = zpeak_dm-zpeak_fp#get_frac(zpeak_dm,zpeak_fp)
    frac_m = get_frac(mpeak_dm,mpeak_fp)
    frac_dm = get_frac(dmpeak_dm,dmpeak_fp)
    proxy_dm = mpeak_dm
elif which == "infall":
    frac_z = zinfall_dm-zinfall_fp#get_frac(zinfall_dm,zinfall_fp)
    frac_m = get_frac(minfall_dm,minfall_fp)
    frac_dm = get_frac(dminfall_dm,dminfall_fp)
    proxy_dm = minfall_dm
    
plt.subplots(1,3,figsize=(3*5.5,4.5))
s = 0.1
xlim = [1.e11,2.e15]

plt.subplot(1,3,1)
plt.plot(proxy_dm,np.zeros(len(proxy_dm)),"k--",alpha=0.5)
plt.scatter(proxy_dm,frac_z,s=s,label="z"+which)
plot_median(proxy_dm,frac_z,np.log10(xlim[0]),np.log10(xlim[1]),n_bins=41)
plt.xlabel("M "+which+" (DM)")
plt.ylabel("(FP-DM)")
plt.xscale('log')
ylim = [-1.,0.5]
plt.ylim(ylim)
plt.xlim(xlim)
plt.legend()

plt.subplot(1,3,2)
plt.plot(proxy_dm,np.zeros(len(proxy_dm)),"k--",alpha=0.5)
plt.scatter(proxy_dm,frac_m,s=s,label="m"+which)
plot_median(proxy_dm,frac_m,np.log10(xlim[0]),np.log10(xlim[1]),n_bins=41)
plt.xlabel("M "+which+" (DM)")
plt.ylabel("(FP-DM)/DM")
plt.xscale('log')
ylim = [-0.5,0.5]
plt.ylim(ylim)
plt.xlim(xlim)
plt.legend()

plt.subplot(1,3,3)
plt.plot(proxy_dm,np.zeros(len(proxy_dm)),"k--",alpha=0.5)
plt.scatter(proxy_dm,frac_dm,s=s,label="dm"+which)
plot_median(proxy_dm,frac_dm,np.log10(xlim[0]),np.log10(xlim[1]),n_bins=41)
plt.xlabel("M "+which+" (DM)")
plt.ylabel("(FP-DM)/DM")
plt.xscale('log')
#ylim = [-2,5]
ylim = [-1,1]
plt.ylim(ylim)
plt.xlim(xlim)
plt.legend()
plt.savefig(which+"stats"+env+".png")
#plt.show()
