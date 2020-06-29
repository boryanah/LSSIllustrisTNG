import matplotlib.pyplot as plt
import numpy as np
import plotparams
plotparams.buba()

proxy = 'vrelax'

bin_cen = np.load("data/scat_bins_"+proxy+".npy")
#hist_fid = np.load("data/scat_"+proxy+"_dm.npy")

wants = ['_matched','']
exts = ['fp','dm']
cs = ['k','orange','dodgerblue']
i = -1
plt.figure(3,figsize=(15,12))
for want_matched in wants:
    for dm_ext in exts:
        if dm_ext == 'fp' and want_matched == '_matched': continue
        i += 1
        hist = np.load("data/scat_med"+want_matched+"_"+proxy+"_"+dm_ext+".npy")
        hist_c = np.load("data/scat_low"+want_matched+"_"+proxy+"_"+dm_ext+".npy")
        hist_s = np.load("data/scat_high"+want_matched+"_"+proxy+"_"+dm_ext+".npy")
        
        plt.plot(bin_cen,hist,linewidth=2.,c=cs[i],linestyle='-',label=dm_ext+' '.join(want_matched.split('_')))
        plt.fill_between(bin_cen,hist_c,hist_s,alpha=0.3)
        #plt.plot(bin_cen,hist_c,linewidth=2.,c=cs[i],linestyle='--')#,label="cens "+dm_ext+' '.join(want_matched.split('_')))
        #plt.plot(bin_cen,hist_s,linewidth=2.,c=cs[i],linestyle='--')#,label="sats "+dm_ext+' '.join(want_matched.split('_')))
#plt.plot(bin_cen,np.zeros(len(hist))-1,linewidth=2.,c='k',linestyle='-',label='median')
#plt.plot(bin_cen,np.zeros(len(hist))-1,linewidth=2.,c='k',linestyle='--',label='68\%')
#plt.plot(bin_cen,np.zeros(len(hist))-1,linewidth=2.,c='k',linestyle='-.',label='high')
plt.ylim([2.e11,2.e15])
plt.xlim([195.,2500.])
plt.yscale('log')
plt.xscale('log')
plt.legend(ncol=1)

plt.xlabel("Vrelax")
plt.ylabel("Subhalo Mass")
plt.savefig("scatter_mass_"+proxy+".png")
plt.show()
