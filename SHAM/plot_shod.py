import matplotlib.pyplot as plt
import numpy as np
import plotparams
plotparams.buba()

proxy = 'vrelax'

bin_cen = np.load("data/hist_bins_"+proxy+".npy")
hist_fid = np.load("data/hist_"+proxy+"_dm.npy")

wants = ['_matched','']
exts = ['fp','dm']
cs = ['gray','orange','dodgerblue']
i = -1
plt.figure(3,figsize=(15,12))
for want_matched in wants:
    for dm_ext in exts:
        if dm_ext == 'fp' and want_matched == '_matched': continue
        i += 1
        hist = np.load("data/hist"+want_matched+"_"+proxy+"_"+dm_ext+".npy")
        hist_c = np.load("data/hist_cen"+want_matched+"_"+proxy+"_"+dm_ext+".npy")
        hist_s = np.load("data/hist_sat"+want_matched+"_"+proxy+"_"+dm_ext+".npy")
        
        plt.plot(bin_cen,hist,linewidth=2.,c=cs[i],linestyle='-',label=dm_ext+' '.join(want_matched.split('_')))
        plt.plot(bin_cen,hist_c,linewidth=2.,c=cs[i],linestyle='--')#,label="cens "+dm_ext+' '.join(want_matched.split('_')))
        plt.plot(bin_cen,hist_s,linewidth=2.,c=cs[i],linestyle='-.')#,label="sats "+dm_ext+' '.join(want_matched.split('_')))
plt.plot(bin_cen,np.zeros(len(hist))-1,linewidth=2.,c='k',linestyle='-',label='total')
plt.plot(bin_cen,np.zeros(len(hist))-1,linewidth=2.,c='k',linestyle='--',label='cens')
plt.plot(bin_cen,np.zeros(len(hist))-1,linewidth=2.,c='k',linestyle='-.',label='sats')
plt.ylim([-0.2,1.2])
plt.xlim([95.,2500.])
#plt.yscale('log')
plt.xscale('log')
plt.legend(ncol=2)

plt.xlabel("Vrelax")
plt.ylabel("SHOD")
plt.savefig("SHOD_"+proxy+".png")
plt.show()
