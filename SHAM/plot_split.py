import numpy as np
import matplotlib.pyplot as plt
import plotparams
#plotparams.default()
plotparams.buba()

proxies = ['vrelax','s2r','mass','vinfall','vpeak','vmax','halfmass','vdisp']
sec_props = ['mstar','g-r','parent','sSFR']
lowm = ''
names = ['low mass (lowm)','high mass',lowm+'blue',lowm+'red',lowm+'satellites',lowm+'centrals',lowm+'SFGs',lowm+'quiescent']
bin_centers = np.load("data_split/bin_centers.npy")
want_matched = '_matched'
proxy = proxies[4]

nprops = len(sec_props)
nrows = 1#2
ncols = nprops
ntot = nrows*ncols
plt.subplots(nrows,ncols,figsize=(ncols*4.4,nrows*5.9))
plot_no = 0
for i in range(nprops):
    for i_bin in range(2):
        sec_prop = sec_props[i]
        label = names[i*2+i_bin]
        ratio = np.load("data_split/SHAM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+".npy")
        error = np.load("data_split/SHAM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+"_error.npy")
        plot_no = i+1

        plt.subplot(nrows,ncols,plot_no)
        plt.plot(bin_centers,np.ones(len(ratio)),'k--',linewidth=2.)
        if i_bin == 0:
            # plot the current proxy
            # blue
            plt.plot(bin_centers,ratio,linewidth=2.,color='#1B2ACC',label=label)#sec_prop+" "+str(i_bin))
            plt.fill_between(bin_centers,ratio+error,ratio-error,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        else:
            # orange ones
            plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c='orange',fmt='o',capsize=4,label=label)#sec_prop+" "+str(i_bin))

        # always plot the fiducial
        #plt.plot(bin_centers,ratio_fid,linewidth=2.,color='#1B2ACC')
        #plt.fill_between(bin_centers,ratio_fid+error_fid,ratio_fid-error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')

        plt.legend()
        plt.ylim([0.7,1.3])
        #origplt.xlim([.7,15])
        plt.xlim([.2,15])
        plt.xscale('log')

        if plot_no >= ntot-ncols+1:
            plt.xlabel(r'$r$ [Mpc/h]')
        if plot_no%ncols == 1:
            plt.ylabel(r'$\xi_{\rm SHAM}/\xi_{\rm TNG300}$')
        else:
            plt.gca().axes.yaxis.set_ticklabels([])

plt.savefig("SHAM"+want_matched+"_split_ratio_all.png")
#plt.show()
