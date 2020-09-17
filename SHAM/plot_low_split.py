import numpy as np
import matplotlib.pyplot as plt
import plotparams
#plotparams.default()
plotparams.buba()

proxies = ['vrelax','s2r','mass','vinfall','vpeak','vmax','halfmass','vdisp']
sec_props = ['mstar','parent','g-r','sSFR','form']
lowm = 'lowm: '
names = ['low mass (lowm)','high mass',lowm+'satellites',lowm+'centrals',lowm+'blue',lowm+'red',lowm+'SFGs',lowm+'quiescent',lowm+'late-forming',lowm+'early-forming']
bin_centers = np.load("data_split/bin_centers.npy")
want_matched = '_matched'
proxy = proxies[4]

import distinct_colours
cs = distinct_colours.get_distinct(4)

nprops = len(sec_props)
nrows = 1#2
ncols = nprops
ntot = nrows*ncols
plt.subplots(nrows,ncols,figsize=(ncols*5.6,nrows*5.9))
plot_no = 0
for i in range(nprops):
    for i_bin in range(2):
        if i == 0:
            color1 = 'dodgerblue'
            color2 = '#CC6677'
            fc1 = 'royalblue'
            
        else:
            color1 = cs[2]#'mediumblue'
            color2 = cs[3]#'mediumorchid'
            fc1 = cs[2]#'navy'#'#089FFF'
            fc2 = cs[3]#'blueviolet'
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
            plt.plot(bin_centers,ratio,linewidth=2.,color=color1,label=label)#sec_prop+" "+str(i_bin))
            plt.fill_between(bin_centers,ratio+error,ratio-error,alpha=0.1, edgecolor=color1, facecolor=fc1)
            #plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c=color1,fmt='o',capsize=4,label=label)#sec_prop+" "+str(i_bin))
        elif i_bin == 1 and i == 0:
            # orange onesc
            plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c=color2,fmt='o',capsize=4,label=label)#sec_prop+" "+str(i_bin))
        else:
            plt.plot(bin_centers,ratio,linewidth=2.,color=color2,label=label)
            plt.fill_between(bin_centers,ratio+error,ratio-error,alpha=0.1, edgecolor=color2, facecolor=fc2)
            
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

plt.savefig("SHAM"+want_matched+"_split_ratio_all.pdf")
#plt.show()
