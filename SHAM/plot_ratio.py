import numpy as np
import matplotlib.pyplot as plt
import plotparams
#plotparams.default()
plotparams.buba()

#proxies = ['vrelax','vpeak','vdisp','vinfall','vmax','s2r','mass','halfmass']
proxies = ['vrelax','vpeak','vmax','vinfall','mpeak','vdisp','minfall','mass']
bin_centers = np.load("data/bin_centers.npy")
want_matched = '_matched'#''#'_matched'
want_fp = ''#'fp'#''

ntot = len(proxies)
nrows = 4
ncols = ntot//nrows
plt.subplots(nrows,ncols,figsize=(ncols*5.3,nrows*4))
for i in range(len(proxies)):
    proxy = proxies[i]
    ratio = np.load("data/SHAM"+want_matched+"_ratio"+want_fp+"_"+proxy+".npy")
    error = np.load("data/SHAM"+want_matched+"_ratio"+want_fp+"_"+proxy+"_error.npy")
    plot_no = i+1
    
    plt.subplot(nrows,ncols,plot_no)
    plt.plot(bin_centers,np.ones(len(ratio)),'k--',linewidth=2.)
    if i != 0:
        # always plot the fiducial
        plt.plot(bin_fid,ratio_fid,linewidth=2.,color='#1B2ACC')
        plt.fill_between(bin_fid,ratio_fid+error_fid,ratio_fid-error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')

        # plot the current proxy
        plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c='orange',fmt='o',capsize=4,label=proxy)

    else:
        # fiducial case is first so record
        #ratio_fid = ratio.copy()
        #error_fid = error.copy()
        bin_fid = bin_centers.copy()
        ratio_fid = np.load("data/SHAM_ratio_shuff.npy")
        error_fid = np.load("data/SHAM_ratio_shuff_error.npy")

        # plot the fiducial
        plt.plot(bin_fid,ratio_fid,linewidth=2.,color='#1B2ACC',label='bHOD')
        plt.fill_between(bin_fid,ratio_fid+error_fid,ratio_fid-error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        
        plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c='orange',fmt='o',capsize=4,label=proxy)

    plt.legend(loc='upper right')
    if want_fp == 'fp':
        if (plot_no-1)//ncols == nrows-1:
            plt.ylim([0.1,1.3])
        elif (plot_no-1)//ncols == nrows-2:
            plt.ylim([0.6,1.3])
        elif (plot_no-1)//ncols == nrows-3:
            plt.ylim([0.8,1.3])
        else:
            plt.ylim([0.8,2.])
    else:#if want_matched == '_matched':
        if (plot_no-1)//ncols == nrows-1:
            plt.ylim([0.1,1.3])
        elif (plot_no-1)//ncols == nrows-2:
            plt.ylim([0.55,1.3])
        elif (plot_no-1)//ncols == nrows-3:
            plt.ylim([0.5,1.3])
        else:
            plt.ylim([0.8,1.3])
    plt.xlim([.2,15])
    plt.xscale('log')

    if plot_no >= ntot-ncols+1:
        plt.xlabel(r'$r$ [Mpc/h]')
    if plot_no%ncols == 1:
        plt.ylabel(r'$\xi_{\rm SHAM}/\xi_{\rm TNG300}$')
    else:
        plt.gca().axes.yaxis.set_ticklabels([])
        
plt.savefig("SHAM"+want_matched+"_ratio"+want_fp+"_all.png")
#plt.show()
