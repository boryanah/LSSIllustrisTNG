import numpy as np
import matplotlib.pyplot as plt
import plotparams
plotparams.buba()

proxies = ['vrelax','vpeak','vmax','vinfall','mpeak','minfall','mass','mmax']#vdisp
lab_proxies = [r'$V_\textrm{relax}$',r'$V_\textrm{peak}$',r'$V_\textrm{circ,max}$',r'$V_\textrm{infall}$',r'$M_\textrm{peak}$',r'$M_\textrm{infall}$',r'$M_{\rm SUBFIND}$',r'$M_\textrm{circ,max}$']
bin_centers = np.load("data/bin_centers.npy")
want_matched = ''#'_matched'#''#'_matched'
want_fps = ['','fp']
colors = ['#CC6677','#DDCC77']

ntot = len(proxies)
nrows = 4
ncols = ntot//nrows
plt.subplots(nrows,ncols,figsize=(ncols*5.3,nrows*4))
for i in range(len(proxies)):
    for j in range(len(want_fps)):
        want_fp = want_fps[j]
        color = colors[j]
        
        proxy = proxies[i]
        lab_proxy = lab_proxies[i]
        ratio = np.load("data/SHAM"+want_matched+"_ratio"+want_fp+"_"+proxy+".npy")
        error = np.load("data/SHAM"+want_matched+"_ratio"+want_fp+"_"+proxy+"_error.npy")
        plot_no = i+1


        plt.subplot(nrows,ncols,plot_no)
        plt.plot(bin_centers,np.ones(len(ratio)),'k--',alpha=0.3,linewidth=2.)
        
        if i != 0:
            # always plot the fiducial
            plt.plot(bin_fid,ratio_fid,linewidth=2.,color='#1B2ACC')
            plt.fill_between(bin_fid,ratio_fid+error_fid,ratio_fid-error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')

            # plot the current proxy
            plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c=color,fmt='o',capsize=4)

        else:
            # fiducial case is first so record
            #ratio_fid = ratio.copy()
            #error_fid = error.copy()
            bin_fid = bin_centers.copy()
            ratio_fid = np.load("data/SHAM_ratio_shuff.npy")
            error_fid = np.load("data/SHAM_ratio_shuff_error.npy")

            # plot the fiducial
            plt.plot(bin_fid,ratio_fid,linewidth=2.,color='#1B2ACC',label='basic HOD')
            plt.fill_between(bin_fid,ratio_fid+error_fid,ratio_fid-error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')

            plt.errorbar(bin_centers,ratio,yerr=error,ls='-',c=color,fmt='o',capsize=4)

        if i == 0 and j == 0:
            plt.errorbar(bin_centers,np.zeros(len(ratio)),yerr=error,ls='-',c=colors[0],fmt='o',capsize=4,label='SHAM (DM)')
            plt.errorbar(bin_centers,np.zeros(len(ratio)),yerr=error,ls='-',c=colors[1],fmt='o',capsize=4,label='SHAM (FP)')
            plt.legend(loc='upper left',fontsize=18)
        if want_fp == 'fp':
            text_x = 4.
            if (plot_no-1)//ncols == nrows-1:
                plt.ylim([0.1,1.3])
                plt.text(text_x,1.1,lab_proxy,fontsize=18)
            elif (plot_no-1)//ncols == nrows-2:
                plt.ylim([0.55,1.3])#([0.text_x,1.3])
                plt.text(text_x,1.1,lab_proxy,fontsize=18)
            elif (plot_no-1)//ncols == nrows-3:
                plt.ylim([0.55,1.3])#([0.8,1.3])
                plt.text(text_x,1.1,lab_proxy,fontsize=18)
            else:
                plt.ylim([0.75,1.5])
                plt.text(text_x,1.3,lab_proxy,fontsize=18)
        else:
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
            plt.ylabel(r'$\xi_{\rm model}/\xi_{\rm TNG300}$')
        else:
            plt.gca().axes.yaxis.set_ticklabels([])
        
plt.savefig("SHAM"+want_matched+"_ratio_all.pdf")
plt.show()
