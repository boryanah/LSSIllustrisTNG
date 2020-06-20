import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Corrfunc

import matplotlib.pyplot as plt
import matplotlib.ticker

import plotparams
#plotparams.default()
plotparams.buba()

proxy = np.array(['m200c', 'm200m', 'mfof', 'vmax', 'vpeak', 'm500c', 's2r'])
#opt = np.array(['shuff', 'env_mass', 'vani', 'vdisp', 'moverr', 'conc', 'form', 'spin'])
opt = np.array(['shuff','partial_fenv', 'partial_s2r','partial_vani'])
#option = np.array(['hydrosim.', 'local env.', 'vel. anis.','vel. disp.',  r'$M_{\rm cent}/R_{\rm halfmass}$', 'concentration','form. epoch', 'halo spin'])
option = np.array(['hydrosim.',r'$N_{\rm gal} \times$ local env.',  r'$N_{\rm gal} \times \sigma^2 R_{\rm halfmass}$', r'$N_{\rm gal} \times$ vel. anis.'])


N = len(opt)

plt.subplots(1,N)

for i in range(N):
    #i = 0

    figname = 'figs/paper/Corr_'+opt[i]+'_'+proxy[1]+'.png'

    plt.subplot(1,N,i+1)
    if i == 0:
        bin_centers, a,b,c,d,xi_fid_shuff,xi_fid_shuff_err,rat_fid_shuff,rat_fid_shuff_err = np.loadtxt('/home/boryanah/lars/test/'+figname[:-4]+'.txt',unpack=True)

    bin_centers, a,b,c,d,Corr_opt_mean,Corr_opt_err,Rat_DM_mean,Rat_DM_err = np.loadtxt('/home/boryanah/lars/test/'+figname[:-4]+'.txt',unpack=True)


    plt.plot(bin_centers,np.ones(len(bin_centers)),'k--',alpha=0.2)

    if i > 0:
        print(1)
        plt.plot(bin_centers,rat_fid_shuff,linewidth=1.,color='#1B2ACC',alpha=0.5)
        plt.fill_between(bin_centers, rat_fid_shuff+rat_fid_shuff_err, rat_fid_shuff-rat_fid_shuff_err,alpha=0.05, edgecolor='#1B2ACC', facecolor='#089FFF')
        plt.plot(bin_centers[bin_centers>1.],rat_fid_shuff[bin_centers>1.],linewidth=1.,color='#1B2ACC')
        plt.fill_between(bin_centers[bin_centers>1.], (rat_fid_shuff+rat_fid_shuff_err)[bin_centers>1.],(rat_fid_shuff-rat_fid_shuff_err)[bin_centers>1.],alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
    else:
        print(2)
        plt.plot(bin_centers,rat_fid_shuff,linewidth=1.,color='#1B2ACC',alpha=0.5)
        plt.fill_between(bin_centers, rat_fid_shuff+rat_fid_shuff_err, rat_fid_shuff-rat_fid_shuff_err,alpha=0.05, edgecolor='#1B2ACC', facecolor='#089FFF')
        plt.plot(bin_centers[bin_centers>1.],rat_fid_shuff[bin_centers>1.],linewidth=1.,color='#1B2ACC',label=''+option[i])
        plt.fill_between(bin_centers[bin_centers>1.], (rat_fid_shuff+rat_fid_shuff_err)[bin_centers>1.],(rat_fid_shuff-rat_fid_shuff_err)[bin_centers>1.],alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        
    if i > 0:
        plt.errorbar(bin_centers,Rat_DM_mean,yerr=Rat_DM_err,ls='-',c='orange',alpha=0.5,fmt='o',capsize=4)
        plt.errorbar(bin_centers[bin_centers>1.],Rat_DM_mean[bin_centers>1.],yerr=Rat_DM_err[bin_centers>1.],ls='-',c='orange',fmt='o',capsize=4,label=''+option[i])
    
    plt.legend(loc='upper left',ncol=2)
    plt.xscale('log')
    if i == 0:
        plt.ylabel(r'$\xi(r)_{\rm sec. prop.} / \xi(r)_{\rm bHOD}$') 
    else:
        plt.gca().axes.yaxis.set_ticklabels([])

    if i > N:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.xlabel(r'$r$ [Mpc/h]')

    tick_spacing = 0.2
    plt.gca().axes.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    locmin = matplotlib.ticker.LogLocator(base=10,subs=(0.1,.2,.3,0.4,.5,0.6,.7,0.8,.9),numticks=10)
    plt.gca().axes.xaxis.set_minor_locator(locmin)
    plt.gca().axes.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.gca().axes.xaxis.set_major_locator(ticker.LogLocator(base=10))
    plt.gca().axes.xaxis.set_major_locator(ticker.FixedLocator([0.1,1,10]))
    plt.ylim([0.65,1.5])
    plt.xlim([0.6,16]) 
    
    if i == 1: plt.text(3, 0.7, r'$r = 0.2$')
    if i == 2: plt.text(3, 0.7, r'$r = 0.45$')
    if i == 3: plt.text(3, 0.7, r'$r = 0.45$')
    # tot pot uses 0.05; min pot uses 0.041
    #plt.text(0.1, 1.3, 'option = '+option[i])
plt.savefig('/home/boryanah/m200m_all_secondary.pdf')
plt.show()
plt.close()

quit()
plt.errorbar(bin_centers, Corr_mean_FP*bin_centers**2,yerr=Corr_err_FP*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label=str(res[i])+ ' FP')
plt.errorbar(bin_centers, Corr_mean_DM*bin_centers**2,yerr=Corr_err_DM*bin_centers**2,lw=2.,ls='-',c='orange',fmt='o',capsize=2,label=str('TNG DM gals'))
plt.errorbar(bin_centers, Corr_mean_DM_shuff*bin_centers**2,yerr=Corr_err_DM_shuff*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label='TNG DM env')
