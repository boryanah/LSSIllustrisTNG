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
#opt = np.array(['shuff', 'env_mass', 's2r', 'vani', 'vdisp', 'conc', 'form', 'spin'])
# TESTING
opt = np.array(['shuff', 'env_mass', 's2r', 'vani', 'vdisp', 'conc', 'min_pot', 'tot_pot'])

#option = np.array(['hydrosim.', 'local env.', r'$\sigma^2 R_{\rm halfmass}$', 'vel. anis.','vel. disp.', 'concentration','form. epoch', 'halo spin'])
# TESTING
option = np.array(['hydrosim.', 'local env.', r'$\sigma^2 R_{\rm halfmass}$', 'vel. anis.','vel. disp.', 'concentration','min. pot.', 'total pot.'])

N = len(opt)

#fs = (30,30)#(8.2/1.2,6.5/1.2)
#fig =plt.figure(figsize=fs)
plt.subplots(2,N//2)#,figsize=fs)

for i in range(N):
    #i = 0

    figname = 'figs/paper/Corr_'+opt[i]+'_'+proxy[1]+'.png'

    plt.subplot(2,N//2,i+1)
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
        '''
          if i> 0:
        plt.plot(bin_centers,rat_fid_shuff,linewidth=1.,color='#1B2ACC')
        plt.fill_between(bin_centers, rat_fid_shuff+rat_fid_shuff_err, rat_fid_shuff-rat_fid_shuff_err,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        else:
        plt.plot(bin_centers,rat_fid_shuff,linewidth=1.,color='#1B2ACC',label=''+option[i])
        plt.fill_between(bin_centers, rat_fid_shuff+rat_fid_shuff_err, rat_fid_shuff-rat_fid_shuff_err,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        '''

    if i > 0:
        plt.errorbar(bin_centers,Rat_DM_mean,yerr=Rat_DM_err,ls='-',c='orange',alpha=0.5,fmt='o',capsize=4)
        plt.errorbar(bin_centers[bin_centers>1.],Rat_DM_mean[bin_centers>1.],yerr=Rat_DM_err[bin_centers>1.],ls='-',c='orange',fmt='o',capsize=4,label=''+option[i])

        '''
        if i > 0:
        plt.errorbar(bin_centers,Rat_DM_mean,yerr=Rat_DM_err,ls='-',c='orange',fmt='o',capsize=4,label=''+option[i])
        '''
    plt.legend(loc='upper left',ncol=2)
    plt.xscale('log')
    if i == 0 or i==N//2:
        #plt.ylabel(r'$\xi(r)_{\rm TNG} / \xi(r)_{\rm option}$') # inv
        #plt.ylabel(r'$\xi(r)_{\rm option} / \xi(r)_{\rm TNG}$')
        plt.ylabel(r'$\xi(r)_{\rm sec. prop.} / \xi(r)_{\rm bHOD}$') # Daniel
    else:
        plt.gca().axes.yaxis.set_ticklabels([])#.set_visible(False)
        #plt.ylabel(r'$\xi(r)_{\rm TNG, opt} / \xi(r)_{\rm TNG}$')

    if i < N//2:
        #plt.gca().axes.get_xticklabels([''])#.set_visible(False)
        plt.gca().axes.xaxis.set_ticklabels([])#.set_visible(False)
    plt.xlabel(r'$r$ [Mpc/h]')

    tick_spacing = 0.2
    plt.gca().axes.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    locmin = matplotlib.ticker.LogLocator(base=10,subs=(0.1,.2,.3,0.4,.5,0.6,.7,0.8,.9),numticks=10)
    plt.gca().axes.xaxis.set_minor_locator(locmin)
    plt.gca().axes.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.gca().axes.xaxis.set_major_locator(ticker.LogLocator(base=10))
    plt.gca().axes.xaxis.set_major_locator(ticker.FixedLocator([0.1,1,10]))
    #plt.ylim([0.38,1.44]) # inv
    #plt.ylim([0.65,1.9])
    plt.ylim([0.65,4])# TESTING

    #plt.xlim([0.075,16])
    plt.xlim([0.6,16]) 
    #plt.ylim([0.95,1.44]) # Daniel
    
    if i == 1: plt.text(8., 0.7, r'I')
    if i == 2: plt.text(8., 0.7, r'II')
    if i == 3: plt.text(8., 0.7, r'III')
    #if i == 4: plt.text(4., 0.7, r'II')
    if i == 4: plt.text(8., 0.7, r'IV')
    if i == 5: plt.text(8., 0.7, r'V')
    if i == 6: plt.text(8., 0.7, r'VI')
    if i == 7: plt.text(8., 0.7, r'VII')
    #plt.text(0.1, 1.3, option = '+option[i])
plt.savefig('/home/boryanah/m200m_all_secondary.pdf')
plt.show()
plt.close()

quit()
plt.errorbar(bin_centers, Corr_mean_FP*bin_centers**2,yerr=Corr_err_FP*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label=str(res[i])+ ' FP')
plt.errorbar(bin_centers, Corr_mean_DM*bin_centers**2,yerr=Corr_err_DM*bin_centers**2,lw=2.,ls='-',c='orange',fmt='o',capsize=2,label=str('TNG DM gals'))
plt.errorbar(bin_centers, Corr_mean_DM_shuff*bin_centers**2,yerr=Corr_err_DM_shuff*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label='TNG DM env')
