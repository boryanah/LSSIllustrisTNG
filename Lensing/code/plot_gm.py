import numpy as np
import matplotlib.pyplot as plt
import plotparams
#plotparams.default()
plotparams.buba()

#opts = np.array(['shuff', 'env_mass', 'min_pot', 'partial_vani', 'partial_env_cw', 'partial_tot_pot', 'conc', 'spin'])
#opts = np.array(['shuff', 'env_mass', 'min_pot', 'partial_vani', 'partial_env_cw', 'partial_tot_pot', 'partial_s2r', 'partial_fenv'])
opts = np.array(['shuff','partial_vani', 'partial_env_cw', 'partial_tot_pot'])
types = ['bias','corr. coeff.']

bin_centers = np.load("../data_gm/bin_fid.npy")
bias_fid = np.load("../data_gm/bias_mean_fid.npy")
corr_coeff_fid = np.load("../data_gm/corr_coeff_mean_fid.npy")
# TESTING
#bias_fid = np.load("../data_gm/corr_gg_mean_fid.npy")*bin_centers**2
#bias_error_fid = np.load("../data_gm/corr_gg_error_fid.npy")*bin_centers**2

bias_error_fid = np.load("../data_gm/bias_error_fid.npy")
corr_coeff_error_fid = np.load("../data_gm/corr_coeff_error_fid.npy")

nprops = len(opts)
nrows = 2
ncols = nprops
ntot = nrows*ncols
plt.subplots(nrows,ncols,figsize=(ncols*4.8,nrows*5))
plot_no = 0
for i in range(nprops):
    for i_type in range(2):
        opt = opts[i]
        print(opt)
        
        if i_type == 0:
            bias = np.load("../data_gm/bias_mean_"+opt+".npy")
            bias_error = np.load("../data_gm/bias_error_"+opt+".npy")
        else:    
            bias = np.load("../data_gm/corr_coeff_mean_"+opt+".npy")
            bias_error = np.load("../data_gm/corr_coeff_error_"+opt+".npy")
        
        plot_no = i_type*ncols+i+1

        plt.subplot(nrows,ncols,plot_no)
        plt.plot(bin_centers,np.ones(len(bias)),'k--',linewidth=2.)
        
        # orange ones
        plt.errorbar(bin_centers,bias,yerr=bias_error,ls='-',c='orange',fmt='o',capsize=4,label="-".join(opt.split('_'))+" "+types[i_type])

        if plot_no == 1:
            lab_fid = 'true'
        else:
            lab_fid = ''
        
        # always plot the fiducial
        if i_type == 0:
            plt.plot(bin_centers,bias_fid,linewidth=2.,color='#1B2ACC',label=lab_fid)
            plt.fill_between(bin_centers,bias_fid+bias_error_fid,bias_fid-bias_error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        else:
            plt.plot(bin_centers,corr_coeff_fid,linewidth=2.,color='#1B2ACC',label=lab_fid)
            plt.fill_between(bin_centers,corr_coeff_fid+corr_coeff_error_fid,corr_coeff_fid-corr_coeff_error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
            
        plt.legend(loc='upper left')
        if i_type == 0:
            plt.ylim([0.7,1.5])
        if i_type == 1:
            plt.ylim([0.5,1.32])
        #origplt.xlim([.7,15])
        plt.xlim([.7,15])
        plt.xscale('log')

        if plot_no >= ntot-ncols+1:
            plt.xlabel(r'$r$ [Mpc/h]')
        if plot_no%ncols == 1 and i_type == 0:
            plt.ylabel(r'$(\xi_{\rm gg}/\xi_{\rm mm})^{1/2}$')
        elif plot_no%ncols == 1 and i_type == 1:
            plt.ylabel(r'$\xi_{\rm gm}/(\xi_{\rm gg} \xi_{\rm mm})^{1/2}$')
        else:
            plt.gca().axes.yaxis.set_ticklabels([])

plt.savefig("bias_corr_coeff.png")
plt.show()
