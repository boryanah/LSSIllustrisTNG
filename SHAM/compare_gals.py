import numpy as np
import matplotlib.pyplot as plt
import Corrfunc
import plotparams
plotparams.default()
import sys
import h5py
from utils import *
plt.rc('text', usetex=True)

want_abundance = int(sys.argv[1])
want_environment = int(sys.argv[2])
want_satellites = 0
comp_sham_true = 1 # abund vs true dm from fp
if want_environment:
    lohi = sys.argv[3]#_lo _hi
else:
    lohi = ''
if want_satellites:
    sat_str = '_sats'
else:
    sat_str = ''
if lohi == '_lo': env_text = ' low env'
elif lohi == '_hi': env_text = ' high env'
elif lohi == '': env_text = ''

if want_abundance: abundance_matched = '_abund'
else: abundance_matched = ''

print("low/high environment, abundance matched = ",lohi,abundance_matched)

# simulation parameters
root = '/home/boryanah/lars/data_2500/'
res = 2500
Lbox = 205.

# matching subhalos in fp and dm
fp_dmo_inds = np.load(root+'sub_match_fp_dmo_'+str(res)+'.npy')
fp_inds = fp_dmo_inds[0]
dmo_inds = fp_dmo_inds[1]

# load data
SubhaloPos_fp = np.load(root+'SubhaloPos_fp.npy')/1000.
SubhaloMassType_fp = np.load(root+'SubhaloMassType_fp.npy')*1.e10
sub_inds_fp = np.arange(len(SubhaloMassType_fp)).astype(int)
SubhaloPos_dm = np.load(root+'SubhaloPos_dm.npy')/1000.
SubhaloMassType_dm = np.load(root+'SubhaloMassType_dm.npy')*1.e10
sub_inds_dm = np.arange(len(SubhaloMassType_dm)).astype(int)

# load other mass proxies
SubhaloVmax_fp = np.load(root+'SubhaloVmax_fp.npy')
SubhaloVelDisp_fp = np.load(root+'SubhaloVelDisp_fp.npy')
SubhaloVpeak_fp = np.load(root+'SubhaloVpeak_fp.npy')
SubhaloMpeak_fp = np.load(root+'SubhaloMpeak_fp.npy')*1.e10
SubhaloVpeak_fp[SubhaloVpeak_fp==0.] = SubhaloVmax_fp[SubhaloVpeak_fp==0.]#*20.
SubhaloVinfall_fp = np.load(root+'SubhaloVinfall_fp.npy')
SubhaloMinfall_fp = np.load(root+'SubhaloMinfall_fp.npy')*1.e10
SubhaloVrelax_fp = np.load(root+'SubhaloVrelax_fp.npy')

SubhaloVmax_dm = np.load(root+'SubhaloVmax_dm.npy')
SubhaloVelDisp_dm = np.load(root+'SubhaloVelDisp_dm.npy')
SubhaloVpeak_dm = np.load(root+'SubhaloVpeak_dm.npy')
SubhaloVpeak_dm[SubhaloVpeak_dm==0.] = SubhaloVmax_dm[SubhaloVpeak_dm==0.]#*20.
SubhaloMpeak_dm = np.load(root+'SubhaloMpeak_dm.npy')*1.e10
SubhaloVinfall_dm = np.load(root+'SubhaloVinfall_dm.npy')
SubhaloMinfall_dm = np.load(root+'SubhaloMinfall_dm.npy')*1.e10
SubhaloVrelax_dm = np.load(root+'SubhaloVrelax_dm.npy')

if want_environment:
    filename = '/home/boryanah/lars/test/CosmicWeb/WEB_CIC_256_DM_TNG300-2.hdf5'
    f = h5py.File(filename, 'r')
    d_smooth = f['density_smooth'][:,:,:] 

    # finding who belongs where in the cosmic web
    N_dim = 256
    gr_size = Lbox/N_dim
    if abundance_matched == '':
        halo_x = SubhaloPos_fp[:,0];halo_y = SubhaloPos_fp[:,1];halo_z = SubhaloPos_fp[:,2]
    else:
        halo_x = SubhaloPos_dm[:,0];halo_y = SubhaloPos_dm[:,1];halo_z = SubhaloPos_dm[:,2]
    i_cw = (halo_x/gr_size).astype(int)
    j_cw = (halo_y/gr_size).astype(int)
    k_cw = (halo_z/gr_size).astype(int)
    j_cw[j_cw == N_dim] = N_dim - 1 # fixing floating point issue
    
    # Environment definition
    env_cw = d_smooth[i_cw,j_cw,k_cw]

if want_satellites:
    i_before = np.arange(len(SubhaloVrelax_fp),dtype=int)
    GroupFirstSub_fp = np.load(root+'GroupFirstSub_fp.npy')
    SubhaloUniqueParent_fp = np.ones(len(SubhaloVrelax_fp),dtype=int)
    #SubhaloUniqueParent_fp[GroupFirstSub_fp[:1000000]] = 0
    SubhaloUniqueParent_fp[GroupFirstSub_fp] = 0
    sats = SubhaloUniqueParent_fp > 0
    ind_sats = i_before[sats]
    
# stellar mass proxy
lum_proxy = SubhaloMassType_fp[:,4]
i_lum_sort = np.argsort(lum_proxy)[::-1]
n_top = 13143#72000#13143#12000
i_lum_sort = i_lum_sort[:n_top]
if want_satellites:
    i_lum_sort, c1, c2_prop = np.intersect1d(ind_sats,i_lum_sort,return_indices=True)
lum_proxy_sorted = lum_proxy[i_lum_sort]
pos_fp_lsorted = SubhaloPos_fp[i_lum_sort]

# match the subhalos
if want_environment and abundance_matched == '':
    env_cw_sorted = env_cw[i_lum_sort]
    if lohi == '_hi':
        hi_quart = np.quantile(env_cw_sorted, .75)
        i_lum_sort = i_lum_sort[env_cw_sorted > hi_quart]
        print(len(i_lum_sort))
    elif lohi == '_lo':
        lo_quart = np.quantile(env_cw_sorted, .25)
        i_lum_sort = i_lum_sort[env_cw_sorted < lo_quart]
        print(len(i_lum_sort))
        
# get the intersected indices
inter_fp, c1, c2 = np.intersect1d(i_lum_sort,fp_inds,return_indices=True)
inter_dm = dmo_inds[c2]

#print(np.sum(SubhaloMpeak_dm[inter_dm]>0.))
#print(np.sum(SubhaloMpeak_fp[inter_fp]>0.))
#np.save("dm_inds"+lohi+sat_str+".npy",inter_dm)
#np.save("fp_inds"+lohi+sat_str+".npy",inter_fp)


print("Matched ",len(inter_dm)," out of ",n_top)

proxies = ['vrelax','vpeak','mpeak','mass','vmax','vinfall','minfall','vdisp']

s = 0.5
al = 0.5
ntot = len(proxies)
nrows = 2
ncols = ntot//nrows
plt.subplots(nrows,ncols,figsize=(ncols*5.3,nrows*4))
for i in range(len(proxies)):
    proxy = proxies[i]
    plot_no = i+1
    xlab = 'true value (DM)'
    ylab = '(true(FP)-true(DM))/true(DM)'
    
    if proxy == 'mass': proxy_dm = SubhaloMassType_dm[:,1]; proxy_fp = SubhaloMassType_fp[:,1]*(4.73/3.98)
    if proxy == 'vpeak': proxy_dm = SubhaloVpeak_dm; proxy_fp = SubhaloVpeak_fp
    if proxy == 'vmax': proxy_dm = SubhaloVmax_dm; proxy_fp = SubhaloVmax_fp
    if proxy == 'halfmass': proxy_dm = SubhaloMassInHalfRad_dm; proxy_fp = SubhaloMassInHalfRad_fp
    if proxy == 'vdisp': proxy_dm = SubhaloVelDisp_dm; proxy_fp = SubhaloVelDisp_fp
    if proxy == 's2r': proxy_dm = SubhaloVelDisp_dm**2*SubhaloHalfmassRad_dm; proxy_fp = SubhaloVelDisp_fp**2*SubhaloHalfmassRad_fp
    if proxy == 'vinfall': proxy_dm = SubhaloVinfall_dm; proxy_fp = SubhaloVinfall_fp
    if proxy == 'minfall': proxy_dm = SubhaloMinfall_dm; proxy_fp = SubhaloMinfall_fp
    if proxy == 'mpeak': proxy_dm = SubhaloMpeak_dm; proxy_fp = SubhaloMpeak_fp
    if proxy == 'vrelax': proxy_dm = SubhaloVrelax_dm; proxy_fp = SubhaloVrelax_fp


    if abundance_matched == '_abund':
        # indices in matched subhalo array ordered by proxy and number determined by fp; basically corresponding abundance matched inds
        inter_dm_abund = (np.argsort(proxy_dm[dmo_inds])[::-1])[:len(inter_fp)]

        # ordering the proxy values
        proxy_dm_abund = (proxy_dm[dmo_inds])[inter_dm_abund]        
        proxy_fp_abund = np.sort(proxy_fp[inter_fp])[::-1]

        if want_environment:
            # ordering the environment values
            #TESTINGenv_cw_sorted = (env_cw[dmo_inds])[inter_dm_abund]
            # how do you want to compute env
            env_cw_sorted = (env_cw[dmo_inds])[inter_dm]
            
            # choose low or high densities
            if lohi == '_hi':
                hi_quart = np.quantile(env_cw_sorted, .75)
                choice = env_cw_sorted > hi_quart
            elif lohi == '_lo':
                lo_quart = np.quantile(env_cw_sorted, .25)
                choice = env_cw_sorted < lo_quart

            # select the hi/lo environment values (proxy_fp/dm changed every iteration)
            proxy_dm_abund = proxy_dm_abund[choice]
            proxy_fp_abund = proxy_fp_abund[choice]

        if comp_sham_true:
            # sort the proxies in terms of the luminosity
            i_sort = np.argsort(SubhaloMassType_fp[inter_fp,4])[::-1]
            proxy_dm = (proxy_dm[inter_dm])[i_sort]
            #proxy_fp = (proxy_fp[inter_fp])[i_sort]
            xlab = 'true value (DM)'
            ylab = '(SHAM(DM)-true(DM))/true(DM)'
            if want_environment:
                proxy_dm = proxy_dm[choice]
                #proxy_fp = proxy_fp[choice]
            
            frac = (proxy_dm-proxy_dm_abund)/proxy_dm
        else:
            frac = (proxy_fp_abund-proxy_dm_abund)/proxy_dm_abund
    else:
        proxy_dm = proxy_dm[inter_dm]
        proxy_fp = proxy_fp[inter_fp]

        frac = (proxy_fp-proxy_dm)/proxy_dm
    
    print(proxy," has over 0 = ",np.sum(proxy_fp>0.)," FP subhalos and ",np.sum(proxy_dm>0.)," DM subhalos")
    
    plt.subplot(nrows,ncols,plot_no)
    line = np.linspace(70.,2.e15,1000)
    plt.plot(line,np.zeros(len(line)),'k--',alpha=0.3,linewidth=2.,label=proxy+env_text)

    
    
    plt.legend()
    ylim = [-0.3,0.3]#[-1,0.8]
    plt.ylim(ylim)
    plt.xscale('log')
    if proxy[0] == 'v':
        xlim = [70,2000.]
        plt.xlim(xlim)
        yscale = 'linear'
        grs = 50
    elif proxy[0] == 'm':
        xlim = [1.e11,2.e15]
        plt.xlim(xlim)
        yscale = 'linear'
        grs = 50

    #plt.scatter(proxy_dm,frac,label=proxy,s=s,alpha=al)
    plt.hexbin(proxy_dm+1.e-3, frac, gridsize=grs, bins='log', extent=(np.log10(xlim[0]),np.log10(xlim[1]),ylim[0],ylim[1]),xscale='log', yscale=yscale, cmap='Greys')

    plot_median(proxy_dm,frac,np.log10(xlim[0]),np.log10(xlim[1]),n_bins=41)
    
    if plot_no >= ntot-ncols+1:
        plt.xlabel(xlab)
    else:
        plt.gca().axes.xaxis.set_ticklabels([])
    if plot_no%ncols == 1:
        plt.ylabel(ylab,fontsize=16)
    else:
        plt.gca().axes.yaxis.set_ticklabels([])

plt.savefig("SHAM"+abundance_matched+"_frac"+lohi+".png")
#plt.show()
quit()

# get the positions of the intersection
xyz_fp = SubhaloPos_fp[inter_fp]
xyz_dm = SubhaloPos_dm[inter_dm]
xyz_fp_all = pos_fp_lsorted

def get_corr(xyz,wei=None,want_plot=None):
    # split into x, y and z
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    # bins
    N_bin = 31
    bins = np.logspace(-1,1.5,N_bin)
    bin_centers = (bins[:-1] + bins[1:])/2.
        
    # compute corrfunc
    results_DM = Corrfunc.theory.xi(X=x,Y=y,Z=z,boxsize=Lbox,nthreads=16,binfile=bins, weights=wei)
    Corr = results_DM['xi']
    
    if want_plot is not None:
        plt.plot(bin_centers, Corr*bin_centers**2,lw=3.,ls='-',label=want_plot)#,c='silver')
        plt.xscale('log')
        plt.ylabel(r'$\xi(r) r^2$')
        plt.xlabel(r'$r$ [Mpc/h]')
        plt.ylim([0,100])
        plt.legend(loc='lower right')

    return Corr

fs = (12,8)
fig = plt.figure(figsize=fs)
corr_dm = get_corr(xyz_dm,want_plot='DM matched')
corr_fp = get_corr(xyz_fp,want_plot='FP matched')
corr_fp_all = get_corr(xyz_fp_all,want_plot='FP all')
plt.legend()
plt.savefig('figs/corr_matched.png')
plt.show()
plt.close()