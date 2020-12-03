import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Corrfunc
import plotparams
plotparams.default()
import sys
from mpi4py import MPI
from utils import *

# parameter choices
proxy = sys.argv[1]
want_matched = ''#'_matched'#'_matched'#'_matched'#'_matched'#''#'_matched'
exclude_most_mass = 0#1 # skip the first galaxies
if exclude_most_mass:
    first_gal = 23450#869
else:
    first_gal = 0
N_gal_choice = int(sys.argv[2]) #12000 # 13173 if not exclude #14126 if exclude
N_div = int(sys.argv[3]) # how many subpops are you dividing it into
isolate_low_mass = 1#1 is only for splits and not for mstar split
if len(sys.argv) > 4:
    sec_prop = sys.argv[4]
    want_matched = '_matched'
    print("It is required you work with the matched samples")
else:
    sec_prop = 'full';
root = '/home/boryanah/lars/data_2500/'
light_dir = "/mnt/gosling1/boryanah/Illustris_light/TNG300/"
res = 2500
Lbox = 205.
i_bin = MPI.COMM_WORLD.Get_rank()

# this is used if we want to add scatter in the Mstar-Vpeak/relax relationship
want_scatter = 0

# START HERE

# matching subhalos in fp and dm
fp_dmo_inds = np.load(root+'sub_match_fp_dmo_'+str(res)+'.npy')
fp_inds = fp_dmo_inds[0]
dmo_inds = fp_dmo_inds[1]


# load the DM data; set this to fp to test SHAM within hydro
dm_ext = 'dm'#'dm'#'fp'
SubhaloVmax_dm = np.load(root+'SubhaloVmax_'+dm_ext+'.npy')
new_type = '_new'# newest is cut at vmax = 100; new at vmax = 50; null is old format/original
SubhaloVpeak_dm = np.load(root+'SubhaloVpeak'+new_type+'_'+dm_ext+'.npy')
SubhaloVpeak_dm[SubhaloVpeak_dm==0.] = SubhaloVmax_dm[SubhaloVpeak_dm==0.]
# og
SubhaloVinfall_dm = np.load(root+'SubhaloVinfall'+new_type+'_'+dm_ext+'.npy')
# TESTING
#SubhaloVinfall_dm = np.load(root+'SubhaloVpeakinfall'+new_type+'_'+dm_ext+'.npy')
SubhaloVinfall_dm[SubhaloVinfall_dm==0.] = SubhaloVmax_dm[SubhaloVinfall_dm==0.]
SubhaloVrelax_dm = np.load(root+'SubhaloVrelax'+new_type+'_'+dm_ext+'.npy')
SubhaloPos_dm = np.load(root+'SubhaloPos_'+dm_ext+'.npy')/1000.
#SubhaloHalfmassRad_dm = np.load(root+'SubhaloHalfmassRad_'+dm_ext+'.npy')
#SubhaloMassInHalfRad_dm = np.load(root+'SubhaloMassInHalfRad_'+dm_ext+'.npy')
#SubhaloVmaxRad_dm = np.load(root+'SubhaloVmaxRad_'+dm_ext+'.npy')
SubhaloMassInMaxRad_dm = np.load(root+'SubhaloMassInMaxRad_'+dm_ext+'.npy')*1.e10
SubhaloMassType_dm = np.load(root+'SubhaloMassType_'+dm_ext+'.npy')
SubhaloVelDisp_dm = np.load(root+'SubhaloVelDisp_'+dm_ext+'.npy')
SubhaloMpeak_dm = np.load(root+'SubhaloMpeak'+new_type+'_'+dm_ext+'.npy')
SubhaloMpeak_dm[SubhaloMpeak_dm==0.] = SubhaloMassType_dm[SubhaloMpeak_dm==0.,1]
# og
SubhaloMinfall_dm = np.load(root+'SubhaloMinfall'+new_type+'_'+dm_ext+'.npy')
# TESTING
#SubhaloMinfall_dm = np.load(root+'SubhaloMpeakinfall'+new_type+'_'+dm_ext+'.npy')
SubhaloMinfall_dm[SubhaloMinfall_dm==0.] = SubhaloMassType_dm[SubhaloMinfall_dm==0.,1]

# load the FP data
fp_ext = 'fp'
GroupFirstSub_fp = np.load(root+'GroupFirstSub_'+fp_ext+'.npy')
SubhaloPos_fp = np.load(root+'SubhaloPos_'+fp_ext+'.npy')/1000.
SubhaloMassType_fp = np.load(root+'SubhaloMassType_'+fp_ext+'.npy')
SubhaloStellarPhotometrics_fp = np.load(root+'SubhaloStellarPhotometrics_'+fp_ext+'.npy')
SubhaloSFR_fp = np.load(root+'SubhaloSFR_'+fp_ext+'.npy')
SubhalosSFR_fp = SubhaloSFR_fp/(SubhaloMassType_fp[:,4]*1.e10)
SubhaloUniqueParent_fp = np.ones(len(SubhaloSFR_fp),dtype=int)
SubhaloUniqueParent_fp[GroupFirstSub_fp[:1000000]] = 0
PeakStats_fp = np.load(root+'peakstats_fp.npy')
SubhaloZpeak_fp = PeakStats_fp[:,0]
SubhaloZform_fp = PeakStats_fp[:,1]
SubhaloGrowth_fp = 1./PeakStats_fp[:,6]
u_fp = np.load(light_dir+"sdss_u_99.npy")
g_fp = np.load(light_dir+"sdss_g_99.npy")
r_fp = np.load(light_dir+"sdss_r_99.npy")
B_fp = SubhaloStellarPhotometrics_fp[:,1]
V_fp = SubhaloStellarPhotometrics_fp[:,2]
gmr_fp = g_fp-r_fp
umg_fp = u_fp-g_fp
BmV_fp = B_fp-V_fp
umr_fp = u_fp-r_fp

# Total number of subhalos
N_sub_dm = SubhaloPos_dm.shape[0]
N_sub_fp = SubhaloPos_fp.shape[0]

# Mass proxy for subhalos
if proxy == 'mass': m_proxy = SubhaloMassType_dm[:,1]*1.e10
if proxy == 'vpeak': m_proxy = SubhaloVpeak_dm
if proxy == 'vmax': m_proxy = SubhaloVmax_dm
if proxy == 'halfmass': m_proxy = SubhaloMassInHalfRad_dm
if proxy == 'mmax': m_proxy = SubhaloMassInMaxRad_dm
if proxy == 'vdisp': m_proxy = SubhaloVelDisp_dm
if proxy == 's2r': m_proxy = SubhaloVelDisp_dm**2*SubhaloHalfmassRad_dm#SubhaloVmaxRad_dm
if proxy == 'vinfall': m_proxy = SubhaloVinfall_dm
if proxy == 'minfall': m_proxy = SubhaloMinfall_dm*1.e10
if proxy == 'mpeak': m_proxy = SubhaloMpeak_dm*1.e10
if proxy == 'vrelax': m_proxy = SubhaloVrelax_dm

i_before = np.arange(len(SubhaloSFR_fp),dtype=int)


if isolate_low_mass:
    M_lim = 5. #10^10 Msun/h
    prop_choice = (SubhaloMassType_fp[:,4] < M_lim) & (SubhaloMassType_fp[:,4] > 1.);
    ind_prop = i_before[prop_choice]
    fp_inds, c1, c2_prop = np.intersect1d(ind_prop,fp_inds,return_indices=True)
    dmo_inds = dmo_inds[c2_prop]

if sec_prop != 'full':
    
    if sec_prop == 'mstar':
        prop_fp = (SubhaloMassType_fp[:,4])
        M_lim = 5.#3. #10^10 Msun/h
        if i_bin == 0: text = "low stellar mass"; prop_choice = (prop_fp < M_lim) & (prop_fp > 1.);
        if i_bin == 1: text = "high stellar mass"; prop_choice = (prop_fp > M_lim);
    elif sec_prop == 'g-r':
        prop_fp = gmr_fp
        gr_lim = .6#.73#0.7
        if i_bin == 0: text = "blue"; prop_choice = (prop_fp < gr_lim)
        if i_bin == 1: text = "red"; prop_choice = (prop_fp > gr_lim)
    elif sec_prop == "parent":
        prop_fp = SubhaloUniqueParent_fp
        if i_bin == 0: text = "satellites"; prop_choice = prop_fp > 0;
        if i_bin == 1: text = "centrals"; prop_choice = prop_fp == 0;
    elif sec_prop == 'sSFR':
        prop_fp = SubhalosSFR_fp
        sfr_lim = -10.4#-10.75
        if i_bin == 0: text = "SFGs"; prop_choice = np.log10(prop_fp) > sfr_lim
        if i_bin == 1: text = "quiescent"; prop_choice = np.log10(prop_fp) < sfr_lim
    elif sec_prop == 'form':
        mode = 'growth'#'peak'#'form'#'growth'
        if mode == 'peak':
            prop_fp = SubhaloZpeak_fp
            form_lim = 0.01
        elif mode == 'form':
            prop_fp = SubhaloZform_fp
            form_lim = 1.35
        elif mode == 'growth':
            prop_fp = SubhaloGrowth_fp
            form_lim = 0.9
                    
        if i_bin == 0: text = "late-forming"; prop_choice = (prop_fp) <= form_lim
        if i_bin == 1: text = "early-forming"; prop_choice = (prop_fp) > form_lim

    ind_prop = i_before[prop_choice]
    fp_inds, c1, c2_prop = np.intersect1d(ind_prop,fp_inds,return_indices=True)
    dmo_inds = dmo_inds[c2_prop]

# mass proxy sorting
chosen_inds_dm = (np.argsort(m_proxy)[::-1])[first_gal:N_gal_choice]

# take only those that have matches
if want_matched == '_matched':
    # note that instersect does a sorting!
    chosen_inds_dm, c1_dm, c2_dm = np.intersect1d(chosen_inds_dm,dmo_inds,return_indices=True)

# number of dm subhalos of interest
n_dm = len(chosen_inds_dm)
print("Matched gals DM = ",len(chosen_inds_dm))

# get the positions and masses of the chosen subhalos
m_proxy_sorted = m_proxy[chosen_inds_dm]
pos_dm_msorted = SubhaloPos_dm[chosen_inds_dm]
m_dm_msorted = SubhaloMassType_dm[chosen_inds_dm,1]*1.e10

# luminosity proxy sorting
lum_proxy = SubhaloMassType_fp[:,4]*1.e10
chosen_inds_fp = (np.argsort(lum_proxy)[::-1])[first_gal:N_gal_choice]

# take only those that have matches
if want_matched == '_matched':
    # note that instersect does a sorting!
    chosen_inds_fp, c1_fp, c2_fp = np.intersect1d(chosen_inds_fp,fp_inds,return_indices=True)
    #dmo_inds[c2_fp] to access equiv in dm

# number of fp subhalos of interest
n_fp = len(chosen_inds_fp)
print("Matched gals FP = ",n_fp)

# smallest common number of gals
n_gal = np.min([n_fp,n_dm])

# conserve the number of galaxies plotted
if want_matched == '_matched':
    i_sort_cons_dm = (np.argsort(m_proxy_sorted)[::-1])[:n_gal]
    chosen_inds_dm = chosen_inds_dm[i_sort_cons_dm]
    pos_dm_msorted = SubhaloPos_dm[chosen_inds_dm]
    m_proxy_sorted = m_proxy[chosen_inds_dm]
    m_dm_msorted = SubhaloMassType_dm[chosen_inds_dm,1]*1.e10
    
    i_sort_cons_fp = (np.argsort(lum_proxy[chosen_inds_fp])[::-1])[:n_gal]
    chosen_inds_fp = chosen_inds_fp[i_sort_cons_fp]

# sorting mass and position and the extra properties
lum_proxy_sorted = lum_proxy[chosen_inds_fp]
if want_scatter:
    dm = 0.125
    lum_proxy = 10.**np.random.normal(np.log10(lum_proxy),dm)
    chosen_inds_fp = (np.argsort(lum_proxy)[::-1])[first_gal:N_gal_choice]
    lum_proxy_sorted = lum_proxy[chosen_inds_fp]
    '''
    #OPTION 2
    min_mass = np.log10(np.min(lum_proxy_sorted))
    max_mass = np.log10(np.max(lum_proxy_sorted))
    dm = 0.125 #dex
    bins = 10.**np.arange(min_mass,max_mass,dm)
    print(len(bins))
    for i in range(len(bins)-1):
        bin_choice = (lum_proxy_sorted > bins[i]) & (lum_proxy_sorted <= bins[i+1])
        idx = chosen_inds_fp[bin_choice]
        np.random.shuffle(idx)
        chosen_inds_fp[bin_choice] = idx
    '''
    
pos_fp_lsorted = SubhaloPos_fp[chosen_inds_fp]
SFR_fp_lsorted = SubhaloSFR_fp[chosen_inds_fp]
gmr_fp_lsorted = gmr_fp[chosen_inds_fp]
umg_fp_lsorted = umg_fp[chosen_inds_fp]
umr_fp_lsorted = umr_fp[chosen_inds_fp]
BmV_fp_lsorted = BmV_fp[chosen_inds_fp]
mstar_fp_lsorted = (SubhaloMassType_fp[:,4]*1.e10)[chosen_inds_fp]
parent_fp_lsorted = SubhaloUniqueParent_fp[chosen_inds_fp]
sSFR_fp_lsorted = SFR_fp_lsorted/mstar_fp_lsorted
zpeak_fp_lsorted = SubhaloZpeak_fp[chosen_inds_fp]

if sec_prop == 'full':
    M_hi = format(np.log10(lum_proxy_sorted[0]),'.1f')
    M_lo = format(np.log10(lum_proxy_sorted[-1]),'.1f')
    text = r'$\log M_\ast = $'+ (str(M_lo)+"-"+str(M_hi))


# END HERE

'''
nbins = 31
m_bins = np.logspace(np.log10(70.),np.log10(2000.),nbins)
bin_cents = .5*(m_bins[1:]+m_bins[:-1])
r_median,r1_low,r1_high,r2_low,r2_high = compute_binned_stats(m_proxy_sorted,m_dm_msorted,m_bins)
np.save("data/scat_med"+want_matched+"_"+proxy+"_"+dm_ext+".npy",r_median)
np.save("data/scat_low"+want_matched+"_"+proxy+"_"+dm_ext+".npy",r1_low)
np.save("data/scat_high"+want_matched+"_"+proxy+"_"+dm_ext+".npy",r1_high)
np.save("data/scat_bins_"+proxy+".npy",bin_cents)
quit()
'''

'''
nbins = 31
bin_edges = np.logspace(np.log10(70.),np.log10(2000.),nbins)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
hist, edges = np.histogram(m_proxy_sorted,bins=bin_edges)
hist_c, edges = np.histogram(m_proxy_sorted[parent_fp_lsorted == 0],bins=bin_edges)
hist_s, edges = np.histogram(m_proxy_sorted[parent_fp_lsorted > 0],bins=bin_edges)
hist_all, edges = np.histogram(m_proxy,bins=bin_edges)

bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
hist_rat = hist/hist_all
hist_rat_c = hist_c/hist_all
hist_rat_s = hist_s/hist_all
np.save("data/hist"+want_matched+"_"+proxy+"_"+dm_ext+".npy",hist_rat)
np.save("data/hist_cen"+want_matched+"_"+proxy+"_"+dm_ext+".npy",hist_rat_c)
np.save("data/hist_sat"+want_matched+"_"+proxy+"_"+dm_ext+".npy",hist_rat_s)
np.save("data/hist_bins_"+proxy+".npy",bin_cen)
quit()
'''

'''
s = 0.01
plt.figure()
plt.scatter(mstar_fp_lsorted,np.log10(sSFR_fp_lsorted),s=s)
plt.xscale('log')

plt.figure()
plt.scatter(mstar_fp_lsorted,gmr_fp_lsorted,s=s)
plt.xscale('log')
plt.show()
quit()
'''

# Use this if you want to compare with shuff with start here and end here
#pos_dm_msorted = np.load("../Lensing/data_2dhod_peak/m200m_shuff_gals.npy")
#pos_fp_lsorted = np.load("../Lensing/data_2dhod_pos/true_gals.npy")

# Corrfunc stuff starts here
xyz_DM = pos_dm_msorted
xyz_FP = pos_fp_lsorted
print("FP gals = ",xyz_FP.shape[0],text)
print("DM gals = ",xyz_DM.shape[0],text)

n_gal = xyz_FP.shape[0]

N_dim = 3
N_bin = 16
Rat_SHAM = np.zeros((N_bin-1,N_dim**3))
for i_x in range(N_dim):
    for i_y in range(N_dim):
        for i_z in range(N_dim):
            xyz_FP_jack = xyz_FP.copy()
            xyz_DM_jack = xyz_DM.copy()
            
            xyz = np.array([i_x,i_y,i_z],dtype=int)
            size = Lbox/N_dim

            bool_arr = np.prod((xyz == (xyz_FP/size).astype(int)),axis=1).astype(bool)
            xyz_FP_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_FP_jack = xyz_FP_jack[np.sum(xyz_FP_jack,axis=1)!=0.]

            bool_arr = np.prod((xyz == (xyz_DM/size).astype(int)),axis=1).astype(bool)
            xyz_DM_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_DM_jack = xyz_DM_jack[np.sum(xyz_DM_jack,axis=1)!=0.]
            
            x_DM_jack = xyz_DM_jack[:,0]
            y_DM_jack = xyz_DM_jack[:,1]
            z_DM_jack = xyz_DM_jack[:,2]
            x_FP_jack = xyz_FP_jack[:,0]
            y_FP_jack = xyz_FP_jack[:,1]
            z_FP_jack = xyz_FP_jack[:,2]

            bins = np.logspace(-0.7,1.5,N_bin)
            bin_centers = (bins[:-1] + bins[1:])/2.

            results_DM = Corrfunc.theory.xi(X=x_DM_jack,Y=y_DM_jack,Z=z_DM_jack,
                                            boxsize=Lbox,nthreads=16,
                                            binfile=bins)
            
            results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,
                                            boxsize=Lbox, nthreads=16,
                                            binfile=bins)
            
            Corr_FP = results_FP['xi']
            Corr_DM = results_DM['xi']
            Rat_SHAM[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_DM/Corr_FP

Rat_SHAM_mean = np.mean(Rat_SHAM,axis=1)
Rat_SHAM_err = np.sqrt(N_dim**3-1)*np.std(Rat_SHAM,axis=1)

if N_div == 1:
    if dm_ext == 'dm': dm_ext = ''
    np.save("data/SHAM"+want_matched+"_ratio"+dm_ext+"_"+proxy+".npy",Rat_SHAM_mean)
    np.save("data/SHAM"+want_matched+"_ratio"+dm_ext+"_"+proxy+"_error.npy",Rat_SHAM_err)
    np.save("data/bin_centers",bin_centers)

if N_div >= 2:
    np.save("data_split/SHAM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+".npy",Rat_SHAM_mean)
    np.save("data_split/SHAM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+"_error.npy",Rat_SHAM_err)
    np.save("data_split/bin_centers",bin_centers)

plt.figure(figsize=(12,8))
plt.title("Auto-correlation ratio for "+str(n_gal)+" galaxies")
plt.errorbar(bin_centers,Rat_SHAM_mean,yerr=Rat_SHAM_err,linewidth=2.,label='SHAM w/ '+proxy+' (DM)/ "true" (FP)')
plt.plot(bin_centers,np.ones(len(Corr_DM)),'k--',linewidth=2.)
plt.xscale('log')
plt.ylim([0,1.3])
plt.xlim([.7,15])
plt.legend()
#plt.text(0.8, .2, text, dict(size=15))
plt.ylabel(r'Ratio')
plt.xlabel(r'$r$ [Mpc/h]')

if N_div == 1:
    plt.savefig("figs/SHAM"+want_matched+"_ratio_"+proxy+"_"+str(N_gal_choice)+"_"+sec_prop+".png")
else:
    plt.savefig("figs/SHAM"+want_matched+"_ratio_"+proxy+"_"+str(N_gal_choice)+"_"+sec_prop+"_"+str(i_bin)+".png")
#plt.show()
plt.close()

bins = 31
'''
plt.figure(1,figsize=(9,7))
plt.title("sSFR")
log_min = -12
log_max = -7#0.9
bin_edges = np.linspace(log_min,log_max,bins)
hist, edges = np.histogram(SFR_fp_lsorted/mstar_fp_lsorted,bins=10**bin_edges,normed=True)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.step(10**bin_cen,hist,where='mid',linewidth=2.)
plt.xscale('log')
plt.xlim([edges[0],edges[-1]])
'''
'''
plt.figure(2,figsize=(9,7))
bin_edges = np.linspace(-.1,0.9,bins+20)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.title("g-r")
hist1, edges = np.histogram(gmr_fp_lsorted[parent_fp_lsorted==0],bins=bin_edges)
hist2, edges = np.histogram(gmr_fp_lsorted[parent_fp_lsorted>0],bins=bin_edges)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.step(bin_cen,hist1,where='mid',linewidth=2.,label='centrals')
plt.step(bin_cen,hist2,where='mid',linewidth=2.,label='subhalos')
plt.legend()

plt.figure(3,figsize=(9,7))
bin_edges = np.linspace(-.1,0.9,bins+20)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.title("sSFR")
hist1, edges = np.histogram(sSFR_fp_lsorted[parent_fp_lsorted==0],bins=bin_edges)
hist2, edges = np.histogram(sSFR_fp_lsorted[parent_fp_lsorted>0],bins=bin_edges)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.step(bin_cen,hist1,where='mid',linewidth=2.,label='centrals')
plt.step(bin_cen,hist2,where='mid',linewidth=2.,label='subhalos')
plt.legend()
plt.show()
'''
'''
plt.figure(3,figsize=(9,7))
log_min = 8.
log_max = 12.8
bin_edges = np.linspace(log_min,log_max,bins)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
hist, edges = np.histogram(mstar_fp_lsorted,bins=bin_edges,normed=True)
plt.title("mstar")
plt.step(10**bin_cen,hist,where='mid',linewidth=2.)
plt.xscale('log')
plt.show()
quit()
'''

s = 0.01
'''
plt.figure(1)
plt.scatter(BmV_fp_sham,umr_fp_sham,s=s)
plt.ylabel("u-r")
plt.xlabel("B-V")
'''
'''
plt.figure(2)
plt.scatter(umg_fp_sham,gmr_fp_sham,s=s)
plt.xlabel("u-g")
plt.ylabel("g-r")

plt.figure(3)
plt.scatter(mstar_fp_sham,gmr_fp_sham,s=s)
plt.xscale('log')
#plt.ylim([10**(-3),10**1.9])

plt.figure(4)
plt.scatter(gmr_fp_sham,par_fp_sham,s=s,alpha=0.5)
#plt.yscale('log')
#plt.ylim([10**(-3),10**1.9])

plt.figure(5)
plt.scatter(mstar_fp_sham,par_fp_sham,s=s,alpha=0.5)
plt.xscale('log')
#plt.show()
'''
'''
plt.figure(6,figsize=(9,7))
bin_edges = np.linspace(-.1,0.9,bins+20)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.title("g-r")
hist1, edges = np.histogram(gmr_fp_sham[par_fp_sham==0],bins=bin_edges)#,normed=True)
hist2, edges = np.histogram(gmr_fp_sham[par_fp_sham>0],bins=bin_edges)#,normed=True)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.step(bin_cen,hist1,where='mid',linewidth=2.,label='centrals')
plt.step(bin_cen,hist2,where='mid',linewidth=2.,label='subhalos')
plt.legend()
'''
'''
plt.figure(7,figsize=(9,7))
log_min = -12
log_max = -7
bin_edges = np.linspace(log_min,log_max,bins+20)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.title("sSFR")
hist1, edges = np.histogram(sSFR_fp_sham[par_fp_sham==0],bins=10**bin_edges)#,normed=True)
hist2, edges = np.histogram(sSFR_fp_sham[par_fp_sham>0],bins=10**bin_edges)#,normed=True)
bin_cen = .5*(bin_edges[1:]+bin_edges[:-1])
plt.step(10**bin_cen,hist1,where='mid',linewidth=2.,label='centrals')
plt.step(10**bin_cen,hist2,where='mid',linewidth=2.,label='subhalos')
plt.xscale('log')
plt.legend()
plt.show()
quit()
'''
