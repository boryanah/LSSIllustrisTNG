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
proxy = 'vpeak'#'vrelax'
want_matched = ''#'_matched'
want_total_pop = 1
first_gal = 0
N_gal = 12000#60000#45000
N_div = 2 # subpopulations
sec_prop = sys.argv[2] # splitting by

root = '/home/boryanah/lars/data_2500/'
light_dir = "/mnt/gosling1/boryanah/Illustris_light/TNG300/"
res = 2500
Lbox = 205.
i_bin = MPI.COMM_WORLD.Get_rank()

# START HERE

# matching subhalos in fp and dm
fp_dmo_inds = np.load(root+'sub_match_fp_dmo_'+str(res)+'.npy')
fp_inds = fp_dmo_inds[0]
dmo_inds = fp_dmo_inds[1]


# load the DM data; set this to fp to test SHAM within hydro
dm_ext = 'dm'#'dm'#'fp'
SubhaloVmax_dm = np.load(root+'SubhaloVmax_'+dm_ext+'.npy')
new_type = '_new'# newest is cut at vmax = 100; new at vmax = 50; null is old format
SubhaloVpeak_dm = np.load(root+'SubhaloVpeak'+new_type+'_'+dm_ext+'.npy')
#SubhaloVpeak_dm = np.load(root+'SubhaloVpeak_new_'+dm_ext+'.npy')
SubhaloVpeak_dm[SubhaloVpeak_dm==0.] = SubhaloVmax_dm[SubhaloVpeak_dm==0.]
SubhaloVinfall_dm = np.load(root+'SubhaloVinfall'+new_type+'_'+dm_ext+'.npy')
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
SubhaloMinfall_dm = np.load(root+'SubhaloMinfall'+new_type+'_'+dm_ext+'.npy')
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

# FP population
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

# mass proxy sorting
chosen_inds_dm = (np.argsort(m_proxy)[::-1])[first_gal:N_gal]


# FP population
# Mass proxy for subhalos
lum_proxy = SubhaloMassType_fp[:,4]
# luminosity proxy sorting
chosen_inds_fp = (np.argsort(lum_proxy)[::-1])[first_gal:N_gal]

if want_total_pop:
    pos_total_dm = SubhaloPos_dm[chosen_inds_dm]
    pos_total_fp = SubhaloPos_fp[chosen_inds_fp]

print("_____________________________________________________")
print("mass of last object = ",lum_proxy[chosen_inds_fp[-1]])

if sec_prop == 'mstar':
    prop_fp = SubhaloMassType_fp[:,4]
    lim = 3.#5. #10^10 Msun/h
    if i_bin == 0: name = 'low-mass'
    if i_bin == 1: name = 'high-mass'
elif sec_prop == 'g-r':
    prop_fp = gmr_fp
    lim = .7#.6#.73#0.7
    if i_bin == 0: name = 'blue'
    if i_bin == 1: name = 'red'
elif sec_prop == "parent":
    prop_fp = SubhaloUniqueParent_fp
    lim = 1
    if i_bin == 0: name = 'centrals'
    if i_bin == 1: name = 'satellites'
elif sec_prop == 'sSFR':
    prop_fp = SubhalosSFR_fp
    lim = -10.#-10.75#-10.4#-10.75
    if i_bin == 0: name = 'quiescent'
    if i_bin == 1: name = 'star-forming'
elif sec_prop == 'form':
    mode = 'growth'#'growth'#'peak'#'form'#'growth'
    if mode == 'peak':
        prop_fp = SubhaloZpeak_fp
        lim = 0.01
    elif mode == 'form':
        prop_fp = SubhaloZform_fp
        lim = 1.35
    elif mode == 'growth':
        prop_fp = SubhaloGrowth_fp
        lim = 0.9
    if i_bin == 0: name = 'late-forming'
    if i_bin == 1: name = 'early-forming'
    
prop_chosen = prop_fp[chosen_inds_fp]
if sec_prop in ['sSFR']:
    lim = np.percentile(prop_chosen,60.)
    
elif sec_prop in ['form','mstar','g-r']:
    lim = np.percentile(prop_chosen,50.)
    

# TESTING
'''
if i_bin == 0:
    chosen_inds_dm = chosen_inds_dm[prop_chosen < lim]
    chosen_inds_fp = chosen_inds_fp[prop_chosen < lim]
if i_bin == 1:
    chosen_inds_dm = chosen_inds_dm[prop_chosen >= lim]
    chosen_inds_fp = chosen_inds_fp[prop_chosen >= lim]

pos_dm_msorted = SubhaloPos_dm[chosen_inds_dm]
pos_fp_lsorted = SubhaloPos_fp[chosen_inds_fp]
'''
if i_bin == 1:
    choice = prop_chosen >= lim
    also_cents = choice & (SubhaloUniqueParent_fp[chosen_inds_fp] == 0)
    print("percentage centrals = %.1f"%(np.sum(also_cents)*100./np.sum(choice)))
    pos_fp_lsorted = (SubhaloPos_fp[chosen_inds_fp])[choice]
    pos_dm_msorted = (SubhaloPos_dm[chosen_inds_dm])[choice]
    #pos_dm_msorted = (SubhaloPos_dm[chosen_inds_dm])[:np.sum(choice)]
if i_bin == 0:
    choice = prop_chosen < lim
    also_cents = choice & (SubhaloUniqueParent_fp[chosen_inds_fp] == 0)
    print("percentage centrals = %.1f"%(np.sum(also_cents)*100./np.sum(choice)))
    pos_fp_lsorted = (SubhaloPos_fp[chosen_inds_fp])[choice]
    pos_dm_msorted = (SubhaloPos_dm[chosen_inds_dm])[choice]
    #pos_dm_msorted = (SubhaloPos_dm[chosen_inds_dm])[:np.sum(choice)]

'''
#plt.scatter((SubhaloMassType_fp[chosen_inds_fp,4])[choice],(SubhalosSFR_fp[chosen_inds_fp])[choice],s=0.5)
plt.scatter((SubhaloMassType_fp[chosen_inds_fp,4])[choice],(gmr_fp[chosen_inds_fp])[choice],s=0.5)
plt.xscale('log')
#plt.yscale('log')
plt.show()
quit()
'''
# END HERE
# Use this if you want to compare with shuff with start here and end here
#pos_dm_msorted = np.load("../Lensing/data_2dhod_peak/m200m_shuff_gals.npy")
#pos_fp_lsorted = np.load("../Lensing/data_2dhod_pos/true_gals.npy")

# Corrfunc stuff starts here
xyz_DM = pos_dm_msorted
xyz_FP = pos_fp_lsorted
print("divisor = ",lim)
print("FP gals = ",xyz_FP.shape[0],name)
print("DM gals = ",xyz_DM.shape[0],name)

n_gal = xyz_FP.shape[0]


N_bin = 16
bins = np.logspace(-0.7,1.5,N_bin)
bin_centers = (bins[:-1] + bins[1:])/2.

N_dim = 3
Rat_SHAM = np.zeros((N_bin-1,N_dim**3))
Rat_DM = np.zeros((N_bin-1,N_dim**3))
Rat_FP = np.zeros((N_bin-1,N_dim**3))
for i_x in range(N_dim):
    for i_y in range(N_dim):
        for i_z in range(N_dim):
            xyz_FP_jack = xyz_FP.copy()
            xyz_DM_jack = xyz_DM.copy()

            if want_total_pop:
                pos_DM_jack = pos_total_dm.copy()
                pos_FP_jack = pos_total_fp.copy()
            
            xyz = np.array([i_x,i_y,i_z],dtype=int)
            size = Lbox/N_dim

            bool_arr = np.prod((xyz == (xyz_FP/size).astype(int)),axis=1).astype(bool)
            xyz_FP_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_FP_jack = xyz_FP_jack[np.sum(xyz_FP_jack,axis=1)!=0.]

            bool_arr = np.prod((xyz == (xyz_DM/size).astype(int)),axis=1).astype(bool)
            xyz_DM_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_DM_jack = xyz_DM_jack[np.sum(xyz_DM_jack,axis=1)!=0.]

            if want_total_pop:
                bool_arr = np.prod((xyz == (pos_total_dm/size).astype(int)),axis=1).astype(bool)
                pos_DM_jack[bool_arr] = np.array([0.,0.,0.])
                pos_DM_jack = pos_DM_jack[np.sum(pos_DM_jack,axis=1)!=0.]
                
                bool_arr = np.prod((xyz == (pos_total_fp/size).astype(int)),axis=1).astype(bool)
                pos_FP_jack[bool_arr] = np.array([0.,0.,0.])
                pos_FP_jack = pos_FP_jack[np.sum(pos_FP_jack,axis=1)!=0.]
            
            x_DM_jack = xyz_DM_jack[:,0]
            y_DM_jack = xyz_DM_jack[:,1]
            z_DM_jack = xyz_DM_jack[:,2]
            x_FP_jack = xyz_FP_jack[:,0]
            y_FP_jack = xyz_FP_jack[:,1]
            z_FP_jack = xyz_FP_jack[:,2]

            results_DM = Corrfunc.theory.xi(X=x_DM_jack,Y=y_DM_jack,Z=z_DM_jack,
                                            boxsize=Lbox,nthreads=16,
                                            binfile=bins)
            
            results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,
                                            boxsize=Lbox, nthreads=16,
                                            binfile=bins)

            Corr_FP = results_FP['xi']
            Corr_DM = results_DM['xi']
            Rat_SHAM[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_DM/Corr_FP

            if want_total_pop:
                Corr_tot_FP = Corrfunc.theory.xi(X=pos_FP_jack[:,0],Y=pos_FP_jack[:,1],Z=pos_FP_jack[:,2],
                                                boxsize=Lbox, nthreads=16, binfile=bins)['xi']
                Corr_tot_DM = Corrfunc.theory.xi(X=pos_DM_jack[:,0],Y=pos_DM_jack[:,1],Z=pos_DM_jack[:,2],
                                                boxsize=Lbox, nthreads=16, binfile=bins)['xi']
                Rat_DM[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_DM/Corr_tot_DM
                Rat_FP[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_FP/Corr_tot_FP

Rat_SHAM_mean = np.mean(Rat_SHAM,axis=1)
Rat_SHAM_err = np.sqrt(N_dim**3-1)*np.std(Rat_SHAM,axis=1)
Rat_DM_mean = np.mean(Rat_DM,axis=1)
Rat_DM_err = np.sqrt(N_dim**3-1)*np.std(Rat_DM,axis=1)
Rat_FP_mean = np.mean(Rat_FP,axis=1)
Rat_FP_err = np.sqrt(N_dim**3-1)*np.std(Rat_FP,axis=1)

if N_div == 1:
    dm_ext = ''
    np.save("data/SHAM"+want_matched+"_ratio"+dm_ext+"_"+proxy+".npy",Rat_SHAM_mean)
    np.save("data/SHAM"+want_matched+"_ratio"+dm_ext+"_"+proxy+"_error.npy",Rat_SHAM_err)
    np.save("data/bin_centers",bin_centers)

if N_div >= 2:
    np.save("data_split/SHAM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+".npy",Rat_SHAM_mean)
    np.save("data_split/SHAM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+"_error.npy",Rat_SHAM_err)
    np.save("data_split/bin_centers",bin_centers)

    if want_total_pop:
        np.save("data_split/DM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+".npy",Rat_DM_mean)
        np.save("data_split/DM"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+"_error.npy",Rat_DM_err)
        np.save("data_split/FP"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+".npy",Rat_FP_mean)
        np.save("data_split/FP"+want_matched+"_ratio_"+proxy+"_"+sec_prop+"_"+str(i_bin)+"_error.npy",Rat_FP_err)
        
    
plt.figure(figsize=(12,8))
plt.title("Auto-correlation ratio for "+str(n_gal)+" galaxies")
plt.errorbar(bin_centers,Rat_SHAM_mean,yerr=Rat_SHAM_err,linewidth=2.,label='SHAM w/ '+proxy+' (DM)/ "true" (FP)')
plt.plot(bin_centers,np.ones(len(Corr_DM)),'k--',linewidth=2.)
plt.xscale('log')
plt.ylim([0,1.3])
plt.xlim([.7,15])
plt.legend()
plt.ylabel(r'Ratio')
plt.xlabel(r'$r$ [Mpc/h]')

if N_div == 1:
    plt.savefig("figs/SHAM"+want_matched+"_ratio_"+proxy+"_"+str(N_gal)+"_"+sec_prop+".png")
else:
    plt.savefig("figs/SHAM"+want_matched+"_ratio_"+proxy+"_"+str(N_gal)+"_"+sec_prop+"_"+str(i_bin)+".png")
#plt.show()
plt.close()
