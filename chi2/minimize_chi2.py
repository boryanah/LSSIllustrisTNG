import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Corrfunc
import plotparams
plotparams.default()
import h5py
from scipy.optimize import minimize
import sys

# info about simulation
reshift = 2.22044604925e-16
h = 0.6774
res = 2500
Lbox = 205.
root = '/home/boryanah/lars/data_2500/'
np.random.seed(300)

# parameter choices
want_fp_pos = 0 # do you want to take the galaxy positions from the FP sim
proxy = 'm200m'
opt = sys.argv[1]#'partial_fenv'
percent_thresh = 5./100.
gal_cut = 10000
N_exc_mm = 100#0# for bias and corr plot #100 # for rest of paper
# r-bins to do the matching
i_r_start = -15
#i_r_end = 0#-1
want_minimize = 0 # we are giving the r coefficient as sys.argv
# where to save images
figname = 'Corr_'+proxy+'.png'

# loading data
GroupFirstSub_dm = np.load(root+'GroupFirstSub_dm.npy')
GroupNsubs_dm = np.load(root+'GroupNsubs_dm.npy')
GroupPos_dm = np.load(root+'GroupPos_dm.npy')/1000.
Group_M_Mean200_dm = np.load(root+'Group_M_Mean200_dm.npy')*1.e10 # Msun/h
Group_Vmax_dm = np.load(root+'Group_Vmax_dm.npy')
SubhaloGrNr_dm = np.load(root+'SubhaloGrNr_dm.npy')
SubhaloPos_dm = np.load(root+'SubhaloPos_dm.npy')/1000.
SubhaloMassType_dm = np.load(root+'SubhaloMassType_dm.npy')
SubhaloVelDisp_dm = np.load(root+'SubhaloVelDisp_dm.npy')
SubhaloHalfmassRad_dm = np.load(root+'SubhaloHalfmassRad_dm.npy')

# Buba-made fields
min_pot_halo_dm = -np.load(root+'GroupMinPotential_dm.npy')
tot_pot_halo_dm = -np.load(root+'GroupTotPotential_dm.npy')
env_mass_halo_dm = np.load(root+'env_mass_dm.npy')
env_halo_dm = np.load(root+'env_dm.npy')
GroupVelAnis_dm_part = np.load(root+'GroupVelAnis_dm_300000.npy')
SubhaloVpeak_dm = np.load(root+'SubhaloVpeak_dm.npy')
if opt == "partial_env_cw":
    # Density in Illustris
    fdir = '/mnt/store1/boryanah/IllustrisTNG/CosmicWeb'
    filename = 'WEB_CIC_256_DM_TNG300-2.hdf5'
    f = h5py.File(fdir+filename, 'r')
    d_smooth = f['density_smooth'][:,:,:] 

    # finding who belongs where in the cosmic web
    N_dim = 256
    gr_size = Lbox/N_dim
    halo_x = GroupPos_dm[:,0];halo_y = GroupPos_dm[:,1];halo_z = GroupPos_dm[:,2]
    i_cw = (halo_x/gr_size).astype(int)
    j_cw = (halo_y/gr_size).astype(int)
    k_cw = (halo_z/gr_size).astype(int)

    # Environment definition
    env_cw_dm = d_smooth[i_cw,j_cw,k_cw]

# loading data
GroupFirstSub_fp = np.load(root+'GroupFirstSub_fp.npy')
GroupNsubs_fp = np.load(root+'GroupNsubs_fp.npy')
SubhaloGrNr_fp = np.load(root+'SubhaloGrNr_fp.npy')
SubhaloPos_fp = np.load(root+'SubhaloPos_fp.npy')/1000.
SubhaloMassType_fp = np.load(root+'SubhaloMassType_fp.npy')
SubhaloLenType_fp = np.load(root+'SubhaloLenType_fp.npy')
# wicked testing
Group_M_Crit200_fp = np.load(root+'Group_M_Crit200_fp.npy')


# numbers of halos
N_halo_dm = len(Group_M_Mean200_dm)
N_halo_fp = len(GroupFirstSub_fp)

# number of subs
N_sub_dm = len(SubhaloGrNr_dm)    
N_sub_fp = len(SubhaloGrNr_fp)

# for getting the group vani and group s2r from the first central 
unique_sub_grnr, firsts = np.unique(SubhaloGrNr_dm,return_index=True)
GroupVelAnis_dm = np.zeros(N_halo_dm)
GroupVelAnis_dm[:] = Group_Vmax_dm[:]
GroupVelAnis_dm[:300000] = GroupVelAnis_dm_part
Group_sigmasq_R = np.zeros(N_halo_dm)
Group_sigmasq_R[:] = Group_Vmax_dm[:]
Group_sigmasq_R[unique_sub_grnr] = (SubhaloVelDisp_dm**2*SubhaloHalfmassRad_dm)[firsts]

# indices of subhalos
inds_sub_dm = np.arange(N_sub_dm)
inds_sub_fp = np.arange(N_sub_fp)


# mass and number of star particles for each subhalo in fp
m_sub_stars = SubhaloMassType_fp[:,4]*1.e10 #M_sun/h
l_sub_stars = SubhaloLenType_fp[:,4]
m_sub_dm = SubhaloMassType_dm[:,1]*1.e10 #M_sun/h

# load the halo matchings
fp_inds_halo = np.load(root+'fp_dmo_inds_halo_'+str(res)+'.npy')[:,0]
dmo_inds_halo = np.load(root+'fp_dmo_inds_halo_'+str(res)+'.npy')[:,1]

# what is the index of the largest subhalo belonging to the first fp halo
mag_exc = 2063

# which doesn't have a direct counterpart in DM
sub_exc = GroupFirstSub_fp[fp_inds_halo[-mag_exc:]]
N_exc = GroupNsubs_fp[fp_inds_halo[-mag_exc:]]
for m in range(mag_exc):
    SubhaloPos_fp[sub_exc[m]:N_exc[m]+sub_exc[m]] = np.array([0.,0.,0.])

#inds_sub_this_fp = inds_sub_fp[:]
inds_sub_fp = np.arange(N_sub_fp)
inds_halo_fp = np.arange(N_halo_fp)

# removing the excluded unmatched
inds_sub_fp = inds_sub_fp[l_sub_stars > gal_cut]
m_sub_gals = m_sub_stars[l_sub_stars > gal_cut]
parent_gal_fp = SubhaloGrNr_fp[l_sub_stars > gal_cut]
pos_gal_fp = SubhaloPos_fp[l_sub_stars > gal_cut]

# remove those outside
inds_sub_fp = inds_sub_fp[np.sum(pos_gal_fp,axis=1)!=0.]
m_sub_gals = m_sub_gals[np.sum(pos_gal_fp,axis=1)!=0.]
parent_gal_fp = parent_gal_fp[np.sum(pos_gal_fp,axis=1)!=0.]
pos_gal_fp = pos_gal_fp[np.sum(pos_gal_fp,axis=1)!=0.]

# galaxy numbers
N_gal = len(inds_sub_fp)
n_gal = N_gal/Lbox**3
print('n_gal = ',n_gal)
print('N_gal = ',N_gal)
            
# remove the double matched ones
fp_inds_halo = fp_inds_halo[:-mag_exc]
dmo_inds_halo = dmo_inds_halo[:-mag_exc]
i_sort = np.argsort(Group_M_Mean200_dm[dmo_inds_halo])[::-1]
N_matched = len(i_sort)
print("N_matched = ",N_matched)

# a little confused about the i_cut
i_cut = np.sort(i_sort)
i_sort_rev = np.argsort(i_sort)
fp_inds_halo = fp_inds_halo[i_cut]
dmo_inds_halo = dmo_inds_halo[i_cut]
i_sort = np.argsort(Group_M_Mean200_dm[dmo_inds_halo])[::-1]
N_matched = len(i_sort)
print("N_matched = ",N_matched)
i_sort_rev = np.argsort(i_sort)
Group_M_Mean200_dm_sorted = (Group_M_Mean200_dm[dmo_inds_halo])[i_sort]

# active are the ones after the first N_exc_mm
active = np.ones(N_halo_dm,dtype=bool)
active_sorted = (active[dmo_inds_halo])[i_sort]
active_sorted[:N_exc_mm] = False

# split into mass bins
n_arr = []
n1 = N_exc_mm
for n2 in range(N_exc_mm+1,N_matched):
    mass = Group_M_Mean200_dm_sorted[n1:n2]
    frac = ((mass[0]-mass[-1])/np.mean(mass))
    
    if frac>percent_thresh:
        n1 = n2-1
        n_arr.append(n1)
n_arr = np.array(n_arr)
n_arr = n_arr.astype(int)
n_bins = len(n_arr)


# counts of gals for all halos
count_halo_fp = np.zeros(N_halo_fp,dtype=int)
count_halo_dm = np.zeros(N_halo_dm,dtype=int)

if want_fp_pos:
    # nstart of gals for all halos
    nstart_halo_fp = np.zeros(N_halo_fp,dtype=int)-1
    nstart_halo_dm = np.zeros(N_halo_dm,dtype=int)-1
    
    # for each halo, how many gal point to it (excluding case of 0) and index of halo
    # order the parents they are already sorted
    new_sort = np.argsort(parent_gal_fp)

    # order the positions by parents. index returns indices of first occurrence
    pos_gal_fp_nsorted = pos_gal_fp[new_sort]
    parent_gal_un_fp, index,counts = np.unique(parent_gal_fp[new_sort],return_index=True,return_counts=True)
    # npstart is like the cumulative npout starting from 0
    nstart = np.cumsum(counts)
    nstart[1:] = nstart[:-1]
    nstart[0] = 0

    count_halo_fp[parent_gal_un_fp] = counts
    nstart_halo_fp[parent_gal_un_fp] = nstart
    count_halo_dm[dmo_inds_halo] += count_halo_fp[fp_inds_halo]
    nstart_halo_dm[dmo_inds_halo] = nstart_halo_fp[fp_inds_halo]
else:
    # for each halo, how many gal point to it (excluding case of 0) and index of halo
    parent_gal_un_fp,counts = np.unique(parent_gal_fp,return_counts=True)
                
    count_halo_fp[parent_gal_un_fp] += counts
    # use me when excluding the 2063 objext and the matched halos are 1:1
    count_halo_dm[dmo_inds_halo] += count_halo_fp[fp_inds_halo]

print("sum halo dm = ",np.sum(count_halo_dm))
print("sum halo fp = ",np.sum(count_halo_fp))

# sort by number of galaxies
count_halo_dm_sorted = (count_halo_dm[dmo_inds_halo])[i_sort]

# sort by pot energy two defs
min_pot_halo_dm_sorted = (min_pot_halo_dm[dmo_inds_halo])[i_sort]
tot_pot_halo_dm_sorted = (tot_pot_halo_dm[dmo_inds_halo])[i_sort]

# sort by env two defs
env_mass_halo_dm_sorted = (env_mass_halo_dm[dmo_inds_halo])[i_sort]
env_halo_dm_sorted = (env_halo_dm[dmo_inds_halo])[i_sort]
if opt == 'partial_env_cw':
    env_cw_halo_dm_sorted = (env_cw_dm[dmo_inds_halo])[i_sort]

# sort by s2R
s2r_halo_dm_sorted = ((Group_sigmasq_R)[dmo_inds_halo])[i_sort]

# velocity anisotropy
vani_halo_dm_sorted = (GroupVelAnis_dm[dmo_inds_halo])[i_sort]

# sorted indices of the matched array
dmo_inds_halo_sorted = dmo_inds_halo[i_sort]

# get the counts of the active objects
active[dmo_inds_halo] = active_sorted[i_sort_rev]
count_halo_dm[np.logical_not(active)] = 0.
print("sum halo dm = ",np.sum(count_halo_dm))
print("gals in the "+str(N_exc_mm)+" most massive halos = ",np.sum(count_halo_dm[~active]))

def get_counts(halo_dm_sorted,r_corr,reverse_ordering):

    # copy the original dmo index array
    dmo_inds_halo_shuff_sorted = dmo_inds_halo_sorted.copy()

    # we are starting from the first 100th halo
    N_old = N_exc_mm
    h = 0
    while N_old < N_matched:
        # get the next bin
        if h < n_bins-1:
            N_temp = n_arr[h]
            h += 1

        # in case you run out of bins
        inc = 500
        if h >= n_bins-1:
            N_temp = N_old+inc
        N_diff = N_temp-N_old

        # makes sense slicing because temp automatically set to length arr
        dm_i = dmo_inds_halo_sorted[N_old:N_temp]
        mass = Group_M_Mean200_dm_sorted[N_old:N_temp]

        # what has been touched
        active_sorted[N_old:N_temp] = True

        # the galaxies in this bin
        gal_count = count_halo_dm_sorted[N_old:N_temp]
        i_count_sort = np.argsort(gal_count)[::-1]
        gal_count_sort = gal_count[i_count_sort]
        i_count_sort_rev = np.argsort(i_count_sort)
        dm_i_sort = dm_i[i_count_sort]

        # choose reordering scheme
        opt_sorted = halo_dm_sorted[N_old:N_temp]

        # in space of satellite ordered
        opt_sorted_sort = (opt_sorted)[i_count_sort]

        # do partial correlation
        if r_corr is not None:
            mean = [0, 0]#.45#.2
            cov = [[1, r_corr], [r_corr, 1]];
            
            # rx corr to Nsat; ry to fenv
            rx, ry = np.random.multivariate_normal(mean, cov, len(opt_sorted_sort)).T
            i_rx_s = np.argsort(rx)[::-1]
            i_ry_s = np.argsort(ry)[::-1]
            i_rx_sr = np.argsort(i_rx_s)
            i_ry_sr = np.argsort(i_ry_s)
            bleh = i_ry_sr[np.argsort(i_rx_sr)]

            # how do we rank order the occupations
            if reverse_ordering:
                i_opt_sss = (np.argsort(opt_sorted_sort))[bleh]
            else:
                i_opt_sss = (np.argsort(opt_sorted_sort)[::-1])[bleh]
        else: i_opt_sss = np.argsort(opt_sorted_sort)[::-1]

        # get the reordered indices and order them back in the original count order
        dm_i = (dm_i_sort[i_opt_sss])[i_count_sort_rev]

        # record the reordered galaxies
        dmo_inds_halo_shuff_sorted[N_old:N_temp] = dm_i
        N_old = N_temp

    # back to original mass order 
    dmo_inds_halo_shuff = dmo_inds_halo_shuff_sorted[i_sort_rev]
    
    # get the counts for the reordered galaxies
    count_halo_dm_shuff = np.zeros(N_halo_dm,dtype=int)
    count_halo_dm_shuff[dmo_inds_halo_shuff] += count_halo_fp[fp_inds_halo]
    count_halo_dm_shuff[np.logical_not(active)] = 0.

    if want_fp_pos:
        # start of the subhalos
        nstart_halo_dm_shuff = np.zeros(N_halo_dm,dtype=int)-1
        nstart_halo_dm_shuff[dmo_inds_halo_shuff] = nstart_halo_fp[fp_inds_halo]
        return count_halo_dm_shuff, nstart_halo_dm_shuff
    else:
        return count_halo_dm_shuff

# assigning galaxies to halos
def get_weights(counts):
    weights = np.zeros(N_sub_dm)

    # assign first gals to one with highest peak
    ind_gal_b1 = inds_dm[counts >= 1]
    ind_fsub_gal = GroupFirstSub_dm[ind_gal_b1]
    Nsubs_gal = GroupNsubs_dm[ind_gal_b1]

    # faster alorithm
    for l in range(len(ind_fsub_gal)):
        N_g = counts[ind_gal_b1[l]]
        i_sub = inds_sub_dm[ind_fsub_gal[l]:(ind_fsub_gal[l]+Nsubs_gal[l])]
        sor_sub = np.argsort(SubhaloVpeak_dm[i_sub])[::-1]
        chosens = i_sub[sor_sub]
        weights[chosens[:N_g]] = 1.
    return weights

def get_xyz_fp_pos(nstart_halo,count_halo,want_pruning=True,wicked_testing=False):

    N_hal = np.sum(count_halo > 0.5)
    N_gal = np.sum(count_halo)
    
    xyz = np.zeros((N_gal,3))

    cou = count_halo[count_halo > 0.5]
    nst = nstart_halo[count_halo > 0.5]
    if wicked_testing:
        grfir = GroupFirstSub_fp[count_halo > 0.5]
    else:
        grfir = GroupFirstSub_dm[count_halo > 0.5]
    cum = 0

    for k in range(N_hal):
        
        pos = pos_gal_fp_nsorted[nst[k]:nst[k]+cou[k]]
        pos -= pos_gal_fp_nsorted[nst[k]]

        if wicked_testing:
            xyz[cum:cum+cou[k]] = pos+SubhaloPos_fp[grfir[k]]
        else:
            xyz[cum:cum+cou[k]] = pos+SubhaloPos_dm[grfir[k]]
        cum += cou[k]


    if want_pruning:
        x = xyz[:,0];y = xyz[:,1];z = xyz[:,2]
        #w = np.full_like(x,1.)

        x[x<0.] = Lbox-x[x<0.]
        z[z<0.] = Lbox-z[z<0.]
        y[y<0.] = Lbox-y[y<0.]
        x[x>Lbox] = -Lbox+x[x>Lbox]
        z[z>Lbox] = -Lbox+z[z>Lbox]
        y[y>Lbox] = -Lbox+y[y>Lbox]
        xyz = np.vstack((x,y,z)).T
        
    return xyz

def get_xyz(weights,want_weights=False,want_pruning=False):
    if want_pruning:
        x = SubhaloPos_dm[weights > 1.e-12,0]
        y = SubhaloPos_dm[weights > 1.e-12,1]
        z = SubhaloPos_dm[weights > 1.e-12,2]
        x[x<0.] = Lbox-x[x<0.]
        z[z<0.] = Lbox-z[z<0.]
        y[y<0.] = Lbox-y[y<0.]
        x[x>Lbox] = -Lbox+x[x>Lbox]
        z[z>Lbox] = -Lbox+z[z>Lbox]
        y[y>Lbox] = -Lbox+y[y>Lbox]
        xyz = np.vstack((x,y,z)).T
    else:
        xyz = SubhaloPos_dm[weights > 1.e-12]

    if want_weights:
        weights = weights[weights > 1.e-12]
        wei = np.full_like(x,1.)
        wei[:] = weights[:]
        
    if want_weights:
        return xyz, wei
    else:
        return xyz

# jackknifing
def get_jacked_corrs(xyz_FP,xyz_DM,xyz_shuff,w_DM=None,w=None,want_plots=False):

    # setting up
    N_dim = 3
    N_bin = 31
    Rat_DM = np.zeros((N_bin-1,N_dim**3))
    Rat_HOD = np.zeros((N_bin-1,N_dim**3))
    Corr_f_DM = np.zeros((N_bin-1,N_dim**3))
    Corr_f_FP = np.zeros((N_bin-1,N_dim**3))
    Corr_f_DM_shuff = np.zeros((N_bin-1,N_dim**3))
    for i_x in range(N_dim):
        for i_y in range(N_dim):
            for i_z in range(N_dim):
                # start by copying the original arrays
                xyz_FP_jack = xyz_FP.copy()
                xyz_shuff_jack = xyz_shuff.copy()
                if w is not None: w_shuff_jack = w.copy()
                xyz_DM_jack = xyz_DM.copy()
                if w_DM is not None: w_DM_jack = w_DM.copy()

                # current indices
                xyz = np.array([i_x,i_y,i_z],dtype=int)
                size = Lbox/N_dim

                # do the jackknifing for the full-physics
                bool_arr = np.prod((xyz == (xyz_FP/size).astype(int)),axis=1).astype(bool)
                xyz_FP_jack[bool_arr] = np.array([0.,0.,0.])
                xyz_FP_jack = xyz_FP_jack[np.sum(xyz_FP_jack,axis=1)!=0.]

                # do the jackknifing for the normal order
                bool_arr = np.prod((xyz == (xyz_DM/size).astype(int)),axis=1).astype(bool)
                xyz_DM_jack[bool_arr] = np.array([0.,0.,0.])
                if w_DM is not None: w_DM_jack = w_DM[np.sum(xyz_DM_jack,axis=1)!=0.]
                else: w_DM_jack = None
                xyz_DM_jack = xyz_DM_jack[np.sum(xyz_DM_jack,axis=1)!=0.]

                # do the jackknifing for the reordered
                bool_arr = np.prod((xyz == (xyz_shuff/size).astype(int)),axis=1).astype(bool)
                xyz_shuff_jack[bool_arr] = np.array([0.,0.,0.])
                if w is not None: w_shuff_jack = w[np.sum(xyz_shuff_jack,axis=1)!=0.]
                else: w_shuff_jack = None
                xyz_shuff_jack = xyz_shuff_jack[np.sum(xyz_shuff_jack,axis=1)!=0.]


                # unpack the x, y and z positions
                x_DM_jack = xyz_DM_jack[:,0]
                y_DM_jack = xyz_DM_jack[:,1]
                z_DM_jack = xyz_DM_jack[:,2]
                x_shuff_jack = xyz_shuff_jack[:,0]
                y_shuff_jack = xyz_shuff_jack[:,1]
                z_shuff_jack = xyz_shuff_jack[:,2]
                x_FP_jack = xyz_FP_jack[:,0]
                y_FP_jack = xyz_FP_jack[:,1]
                z_FP_jack = xyz_FP_jack[:,2]


                # bins
                bins = np.logspace(-1,1.5,N_bin)
                bin_centers = (bins[:-1] + bins[1:])/2.

                # compute corrfuncs
                results_DM = Corrfunc.theory.xi(X=x_DM_jack,Y=y_DM_jack,Z=z_DM_jack,boxsize=Lbox,nthreads=16,binfile=bins, weights=w_DM_jack)
                results_shuff = Corrfunc.theory.xi(X=x_shuff_jack,Y=y_shuff_jack,Z=z_shuff_jack,boxsize=Lbox,nthreads=16,binfile=bins, weights=w_shuff_jack)
                results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,boxsize=Lbox, nthreads=16,binfile=bins)

                # record results
                Corr_f_FP[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_FP['xi']
                Corr_f_DM[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_DM['xi']
                Corr_f_DM_shuff[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_shuff['xi']
                Rat_DM[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_shuff['xi']/results_DM['xi']
                Rat_HOD[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_shuff['xi']/(xi_fid_shuff/bin_centers**2)

    Rat_DM_mean = np.mean(Rat_DM,axis=1)
    Rat_DM_err = np.sqrt(N_dim**3-1)*np.std(Rat_DM,axis=1)

    Rat_HOD_mean = np.mean(Rat_HOD,axis=1)
    Rat_HOD_err = np.sqrt(N_dim**3-1)*np.std(Rat_HOD,axis=1)

    Corr_mean_DM = np.mean(Corr_f_DM,axis=1)
    Corr_err_DM = np.sqrt(N_dim**3-1)*np.std(Corr_f_DM,axis=1)

    Corr_mean_DM_shuff = np.mean(Corr_f_DM_shuff,axis=1)
    Corr_err_DM_shuff = np.sqrt(N_dim**3-1)*np.std(Corr_f_DM_shuff,axis=1)

    Corr_mean_FP = np.mean(Corr_f_FP,axis=1)
    Corr_err_FP = np.sqrt(N_dim**3-1)*np.std(Corr_f_FP,axis=1)

    if want_plots:
        fs = (12,8)
        fig = plt.figure(figsize=fs)
        plt.errorbar(bin_centers, Corr_mean_FP*bin_centers**2,yerr=Corr_err_FP*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label=str(res)+ ' FP')
        plt.errorbar(bin_centers, Corr_mean_DM*bin_centers**2,yerr=Corr_err_DM*bin_centers**2,lw=2.,ls='-',c='orange',fmt='o',capsize=2,label=str('TNG DM gals'))
        plt.errorbar(bin_centers, Corr_mean_DM_shuff*bin_centers**2,yerr=Corr_err_DM_shuff*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label='TNG DM env')
        plt.xscale('log')
        plt.ylabel(r'$\xi(r) r^2$')
        plt.xlabel(r'$r$ [Mpc/h]')
        plt.legend()
        plt.ylim([0,100])
        plt.text(0.1, 6., str(int(percent_thresh*100))+'% mass bin', dict(size=15))
        plt.legend(loc='lower right')
        plt.savefig('figs/'+figname)
        plt.close()

        fig = plt.figure(figsize=fs)
        plt.plot(bin_centers,np.ones(len(bin_centers)),lw=1.5,ls='--',alpha=0.2)
        plt.errorbar(bin_centers,Rat_DM_mean,yerr=Rat_DM_err,lw=2.,ls='-',c='orange',fmt='o',capsize=4,label=opt)
        plt.xscale('log')
        plt.ylabel(r'$\xi(r)_{\rm TNG, shuff} / \xi(r)_{\rm TNG}$')
        plt.xlabel(r'$r$ [Mpc/h]')
        plt.ylim([0.7,1.3]);
        plt.text(0.1, 0.7, str(int(percent_thresh*100))+'% mass bin', dict(size=15))
        plt.savefig('figs/'+figname[:-4]+'_rat.png')
        plt.close()

    return bin_centers, Corr_mean_FP, Corr_err_FP, Corr_mean_DM, Corr_err_DM, Corr_mean_DM_shuff, Corr_err_DM_shuff, Rat_DM_mean, Rat_DM_err, Rat_HOD_mean, Rat_HOD_err

def get_corr(xyz,wei=None,want_plot=None):
    # split into x, y and z
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    # bins
    N_bin = 31
    bins = np.logspace(-1,1.5,N_bin)
    # TODO: make nicer
    bins = bins[i_r_start-1:]#i_r_end]
    bin_centers = (bins[:-1] + bins[1:])/2.
        
    # compute corrfunc
    results_DM = Corrfunc.theory.xi(X=x,Y=y,Z=z,boxsize=Lbox,nthreads=16,binfile=bins, weights=wei)
    Corr = results_DM['xi']
    
    if want_plot is not None:
        plt.plot(bin_centers, Corr*bin_centers**2,lw=2.,ls='-',label=want_plot)#,c='silver')
        plt.xscale('log')
        plt.ylabel(r'$\xi(r) r^2$')
        plt.xlabel(r'$r$ [Mpc/h]')
        plt.ylim([0,100])
        plt.text(0.1, 6., str(int(percent_thresh*100))+'% mass bin', dict(size=15))
        plt.legend(loc='lower right')

    return Corr

def get_chi2(r_corr,halo_dm_sorted,rev,verbose=False):
    if r_corr > 1. or r_corr < 0.: return np.inf
    if want_fp_pos:
        count_halo, nstart_halo = get_counts(halo_dm_sorted,r_corr=r_corr,reverse_ordering=rev)
        xyz = get_xyz_fp_pos(nstart_halo,count_halo)
    else:
        count_halo = get_counts(halo_dm_sorted,r_corr=r_corr,reverse_ordering=rev)
        weight_sub = get_weights(count_halo)
        xyz = get_xyz(weight_sub)
    corr = get_corr(xyz,want_plot=None)
    dcorr = corr_dm-corr
    chi2 = np.dot(dcorr,np.dot(prec,dcorr))
    if verbose:
        print("chi2 = ",chi2)
        print("r = ",r_corr[0])
        print("________________________")
    return chi2

# indices of all halos and subhalos
inds_dm = np.arange(N_halo_dm)
inds_sub_dm = np.arange(N_sub_dm)

# corr func of the true galaxy distribution
if want_fp_pos:
    xyz_dm = get_xyz_fp_pos(nstart_halo_dm,count_halo_dm)
else:
    weights_sub = get_weights(count_halo_dm)
    xyz_dm = get_xyz(weights_sub)
corr_dm = get_corr(xyz_dm)

# get the fiducial values
bin_centers, xi_fp_mean,xi_fp_err,xi_dm_mean,xi_dm_err,xi_fid_shuff,xi_fid_shuff_err,rat_fid_shuff,rat_fid_shuff_err = np.loadtxt('/home/boryanah/lars/LSSIllustrisTNG/figs/paper/Corr_shuff_'+proxy+'.txt',unpack=True)
xi_fp_err /= bin_centers**2

print("starting bin = ",bin_centers[i_r_start])
print("finishing bin = ",bin_centers[-1])#[i_r_end-1])
xi_fp_err = xi_fp_err[i_r_start:]#i_r_end]
prec = np.diag(1./xi_fp_err**2)

# make choice for which parameter you want
if opt == 'partial_fenv':     par_halo_dm_sorted = env_mass_halo_dm_sorted; reverse = 0
if opt == 'partial_env_cw':     par_halo_dm_sorted = env_cw_halo_dm_sorted; reverse = 0
if opt == 'partial_env':      par_halo_dm_sorted = env_halo_dm_sorted; reverse = 0
if opt == 'partial_s2r':      par_halo_dm_sorted = s2r_halo_dm_sorted; reverse = 1
if opt == 'partial_vani':     par_halo_dm_sorted = vani_halo_dm_sorted; reverse = 1
if opt == 'partial_tot_pot':  par_halo_dm_sorted = tot_pot_halo_dm_sorted; reverse = 0
if opt == 'partial_min_pot':  par_halo_dm_sorted = min_pot_halo_dm_sorted; reverse = 0

# run minimizer
if want_minimize:
    r = 0.5
    res = minimize(get_chi2, r, args=(par_halo_dm_sorted,reverse), method='powell',options={'xtol': 1.e-4, 'disp': True})
    r = res.x
    print("r =",r)
else:
    r = float(sys.argv[2])


# get counts and weights of the final answer
if want_fp_pos:
    count_halo_dm_shuff, nstart_halo_dm_shuff = get_counts(par_halo_dm_sorted,r_corr=r,reverse_ordering=reverse)
    xyz_sh = get_xyz_fp_pos(nstart_halo_dm_shuff,count_halo_dm_shuff)
else:
    count_halo_dm_shuff = get_counts(par_halo_dm_sorted,r_corr=r,reverse_ordering=reverse)
    weights_sub_shuff = get_weights(count_halo_dm_shuff)
    xyz_sh = get_xyz(weights_sub_shuff)

if want_fp_pos:
    ext = "pos"
else:
    ext = "peak"

# get positions of the full physics galaxies
# TESTING WICKED IDEA
#xyz_fp = pos_gal_fp
'''
i_sort = np.argsort(Group_M_Crit200_fp[fp_inds_halo])[::-1]
i_sort_rev = np.argsort(i_sort)
active = np.ones(N_halo_fp,dtype=bool)
active_sorted = (active[fp_inds_halo])[i_sort]
active_sorted[:N_exc_mm+2] = False
active[fp_inds_halo] = active_sorted[i_sort_rev]
count_halo_fp[np.logical_not(active)] = 0.
xyz_fp = get_xyz_fp_pos(nstart_halo_fp,count_halo_fp,wicked_testing=True)
print(xyz_fp.shape)
'''

# save positions of the reordered galaxies
if N_exc_mm == 0:
    np.save("../Lensing/data_2dhod_"+ext+"/"+proxy+"_"+opt+"_gals_all.npy",xyz_sh)
else:
    # TESTING
    #np.save("../Lensing/data_2dhod_"+ext+"/"+proxy+"_"+opt+"_gals.npy",xyz_sh)
    np.save("../Lensing/data_2dhod_"+ext+"/"+proxy+"_"+opt+"_gals_perfect.npy",xyz_sh)
    
fs = (12,8)
fig = plt.figure(figsize=fs)
#corr_fp = get_corr(xyz_fp,want_plot='Full-physics')
corr_dm = get_corr(xyz_dm,want_plot='Dark Matter')
corr_sh = get_corr(xyz_sh,want_plot='Reordered')
plt.legend()
plt.savefig('figs/Corrs_'+opt+'.png')
#plt.show()
plt.close()
