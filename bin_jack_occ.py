import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Corrfunc
import h5py
import plotparams
plotparams.default()
import sys

st_cut = 10000

want_fp_pos = 1 # do you want to take the galaxy position from the FP galaxies
if len(sys.argv) > 1:
    if sys.argv[1] == '--help':
        print('proxy = m200c m200m mfof vmax vpeak m500c s2r')
        print('opt = env_cw spin vdisp min_pot tot_pot s2r vani vir moverr nsubs form env_mass vpeak concs concl env shuff conc partial_fenv partial_min_pot partial_tot_pot partial_s2r partial_vani submass m200m')
        exit(1)
    
    opt = sys.argv[2]
    if opt[:7] == 'partial':
        r = sys.argv[3]; mean = [0, 0]#.5#.2
        cov = [[1, float(r)], [float(r), 1]]; print(r)                  
    proxy = sys.argv[1]
else:
    opt = 'shuff'#'spin'#'vdisp'#'vani'#'env_mass'#'shuff'#'form'#'env_mass'#'vpeak'#'concs'#'concl'#'env'#'shuff' #'conc'
    proxy = 'm200m'#'m200m'#'mfof'#'vmax'#'vpeak'#'' 
binning = 'range'#'std' #'range'
fix_1bin = True
exclude_most_mass = True
N_exc_mm = 100#100#30#100
show_fid = False#True
if fix_1bin == True:
    skip = 500
    from mpi4py import MPI
    i_bin = MPI.COMM_WORLD.Get_rank()*skip
    offset = 0
    i_bin += offset
    
percent_num = 5


percent_thresh = percent_num/100.
# bro nikva ideq tova kakvo
N1 = 0
N2 = 100
N_div = 1

if fix_1bin == True:
    fname_shuff = 'figs/many/Corr_poiss_shuff_'+str(st_cut)+'_thr5_'+str(i_bin)+'_'+str(i_bin+skip)+'.txt'
    #figname = 'figs/many/Corr_poiss_'+opt+'_'+str(st_cut)+'_thr'+str(percent_num)+'_'+str(i_bin)+'_'+str(i_bin+skip)+'.png'
    if exclude_most_mass == False: figname = 'figs/paper/Corr_'+opt+'_'+proxy+'_withmostmassive.png'
    else: figname = 'figs/paper/Corr_'+opt+'_'+proxy+'.png'
else:
    fname_shuff = 'Corr_poiss_shuff_'+str(st_cut)+'_thr1_N1_'+str(N1)+'_N2_'+str(N2)+'.txt'    
    figname = 'Corr_poiss_'+opt+'_'+str(st_cut)+'_thr'+str(percent_num)+'_N1_'+str(N1)+'_N2_'+str(N2)+'.png'

roots = ['/home/boryanah/lars/data_1820/','/home/boryanah/lars/data_455/','/home/boryanah/lars/data_2500/']
N = len(roots)
res = [1820,455,2500]
#res = ['TNG100','TNG100','TNG300']
alphas = [.005,.002,.01]
m_stars = 100*[9.44e5, 6.04e7, 7.44e6]
Lboxs = [75.,75.,205.]# and bigger!

gal_cuts = [40800,16,st_cut]#[110800,16,st_cut]#
#kmin = 2*np.pi/Lboxs[0]/1000.
kmin = 0.
dk = 0.05/1000.
cs = ['r','g','b']

i = 2; root = roots[i]

# matching subhalos in fp and dm
fp_dmo_inds = np.load(root+'sub_match_fp_dmo_'+str(res[i])+'.npy')
fp_inds = fp_dmo_inds[0]
dmo_inds = fp_dmo_inds[1]


# loading data
SubhaloGrNr_dm = np.load(root+'SubhaloGrNr_dm.npy')
GroupFirstSub_dm = np.load(root+'GroupFirstSub_dm.npy')

#GroupNsubs_dm = np.load(root+'GroupNsubs_200_dm.npy')
#Subhalo_200_bool_dm = np.load(root+'Subhalo_200_bool_dm.npy')
GroupNsubs_dm = np.load(root+'GroupNsubs_dm.npy')
SubhaloVmax_dm = np.load(root+'SubhaloVmax_dm.npy')

SubhaloVpeak_dm = np.load(root+'SubhaloVpeak_dm.npy')
SubhaloVpeak_dm[SubhaloVpeak_dm==0.] = SubhaloVmax_dm[SubhaloVpeak_dm==0.]#*20.
SubhaloPos_dm = np.load(root+'SubhaloPos_dm.npy')/1000.
SubhaloSpin_dm = np.load(root+'SubhaloSpin_dm.npy')
SubhaloSpin_dm = np.sum(SubhaloSpin_dm,axis=1)
SubhaloHalfmassRad_dm = np.load(root+'SubhaloHalfmassRad_dm.npy')
SubhaloMassInHalfRad_dm = np.load(root+'SubhaloMassInHalfRad_dm.npy')

Group_M_Crit200_dm = np.load(root+'Group_M_Crit200_dm.npy')
Group_M_Crit500_dm = np.load(root+'Group_M_Crit500_dm.npy')

Group_M_Mean200_dm = np.load(root+'Group_M_Mean200_dm.npy')
Group_R_Mean200_dm = np.load(root+'Group_R_Mean200_dm.npy')
G_N = 4.302e-3
Group_V_Mean200_dm = np.sqrt(G_N*Group_M_Mean200_dm*1.e10/(Group_R_Mean200_dm*1.e3))
Group_Vmax_dm = np.load(root+'Group_Vmax_dm.npy')
#Group_M_over_R = Group_R_Mean200_dm
#Group_M_over_R = Group_M_Mean200_dm/Group_R_Mean200_dm


GroupVmaxRad_dm = np.load(root+'GroupVmaxRad_dm.npy')
SubhaloVmaxRad_dm = np.load(root+'SubhaloVmaxRad_dm.npy')
Group_V_Crit200_dm = np.load(root+'Group_V_Crit200_dm.npy')
GroupPos_dm = np.load(root+'GroupPos_dm.npy')/1000.
GroupMassType_dm = np.load(root+'GroupMassType_dm.npy')

GroupSnapForm_dm = np.load(root+'GroupSnapForm_dm.npy') #deprecated to first 78962
SubhaloMassType_dm = np.load(root+'SubhaloMassType_dm.npy')
# loading data
SubhaloGrNr_fp = np.load(root+'SubhaloGrNr_fp.npy')
SubhaloVelDisp_dm = np.load(root+'SubhaloVelDisp_dm.npy')
#SubhaloVelAnis_dm = np.load(root+'SubhaloVelAnis_dm_GroupFirst.npy')
GroupVelAnis_dm_part = np.load(root+'GroupVelAnis_dm_300000.npy')
GroupFirstSub_fp = np.load(root+'GroupFirstSub_fp.npy')
GroupNsubs_fp = np.load(root+'GroupNsubs_fp.npy')
SubhaloPos_fp = np.load(root+'SubhaloPos_fp.npy')/1000.
Group_M_Crit200_fp = np.load(root+'Group_M_Crit200_fp.npy')

SubhaloMassType_fp = np.load(root+'SubhaloMassType_fp.npy')
SubhaloLenType_fp = np.load(root+'SubhaloLenType_fp.npy')
    
GroupMassType_fp = np.load(root+'GroupMassType_fp.npy')
GroupPos_fp = np.load(root+'GroupPos_fp.npy')/1000.

#print(np.sum(GroupFirstSub_fp == -1))
#print(np.sum(GroupFirstSub_dm == -1))
#2860643
#499242
        
# redshift of slice and mesh size
reshift = 2.22044604925e-16
Nmesh = 128

# for getting the group vpeak and vdisp and spin
unique_sub_grnr, firsts = np.unique(SubhaloGrNr_dm,return_index=True)
GroupVpeak_dm = np.zeros(len(GroupFirstSub_dm))
GroupVelDisp_dm = np.zeros(len(GroupFirstSub_dm))
GroupVelAnis_dm = np.zeros(len(GroupFirstSub_dm))
Group_M_over_R = np.zeros(len(GroupFirstSub_dm))
GroupSpin_dm = np.zeros(len(GroupFirstSub_dm))
Group_sigmasq_R = np.zeros(len(GroupFirstSub_dm))
Group_virial = np.zeros(len(GroupFirstSub_dm))
GroupSubmass_dm = np.zeros(len(GroupFirstSub_dm))

# do peak and disp and spin and vani and s2r
GroupVelDisp_dm[:] = Group_Vmax_dm[:]
GroupVelDisp_dm[unique_sub_grnr] = SubhaloVelDisp_dm[firsts]
Group_M_over_R[:] = Group_Vmax_dm[:]
Group_M_over_R[unique_sub_grnr] = (SubhaloMassType_dm[:,1]/SubhaloHalfmassRad_dm)[firsts]
#Group_M_over_R[unique_sub_grnr] = (SubhaloMassInHalfRad_dm/SubhaloHalfmassRad_dm)[firsts]
GroupVelAnis_dm[:] = Group_Vmax_dm[:]
GroupVelAnis_dm[:300000] = GroupVelAnis_dm_part
GroupVpeak_dm[:] = Group_Vmax_dm[:]
GroupVpeak_dm[unique_sub_grnr] = SubhaloVpeak_dm[firsts]
GroupSubmass_dm[:] = Group_Vmax_dm[:]
GroupSubmass_dm[unique_sub_grnr] = SubhaloMassType_dm[:,1][firsts]#SubhaloMassInHalfRad_dm[firsts]#SubhaloMassType_dm[:,1][firsts]
GroupSpin_dm[:] = Group_Vmax_dm[:]
GroupSpin_dm[unique_sub_grnr] = SubhaloSpin_dm[firsts]
GroupSpin_dm /= (Group_V_Mean200_dm*Group_M_Mean200_dm*Group_R_Mean200_dm)
Group_sigmasq_R[:] = Group_Vmax_dm[:]
Group_sigmasq_R[unique_sub_grnr] = (SubhaloVelDisp_dm**2*SubhaloHalfmassRad_dm)[firsts]
#GrNr_fp, GrNr_N_fp = np.unique(SubhaloGrNr_fp,return_counts=True)
#GrNr_dm, GrNr_N_dm = np.unique(SubhaloGrNr_dm,return_counts=True)
Group_virial[:] = Group_Vmax_dm[:]
Group_virial[unique_sub_grnr] = (SubhaloMassType_dm[:,1]/SubhaloVelDisp_dm**2*SubhaloHalfmassRad_dm)[firsts]


# hubble
h = 0.704

# mass and pos of halo and subhalo in dm
if proxy == 'm200c':    M_200c_dm = Group_M_Crit200_dm*1.e10 #M_sun/h
if proxy == 'm500c':    M_200c_dm = Group_M_Crit500_dm*1.e10 #M_sun/h
if proxy == 'm200m':    M_200c_dm = Group_M_Mean200_dm*1.e10 #M_sun/h
if proxy == 'vpeak':    M_200c_dm = GroupVpeak_dm
if proxy == 'vmax':     M_200c_dm = Group_Vmax_dm
if proxy == 'mfof':     M_200c_dm = GroupMassType_dm[:,1]*1.e10 #M_sun/h
if proxy == 's2r':      M_200c_dm = Group_sigmasq_R

M_200c_fp = Group_M_Crit200_fp*1.e10 #M_sun/h
halo_pos_dm = GroupPos_dm
halo_pos_fp = GroupPos_fp
subhalo_pos_dm = SubhaloPos_dm

# numbers of halos
N_halo_dm = len(M_200c_dm)
N_halo_fp = len(M_200c_fp)
print('N_halo_fp = ', N_halo_fp)
print('N_halo_dm = ', N_halo_dm)
   
'''
# mass cut for halos
inds_dm = np.arange(N_halo_dm)
inds_cut_dm2 = inds_dm[M_200c_dm > 5.e10]
inds_cut_dm = inds_dm[M_200c_dm > 1.e10]
M_200c_cut_dm = M_200c_dm[M_200c_dm > 1.e10]
N_halo_cut_dm = len(M_200c_cut_dm)
'''

# contains pointer to halo parent for each subhalo
#parent_dm = SubhaloGrNr_dm[Subhalo_200_bool_dm == 1]
parent_dm = SubhaloGrNr_dm
parent_fp = SubhaloGrNr_fp
# number of subs
N_sub_dm = len(SubhaloGrNr_dm)
    
N_sub_fp = len(SubhaloGrNr_fp)
print('N_sub_fp = ', N_sub_fp)
print('N_sub_dm = ', N_sub_dm)
# indices of subhalos
inds_sub_dm = np.arange(N_sub_dm)
inds_sub_fp = np.arange(N_sub_fp)

# mass of the halo parent for each sub
parent_mass_fp = M_200c_fp[parent_fp]

# mass and number of star particles for each subhalo in fp

m_sub_stars = SubhaloMassType_fp[:,4]*1.e10 #M_sun/h
l_sub_stars = SubhaloLenType_fp[:,4]
m_sub_dm = SubhaloMassType_dm[:,1]*1.e10 #M_sun/h
    
# applying cut to obtain index of gals, mass of stars, halo parent and pos for each gal
#inds_sub_fp = inds_sub_fp[ m_sub_stars > m_stars[i]]
m_sub_gals = m_sub_stars[m_sub_stars > m_stars[i]]
#parent_gal_fp = parent_fp[m_sub_stars > m_stars[i]]

'''
# since we're excluding no matched
inds_sub_fp = inds_sub_fp[l_sub_stars > gal_cuts[i]]
m_sub_gals = m_sub_stars[l_sub_stars > gal_cuts[i]]
parent_gal_fp = parent_fp[l_sub_stars > gal_cuts[i]]
pos_gal_fp = SubhaloPos_fp[l_sub_stars > gal_cuts[i]]
'''

# load the halo matchings
fp_inds_halo = np.load(root+'fp_dmo_inds_halo_'+str(res[i])+'.npy')[:,0]
dmo_inds_halo = np.load(root+'fp_dmo_inds_halo_'+str(res[i])+'.npy')[:,1]


mag_exc = 2063#1000#2063

# what is the index of the largest subhalo belonging to the first fp halo
# which doesn't have a direct counterpart in DM
sub_exc = GroupFirstSub_fp[fp_inds_halo[-mag_exc:]]
N_exc = GroupNsubs_fp[fp_inds_halo[-mag_exc:]]
for m in range(mag_exc):
    SubhaloPos_fp[sub_exc[m]:N_exc[m]+sub_exc[m]] = np.array([0.,0.,0.])

fp_inds_halo_tot = fp_inds_halo[:]
dmo_inds_halo_tot = dmo_inds_halo[:]

#inds_sub_this_fp = inds_sub_fp[:]
inds_sub_fp = np.arange(N_sub_fp)
inds_halo_fp = np.arange(N_halo_fp)

# removing the excluded unmatched
inds_sub_fp = inds_sub_fp[l_sub_stars > gal_cuts[i]]
m_sub_gals = m_sub_stars[l_sub_stars > gal_cuts[i]]
parent_gal_fp = parent_fp[l_sub_stars > gal_cuts[i]]
pos_gal_fp = SubhaloPos_fp[l_sub_stars > gal_cuts[i]]

# remove those outside
inds_sub_fp = inds_sub_fp[np.sum(pos_gal_fp,axis=1)!=0.]
m_sub_gals = m_sub_gals[np.sum(pos_gal_fp,axis=1)!=0.]
parent_gal_fp = parent_gal_fp[np.sum(pos_gal_fp,axis=1)!=0.]
pos_gal_fp = pos_gal_fp[np.sum(pos_gal_fp,axis=1)!=0.]
            
# galaxy numbers
N_gal = len(inds_sub_fp)
n_gal = N_gal/Lboxs[i]**3
print('n_gal = ',n_gal)
print('N_gal = ',N_gal)

# start anew
fp_inds_halo = fp_inds_halo_tot[:]
dmo_inds_halo = dmo_inds_halo_tot[:]        
            
# remove the double matched ones
fp_inds_halo = fp_inds_halo[:-mag_exc]
dmo_inds_halo = dmo_inds_halo[:-mag_exc]

i_sub_fp = np.arange(N_sub_fp,dtype=int)
i_sub_dm = np.arange(N_sub_dm,dtype=int)

i_sort = np.argsort(M_200c_dm[dmo_inds_halo])[::-1]
N_matched = len(i_sort)
print("N_matched = ",N_matched)
    
i_sort = i_sort[(N1*N_matched//N_div):(N2*N_matched//N_div)]
i_cut = np.sort(i_sort)
i_sort_rev = np.argsort(i_sort)
fp_inds_halo = fp_inds_halo[i_cut]
dmo_inds_halo = dmo_inds_halo[i_cut]
if proxy == 's2r':
    #M_200c_dm = 1./M_200c_dm
    i_sort = np.argsort(M_200c_dm[dmo_inds_halo])[::-1]
else:
    i_sort = np.argsort(M_200c_dm[dmo_inds_halo])[::-1]
N_matched = len(i_sort)
print("N_matched = ",N_matched)
i_sort_rev = np.argsort(i_sort)
M_200c_dm_sorted = (M_200c_dm[dmo_inds_halo])[i_sort]
if fix_1bin == True:
    active = np.zeros(len(M_200c_dm),dtype=bool)
    active_sorted = (active[dmo_inds_halo])[i_sort]

print("N1 = ",N1)
print("N2 = ",N2)
print("N2-N1 = ",N2-N1)

n_arr = []
n1 = N_exc_mm
for n2 in range(N_exc_mm+1,N_matched):
    mass = M_200c_dm_sorted[n1:n2]
    if binning == 'range':
        frac = ((mass[0]-mass[-1])/np.mean(mass))
    elif binning == 'std':
        frac = (np.std(mass)/np.mean(mass))
    if frac>percent_thresh:
        #print(frac)
        #print("%10.5E %10.5E" %(mass[0],mass[-1]))
        #print("n2-n1",n2-n1)

        n1 = n2-1
        n_arr.append(n1)
n_arr = np.array(n_arr)

# minimum halo mass in the HOD
M_200c_dm_min = np.min(M_200c_dm[dmo_inds_halo])

if want_fp_pos:
    # for each halo, how many gal point to it (excluding case of 0) and index of halo
    # order the parents they are already sorted
    new_sort = np.argsort(parent_gal_fp)
    #new_sort = np.arange(len(new_sort))))

    # order the positions by parents. index returns indices of first occurrence
    pos_gal_fp_nsorted = pos_gal_fp[new_sort]
    parent_gal_un_fp, index,counts = np.unique(parent_gal_fp[new_sort],return_index=True,return_counts=True)
    # npstart is like the cumulative npout starting from 0
    nstart = np.cumsum(counts)
    nstart[1:] = nstart[:-1]
    nstart[0] = 0
    #pos_gal_first_fp = pos_gal_fp_nsorted[index] # the position of the central galaxy
    # and then use this npstart and npout arrays to extract from the sorted positions by parent array

    # counts of gals for all halos
    count_halo_fp = np.zeros(len(M_200c_fp),dtype=int)
    nstart_halo_fp = np.zeros(len(M_200c_fp),dtype=int)-1
    count_halo_dm = np.zeros(len(M_200c_dm),dtype=int)
    nstart_halo_dm = np.zeros(len(M_200c_dm),dtype=int)-1
    count_halo_dm_shuff = np.zeros(len(M_200c_dm),dtype=int)
    nstart_halo_dm_shuff = np.zeros(len(M_200c_dm),dtype=int)-1

    count_halo_fp[parent_gal_un_fp] = counts
    nstart_halo_fp[parent_gal_un_fp] = nstart
    # use me when excluding the 2063 objext and the matched halos are 1:1
    count_halo_dm[dmo_inds_halo] += count_halo_fp[fp_inds_halo]

else:
    # mass of parent halo for each gal
    parent_mass_gal = M_200c_fp[parent_gal_fp]
    # for each halo, how many gal point to it (excluding case of 0) and index of halo
    parent_gal_un_fp,counts = np.unique(parent_gal_fp,return_counts=True)
    #parent_mass_un_gal = parent_mass_gal[parent_gal_un_fp]
    # counts of gals for all halos
    count_halo_fp = np.zeros(len(M_200c_fp),dtype=int)
    count_halo_dm = np.zeros(len(M_200c_dm),dtype=int)
    count_halo_dm_shuff = np.zeros(len(M_200c_dm),dtype=int)
            
    count_halo_fp[parent_gal_un_fp] += counts
    # use me when excluding the 2063 objext and the matched halos are 1:1
    count_halo_dm[dmo_inds_halo] += count_halo_fp[fp_inds_halo]

# count the galaxies in the corresponding DM halos
'''
# USE ME WHEN DOUBLE COUNTING
u,b,c = np.unique(dmo_inds_halo,return_index=True,return_counts=True)
print(len(u))
c1 = c==1
print(np.sum(c1))
count_halo_dm[dmo_inds_halo[b[c1]]] += count_halo_fp[fp_inds_halo[b[c1]]]
rs = u[c>1]
for r in rs:
    count_halo_dm[r] += np.sum(count_halo_fp[fp_inds_halo[dmo_inds_halo == r]])
'''

print("sum halo dm = ",np.sum(count_halo_dm))
print("sum halo fp = ",np.sum(count_halo_fp))

# sort by number of galaxies
count_halo_dm_sorted = (count_halo_dm[dmo_inds_halo])[i_sort]

'''
#parent_mass_gal_un_fp = M_200c_fp[parent_gal_un_fp]
#plt.scatter(parent_mass_gal_un_fp,counts,c='y',marker='+',lw=1,alpha=alphas[i])
plt.scatter(GroupVelAnis_dm[dmo_inds_halo],count_halo_dm[dmo_inds_halo],c='y',marker='+',lw=1,alpha=alphas[i])
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Vel. anisotropy of halo')
plt.ylabel('Number of galaxies')
plt.ylim([-1,51])
#plt.xlim([1.e11,1.1e15])
plt.savefig('Ngal_vani_scatter.png')
plt.show()
# for dm plt.scatter(M_200c_dm[dmo_inds_halo],count_halo_dm[dmo_inds_halo],c='y',marker='+',lw=1,alpha=alphas[i])
'''

# sort by env two defs
min_pot_halo_dm = -np.load(root+'GroupMinPotential_dm.npy')
tot_pot_halo_dm = -np.load(root+'GroupTotPotential_dm.npy')
min_pot_halo_dm_sorted = (min_pot_halo_dm[dmo_inds_halo])[i_sort]
tot_pot_halo_dm_sorted = (tot_pot_halo_dm[dmo_inds_halo])[i_sort]
env_mass_halo_dm = np.load(root+'env_mass_dm.npy')
env_halo_dm = np.load(root+'env_dm.npy')
env_mass_halo_dm_sorted = (env_mass_halo_dm[dmo_inds_halo])[i_sort]
env_halo_dm_sorted = (env_halo_dm[dmo_inds_halo])[i_sort]

if opt == "env_cw":
    # Density in Illustris
    filename = 'CosmicWeb/WEB_CIC_256_DM_TNG300-2.hdf5'
    f = h5py.File(filename, 'r')
    d_smooth = f['density_smooth'][:,:,:] 

    # finding who belongs where in the cosmic web
    N_dim = 256
    gr_size = Lboxs[i]/N_dim
    halo_x = GroupPos_dm[:,0];halo_y = GroupPos_dm[:,1];halo_z = GroupPos_dm[:,2]
    i_cw = (halo_x/gr_size).astype(int)
    j_cw = (halo_y/gr_size).astype(int)
    k_cw = (halo_z/gr_size).astype(int)

    # Environment definition
    env_cw = d_smooth[i_cw,j_cw,k_cw]
    env_cw_halo_dm_sorted = (env_cw[dmo_inds_halo])[i_sort]

# sort Nsubs
nsubs_halo_dm_sorted = ((GroupNsubs_dm)[dmo_inds_halo])[i_sort]

# sort m200m
m200m_halo_dm_sorted = ((Group_M_Mean200_dm)[dmo_inds_halo])[i_sort]

# sort by M over R
moverr_halo_dm_sorted = ((Group_M_over_R)[dmo_inds_halo])[i_sort]

# sort by measure of virialization
vir_halo_dm_sorted = ((Group_virial)[dmo_inds_halo])[i_sort]

# sort by s2R
s2r_halo_dm_sorted = ((Group_sigmasq_R)[dmo_inds_halo])[i_sort]

# sort by conc
conc_VRrat_halo_dm = Group_Vmax_dm/GroupVmaxRad_dm
conc_halo_dm = (Group_Vmax_dm/Group_V_Crit200_dm)
conc_halo_dm_sorted = ((conc_VRrat_halo_dm)[dmo_inds_halo])[i_sort]
#conc_halo_dm_sorted = ((conc_halo_dm)[dmo_inds_halo])[i_sort]

# sort by snapshot of formation
GroupSnapForm_dm_sorted = np.zeros(len(M_200c_dm))
GroupSnapForm_dm_sorted[:len(GroupSnapForm_dm)] = GroupSnapForm_dm
GroupSnapForm_dm_sorted = (GroupSnapForm_dm_sorted[dmo_inds_halo])[i_sort]

# Vdispersion
vdisp_halo_dm_sorted = (GroupVelDisp_dm[dmo_inds_halo])[i_sort]

# Vanisotropy
vani_halo_dm_sorted = (GroupVelAnis_dm[dmo_inds_halo])[i_sort]

# submass
submass_halo_dm_sorted = (GroupSubmass_dm[dmo_inds_halo])[i_sort] 

# V200m
v200m_halo_dm_sorted = (Group_V_Mean200_dm[dmo_inds_halo])[i_sort]

# spin
spin_halo_dm_sorted = (GroupSpin_dm[dmo_inds_halo])[i_sort]

# Vpeak
vmax_halo_dm_sorted = (GroupVpeak_dm[dmo_inds_halo])[i_sort]#(Group_Vmax_dm[dmo_inds_halo])[i_sort]#(GroupMassType_dm[dmo_inds_halo,1])[i_sort]
Group_Vmax_dm_sorted = (Group_Vmax_dm[dmo_inds_halo])[i_sort]

N_hal_bin = np.load('N_hal_bin.npy')
b = 0#11*6
ns = 100+np.sum(N_hal_bin[:b])
print(ns)
print(N_hal_bin[b])
ne = ns+np.sum(N_hal_bin[:(b+40)])
print(ne)

'''
x = vani_halo_dm_sorted[5000:5500]
y = s2r_halo_dm_sorted[5000:5500]#count_halo_dm_sorted[2000:2500]#
'''
'''
np.save('env',x)
np.save('Ngal',y)
print("gals = ", np.sum(y))
xedges = np.linspace(np.min(x),np.max(x),30)
yedges = np.linspace(np.min(y),np.max(y),30)
H, xedges, yedges = np.histogram2d(x, y,normed=True, bins=(xedges, yedges))
H = H.T  # Let each row list bins with common y range.
plt.imshow(np.log10(H+1.e-6), interpolation='nearest', origin='low', aspect='auto',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar()
plt.show()
quit()
'''
'''
#plt.scatter(M_200c_dm_sorted,GroupSnapForm_dm_sorted)
#plt.scatter(M_200c_dm_sorted,vmax_halo_dm_sorted)
#plt.scatter(M_200c_dm_sorted,Group_Vmax_dm_sorted)
#plt.scatter(M_200c_dm_sorted,vdisp_halo_dm_sorted)
#plt.scatter(M_200c_dm_sorted[:100000],vani_halo_dm_sorted[:100000],c='y',marker='+',lw=1,alpha=alphas[i])
#plt.scatter(env_halo_dm_sorted[count_halo_dm_sorted>0][100:],count_halo_dm_sorted[count_halo_dm_sorted>0][100:],c='b',marker='+',lw=1,alpha=alphas[i]*15)
plt.scatter(x,y,c='b',marker='+',lw=1,alpha=alphas[i]*50)
#plt.scatter(M_200c_dm_sorted,spin_halo_dm_sorted
# COOL
#plt.scatter(vdisp_halo_dm_sorted/v200m_halo_dm_sorted,conc_halo_dm_sorted,c='y',marker='+',lw=1,alpha=alphas[i])
#plt.scatter(GroupSnapForm_dm_sorted,conc_halo_dm_sorted,c='y',marker='+',lw=1,alpha=alphas[i])
#plt.scatter(GroupSnapForm_dm_sorted,vdisp_halo_dm_sorted/v200m_halo_dm_sorted,c='y',marker='+',lw=1,alpha=alphas[i])
#plt.scatter(vmax_halo_dm_sorted,Group_Vmax_dm_sorted)
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim([1.,4.e5])
#plt.xlim([1.e10,2.e15])
#plt.xlim([10.,100.])
#plt.savefig('show_vmax_m200c.png')
#plt.savefig('show_m200m_s2r.png')
plt.savefig('show_vani_Ngal.png')
#plt.savefig('show_vdisp_m200c.png')
plt.show()
quit()
'''
dmo_inds_halo_sorted = dmo_inds_halo[i_sort]
    
# fluctuation version
#percent_thresh = 5/100.#.0000000000001/100.#1 5 10 20 50
#n_arr = np.loadtxt("n_arr_"+str(percent_thresh)+".txt")
n_arr = n_arr.astype(int)
N_old = N_exc_mm
if fix_1bin == True and exclude_most_mass == False:
    active_sorted[:N_old] = True
h = 0
n_bins = len(n_arr)
N_diff = 0
nbin = []
k = []

while N_old < N_matched:
    if fix_1bin == True and h < n_bins-1:
        if h < i_bin:
            print("HERO",i_bin,h)
            N_old = n_arr[h]
            h += 1

            N_temp = n_arr[h]
            continue
        elif i_bin+skip-1 < h:
            N_old = 1.e12
            continue
    if h < n_bins-1:
        N_temp = n_arr[h]
        h += 1


    N_diff = N_temp-N_old

    inc = 500
    if h >= n_bins-1: #normal
        #if N_diff > inc or h >= len(n_arr)-1: # avoiding too large bins
        #if N_diff < inc: # testing which parts mess up large scales
        #if N_diff > inc or h >= len(n_arr)-1 or N_diff < inc: # testing which parts mess up large scales constant bins i think
        #if N_diff > inc or h <= len(n_arr)-1: # testing should lead to overfl
        N_temp = N_old+inc
        #print(N_temp)

    N_diff = N_temp-N_old
    nbin.append(N_diff)
    # makes sense slicing because temp automatically set to length arr
    dm_i = dmo_inds_halo_sorted[N_old:N_temp]
    mass = M_200c_dm_sorted[N_old:N_temp]

    if fix_1bin == True:
        active_sorted[N_old:N_temp] = True
    gal_count = count_halo_dm_sorted[N_old:N_temp]#M_200c_dm_sorted[N_old:N_temp]#count_halo_dm_sorted[N_old:N_temp]

    if binning == 'range':
        frac = ((mass[0]-mass[-1])/np.mean(mass))
    elif binning == 'std':
        frac = (np.std(mass)/np.mean(mass))
    #print("frac[%] = ",frac*100.)

    if fix_1bin == True and h < n_bins-2:
        if h-1 == i_bin:
            mass_beg = mass[0]
        mass_end = mass[-1]


    i_count_sort = np.argsort(gal_count)[::-1]
    gal_count_sort = gal_count[i_count_sort]
    i_count_sort_rev = np.argsort(i_count_sort)
    dm_i_sort = dm_i[i_count_sort]

    if np.sum(gal_count_sort) == 0: True;#print("IT'S EMPTY"); 
    else:
        True
        #print("mass range: %10.5E - %10.5E, bin size:" %(mass[0],mass[-1]),N_diff)


    #if opt == 'shuff': np.random.seed(3334*len(mass)); np.random.shuffle(dm_i);
    if opt == 'shuff': np.random.shuffle(dm_i);
    else:
        if opt == 'vpeak':
            opt_sorted = vmax_halo_dm_sorted[N_old:N_temp]
        elif opt == 'vdisp':
            opt_sorted = vdisp_halo_dm_sorted[N_old:N_temp]
        elif opt == 'vani' or opt == 'partial_vani':
            opt_sorted = vani_halo_dm_sorted[N_old:N_temp]
        elif opt == 'moverr':
            opt_sorted = moverr_halo_dm_sorted[N_old:N_temp]
        elif opt == 'nsubs':
            opt_sorted = nsubs_halo_dm_sorted[N_old:N_temp]
        elif opt == 's2r' or opt == 'partial_s2r':
            opt_sorted = s2r_halo_dm_sorted[N_old:N_temp]
        elif opt == 'vir':
            opt_sorted = vir_halo_dm_sorted[N_old:N_temp]
        elif opt == 'spin':
            opt_sorted = spin_halo_dm_sorted[N_old:N_temp]
        elif opt == 'conc':
            opt_sorted = conc_halo_dm_sorted[N_old:N_temp]
        elif opt == 'form':
            opt_sorted = GroupSnapForm_dm_sorted[N_old:N_temp]
        elif opt == 'concl':
            opt_sorted = conc_halo_dm_sorted[N_old:N_temp]
        elif opt == 'concs':
            opt_sorted = conc_halo_dm_sorted[N_old:N_temp]
        elif opt == 'm200m':
            opt_sorted = m200m_halo_dm_sorted[N_old:N_temp]
        elif opt == 'env_mass' or opt == 'partial_fenv':
            opt_sorted = env_mass_halo_dm_sorted[N_old:N_temp]
        elif opt == 'env':
            opt_sorted = env_halo_dm_sorted[N_old:N_temp]
        elif opt == 'env_cw':
            opt_sorted = env_cw_halo_dm_sorted[N_old:N_temp]
        elif opt == 'min_pot' or opt == 'partial_min_pot':
            opt_sorted = min_pot_halo_dm_sorted[N_old:N_temp]
        elif opt == 'tot_pot' or opt == 'partial_tot_pot':
            opt_sorted = tot_pot_halo_dm_sorted[N_old:N_temp]
        elif opt == 'submass':
            opt_sorted = submass_halo_dm_sorted[N_old:N_temp]

        # in space of satellite ordered
        opt_sorted_sort = (opt_sorted)[i_count_sort]
        if opt == 'conc':
            # sorted indices, gals and conc in sorted gal space
            i_opt_sss = np.argsort(opt_sorted_sort)
            opt_sss = opt_sorted_sort[i_opt_sss]

            count_bin = np.sort(gal_count_sort)[::-1]
            i_1_0 = np.argmax(count_bin <= 1)
            i_opt_sss[i_1_0:] = (i_opt_sss[i_1_0:])[np.argsort(opt_sss[i_1_0:])[::-1]]


            '''
            elif opt == 'vani': # TESTING
            # sorted indices, gals and conc in sorted gal space
            i_opt_sss = np.argsort(opt_sorted_sort)
            opt_sss = opt_sorted_sort[i_opt_sss]

            count_bin = np.sort(gal_count_sort)[::-1]
            i_1_0 = np.argmax(count_bin <= 1)
            i_opt_sss[i_1_0:] = (i_opt_sss[i_1_0:])[np.argsort(opt_sss[i_1_0:])[::-1]]
            '''    
        elif opt == 'concs': i_opt_sss = np.argsort(opt_sorted_sort)
        elif opt == 'vdisp': i_opt_sss = np.argsort(opt_sorted_sort)
        elif opt == 'moverr': i_opt_sss = np.argsort(opt_sorted_sort)
        #elif opt == 'nsubs': i_opt_sss = np.argsort(opt_sorted_sort)
        elif opt == 'vani': i_opt_sss = np.argsort(opt_sorted_sort)
        elif opt == 's2r': i_opt_sss = np.argsort(opt_sorted_sort)
        elif opt == 'submass': i_opt_sss = np.argsort(opt_sorted_sort)
        elif opt == 'partial_fenv':
            # rx corr to Nsat; ry to fenv
            rx, ry = np.random.multivariate_normal(mean, cov, len(opt_sorted_sort)).T
            i_rx_s = np.argsort(rx)[::-1]
            i_ry_s = np.argsort(ry)[::-1]
            i_rx_sr = np.argsort(i_rx_s)
            i_ry_sr = np.argsort(i_ry_s)

            bleh = i_ry_sr[np.argsort(i_rx_sr)]

            i_opt_sss = (np.argsort(opt_sorted_sort)[::-1])[bleh]#(dm_i_sort[np.argsort(opt_sorted_sort)[::-1]])
            '''
            if len(opt_sorted_sort)>100:
                print(len(opt_sorted_sort))
                plt.scatter((np.sort(opt_sorted_sort)[::-1])[i_opt_sss],gal_count_sort)# is how the elements are ordered in space of N_sat sorted
                plt.show()
            '''
        elif opt == 'partial_min_pot' or opt == 'partial_tot_pot':
            # rx corr to Nsat; ry to fenv
            rx, ry = np.random.multivariate_normal(mean, cov, len(opt_sorted_sort)).T
            i_rx_s = np.argsort(rx)[::-1]
            i_ry_s = np.argsort(ry)[::-1]
            i_rx_sr = np.argsort(i_rx_s)
            i_ry_sr = np.argsort(i_ry_s)

            bleh = i_ry_sr[np.argsort(i_rx_sr)]

            i_opt_sss = (np.argsort(opt_sorted_sort)[::-1])[bleh]
            
        elif opt == 'partial_vani':
            # rx corr to Nsat; ry to vani
            rx, ry = np.random.multivariate_normal(mean, cov, len(opt_sorted_sort)).T
            i_rx_s = np.argsort(rx)[::-1]
            i_ry_s = np.argsort(ry)[::-1]
            i_rx_sr = np.argsort(i_rx_s)
            i_ry_sr = np.argsort(i_ry_s)

            bleh = i_ry_sr[np.argsort(i_rx_sr)]

            print(bleh)
            i_opt_sss = (np.argsort(opt_sorted_sort))[bleh]
        elif opt == 'partial_s2r':
            # rx corr to Nsat; ry to vani
            rx, ry = np.random.multivariate_normal(mean, cov, len(opt_sorted_sort)).T
            i_rx_s = np.argsort(rx)[::-1]
            i_ry_s = np.argsort(ry)[::-1]
            i_rx_sr = np.argsort(i_rx_s)
            i_ry_sr = np.argsort(i_ry_s)

            bleh = i_ry_sr[np.argsort(i_rx_sr)]

            i_opt_sss = (np.argsort(opt_sorted_sort))[bleh]
        else: i_opt_sss = np.argsort(opt_sorted_sort)[::-1]

        dm_i = (dm_i_sort[i_opt_sss])[i_count_sort_rev]
        #print(dm_i)
    dmo_inds_halo_sorted[N_old:N_temp] = dm_i
    N_old = N_temp

nbin = np.array(nbin)
dmo_inds_halo_shuff = dmo_inds_halo_sorted[i_sort_rev]
#print(N_temp)
#print(n_arr[-1])

if fix_1bin == True:
    #active = active_sorted[i_sort_rev]
    active[dmo_inds_halo] = active_sorted[i_sort_rev]
    #dmo_inds_halo_shuff = dmo_inds_halo_shuff[active]

    count_halo_dm = np.zeros(len(M_200c_dm),dtype=int)
    #dmo_inds_halo = dmo_inds_halo[active]

    # count the galaxies in the corresponding DM halos
    u,b,c = np.unique(dmo_inds_halo,return_index=True,return_counts=True)
    print(len(u))
    c1 = c==1
    print(np.sum(c1))
    count_halo_dm[dmo_inds_halo[b[c1]]] += count_halo_fp[fp_inds_halo[b[c1]]]
    rs = u[c>1]
    for r in rs:
        count_halo_dm[r] += np.sum(count_halo_fp[fp_inds_halo[dmo_inds_halo == r]])
    # !!!
    count_halo_dm[np.logical_not(active)] = 0.
    print("sum halo dm = ",np.sum(count_halo_dm))
    print("sum halo fp = ",np.sum(count_halo_fp))


# BEG SHUFF
'''
# USE ME WHEN DOUBLE COUNTING
# counting galaxies
u_shuff,b_shuff,c_shuff = np.unique(dmo_inds_halo_shuff,return_index=True,return_counts=True)

print(len(u))
c1 = c==1
print(np.sum(c1))
c1_shuff = c_shuff==1
print(np.sum(c1))
c2 = c==2
print(np.sum(c2))
c2 = c==3
print(np.sum(c2))
c2 = c==4
print(np.sum(c2))
c2 = c==5
print(np.sum(c2))

count_halo_dm_shuff[dmo_inds_halo_shuff[b_shuff[c1_shuff]]] += count_halo_fp[fp_inds_halo[b_shuff[c1_shuff]]]
rs = u_shuff[c_shuff>1]
for r in rs:
    count_halo_dm_shuff[r] += np.sum(count_halo_fp[fp_inds_halo[dmo_inds_halo_shuff == r]])

'''
# use me when excluding the 2063 objext
count_halo_dm_shuff[dmo_inds_halo_shuff] += count_halo_fp[fp_inds_halo]

if want_fp_pos:
    # start of the subhalos
    nstart_halo_dm[dmo_inds_halo] = nstart_halo_fp[fp_inds_halo]
    nstart_halo_dm_shuff[dmo_inds_halo_shuff] = nstart_halo_fp[fp_inds_halo]

# !!!
if fix_1bin == True:
    count_halo_dm_shuff[np.logical_not(active)] = 0.

if want_fp_pos:

    N_hal = np.sum(count_halo_dm > 0.5)
    N_gal = np.sum(count_halo_dm)

    print("N_gal = ",N_gal)
    xyz_shuff = np.zeros((N_gal,3))
    xyz_DM = np.zeros((N_gal,3))

    cou = count_halo_dm[count_halo_dm > 0.5]
    nst = nstart_halo_dm[count_halo_dm > 0.5]
    grfir = GroupFirstSub_dm[count_halo_dm > 0.5]
    cum = 0
    cou_shuff = count_halo_dm_shuff[count_halo_dm_shuff > 0.5]
    nst_shuff = nstart_halo_dm_shuff[count_halo_dm_shuff > 0.5]
    print(np.sum(nst_shuff==-1))
    print(np.sum(nst==-1))
    print(np.sum(cou_shuff))

    grfir_shuff = GroupFirstSub_dm[count_halo_dm_shuff > 0.5]
    cum_shuff = 0
    for k in range(N_hal):
        
        pos = pos_gal_fp_nsorted[nst[k]:nst[k]+cou[k]]
        pos -= pos_gal_fp_nsorted[nst[k]]#pos_gal_first_fp[k]
        
        xyz_DM[cum:cum+cou[k]] = pos+SubhaloPos_dm[grfir[k]]
        cum += cou[k]


        #if k == 23: quit()
        pos_shuff = pos_gal_fp_nsorted[nst_shuff[k]:nst_shuff[k]+cou_shuff[k]]
        pos_shuff -= pos_gal_fp_nsorted[nst_shuff[k]]#pos_gal_first_fp[k]
        xyz_shuff[cum_shuff:cum_shuff+cou_shuff[k]] = pos_shuff+SubhaloPos_dm[grfir_shuff[k]]
        cum_shuff += cou_shuff[k]
    # subtract the central galaxy tuki

    # BSING TESTING TUKSI
    x_DM = xyz_DM[:,0];y_DM = xyz_DM[:,1];z_DM = xyz_DM[:,2]
    w_DM = np.full_like(x_DM,1.)
    x = xyz_shuff[:,0];y = xyz_shuff[:,1];z = xyz_shuff[:,2]
    w = np.full_like(x,1.)

    x_DM[x_DM<0.] = Lboxs[i]-x_DM[x_DM<0.]
    z_DM[z_DM<0.] = Lboxs[i]-z_DM[z_DM<0.]
    y_DM[y_DM<0.] = Lboxs[i]-y_DM[y_DM<0.]
    x_DM[x_DM>Lboxs[i]] = -Lboxs[i]+x_DM[x_DM>Lboxs[i]]
    z_DM[z_DM>Lboxs[i]] = -Lboxs[i]+z_DM[z_DM>Lboxs[i]]
    y_DM[y_DM>Lboxs[i]] = -Lboxs[i]+y_DM[y_DM>Lboxs[i]]

    x[x<0.] = Lboxs[i]-x[x<0.]
    z[z<0.] = Lboxs[i]-z[z<0.]
    y[y<0.] = Lboxs[i]-y[y<0.]
    x[x>Lboxs[i]] = -Lboxs[i]+x[x>Lboxs[i]]
    z[z>Lboxs[i]] = -Lboxs[i]+z[z>Lboxs[i]]
    y[y>Lboxs[i]] = -Lboxs[i]+y[y>Lboxs[i]]

    xyz_DM = np.vstack((x_DM,y_DM,z_DM)).T
    xyz_shuff = np.vstack((x,y,z)).T
    
else:
    weights_sub = np.zeros(N_sub_dm)

    #SubhaloNsubs_dm = GroupNsubs_dm[parent_dm]
    # assign instead first gals to one with highest peak and not first
    inds_dm = np.arange(N_halo_dm)
    inds_sub_dm = np.arange(N_sub_dm)
    ind_gal_b1 = inds_dm[count_halo_dm >= 1]
    ind_fsub_gal = GroupFirstSub_dm[ind_gal_b1]
    Nsubs_gal = GroupNsubs_dm[ind_gal_b1]
    print(len(ind_gal_b1))
    # faster
    for l in range(len(ind_fsub_gal)):
        N_g = count_halo_dm[ind_gal_b1[l]]
        i_sub = inds_sub_dm[ind_fsub_gal[l]:(ind_fsub_gal[l]+Nsubs_gal[l])]
        #sor_sub = np.argsort(SubhaloVmax_dm[parent_dm == dmh])[::-1]
        sor_sub = np.argsort(SubhaloVpeak_dm[i_sub])[::-1]
        #sor_sub = np.argsort(m_sub_dm[parent_dm == dmh])[::-1]
        #sor_sub = np.argsort(SubhaloGrNr_dm[parent_dm == dmh])[::-1]

        chosens = i_sub[sor_sub]
        weights_sub[chosens[:N_g]] = 1.

    print("weight sum = ",np.sum(weights_sub))

    # BEG SHUFF
    weights_sub_shuff = np.zeros(N_sub_dm)

    inds_dm = np.arange(N_halo_dm)
    inds_sub_dm = np.arange(N_sub_dm)
    ind_gal_b1 = inds_dm[count_halo_dm_shuff >= 1]
    ind_fsub_gal = GroupFirstSub_dm[ind_gal_b1]
    Nsubs_gal = GroupNsubs_dm[ind_gal_b1]
    print(len(ind_gal_b1))

    # faster

    for l in range(len(ind_fsub_gal)):
        N_g = count_halo_dm_shuff[ind_gal_b1[l]]
        i_sub = inds_sub_dm[ind_fsub_gal[l]:(ind_fsub_gal[l]+Nsubs_gal[l])]
        if len(i_sub) <= N_g:
            weights_sub_shuff[i_sub] = N_g*1./len(i_sub)
            continue
        #sor_sub = np.argsort(SubhaloVmax_dm[parent_dm == dmh])[::-1]
        sor_sub = np.argsort(SubhaloVpeak_dm[i_sub])[::-1]
        #sor_sub = np.argsort(m_sub_dm[parent_dm == dmh])[::-1]
        #sor_sub = np.argsort(SubhaloGrNr_dm[parent_dm == dmh])[::-1]
        chosens = i_sub[sor_sub]
        weights_sub_shuff[chosens[:N_g]] = 1.

    print("weight sum shuff = ",np.sum(weights_sub_shuff))
    print("number of gals in "+str(N_exc_mm)+" most massive halos = ",np.sum(count_halo_dm_sorted[:N_exc_mm]))

    # This is where the indentation used to end
    # CORRFUNC
    x_DM = subhalo_pos_dm[:,0]
    y_DM = subhalo_pos_dm[:,1]
    #print("unique y's = ", len(np.unique(y_DM)))
    z_DM = subhalo_pos_dm[:,2]
    #print("len z = ",len(z_DM))

    # cut for only positive weights

    x_DM = x_DM[weights_sub > 1.e-12]
    y_DM = y_DM[weights_sub > 1.e-12]
    z_DM = z_DM[weights_sub > 1.e-12]

    x_DM[x_DM<0.] = Lboxs[i]-x_DM[x_DM<0.]
    z_DM[z_DM<0.] = Lboxs[i]-z_DM[z_DM<0.]
    y_DM[y_DM<0.] = Lboxs[i]-y_DM[y_DM<0.]
    x_DM[x_DM>Lboxs[i]] = -Lboxs[i]+x_DM[x_DM>Lboxs[i]]
    z_DM[z_DM>Lboxs[i]] = -Lboxs[i]+z_DM[z_DM>Lboxs[i]]
    y_DM[y_DM>Lboxs[i]] = -Lboxs[i]+y_DM[y_DM>Lboxs[i]]


    weights_sub = weights_sub[weights_sub > 1.e-12]

    #print("unique z's = ",len(np.unique(z_DM)))
    #print("unique y's = ",len(np.unique(y_DM)))
    #print("unique x's = ",len(np.unique(x_DM)))
    w_DM = np.full_like(x_DM,1.)
    w_DM[:] = weights_sub[:]
    #print("max weight = ",np.max(w_DM))

    x = subhalo_pos_dm[:,0]
    y = subhalo_pos_dm[:,1]
    #print("unique y's = ", len(np.unique(y)))
    z = subhalo_pos_dm[:,2]
    #print("len z = ",len(z))

    x = x[weights_sub_shuff > 1.e-12]
    y = y[weights_sub_shuff > 1.e-12]
    z = z[weights_sub_shuff > 1.e-12]

    x[x<0.] = Lboxs[i]-x[x<0.]
    z[z<0.] = Lboxs[i]-z[z<0.]
    y[y<0.] = Lboxs[i]-y[y<0.]
    x[x>Lboxs[i]] = -Lboxs[i]+x[x>Lboxs[i]]
    z[z>Lboxs[i]] = -Lboxs[i]+z[z>Lboxs[i]]
    y[y>Lboxs[i]] = -Lboxs[i]+y[y>Lboxs[i]]
    weights_sub_shuff = weights_sub_shuff[weights_sub_shuff > 1.e-12]

    #print("unique z's = ",len(np.unique(z)))
    #print("unique y's = ",len(np.unique(y)))
    #print("unique x's = ",len(np.unique(x)))

    w = np.full_like(x,1.)
    w[:] = weights_sub_shuff[:]
    print("max weight = ",np.max(w))

x_FP = pos_gal_fp[:,0]
y_FP = pos_gal_fp[:,1]
z_FP = pos_gal_fp[:,2]

xyz_FP = np.vstack((x_FP,y_FP,z_FP)).T
xyz_shuff = np.vstack((x,y,z)).T
xyz_DM = np.vstack((x_DM,y_DM,z_DM)).T

if want_fp_pos:
    ext = 'pos'
else:
    ext = 'peak'
np.save("Lensing/data_2dhod_"+ext+"/"+proxy+"_"+opt+"_gals.npy",xyz_shuff)
np.save("Lensing/data_2dhod_"+ext+"/true_gals.npy",xyz_DM)
#quit()

# REMOVE for new proxies
bin_centers, a,b,c,d,xi_fid_shuff,xi_fid_shuff_err,rat_fid_shuff,rat_fid_shuff_err = np.loadtxt('/home/boryanah/lars/LSSIllustrisTNG/figs/paper/Corr_shuff_'+proxy+'.txt',unpack=True)

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
            xyz_FP_jack = xyz_FP.copy()
            xyz_shuff_jack = xyz_shuff.copy()
            w_shuff_jack = w.copy()
            xyz_DM_jack = xyz_DM.copy()
            w_DM_jack = w_DM.copy()
            
            xyz = np.array([i_x,i_y,i_z],dtype=int)
            size = Lboxs[i]/N_dim

            bool_arr = np.prod((xyz == (xyz_FP/size).astype(int)),axis=1).astype(bool)
            xyz_FP_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_FP_jack = xyz_FP_jack[np.sum(xyz_FP_jack,axis=1)!=0.]

            bool_arr = np.prod((xyz == (xyz_DM/size).astype(int)),axis=1).astype(bool)
            xyz_DM_jack[bool_arr] = np.array([0.,0.,0.])
            w_DM_jack = w_DM[np.sum(xyz_DM_jack,axis=1)!=0.]
            xyz_DM_jack = xyz_DM_jack[np.sum(xyz_DM_jack,axis=1)!=0.]
            

            bool_arr = np.prod((xyz == (xyz_shuff/size).astype(int)),axis=1).astype(bool)
            xyz_shuff_jack[bool_arr] = np.array([0.,0.,0.])
            w_shuff_jack = w[np.sum(xyz_shuff_jack,axis=1)!=0.]
            xyz_shuff_jack = xyz_shuff_jack[np.sum(xyz_shuff_jack,axis=1)!=0.]
            
            
            x_DM_jack = xyz_DM_jack[:,0]
            y_DM_jack = xyz_DM_jack[:,1]
            z_DM_jack = xyz_DM_jack[:,2]
            x_shuff_jack = xyz_shuff_jack[:,0]
            y_shuff_jack = xyz_shuff_jack[:,1]
            z_shuff_jack = xyz_shuff_jack[:,2]
            x_FP_jack = xyz_FP_jack[:,0]
            y_FP_jack = xyz_FP_jack[:,1]
            z_FP_jack = xyz_FP_jack[:,2]

            #print('xs = ',len(x_FP_jack))
            #print('xs = ',len(x_shuff_jack))
            
            #bins = np.logspace(-1,1,N_bin)
            #bins = np.logspace(-1,1.3,N_bin)
            bins = np.logspace(-1,1.5,N_bin)
            bin_centers = (bins[:-1] + bins[1:])/2.
            
            results_DM = Corrfunc.theory.xi(X=x_DM_jack,Y=y_DM_jack,Z=z_DM_jack,boxsize=Lboxs[i],nthreads=16,binfile=bins, weights=w_DM_jack, weight_type='pair_product')
            
            #plt.plot(bin_centers, results['xi']*bin_centers**2,cs[i]+'-',linewidth=2.,label=str(res[i])+ ' DM')
            Corr_f_DM[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_DM['xi']
    
            results_shuff = Corrfunc.theory.xi(X=x_shuff_jack,Y=y_shuff_jack,Z=z_shuff_jack,boxsize=Lboxs[i],nthreads=16,binfile=bins, weights=w_shuff_jack, weight_type='pair_product')

            rat_xi = results_shuff['xi']/results_DM['xi']
            #rat_xi = results_shuff['xi']/(xi_fid_shuff/bin_centers**2)
            #rat_xi_shuff = (xi_fid_shuff/bin_centers**2)/results_DM['xi'] 
            plt.plot(bin_centers,rat_xi,linewidth=2.,label=opt)
            
            
            Rat_DM[:,i_x+N_dim*i_y+N_dim**2*i_z] = rat_xi
            if opt == 'shuff':
                # REMOVE for new proxies
                #Rat_HOD[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_DM['xi']/results_shuff['xi']
                Rat_HOD[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_DM['xi']/(xi_fid_shuff/bin_centers**2)
            else:
                Rat_HOD[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_shuff['xi']/(xi_fid_shuff/bin_centers**2)
            Corr_f_DM_shuff[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_shuff['xi']
            
            results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,
                                         boxsize=Lboxs[i], nthreads=16,
                                         binfile=bins)
            #plt.plot(bin_centers, results['xi']*bin_centers**2,'r',linewidth=2.,label=str(res[i])+ ' FP')
            Corr_f_FP[:,i_x+N_dim*i_y+N_dim**2*i_z] = results_FP['xi']

            #if True:
            if np.product(xyz == np.array([N_dim-1,N_dim-1,N_dim-1])).astype(bool):
            
                plt.xscale('log')
                plt.ylabel(r'$\xi(r) r^2$')
                plt.xlabel(r'$r$ [Mpc/h]')
                plt.legend()#(loc='best')
                plt.ylim([0.8,1.5])
                #plt.ylim([0,80])
                if binning == 'range':
                    plt.text(0.1, .85, str(int(percent_thresh*100))+'% mass bin', dict(size=15))
                elif binning == 'std':
                    plt.text(0.1, .85, str(int(percent_thresh*100))+'% mass fluct', dict(size=15))
                if fix_1bin == True:
                    plt.legend(loc='lower right')
                    plt.text(0.1, 1.45, "mass range: %10.5E - %10.5E, bin size: %4i" %(mass_beg,mass_end,N_diff), dict(size=12))
                    plt.close()#plt.savefig('/home/boryanah/lars/LSSIllustrisTNG/figs/many/'+figname[:-4]+'_'+str(i_bin)+'_all.png')
                else:
                    plt.savefig('/home/boryanah/lars/LSSIllustrisTNG/figs/'+figname[:-4]+'_all.png')
                plt.show()

if show_fid == True:
    bin_centers, a,b,c,d,xi_fid_shuff,xi_fid_shuff_err,rat_fid_shuff,rat_fid_shuff_err = np.loadtxt(fname_shuff,unpack=True)

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

'''
# REMOVE for new proxies
np.savetxt('/home/boryanah/lars/LSSIllustrisTNG/figs/paper/Corr_shuff_'+proxy+'.txt',(np.vstack((bin_centers, Corr_mean_FP*bin_centers**2,Corr_err_FP*bin_centers**2,Corr_mean_DM*bin_centers**2,Corr_err_DM*bin_centers**2,Corr_mean_DM_shuff*bin_centers**2,Corr_err_DM_shuff*bin_centers**2,Rat_HOD_mean,Rat_HOD_err)).T))
quit()
'''

np.savetxt('/home/boryanah/lars/LSSIllustrisTNG/'+figname[:-4]+'.txt',(np.vstack((bin_centers, Corr_mean_FP*bin_centers**2,Corr_err_FP*bin_centers**2,Corr_mean_DM*bin_centers**2,Corr_err_DM*bin_centers**2,Corr_mean_DM_shuff*bin_centers**2,Corr_err_DM_shuff*bin_centers**2,Rat_HOD_mean,Rat_HOD_err)).T))
fs = (12,8)#(8.2/1.4,6.5/1.4)#(8.2,6.5)
fig = plt.figure(figsize=fs)
plt.errorbar(bin_centers, Corr_mean_FP*bin_centers**2,yerr=Corr_err_FP*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label=str(res[i])+ ' FP')
#tukplt.errorbar(bin_centers, Corr_mean_DM*bin_centers**2,yerr=Corr_err_DM*bin_centers**2,lw=2.,ls='-',c='orange',fmt='o',capsize=2,label=str(res[i])+ ' DM')
#tukplt.errorbar(bin_centers, Corr_mean_FP*bin_centers**2,yerr=Corr_err_FP*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label=str('TNG FP gals'))
plt.errorbar(bin_centers, Corr_mean_DM*bin_centers**2,yerr=Corr_err_DM*bin_centers**2,lw=2.,ls='-',c='orange',fmt='o',capsize=2,label=str('TNG DM gals'))
plt.errorbar(bin_centers, Corr_mean_DM_shuff*bin_centers**2,yerr=Corr_err_DM_shuff*bin_centers**2,lw=2.,ls='-',c='silver',fmt='o',capsize=2,label='TNG DM env')#pdf

if show_fid == True:
    print('bla')
    #tukplt.plot(bin_centers,xi_fid_shuff,linewidth=1.,color='#1B2ACC',label='fid shuff')
    #tukplt.fill_between(bin_centers, xi_fid_shuff+xi_fid_shuff_err, xi_fid_shuff-xi_fid_shuff_err,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
#tukplt.title("Bin numbers: "+str(i_bin)+"-"+str(i_bin+skip), dict(size=15))

plt.xscale('log')
plt.ylabel(r'$\xi(r) r^2$')
plt.xlabel(r'$r$ [Mpc/h]')
plt.legend()
#tukplt.ylim([0,80])
plt.ylim([0,100])
if binning == 'range':
    plt.text(0.1, 6., str(int(percent_thresh*100))+'% mass bin', dict(size=15))
elif binning == 'std':
    plt.text(0.1, 6., str(int(percent_thresh*100))+'% mass fluct', dict(size=15))
if fix_1bin == True:
    plt.legend(loc='lower right')
    #tukplt.text(0.1, 74, "mass range: %10.5E - %10.5E" %(mass_beg,mass_end), dict(size=15))
    #tukplt.savefig('/home/boryanah/lars/LSSIllustrisTNG/'+figname[:-4]+'.pdf')
    #plt.savefig('/home/boryanah/lars/LSSIllustrisTNG/'+figname[:-4]+'.png')
else:
    plt.savefig('/home/boryanah/lars/LSSIllustrisTNG/figs/'+figname)
#plt.show()
plt.close()

fig = plt.figure(figsize=fs)#fig = plt.figure(figsize=(7,5.))
plt.plot(bin_centers,np.ones(len(bin_centers)),lw=1.5,ls='--',alpha=0.2)
#if fix_1bin == False:
if show_fid == True:
    print('bla')
    #tukplt.plot(bin_centers,rat_fid_shuff,linewidth=1.,color='#1B2ACC',label='fid shuff')
    #tukplt.fill_between(bin_centers, rat_fid_shuff+rat_fid_shuff_err, rat_fid_shuff-rat_fid_shuff_err,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')


plt.errorbar(bin_centers,Rat_DM_mean,yerr=Rat_DM_err,lw=2.,ls='-',c='orange',fmt='o',capsize=4,label=opt)
plt.xscale('log')
#tukplt.ylabel(r'$\xi(r)_{\rm opt} / \xi(r)_{\rm DM}$',dict(size=15))
plt.ylabel(r'$\xi(r)_{\rm TNG, shuff} / \xi(r)_{\rm TNG}$')
#pdf
plt.xlabel(r'$r$ [Mpc/h]')
#tukplt.legend()

if exclude_most_mass == False:
    #plt.ylim([0.8,1.5])
    plt.ylim([0.8,1.5])
    if binning == 'range':
        plt.text(0.1, .85, str(int(percent_thresh*100))+'% mass bin', dict(size=15))
    elif binning == 'std':
        plt.text(0.1, .85, str(int(percent_thresh*100))+'% mass fluct', dict(size=15))

if fix_1bin == True:
    #plt.title("Bin: "+str(i_bin)+"-"+str(i_bin+skip)+", size: "+ str(nbin), dict(size=8))
    #tukplt.title("Bin numbers: "+str(i_bin)+"-"+str(i_bin+skip), dict(size=15))
    print(repr(np.array(nbin)))
    #tukplt.legend(loc='upper right')
    if exclude_most_mass == True:
        if i_bin == 50:
            plt.ylim([0.15,3.]);
            plt.text(0.1, 0.23, str(int(percent_thresh*100))+'% mass bin', dict(size=15))
            #tukplt.text(0.1, 2.87, ("mass range: %8.3E - %8.3E" %(mass_beg,mass_end)), dict(size=15))
        else:
            # EI TUKA
            plt.ylim([0.7,1.3]);
            
            #pdf
            plt.text(0.1, 0.7, str(int(percent_thresh*100))+'% mass bin', dict(size=15))
            
            #pdf
            #tukplt.text(0.1, 1.55, ("mass range: %8.3E - %8.3E" %(mass_beg,mass_end)), dict(size=15))
    else:
        print('bla')
        #tukplt.text(0.1, 1.45, ("mass range: %8.3E - %8.3E" %(mass_beg,mass_end)), dict(size=15))
    #tukplt.savefig('/home/boryanah/lars/LSSIllustrisTNG/'+figname[:-4]+'_rat.pdf')
    plt.savefig('/home/boryanah/lars/LSSIllustrisTNG/'+figname[:-4]+'_rat.png')
    #plt.show()
    plt.close()
else:
    plt.savefig('/home/boryanah/lars/LSSIllustrisTNG/figs/'+figname[:-4]+'_rat.png')
    #plt.show()
    plt.close()

print('________________________________________')
print('              STATISTICS')
print('________________________________________')

print(bin_centers)

print('mass proxy = ',proxy)
print('secondary parameter = ',opt)
N_beg = int(N_bin/2-3)+5
N_end = N_bin-4-2
print('distance range [Mpc/h] = ',bin_centers[N_beg],'-',bin_centers[N_end])

# it does make sense for shuff cause it is actually hydrosim and not shuff and therefore constant
# only for that case do you divide by a different factor

print('percentage difference = ',(1-np.mean(Rat_DM_mean[N_beg:N_end]))*100.,'%')
print('cumulative error = ',np.sqrt(np.sum((Rat_DM_err[N_beg:N_end])**2)))
print('percentage difference wrt HOD = ',(1-np.mean(Rat_HOD_mean[N_beg:N_end]))*100.,'%')
print('cumulative error wrt HOD = ',np.sqrt(np.sum((Rat_HOD_err[N_beg:N_end])**2)))
print('________________________________________')
quit()

# bin centers
bins = np.logspace(-1,1,20)
bin_centers = (bins[:-1] + bins[1:])/2.

# positions of centrals
inds_fp = np.arange(N_halo_fp)
inds_sub_fp = np.arange(N_sub_fp)
ind_gal_b1 = inds_fp[count_halo_fp >= 1]
ind_fsub_gal = GroupFirstSub_fp[ind_gal_b1]
xyz_central = SubhaloPos_fp[ind_fsub_gal]


# positions of all galaxies
xyz_FP = np.vstack((x_FP,y_FP,z_FP)).T
xyz_FP_jack = xyz_FP.copy()
x_FP_jack = xyz_FP_jack[:,0]
y_FP_jack = xyz_FP_jack[:,1]
z_FP_jack = xyz_FP_jack[:,2]

# for comparing with and without satellites
results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,
                                boxsize=Lboxs[i], nthreads=16,
                                binfile=bins)

plt.plot(bin_centers, results_FP['xi']*bin_centers**2,'r',linewidth=2.,label='all FP '+str(st_cut))

x_FP = pos_gal_fp_2[:,0]
y_FP = pos_gal_fp_2[:,1]
z_FP = pos_gal_fp_2[:,2]

# positions of all galaxies
xyz_FP = np.vstack((x_FP,y_FP,z_FP)).T
xyz_FP_jack = xyz_FP.copy()
x_FP_jack = xyz_FP_jack[:,0]
y_FP_jack = xyz_FP_jack[:,1]
z_FP_jack = xyz_FP_jack[:,2]

# for comparing with and without satellites
results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,
                                boxsize=Lboxs[i], nthreads=16,
                                binfile=bins)

plt.plot(bin_centers, results_FP['xi']*bin_centers**2,'r',linewidth=2.,label='all FP '+str(2*st_cut))

# positions of all central galaxies
xyz_FP_jack = xyz_central.copy()
x_FP_jack = xyz_FP_jack[:,0]
y_FP_jack = xyz_FP_jack[:,1]
z_FP_jack = xyz_FP_jack[:,2]
print('N central = ',len(x_FP_jack))

# for comparing with and without satellites
results_FP = Corrfunc.theory.xi(X=x_FP_jack,Y=y_FP_jack,Z=z_FP_jack,
                                boxsize=Lboxs[i], nthreads=16,
                                binfile=bins)

plt.plot(bin_centers, results_FP['xi']*bin_centers**2,'b',linewidth=2.,label='centrals FP')
plt.xscale('log')
plt.ylabel(r'$\xi(r) r^2$')
plt.xlabel(r'$r$ [Mpc/h]')
plt.legend()
plt.ylim([0,80])
plt.show()


if False:
    if False:
        if False:
            if False:
                if False:
                    '''
                    if len(opt_sorted_sort)>100:
                    
                    plt.scatter(opt_sorted_sort[i_opt_sss],gal_count_sort)#  is how the elements are ordered in space of N_sat sorted
                    plt.show()
                    '''



                '''
                # for testing of whether preferentially put into lighter halos i_dm_i = np.arange(len(dm_i)) and shuffle that 
                # and after EMPTY print(np.mean(np.array(k)));
                m = (np.mean(gal_count*np.arange(len(mass))[::-1])/np.mean(gal_count[i_dm_i]*np.arange(len(mass[i_dm_i]))[::-1]))
                if np.isnan(m) != True:
                    k.append(m)
                print("expect larger than 1 median mean",k[-1]) # YES WE DO
                '''
