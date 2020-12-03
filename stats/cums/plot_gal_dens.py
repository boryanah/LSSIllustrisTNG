import numpy as np
import scipy.spatial as spatial
import itertools
import sys
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import h5py

import plotparams
plotparams.buba()

# do the computation for the true galaxies?
want_true = 0

# what are we displaying
proxy = "m200m"
opt = 'shuff'#sys.argv[1]#"partial_fenv"#"partial_s2r"#"shuff"#"spin"#"shuff"#"conc"#"vani"#"shuff"#"env"#"partial_vani""partial_s2r""vdisp"

# parameter choices
Lbox = 205. # in Mpc/h
bounds = np.array([Lbox,Lbox,Lbox])
N_dim = 256#int(sys.argv[2])#32*2
gr_size = Lbox/N_dim
R = 3.#float(sys.argv[3])#4

# load them galaxies
ext2 = 'data_2dhod_peak'
ext1 = 'data_2dhod_pos'
gal_dir = "/home/boryanah/lars/LSSIllustrisTNG/Lensing/"
test_name = '-'.join(opt.split('_'));print(test_name)
pos_g = np.load(gal_dir+ext1+"/"+"true_gals.npy")
pos_g_opt = np.load(gal_dir+ext2+"/"+proxy+"_"+opt+"_gals.npy")
pos_m = np.load(gal_dir+"pos_parts_down_61035156.npy")

# how mangy galaxies
N_g = pos_g.shape[0]
N_g_opt = pos_g_opt.shape[0]
print("Number of gals = ",N_g_opt)

def get_density(pos):
    g_x = pos[:,0]
    g_y = pos[:,1]
    g_z = pos[:,2]
    D, edges = np.histogramdd(np.transpose([g_x,g_y,g_z]),bins=N_dim,range=[[0,Lbox],[0,Lbox],[0,Lbox]])
    D_avg = N_g*1./N_dim**3
    D /= D_avg
    D -= 1.
    return D

D_g = get_density(pos_g)
D_g_opt = get_density(pos_g_opt)
D_m = get_density(pos_m)

def Wg(k2, r):
    return np.exp(-k2*r*r/2.)


def get_smooth_density(D,rsmo=4.):
    karr = np.fft.fftfreq(N_dim, d=Lbox/(2*np.pi*N_dim))
    dfour = np.fft.fftn(D)
    dksmo = np.zeros(shape=(N_dim, N_dim, N_dim),dtype=complex)
    ksq = np.zeros(shape=(N_dim, N_dim, N_dim),dtype=complex)
    ksq[:,:,:] = karr[None,None,:]**2+karr[None,:,None]**2+karr[:,None,None]**2
    dksmo[:,:,:] = Wg(ksq,rsmo)*dfour
    drsmo = np.real(np.fft.ifftn(dksmo))
    return drsmo

# slow
#D_g_smo = get_smooth_density(D_g,rsmo=R)
# fast
D_g_smo = gaussian_filter(D_g,R)
# slow
#D_g_opt_smo = get_smooth_density(D_g_opt,rsmo=R)
# fast
D_g_opt_smo = gaussian_filter(D_g_opt,R)
D_m_smo = gaussian_filter(D_m,1.)

slice_id = 200#20
Dens = D_g_smo[:,slice_id,:]
Dopt = D_g_opt_smo[:,slice_id,:]
Dens_m = D_m_smo[:,slice_id,:]
Diff = Dopt-Dens
Dpos = Diff.copy()
Dpos[Diff<0.] = 0.
Dneg = Diff.copy()
Dneg[Diff>0.] = 0.


cm_front = 'gist_stern'
cm_back = 'terrain'

al = 0.4

plt.imshow(np.log10(Dens_m+1),cmap='Greys')
plt.imshow(Dpos,cmap=cm_front,alpha=al)#'coolwarm'
#plt.colorbar()
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
plt.xlabel(r'$X \ [{\rm Mpc}/h]$')
plt.ylabel(r'$Y \ [{\rm Mpc}/h]$')

plt.savefig("cums_"+proxy+"_"+opt+".png")
plt.show()
plt.close()
quit()
nrows = 1
ncols = 2
plt.subplots(nrows,ncols,figsize=(15,5))
plt.subplot(nrows,ncols,1)
plt.imshow(Dens,cmap=cm_back)
plt.imshow(Dpos,cmap=cm_front,alpha=al)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
plt.xlabel(r'$X \ [{\rm Mpc}/h]$')
plt.ylabel(r'$Y \ [{\rm Mpc}/h]$')

plt.subplot(nrows,ncols,2)
plt.imshow(np.log10(Dens_m+1),cmap='Greys')
plt.imshow(Dpos,cmap=cm_front,alpha=al)#'coolwarm'
#plt.colorbar()
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
plt.xlabel(r'$X \ [{\rm Mpc}/h]$')
plt.ylabel(r'$Y \ [{\rm Mpc}/h]$')

plt.savefig("cums_"+proxy+"_"+opt+".pdf")
plt.show()
plt.close()
