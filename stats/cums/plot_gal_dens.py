import numpy as np
import scipy.spatial as spatial
import itertools
import sys
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

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
gal_dir = "/home/boryanah/lars/test/Lensing/"
test_name = '-'.join(opt.split('_'));print(test_name)
pos_g = np.load(gal_dir+"true_gals.npy")
pos_g_opt = np.load(gal_dir+proxy+"_"+opt+"_gals.npy")
    
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

slice_id = 20
Dens = D_g_smo[:,slice_id,:]
Dopt = D_g_opt_smo[:,slice_id,:]
Diff = Dens-Dopt
Dpos = Diff.copy()
Dpos[Diff<0.] = 0.
Dneg = Diff.copy()
Dneg[Diff>0.] = 0.
Dneg = np.abs(Dneg)

cm = 'gist_stern'#'prism'#'jet'#'gist_stern'#'jet'#'gist_stern'
cm2 = 'terrain'#'Spectral'#'gist_rainbow'
al = 0.4
plt.subplots(1,3,figsize=(15,5))
plt.subplot(1,3,1)
plt.imshow(Dens,cmap=cm2)
plt.imshow(Dpos,cmap=cm,alpha=al)
plt.colorbar()
plt.subplot(1,3,2)
plt.imshow(Dens,cmap=cm2)
plt.colorbar()
plt.subplot(1,3,3)
plt.imshow(Dens,cmap=cm2)
plt.imshow(Dneg,cmap=cm,alpha=al)
plt.colorbar()
#plt.savefig("figs/cums_"+proxy+"_"+opt+".png")
#plt.show()
plt.close()

nrows = 1
ncols = 2
plt.subplots(nrows,ncols,figsize=(15,5))
plt.subplot(nrows,ncols,1)
plt.imshow(Dens,cmap=cm2)
#plt.imshow(Dpos,cmap=cm,alpha=al)
plt.colorbar()
plt.subplot(nrows,ncols,2)
plt.imshow(Dens,cmap=cm2)
plt.imshow(Dneg,cmap=cm,alpha=al)
plt.colorbar()
plt.savefig("figs/cums_"+proxy+"_"+opt+".png")
plt.show()
plt.close()
