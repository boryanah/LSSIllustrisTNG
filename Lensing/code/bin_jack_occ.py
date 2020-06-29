import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Corrfunc
from Corrfunc.theory.DD import DD
import plotparams
plotparams.default()
import sys
#from mpi4py import MPI
from halotools.utils import randomly_downsample_data

plt.rc('text', usetex=True)

n_thread = 16
periodic = True
dir_part = "/mnt/gosling1/boryanah/TNG300/"
part_fn = "parts_position_tng300-3_99.npy"
proxy = "m200m"
opt = sys.argv[1]#"partial_fenv"#"partial_s2r"#"shuff"#"spin"#"shuff"#"conc"#"vani"#"shuff"#"env"#"partial_vani""partial_s2r""vdisp"

Lbox = 205. # in Mpc/h
pcle_mass = 3.03e9 #M_solar/h
# fiducial means are you looking at true galaxies or not
if len(sys.argv) > 2:
    fiducial = int(sys.argv[2])#True
else:
    fiducial = False
    
test = False
if test: test_name = "test"
else: test_name = '-'.join(opt.split('_'));print(test_name)#proxy+"_"+opt
if test:
    Npts = int(1e5)
    x = np.random.uniform(0, Lbox, Npts)
    y = np.random.uniform(0, Lbox, Npts)
    z = np.random.uniform(0, Lbox, Npts)
    pos_m = np.vstack((x, y, z)).T
else:
    pos_m = np.load(dir_part+part_fn)/1000. # in Mpc/h

if test:
    Npts = int(1.2e4)
    x = np.random.uniform(0, Lbox, Npts)
    y = np.random.uniform(0, Lbox, Npts)
    z = np.random.uniform(0, Lbox, Npts)
    pos_g = np.vstack((x, y, z)).T
else:
    pos_g = np.load("../true_gals.npy")

if test:
    Npts = int(1.2e4)
    x = np.random.uniform(0, Lbox, Npts)
    y = np.random.uniform(0, Lbox, Npts)
    z = np.random.uniform(0, Lbox, Npts)
    pos_g_opt = np.vstack((x, y, z)).T
else:
    pos_g_opt = np.load("../"+proxy+"_"+opt+"_gals.npy")

if fiducial == False:
    pos_g = pos_g_opt.copy()

N_m = pos_m.shape[0]
down = 5000
N_m_down = N_m//down
print("Number of pcles = ",N_m_down)
print("Downsampling...")
try:
    pos_m = np.load("../pos_m_down_"+str(down)+".npy")
except:
    #x_m = pos_m[::down,0]
    #y_m = pos_m[::down,1]
    #z_m = pos_m[::down,2]
    #pos_m = np.vstack((x_m,y_m,z_m)).T

    pos_m = randomly_downsample_data(pos_m, N_m_down)

    print(pos_m.shape[0])
    np.save("../pos_m_down_"+str(down)+".npy",pos_m)
    plt.scatter(pos_m[:,0],pos_m[:,1],s=0.01)
    plt.show()
print("Downsampled O.o!")

down_fac = N_m/N_m_down
pcle_mass *= down_fac

N_g = pos_g.shape[0]
N_g_opt = pos_g.shape[0]
print("Number of gals = ",N_g)

pos_m = pos_m.astype(np.float32)
X_m = pos_m[:,0]
Y_m = pos_m[:,1]
Z_m = pos_m[:,2]

X_g = pos_g[:,0]
Y_g = pos_g[:,1]
Z_g = pos_g[:,2]

# Random points
print("Getting random data...")
try:
    X_r_m = np.load("../X_r_m.npy")
    Y_r_m = np.load("../Y_r_m.npy")
    Z_r_m = np.load("../Z_r_m.npy")
    X_r_g = np.load("../X_r_g.npy")
    Y_r_g = np.load("../Y_r_g.npy")
    Z_r_g = np.load("../Z_r_g.npy")
    N_r_m = len(X_r_m)
    N_r_g = len(X_r_g)
    print("N_r_m = ",N_r_m)
except:
    N_r_m = N_m_down*35#5 # 5 suff on large scales
    N_r_g = N_g*35#5 #35 # 5 suff on large scales
    X_r_m = np.random.uniform(0.,Lbox,N_r_m).astype(np.float32)
    Y_r_m = np.random.uniform(0.,Lbox,N_r_m).astype(np.float32)
    Z_r_m = np.random.uniform(0.,Lbox,N_r_m).astype(np.float32)
    X_r_g = np.random.uniform(0.,Lbox,N_r_g).astype(np.float32)
    Y_r_g = np.random.uniform(0.,Lbox,N_r_g).astype(np.float32)
    Z_r_g = np.random.uniform(0.,Lbox,N_r_g).astype(np.float32)
    np.save("../X_r_m.npy",X_r_m)
    np.save("../Y_r_m.npy",Y_r_m)
    np.save("../Z_r_m.npy",Z_r_m)
    np.save("../X_r_g.npy",X_r_g)
    np.save("../Y_r_g.npy",Y_r_g)
    np.save("../Z_r_g.npy",Z_r_g)
    
print("Got random data o.O!")


test_g = False # Has been tested and works!
if test_g:
    X_r_m = X_r_g.copy()
    Y_r_m = Y_r_g.copy()
    Z_r_m = Z_r_g.copy()
    X_m = X_g.copy()
    Y_m = Y_g.copy()
    Z_m = Z_g.copy()
test_m = False
if test_m:
    X_r_g = X_r_m.copy()
    Y_r_g = Y_r_m.copy()
    Z_r_g = Z_r_m.copy()
    X_g = X_m.copy()
    Y_g = Y_m.copy()
    Z_g = Z_m.copy()

N_bin = 16
N_dim = 3
lb_min = -0.7
lb_max = 1.2
bins = np.logspace(lb_min,lb_max,N_bin)
x_lim_min = 10**(lb_min)
x_lim_max = 10**(lb_max)
bin_centers = (bins[:-1] + bins[1:])/2.
eps = bin_centers/50.
eps_m = -bin_centers/50.

xyz_m = np.vstack((X_m,Y_m,Z_m)).T
xyz_g = np.vstack((X_g,Y_g,Z_g)).T
xyz_r_m = np.vstack((X_r_m,Y_r_m,Z_r_m)).T
xyz_r_g = np.vstack((X_r_g,Y_r_g,Z_r_g)).T
bias = np.zeros((N_bin-1,N_dim**3))
corr_coeff = np.zeros((N_bin-1,N_dim**3))
corr_g = np.zeros((N_bin-1,N_dim**3))
corr_m = np.zeros((N_bin-1,N_dim**3))
corr_gm = np.zeros((N_bin-1,N_dim**3))
# JACKKNIFE ERROR ESTIMATION
for i_x in range(N_dim):
    for i_y in range(N_dim):
        for i_z in range(N_dim):
            xyz_g_jack = xyz_g.copy()
            xyz_m_jack = xyz_m.copy()
            xyz_r_g_jack = xyz_r_g.copy()
            xyz_r_m_jack = xyz_r_m.copy()
            
            xyz = np.array([i_x,i_y,i_z],dtype=int)
            size = Lbox/N_dim

            bool_arr = np.prod((xyz == (xyz_g/size).astype(int)),axis=1).astype(bool)
            xyz_g_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_g_jack = xyz_g_jack[np.sum(xyz_g_jack,axis=1)!=0.]

            bool_arr = np.prod((xyz == (xyz_m/size).astype(int)),axis=1).astype(bool)
            xyz_m_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_m_jack = xyz_m_jack[np.sum(xyz_m_jack,axis=1)!=0.]

            bool_arr = np.prod((xyz == (xyz_r_g/size).astype(int)),axis=1).astype(bool)
            xyz_r_g_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_r_g_jack = xyz_r_g_jack[np.sum(xyz_r_g_jack,axis=1)!=0.]

            bool_arr = np.prod((xyz == (xyz_r_m/size).astype(int)),axis=1).astype(bool)
            xyz_r_m_jack[bool_arr] = np.array([0.,0.,0.])
            xyz_r_m_jack = xyz_r_m_jack[np.sum(xyz_r_m_jack,axis=1)!=0.]
            
            
            X_jack_m = xyz_m_jack[:,0]
            Y_jack_m = xyz_m_jack[:,1]
            Z_jack_m = xyz_m_jack[:,2]
            
            X_jack_g = xyz_g_jack[:,0]
            Y_jack_g = xyz_g_jack[:,1]
            Z_jack_g = xyz_g_jack[:,2]

            X_jack_r_m = xyz_r_m_jack[:,0]
            Y_jack_r_m = xyz_r_m_jack[:,1]
            Z_jack_r_m = xyz_r_m_jack[:,2]

            X_jack_r_g = xyz_r_g_jack[:,0]
            Y_jack_r_g = xyz_r_g_jack[:,1]
            Z_jack_r_g = xyz_r_g_jack[:,2]


            # Cross-correlation for gm
            print("Nightmare is starting")
            autocorr = 0
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_g, Y1=Y_jack_g, Z1=Z_jack_g,
                         X2=X_jack_m, Y2=Y_jack_m, Z2=Z_jack_m,
                         boxsize=Lbox,periodic=periodic)

            DD_gm = results['npairs'].astype(float)
            DD_gm /= (N_g*1.*N_m_down)



            autocorr = 0
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_g, Y1=Y_jack_g, Z1=Z_jack_g,
                         X2=X_jack_r_m, Y2=Y_jack_r_m, Z2=Z_jack_r_m,
                         boxsize=Lbox,periodic=periodic)


            DR_gm = results['npairs'].astype(float)
            DR_gm /= (N_g*1.*N_r_m)

            autocorr = 0
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_r_g, Y1=Y_jack_r_g, Z1=Z_jack_r_g,
                         X2=X_jack_m, Y2=Y_jack_m, Z2=Z_jack_m,
                         boxsize=Lbox,periodic=periodic)

            RD_gm = results['npairs'].astype(float)
            RD_gm /= (N_r_g*1.*N_m_down)

            autocorr = 0
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_r_g, Y1=Y_jack_r_g, Z1=Z_jack_r_g,
                         X2=X_jack_r_m, Y2=Y_jack_r_m, Z2=Z_jack_r_m,
                         boxsize=Lbox,periodic=periodic)


            RR_gm = results['npairs'].astype(float)
            RR_gm /= (N_r_g*1.*N_r_m)

            Corr_gm = (DD_gm-DR_gm-RD_gm+RR_gm)/RR_gm

            # Corr_g
            autocorr = 1
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_g, Y1=Y_jack_g, Z1=Z_jack_g,
                         boxsize=Lbox,periodic=periodic)

            DD_gg = results['npairs'].astype(float)
            DD_gg /= (N_g*1.*N_g)
            
            autocorr = 1
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_r_g, Y1=Y_jack_r_g, Z1=Z_jack_r_g,
                         boxsize=Lbox,periodic=periodic)

            RR_gg = results['npairs'].astype(float)
            RR_gg /= (N_r_g*1.*N_r_g)

            Corr_g = DD_gg/RR_gg-1.
            
            # Corr_m
            autocorr = 1
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_m, Y1=Y_jack_m, Z1=Z_jack_m,
                         boxsize=Lbox,periodic=periodic)

            DD_mm = results['npairs'].astype(float)
            DD_mm /= (N_m_down*1.*N_m_down)

            autocorr = 1
            results = DD(autocorr,nthreads=n_thread,binfile=bins,
                         X1=X_jack_r_m, Y1=Y_jack_r_m, Z1=Z_jack_r_m,
                         boxsize=Lbox,periodic=periodic)

            RR_mm = results['npairs'].astype(float)
            RR_mm /= (N_r_m*1.*N_r_m)
            
            Corr_m = DD_mm/RR_mm-1.
            
            bias[:,i_x+N_dim*i_y+N_dim**2*i_z] = np.sqrt(Corr_g/Corr_m)
            corr_coeff[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_gm/np.sqrt(Corr_g*Corr_m) 
            corr_g[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_g
            corr_m[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_m
            corr_gm[:,i_x+N_dim*i_y+N_dim**2*i_z] = Corr_gm


bias_mean = np.mean(bias,axis=1)
bias_error = np.sqrt(N_dim**3-1)*np.std(bias,axis=1)
corr_coeff_mean = np.mean(corr_coeff,axis=1)
corr_coeff_error = np.sqrt(N_dim**3-1)*np.std(corr_coeff,axis=1)
Corr_gm_mean = np.mean(corr_gm,axis=1)
Corr_gm_error = np.sqrt(N_dim**3-1)*np.std(corr_gm,axis=1)
Corr_g_mean = np.mean(corr_g,axis=1)
Corr_g_error = np.sqrt(N_dim**3-1)*np.std(corr_g,axis=1)
Corr_m_mean = np.mean(corr_m,axis=1)
Corr_m_error = np.sqrt(N_dim**3-1)*np.std(corr_m,axis=1)


if fiducial:
    np.save("../data_gm/bin_fid.npy",bin_centers)
    np.save("../data_gm/bias_mean_fid.npy",bias_mean)
    np.save("../data_gm/corr_coeff_mean_fid.npy",corr_coeff_mean)
    np.save("../data_gm/corr_gm_mean_fid.npy",Corr_gm_mean)
    np.save("../data_gm/corr_gg_mean_fid.npy",Corr_g_mean)
    np.save("../data_gm/bias_error_fid.npy",bias_error)
    np.save("../data_gm/corr_coeff_error_fid.npy",corr_coeff_error)
    np.save("../data_gm/corr_gm_error_fid.npy",Corr_gm_error)
    np.save("../data_gm/corr_gg_error_fid.npy",Corr_g_error)
else:

    np.save("../data_gm/bin_cents.npy",bin_centers)
    np.save("../data_gm/bias_mean_"+opt+".npy",bias_mean)
    np.save("../data_gm/corr_coeff_mean_"+opt+".npy",corr_coeff_mean)
    np.save("../data_gm/corr_gm_mean_"+opt+".npy",Corr_gm_mean)
    np.save("../data_gm/corr_gg_mean_"+opt+".npy",Corr_g_mean)
    np.save("../data_gm/bias_error_"+opt+".npy",bias_error)
    np.save("../data_gm/corr_coeff_error_"+opt+".npy",corr_coeff_error)
    np.save("../data_gm/corr_gm_error_"+opt+".npy",Corr_gm_error)
    np.save("../data_gm/corr_gg_error_"+opt+".npy",Corr_g_error)
    
    Bin_fid = np.load("../data_gm/bin_fid.npy")
    Bias_mean_fid = np.load("../data_gm/bias_mean_fid.npy")
    Corr_coeff_mean_fid = np.load("../data_gm/corr_coeff_mean_fid.npy")
    Corr_gm_mean_fid = np.load("../data_gm/corr_gm_mean_fid.npy")
    Corr_g_mean_fid = np.load("../data_gm/corr_gg_mean_fid.npy")
    Bias_error_fid = np.load("../data_gm/bias_error_fid.npy")
    Corr_coeff_error_fid = np.load("../data_gm/corr_coeff_error_fid.npy")
    Corr_gm_error_fid = np.load("../data_gm/corr_gm_error_fid.npy")
    Corr_g_error_fid = np.load("../data_gm/corr_gg_error_fid.npy")


quit()
size = (12,8)
plt.figure(1,figsize=size)
plt.title(r'Bias')
if fiducial:
    plt.errorbar(bin_centers,bias_mean,yerr=bias_error,linewidth=2.,fmt='o',capsize=2,label="gg/mm")
else:
    eb_f = plt.errorbar(Bin_fid,Bias_mean_fid,yerr=Bias_error_fid,linewidth=2.,fmt='d',capsize=2,color='dodgerblue',ls='--',label='fiducial')
    eb_f[-1][0].set_linestyle('--')
    eb = plt.errorbar(bin_centers+eps,bias_mean,yerr=bias_error,fmt='o',capsize=2,color='dodgerblue',linewidth=2.,ls='-',label=test_name)
    #eb[-1][0].set_linestyle('-')
plt.plot(bin_centers,np.ones(N_bin-1),'k--',linewidth=2.)
plt.xscale('log')
plt.ylim([0,1.8])
plt.xlim([x_lim_min,x_lim_max])
plt.legend()
plt.ylabel(r'$b(r)=(\xi_{gg}(r)/\xi_{mm}(r))^{\frac{1}{2}}$')
plt.xlabel(r'$r$ [Mpc/h]')
plt.savefig("figs/gg_to_mm_ratio.png")

plt.figure(3,figsize=size)
plt.title(r'Correlation coeffcient')
if fiducial:
    plt.errorbar(bin_centers,corr_coeff_mean,yerr=corr_coeff_error,fmt='o',capsize=2,linewidth=2.,label=r'gm*gm/(gg*mm)')
else:
    eb_f = plt.errorbar(Bin_fid,Corr_coeff_mean_fid,yerr=Corr_coeff_error_fid,fmt='d',capsize=2,linewidth=2.,color='dodgerblue',ls='--',label='fiducial')
    eb_f[-1][0].set_linestyle('--')
    eb = plt.errorbar(bin_centers+eps,corr_coeff_mean,yerr=corr_coeff_error,fmt='o',capsize=2,color='dodgerblue',linewidth=2.,ls='-',label=test_name)
    #eb[-1][0].set_linestyle('-')
plt.plot(bin_centers,np.ones(N_bin-1),'k--',linewidth=2.)
plt.xscale('log')
plt.ylim([0.5,1.2])
plt.xlim([x_lim_min,x_lim_max])
plt.legend()
plt.ylabel(r'$r(r)=\xi_{gm}(r)/(\xi_{gg}(r) \xi_{mm}(r))^{\frac{1}{2}}$')
plt.xlabel(r'$r$ [Mpc/h]')
plt.savefig("figs/corr_coeff.png")


plt.figure(2,figsize=size)
plt.title("Correlation function")
if fiducial:
    plt.errorbar(bin_centers,Corr_m_mean*bin_centers**2,yerr=Corr_m_error*bin_centers**2,linewidth=2.,fmt='o',capsize=2,label='mm')
    plt.errorbar(bin_centers,Corr_g_mean*bin_centers**2,yerr=Corr_g_error*bin_centers**2,linewidth=2.,fmt='o',capsize=2,label='gg')
    plt.errorbar(bin_centers,Corr_gm_mean*bin_centers**2,yerr=Corr_gm_error*bin_centers**2,linewidth=2.,fmt='o',capsize=2,label='gm')
else:
    eb_f = plt.errorbar(Bin_fid,Corr_gm_mean_fid*Bin_fid**2,yerr=Corr_gm_error_fid*bin_centers**2,linewidth=2.,color='darkviolet',ls='--',fmt='d',capsize=2,label='gm fiducial')
    eb_f[-1][0].set_linestyle('--')
    eb_f = plt.errorbar(Bin_fid,Corr_g_mean_fid*Bin_fid**2,yerr=Corr_g_error_fid*bin_centers**2,linewidth=2.,color='dodgerblue',ls='--',fmt='d',capsize=2,label='gg fiducial')
    eb_f[-1][0].set_linestyle('--')
    eb = plt.errorbar(bin_centers+eps,Corr_g_mean*bin_centers**2,yerr=Corr_g_error*bin_centers**2,linewidth=2.,color='dodgerblue',fmt='o',capsize=2,ls='-',label="gg "+test_name)
    #eb[-1][0].set_linestyle('-')
    eb = plt.errorbar(bin_centers+eps,Corr_gm_mean*bin_centers**2,yerr=Corr_gm_error*bin_centers**2,linewidth=2.,color='darkviolet',fmt='o',capsize=2,ls='-',label="gm "+test_name)
    eb = plt.errorbar(bin_centers+eps_m,Corr_m_mean*bin_centers**2,yerr=Corr_m_error*bin_centers**2,linewidth=2.,color='orange',fmt='o',capsize=2,ls='-',label='mm')
plt.xscale('log')
#plt.ylim([-1,15])
plt.xlim([x_lim_min,x_lim_max])
plt.ylabel(r'$\xi(r) r^2$')
plt.xlabel(r'$r$ [Mpc/h]')
plt.legend()
plt.savefig("figs/corr_funcs.png")
plt.show()

if True:
    '''
    # Auto-correlation of gg
    results_g = Corrfunc.theory.xi(X=X_jack_g,Y=Y_jack_g,Z=Z_jack_g,
    boxsize=Lbox, nthreads=n_thread,
    binfile=bins)
    
    Corr_g = results_g['xi']
    
    # Auto-correlation of mm
    results_m = Corrfunc.theory.xi(X=X_jack_m,Y=Y_jack_m,Z=Z_jack_m,
    boxsize=Lbox,nthreads=n_thread,
    binfile=bins)
    
    Corr_m = results_m['xi']
    '''
