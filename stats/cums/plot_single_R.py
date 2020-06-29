import numpy as np
import matplotlib.pyplot as plt
import plotparams
plotparams.buba()

proxy = 'm200m'

types = ['second moment','third moment']

bin_centers = np.load("data_many/rs.npy")

# if not ratios
second_fid = np.load("data_many/second_true_mean.npy")
second_error_fid = np.load("data_many/second_true_err.npy")
third_fid = np.load("data_many/third_true_mean.npy")
third_error_fid = np.load("data_many/third_true_err.npy")

fig,axs = plt.subplots(2,1,figsize=(6.5,11))
for i in range(2):
    lab_fid = "``true'' galaxies"
    plt.subplot(2,1,i+1)
    if i == 0:
        plt.plot(bin_centers,second_fid,linewidth=2.,color='#1B2ACC',label=lab_fid)
        plt.fill_between(bin_centers,second_fid+second_error_fid,second_fid-second_error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        plt.ylabel(r'$\langle \delta_{R,\ \textrm{true}}^2 \rangle$')
    else:
        plt.plot(bin_centers,third_fid,linewidth=2.,color='#1B2ACC',label=lab_fid)
        plt.fill_between(bin_centers,third_fid+third_error_fid,third_fid-third_error_fid,alpha=0.1, edgecolor='#1B2ACC', facecolor='#089FFF')
        plt.ylabel(r'$\langle \delta_{R,\ \textrm{true}}^3 \rangle$')

    plt.yscale('log')
    if i == 0: plt.legend()#loc='upper left')
    plt.xlim([.5,8.5])        
    plt.xticks(np.arange(int(min(bin_centers)), int(max(bin_centers))+1, 1.0))
    plt.text(0.5, 0.1,types[i], ha='center', va='center', transform=axs[i].transAxes)
    plt.xlabel(r'$R$ [Mpc/h]')

plt.savefig("single_second_third.png")
plt.show()
