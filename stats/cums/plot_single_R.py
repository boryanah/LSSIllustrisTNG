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
col_fid = 'black'
second_shuff = np.load("data_many/second_m200m_shuff_mean.npy")
second_error_shuff = np.load("data_many/second_m200m_shuff_err.npy")
third_shuff = np.load("data_many/third_m200m_shuff_mean.npy")
third_error_shuff = np.load("data_many/third_m200m_shuff_err.npy")
col_shuff = '#1B2ACC'#'dodgerblue'#
col_shuff2 = 'dodgerblue'#'#089FF'

fig,axs = plt.subplots(2,1,figsize=(8.5,11))
for i in range(2):
    lab_fid = "TNG300"
    lab_shuff = "basic HOD"
    plt.subplot(2,1,i+1)
    if i == 0:
        plt.plot(bin_centers,second_fid,linewidth=2.,color=col_fid,label=lab_fid)
        plt.fill_between(bin_centers,second_fid+second_error_fid,second_fid-second_error_fid,alpha=0.1, edgecolor=col_fid, facecolor=col_fid)
        plt.plot(bin_centers,second_shuff,linewidth=2.,color=col_shuff,label=lab_shuff)
        plt.fill_between(bin_centers,second_shuff+second_error_shuff,second_shuff-second_error_shuff,alpha=0.1, edgecolor=col_shuff, facecolor=col_shuff2)
        plt.ylabel(r'$\langle \delta_{R,\ \textrm{TNG300}}^2 \rangle$')
    else:
        plt.plot(bin_centers,third_fid,linewidth=2.,color=col_fid,label=lab_fid)
        plt.fill_between(bin_centers,third_fid+third_error_fid,third_fid-third_error_fid,alpha=0.1, edgecolor=col_fid, facecolor=col_fid)
        plt.plot(bin_centers,third_shuff,linewidth=2.,color=col_shuff,label=lab_shuff)
        plt.fill_between(bin_centers,third_shuff+third_error_shuff,third_shuff-third_error_shuff,alpha=0.1, edgecolor=col_shuff, facecolor=col_shuff2)
        plt.ylabel(r'$\langle \delta_{R,\ \textrm{TNG300}}^3 \rangle$')

    plt.yscale('log')
    if i == 0: plt.legend()#loc='upper left')
    plt.xlim([2.5,8.5])        
    plt.xticks(np.arange(int(min(bin_centers)), int(max(bin_centers))+1, 1.0))
    plt.text(0.5, 0.1,types[i], ha='center', va='center', transform=axs[i].transAxes)
    plt.xlabel(r'$R$ [Mpc/h]')

plt.savefig("single_second_third.pdf")
plt.show()
