import scipy.stats as scist
import numpy as np
import matplotlib.pyplot as plt

def compute_binned_stats(mass,r,bins):
    # compute 2d histograms
    r_median, bins, inds = scist.binned_statistic(mass,r,statistic='median',bins=bins)
    r1_low, bins, inds = scist.binned_statistic(mass,r,statistic=lambda y: np.percentile(y,q1_low),bins=bins)
    r1_high, bins, inds = scist.binned_statistic(mass,r,statistic=lambda y: np.percentile(y,q1_high),bins=bins)
    r2_low, bins, inds = scist.binned_statistic(mass,r,statistic=lambda y: np.percentile(y,q2_low),bins=bins)
    r2_high, bins, inds = scist.binned_statistic(mass,r,statistic=lambda y: np.percentile(y,q2_high),bins=bins)
    return r_median,r1_low,r1_high,r2_low,r2_high

def plot_median(M,r,logm_min,logm_max,n_bins=41):
    m_bins = np.logspace(logm_min,logm_max,n_bins)
    bin_cents = .5*(m_bins[1:]+m_bins[:-1])
    r_median,r1_low,r1_high,r2_low,r2_high = compute_binned_stats(M,r,m_bins)

    plt.plot(bin_cents,r_median,'orange')
    #plt.plot(bin_cents,r2_high,'r--',label='95%')
    #plt.plot(bin_cents,r2_low,'r--')
    plt.plot(bin_cents,r1_high,'orange',ls='--',lw=1.5,label='68%')
    plt.plot(bin_cents,r1_low,'orange',ls='--',lw=1.5)
    return

lim1 = 68
q1_low = (100-lim1)//2
q1_high = 100-q1_low
lim2 = 95
q2_low = (100-lim2)//2
q2_high = 100-q2_low
