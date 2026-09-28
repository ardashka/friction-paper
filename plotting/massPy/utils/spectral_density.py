import scipy.stats as stats
import numpy as np

def power_spectrum(Fx, Fy=0.):
    # fast fourier transform (complex field)
    psi_fft = np.fft.fft2( Fx + 1J*Fy )
    # unit-wave vector
    k_hat = (     np.fft.fftfreq(psi_fft.shape[0])[:,None] * psi_fft.shape[0]
             + 1J*np.fft.fftfreq(psi_fft.shape[1])[None,:] * psi_fft.shape[1] )
    k_nrm = np.abs(k_hat)
    np.divide(k_hat, k_nrm, out=k_hat, where=k_nrm!=0)
    k_nrm = k_nrm.flatten()
    
    k_max = min(Fx.shape)//2 + 1
    k_bins = np.arange(.5, k_max, 1.)
    k_vals = .5 * (k_bins[1:] + k_bins[:-1])
    
    psi2_fft = (.5 * np.abs(psi_fft)**2).flatten()
    Pk, _, _ = stats.binned_statistic(k_nrm, psi2_fft, statistic="mean", bins=k_bins)
    
    return Pk#, k_vals
