import numpy as np

def decomposition_fft(Fx, Fy):
    """ Hodge-Helmholtz decomposition """
    # fast fourier transform
    psi_fft = np.fft.fft2(Fx + 1J*Fy)
    # unit-wave vector
    k_hat = (     np.fft.fftfreq(psi_fft.shape[0])[:,None] * psi_fft.shape[0]
             + 1J*np.fft.fftfreq(psi_fft.shape[1])[None,:] * psi_fft.shape[1] )
    k_nrm = np.abs(k_hat)
    np.divide(k_hat, k_nrm, out=k_hat, where=k_nrm!=0)
    # longitudinal/irrotational component
    Fl_fft = np.real(k_hat * np.conj(psi_fft)) * k_hat
    # transverse (solenoidal) component
    Ft_fft = psi_fft - Fl_fft
    return Fl_fft, Ft_fft, k_nrm

def decomposition(Fx, Fy):
    """ Hodge-Helmholtz decomposition """
    Fl_fft, Ft_fft = decomposition_fft(Fx, Fy)
    return np.fft.irfft2(Fl_fft), np.fft.irfft2(Ft_fft)
