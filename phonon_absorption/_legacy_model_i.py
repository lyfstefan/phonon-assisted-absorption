"""Pure functions extracted verbatim from author research scripts.
Historical reference only; see docs/provenance.md for limitations.
"""
import numpy as np
import scipy.special


def n_B(Omega, T):
    return 1.0 / (np.exp(Omega/T) - 1)

def alpha(omega,g_cv,Omega,E_0,gamma,wf,T):
    term1 = np.exp(-g_cv**2 * (2 * n_B(Omega, T) + 1))
    summation = 0
    for m in range(-50, 50):
        numerator = np.exp(0.5 * m  * Omega/T)  * gamma
        denominator = (omega - m*Omega - E_0)**2 + gamma**2
        factor = numerator/denominator
        bessel_term = scipy.special.iv(m,g_cv**2 / np.sinh(Omega/2/T))
        summation += factor * bessel_term
    return term1 * summation * wf**2

def Alpha_fit(omega, Omega, gamma, T, g1, E_1, wf1, g2, E_2, wf2):
    return (alpha(omega, g1, Omega, E_1, gamma, wf1, T)+
            alpha(omega, g2, Omega, E_2, gamma, wf2, T))
