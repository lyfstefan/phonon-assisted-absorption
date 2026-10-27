"""Extracted candidate independent-frequency formulas; not a confirmed final fit.
No original arrays, plotting or top-level optimization included.
"""
import numpy as np
from scipy.special import iv
from math import sinh

def nB(Omega, T):
    return 1 / (np.exp(Omega / T) - 1)

def tildeG(Omega, T, g):
    return (2 * nB(Omega, T) + 1) * g**2

def LorenzianPeak(gamma, omega):
    return gamma / (omega**2 + gamma**2)

def W0(theta, g, Omega, T):
    term1 = (np.cos(theta / 2)**4 * np.sin(theta / 2)**4)
    term2 = ((np.cos(theta / 2)**6 * np.sin(theta / 2)**2) +
             (np.cos(theta / 2)**2 * np.sin(theta / 2)**6)) * np.exp(-2 * tildeG(Omega, T, g))
    return 4 * (term1 + term2)

def W1(theta):
    return 2 * (np.cos(theta / 2)**4 * np.sin(theta / 2)**4)

def W2(theta):
    return (np.cos(theta / 2)**8 + np.sin(theta / 2)**8)

def Alpha0(omega, theta, g, Omega, E_0, gamma, wf, T):
    return W0(theta, g, Omega, T) * LorenzianPeak(gamma, omega - E_0) * wf**2

def Alpha1(omega, theta, g, Omega, E_0, gamma, wf, T):
    sum_result = 0
    for m in range(-20, 20):
        prefactor = (-1)**m * np.exp(m * Omega / (2 * T))
        lorentz = LorenzianPeak(gamma, omega - m * Omega - E_0)
        bessel = iv(m, 4 * g**2 / sinh(Omega / (2 * T)))
        sum_result += prefactor * lorentz * bessel
    return W1(theta) * np.exp(-4 * tildeG(Omega, T, g)) * sum_result * wf**2

def Alpha2(omega, theta, g, Omega, E_0, gamma, wf, T):
    sum_result = 0
    for m in range(-20, 20):
        prefactor = np.exp(m * Omega / (2 * T))
        lorentz = LorenzianPeak(gamma, omega - m * Omega - E_0)
        bessel = iv(m, 4 * g**2 / sinh(Omega / (2 * T)))
        sum_result += prefactor * lorentz * bessel
    return W2(theta) * np.exp(-4 * tildeG(Omega, T, g)) * sum_result * wf**2

def alpha(omega, theta, g, Omega, E_0, gamma, wf, T):
    return (Alpha0(omega, theta, g, Omega, E_0, gamma, wf, T) +
            Alpha1(omega, theta, g, Omega, E_0, gamma, wf, T) +
            Alpha2(omega, theta, g, Omega, E_0, gamma, wf, T))

def Alpha_fit(omega, Omega1, Omega2, gamma, T, theta1, g1, E_1, wf1, theta2, g2, E_2, wf2):
    return (alpha(omega, theta1, g1, Omega1, E_1, gamma, wf1, T)+
            alpha(omega, theta2, g2, Omega2, E_2, gamma, wf2, T))


def Alpha_fit_gamma12(omega, Omega1, Omega2, gamma1, gamma2, T, theta1, g1, E_1, wf1, theta2, g2, E_2, wf2):
        return (alpha(omega, theta1, g1, Omega1, E_1, gamma1, wf1, T)+
                alpha(omega, theta2, g2, Omega2, E_2, gamma2, wf2, T))
