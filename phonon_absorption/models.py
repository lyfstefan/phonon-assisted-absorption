"""Single-mode finite-temperature polaron spectra. All energies are in eV.

Public amplitudes are |psi|^2, whereas the recovered source uses wf**2.
The thermal Bessel weights are evaluated as a Skellam distribution to avoid
the original product of separately overflowing and underflowing factors.
"""
from dataclasses import dataclass
import numpy as np
from scipy.stats import skellam, poisson


@dataclass(frozen=True)
class Channel:
    exciton_energy_eV: float
    coupling: float
    amplitude: float
    mixing_angle_rad: float = 0.0

    def __post_init__(self):
        if not np.all(np.isfinite([self.exciton_energy_eV,self.coupling,self.amplitude,self.mixing_angle_rad])):
            raise ValueError("channel parameters must be finite")
        if self.coupling < 0 or self.amplitude < 0:
            raise ValueError("coupling magnitude and amplitude must be nonnegative")


def _thermal(phonon_energy_eV,kBT_eV):
    if not np.isfinite(phonon_energy_eV) or phonon_energy_eV <= 0:
        raise ValueError("phonon energy must be finite and positive")
    if not np.isfinite(kBT_eV) or kBT_eV < 0:
        raise ValueError("kBT must be finite and nonnegative (eV, not kelvin)")
    if kBT_eV == 0:
        return 0.0
    ratio=phonon_energy_eV/kBT_eV
    return 0.0 if ratio>700 else 1/np.expm1(ratio)


def sideband_weights(huang_rhys,phonon_energy_eV,kBT_eV,cutoff=80):
    """Return symmetric m=-cutoff...cutoff and UNRENORMALIZED thermal weights.

    P_m = exp[-S(2nB+1)] exp[m Omega/(2 kBT)]
          I_m(2S sqrt(nB(nB+1))).
    This is the difference of independent Poisson variables with means
    S(nB+1) and S*nB. At T=0 it is a one-sided Poisson law. Truncated
    weights are not rescaled: callers must check cutoff convergence.
    """
    if not isinstance(cutoff,int) or cutoff<1:
        raise ValueError("cutoff must be an integer >= 1")
    if not np.isfinite(huang_rhys) or huang_rhys < 0:
        raise ValueError("Huang-Rhys factor must be finite and nonnegative")
    n=_thermal(phonon_energy_eV,kBT_eV)
    m=np.arange(-cutoff,cutoff+1)
    if huang_rhys == 0:
        weights=(m==0).astype(float)
    elif n == 0:
        weights=poisson.pmf(m,huang_rhys)
    else:
        weights=skellam.pmf(m,huang_rhys*(n+1),huang_rhys*n)
    if not np.isfinite(weights).all():
        raise FloatingPointError("thermal weights are nonfinite at these parameters")
    return m,weights


def _lorentz(energy,centers,linewidth):
    if not np.isfinite(linewidth) or linewidth<=0:
        raise ValueError("linewidth (HWHM) must be finite and positive")
    energy=np.asarray(energy,dtype=float)
    if not np.isfinite(energy).all():
        raise ValueError("energy coordinates must be finite")
    return linewidth/((energy[...,None]-np.asarray(centers))**2+linewidth**2)


def model_i(energy_eV,channel,phonon_energy_eV,linewidth_eV,kBT_eV,cutoff=80):
    """Intraband-only model: channel.coupling means g_cv (not g_norm)."""
    m,w=sideband_weights(channel.coupling**2,phonon_energy_eV,kBT_eV,cutoff)
    l=_lorentz(energy_eV,channel.exciton_energy_eV+m*phonon_energy_eV,linewidth_eV)
    return channel.amplitude*np.sum(l*w,axis=-1)


def model_ii(energy_eV,channel,phonon_energy_eV,linewidth_eV,kBT_eV,cutoff=80):
    """Intra/interband model: channel.coupling means sqrt(g_inter^2+g_diag^2)."""
    g,theta=channel.coupling,channel.mixing_angle_rad
    n=_thermal(phonon_energy_eV,kBT_eV)
    c,s=np.cos(theta/2),np.sin(theta/2)
    w0=4*(c**4*s**4+(c**6*s**2+c**2*s**6)*np.exp(-2*g*g*(2*n+1)))
    w1=2*c**4*s**4
    w2=c**8+s**8
    m,thermal=sideband_weights(4*g*g,phonon_energy_eV,kBT_eV,cutoff)
    l=_lorentz(energy_eV,channel.exciton_energy_eV+m*phonon_energy_eV,linewidth_eV)
    direct=_lorentz(energy_eV,[channel.exciton_energy_eV],linewidth_eV)[...,0]
    return channel.amplitude*(w0*direct+np.sum(l*thermal*(w2+(-1.)**m*w1),axis=-1))


def spectrum(energy_eV,channels,model,phonon_energy_eV,linewidth_eV,kBT_eV,cutoff=80):
    """Independent exciton channels sharing phonon energy, linewidth and kBT."""
    if model not in ("I","II"):
        raise ValueError("model must be I or II")
    fn=model_i if model=="I" else model_ii
    return sum((fn(energy_eV,c,phonon_energy_eV,linewidth_eV,kBT_eV,cutoff) for c in channels),
               start=np.zeros_like(np.asarray(energy_eV,dtype=float)))
