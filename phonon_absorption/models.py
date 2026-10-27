"""Single-mode polaron spectra adapted from author source formulas; energies in eV.

Public amplitude is |psi|^2; the original scripts use wf**2. Channel-dependent
phonon energies describe independent excitons, not one exciton's multimode sum.
"""

from dataclasses import dataclass
import numpy as np
from .thermal import bose_occupation, distribution


@dataclass(frozen=True)
class Channel:
    exciton_energy_eV: float
    coupling: float
    amplitude: float
    mixing_angle_rad: float = 0.0

    def __post_init__(self):
        values = [
            self.exciton_energy_eV,
            self.coupling,
            self.amplitude,
            self.mixing_angle_rad,
        ]
        if not np.isfinite(values).all():
            raise ValueError("channel parameters must be finite")
        if self.coupling < 0 or self.amplitude < 0:
            raise ValueError("coupling magnitude and amplitude must be nonnegative")


@dataclass(frozen=True)
class VibronicChannel:
    channel: Channel
    phonon_energy_eV: float
    linewidth_eV: float

    def __post_init__(self):
        if not isinstance(self.channel, Channel):
            raise ValueError("channel must be a Channel")
        values = [self.phonon_energy_eV, self.linewidth_eV]
        if not np.isfinite(values).all() or min(values) <= 0:
            raise ValueError("phonon energy and linewidth must be finite and positive")


def sideband_weights(huang_rhys, phonon_energy_eV, kBT_eV, cutoff=80):
    """Compatible (m, weight) API; truncated weights are not renormalized."""
    result = distribution(huang_rhys, phonon_energy_eV, kBT_eV, cutoff)
    return result.indices, result.weights


def _model_ii_factors(theta, g, n):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    w0 = 4 * (
        c**4 * s**4 + (c**6 * s**2 + c**2 * s**6) * np.exp(-2 * g * g * (2 * n + 1))
    )
    return w0, 2 * c**4 * s**4, c**8 + s**8


def _weighted_lorentzians(energy, centers, weights, linewidth, max_working_mb=32.0):
    """Block large energy grids; budget estimates kernel temporaries, not total RAM."""
    energy = np.asarray(energy, dtype=float)
    if not np.isfinite(energy).all():
        raise ValueError("energy coordinates must be finite")
    if not np.isfinite(max_working_mb) or max_working_mb <= 0:
        raise ValueError("working-array budget must be positive")
    flat = energy.ravel()
    out = np.empty_like(flat)
    block = max(1, int(max_working_mb * 1024**2 / (8 * len(centers) * 3)))
    for start in range(0, len(flat), block):
        delta = flat[start : start + block, None] - centers
        out[start : start + block] = np.sum(
            linewidth / (delta * delta + linewidth * linewidth) * weights, axis=1
        )
    return out.reshape(energy.shape)


@dataclass(frozen=True)
class PreparedSpectrum:
    """Fixed-parameter thermal weights and centers, reusable across energy grids."""

    kernels: tuple
    diagnostics: tuple

    def evaluate(self, energy_eV, max_working_mb=32.0, components=False):
        values = [
            _weighted_lorentzians(energy_eV, centers, weights, line, max_working_mb)
            for centers, weights, line in self.kernels
        ]
        return np.asarray(values) if components else np.sum(values, axis=0)


def prepare_spectrum(channels, model, kBT_eV, cutoff=80, tail_tolerance=1e-12):
    """Prepare independent kernels. Explicit cutoff=None enables adaptive support."""
    if model not in ("I", "II"):
        raise ValueError("model must be I or II")
    channels = tuple(channels)
    if not channels:
        raise ValueError("provide at least one vibronic channel")
    kernels, diagnostics = [], []
    for item in channels:
        if not isinstance(item, VibronicChannel):
            raise ValueError("expected VibronicChannel objects")
        ch, omega, line = item.channel, item.phonon_energy_eV, item.linewidth_eV
        s = ch.coupling**2 * (1 if model == "I" else 4)
        thermal = distribution(s, omega, kBT_eV, cutoff, tail_tolerance)
        m = thermal.indices
        weights = thermal.weights.copy()
        centers = ch.exciton_energy_eV + m * omega
        bound = ch.amplitude * thermal.omitted_mass / line
        if model == "II":
            w0, w1, w2 = _model_ii_factors(
                ch.mixing_angle_rad, ch.coupling, bose_occupation(omega, kBT_eV)
            )
            weights *= w2 + (-1.0) ** m * w1
            centers = np.r_[centers, ch.exciton_energy_eV]
            weights = np.r_[weights, w0]
            bound *= w1 + w2
        weights *= ch.amplitude
        centers.setflags(write=False)
        weights.setflags(write=False)
        kernels.append((centers, weights, line))
        diagnostics.append(
            {
                **thermal.diagnostics(),
                "exciton_energy_eV": ch.exciton_energy_eV,
                "spectral_absolute_truncation_bound_estimate": float(bound),
            }
        )
    return PreparedSpectrum(tuple(kernels), tuple(diagnostics))


def spectrum_independent(
    energy_eV,
    channels,
    model,
    kBT_eV,
    cutoff=80,
    tail_tolerance=1e-12,
    max_working_mb=32.0,
):
    return prepare_spectrum(channels, model, kBT_eV, cutoff, tail_tolerance).evaluate(
        energy_eV, max_working_mb
    )


def spectrum(
    energy_eV,
    channels,
    model,
    phonon_energy_eV,
    linewidth_eV,
    kBT_eV,
    cutoff=80,
    tail_tolerance=1e-12,
    max_working_mb=32.0,
):
    """Compatible shared-frequency/linewidth API."""
    modes = [VibronicChannel(c, phonon_energy_eV, linewidth_eV) for c in channels]
    return spectrum_independent(
        energy_eV, modes, model, kBT_eV, cutoff, tail_tolerance, max_working_mb
    )


def model_i(
    energy_eV, channel, phonon_energy_eV, linewidth_eV, kBT_eV, cutoff=80, **kwargs
):
    """Intraband-only model: channel.coupling is g_cv."""
    return spectrum(
        energy_eV,
        [channel],
        "I",
        phonon_energy_eV,
        linewidth_eV,
        kBT_eV,
        cutoff,
        **kwargs,
    )


def model_ii(
    energy_eV, channel, phonon_energy_eV, linewidth_eV, kBT_eV, cutoff=80, **kwargs
):
    """Intra/interband model: channel.coupling is sqrt(g_inter**2+g_diag**2)."""
    return spectrum(
        energy_eV,
        [channel],
        "II",
        phonon_energy_eV,
        linewidth_eV,
        kBT_eV,
        cutoff,
        **kwargs,
    )


def model_ii_components(
    energy_eV,
    channel,
    phonon_energy_eV,
    linewidth_eV,
    kBT_eV,
    cutoff=80,
    tail_tolerance=1e-12,
):
    """Original alpha0/alpha1/alpha2 terms; alpha1 is a signed contribution."""
    item = VibronicChannel(channel, phonon_energy_eV, linewidth_eV)
    thermal = distribution(
        4 * channel.coupling**2, phonon_energy_eV, kBT_eV, cutoff, tail_tolerance
    )
    w0, w1, w2 = _model_ii_factors(
        channel.mixing_angle_rad,
        channel.coupling,
        bose_occupation(phonon_energy_eV, kBT_eV),
    )
    centers = channel.exciton_energy_eV + thermal.indices * phonon_energy_eV

    def evaluate(c, w):
        return _weighted_lorentzians(
            energy_eV,
            np.asarray(c),
            np.asarray(w) * channel.amplitude,
            item.linewidth_eV,
        )

    return {
        "alpha0": evaluate([channel.exciton_energy_eV], [w0]),
        "alpha1": evaluate(centers, thermal.weights * (-1.0) ** thermal.indices * w1),
        "alpha2": evaluate(centers, thermal.weights * w2),
    }
