"""Stable sideband probabilities with explicit truncation-tail estimates."""

from dataclasses import dataclass
import numpy as np
from scipy.stats import skellam, poisson
from scipy.special import gammaln, hyp0f1, xlogy, ive


def bose_occupation(phonon_energy_eV, kBT_eV):
    if not np.isfinite(phonon_energy_eV) or phonon_energy_eV <= 0:
        raise ValueError("phonon energy must be finite and positive")
    if not np.isfinite(kBT_eV) or kBT_eV < 0:
        raise ValueError("kBT must be finite and nonnegative, in eV")
    if kBT_eV == 0:
        return 0.0
    ratio = phonon_energy_eV / kBT_eV
    return 0.0 if ratio > 700 else float(1 / np.expm1(ratio))


@dataclass(frozen=True)
class Sidebands:
    indices: np.ndarray
    weights: np.ndarray
    omitted_mass: float
    adaptive: bool
    tail_method: str = "scipy_cdf_sf"
    weight_method: str = "scipy_pmf"

    def diagnostics(self):
        return {
            "first_sideband": int(self.indices[0]),
            "last_sideband": int(self.indices[-1]),
            "terms": len(self.indices),
            "retained_weight": float(self.weights.sum()),
            "omitted_thermal_mass_estimate": self.omitted_mass,
            "adaptive": self.adaptive,
            "tail_method": self.tail_method,
            "weight_method": self.weight_method,
        }


def _tails(rv, lo, hi, mu_plus, mu_minus):
    try:
        value = float(rv.cdf(lo - 1) + rv.sf(hi))
        if np.isfinite(value):
            return value, "scipy_cdf_sf"
    except (OverflowError, FloatingPointError, RuntimeError):
        pass
    # D=X-Y, X,Y independent Poisson: D>hi implies X>hi;
    # D<lo implies Y>-lo for lo<=0. This union bound remains
    # conservative when lo>0 (then its left term is simply one).
    bound = float(poisson.sf(-lo, mu_minus) + poisson.sf(hi, mu_plus))
    return min(1.0, bound), "independent_poisson_tail_bound"


def _weights(rv, m, mu_plus, mu_minus):
    try:
        values = rv.pmf(m)
        if np.isfinite(values).all():
            return values, "scipy_pmf"
    except (OverflowError, FloatingPointError, RuntimeError):
        pass
    order = abs(m)
    product = mu_plus * mu_minus
    if product <= 100:
        # Exact positive Poisson-convolution / hypergeometric expression,
        # evaluated in log space; avoids high-order Boost gamma overflow.
        mean = np.where(m >= 0, mu_plus, mu_minus)
        logp = (
            -mu_plus
            - mu_minus
            + xlogy(order, mean)
            - gammaln(order + 1)
            + np.log(hyp0f1(order + 1, product))
        )
        return np.exp(logp), "log_poisson_hypergeometric"
    argument = 2 * np.sqrt(product)
    with np.errstate(divide="ignore"):
        logp = (
            -mu_plus
            - mu_minus
            + argument
            + np.log(ive(order, argument))
            + m * 0.5 * (np.log(mu_plus) - np.log(mu_minus))
        )
    return np.exp(logp), "scaled_bessel_log"


def distribution(
    huang_rhys,
    phonon_energy_eV,
    kBT_eV,
    cutoff=80,
    tail_tolerance=1e-12,
    max_terms=16385,
):
    """Unrenormalized Skellam weights; Poisson at T=0, delta at zero coupling.

    cutoff=None expands an initial mean/variance-based bracket until CDF/SF
    tails meet tolerance. Weights never use a normal approximation. Tail
    estimates use SciPy distribution functions, not certified arithmetic.
    """
    if not np.isfinite(huang_rhys) or huang_rhys < 0:
        raise ValueError("S must be finite and nonnegative")
    if not np.isfinite(tail_tolerance) or not 0 < tail_tolerance < 1:
        raise ValueError("tail tolerance must lie strictly between zero and one")
    if not isinstance(max_terms, int) or max_terms < 3:
        raise ValueError("max_terms must be an integer >= 3")
    if cutoff is not None and (
        not isinstance(cutoff, int) or isinstance(cutoff, bool) or cutoff < 1
    ):
        raise ValueError("cutoff must be None or an integer >= 1")
    n = bose_occupation(phonon_energy_eV, kBT_eV)
    if huang_rhys == 0:
        if cutoff is not None and 2 * cutoff + 1 > max_terms:
            raise ValueError("requested support exceeds max_terms")
        m = np.array([0]) if cutoff is None else np.arange(-cutoff, cutoff + 1)
        return Sidebands(m, (m == 0).astype(float), 0.0, cutoff is None)
    mu_plus, mu_minus = huang_rhys * (n + 1), huang_rhys * n
    rv = poisson(huang_rhys) if n == 0 else skellam(mu_plus, mu_minus)
    if cutoff is None:
        variance = huang_rhys * (2 * n + 1)
        if not np.isfinite(variance):
            raise FloatingPointError("nonfinite sideband variance")
        radius = 8 * np.sqrt(variance + 1)
        for _ in range(20):
            lo = int(np.floor(huang_rhys - radius))
            hi = int(np.ceil(huang_rhys + radius))
            if n == 0:
                lo = max(0, lo)
            if hi - lo + 1 > max_terms:
                raise ValueError("adaptive support exceeds max_terms")
            omitted, tail_method = _tails(rv, lo, hi, mu_plus, mu_minus)
            if np.isfinite(omitted) and omitted <= tail_tolerance:
                break
            radius *= 1.5
        else:
            raise FloatingPointError("sideband tail search failed")
    else:
        lo, hi = -cutoff, cutoff
        if hi - lo + 1 > max_terms:
            raise ValueError("requested support exceeds max_terms")
        omitted, tail_method = _tails(rv, lo, hi, mu_plus, mu_minus)
    m = np.arange(lo, hi + 1)
    weights, weight_method = _weights(rv, m, mu_plus, mu_minus)
    if not np.isfinite(weights).all() or not np.isfinite(omitted):
        raise FloatingPointError("nonfinite thermal distribution")
    return Sidebands(
        m,
        weights,
        float(np.clip(omitted, 0, 1)),
        cutoff is None,
        tail_method,
        weight_method,
    )
