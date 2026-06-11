"""
Helper functions to compute the per-binary formation efficiency of double
compact objects (DCOs) in *legacy* COMPAS output files (e.g. the models from
Broekgaarden et al. 2022, https://arxiv.org/abs/2103.02608).

The "formation efficiency" of a binary is its sampling weight divided by the
total stellar mass that the COMPAS simulation represents at that binary's
metallicity:

    formation_efficiency = weight / total_mass_evolved_at_that_metallicity   [Msun^-1]

These functions read directly from the COMPAS HDF5 output with h5py. Unlike
``legacy_compas_scripts/ClassCOMPAS.py``, they do not load the full DCO table
into a ``COMPASData`` object (which keeps dozens of unused 17M+ element
arrays in memory) -- only the handful of columns needed for the formation
efficiency calculation are read.
"""

from __future__ import annotations

from typing import Optional, Sequence

import h5py as h5
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Reference metallicity grid used in Broekgaarden+22 (53 points).
#
# This is provided for convenience (e.g. axis limits in plots). The actual
# metallicity grid of a given simulation is always derived directly from its
# HDF5 file by `get_total_mass_evolved_per_metallicity` /
# `calculate_formation_efficiency_per_binary` below, so you do not need to
# pass this in for the calculation itself.
# ---------------------------------------------------------------------------
METALLICITY_GRID_BROEKGAARDEN22 = np.array([
    0.0001, 0.00011, 0.00012, 0.00014, 0.00016, 0.00017,
    0.00019, 0.00022, 0.00024, 0.00027, 0.0003, 0.00034,
    0.00037, 0.00042, 0.00047, 0.00052, 0.00058, 0.00065,
    0.00073, 0.00081, 0.0009, 0.00101, 0.00113, 0.00126,
    0.0014, 0.00157, 0.00175, 0.00195, 0.00218, 0.00243,
    0.00272, 0.00303, 0.00339, 0.00378, 0.00422, 0.00471,
    0.00526, 0.00587, 0.00655, 0.00732, 0.00817, 0.00912,
    0.01018, 0.01137, 0.01269, 0.01416, 0.01581, 0.01765,
    0.01971, 0.022, 0.0244, 0.02705, 0.03,
])


# ---------------------------------------------------------------------------
# The 20 binary-population-synthesis (BPS) variations from Broekgaarden+22
# (https://arxiv.org/abs/2103.02608), labelled A-T. Each entry gives the
# sub-directory (relative to a top-level data directory containing one
# COMPASOutput.h5 per variation), whether the "optimistic" common-envelope
# assumption should be used (see `get_dco_mask`), and a short label for
# plotting/legends.
# ---------------------------------------------------------------------------
BROEKGAARDEN22_MODELS = [
    {'name': 'A', 'directory': 'fiducial', 'optimistic': False, 'label': 'fiducial'},
    {'name': 'B', 'directory': 'massTransferEfficiencyFixed_0_25', 'optimistic': False, 'label': r'$\beta=0.25$'},
    {'name': 'C', 'directory': 'massTransferEfficiencyFixed_0_5', 'optimistic': False, 'label': r'$\beta=0.5$'},
    {'name': 'D', 'directory': 'massTransferEfficiencyFixed_0_75', 'optimistic': False, 'label': r'$\beta=0.75$'},
    {'name': 'E', 'directory': 'unstableCaseBB', 'optimistic': False, 'label': 'unstable case BB'},
    {'name': 'F', 'directory': 'unstableCaseBB', 'optimistic': True, 'label': 'unstable case BB + optimistic CE'},
    {'name': 'G', 'directory': 'alpha0_1', 'optimistic': False, 'label': r'$\alpha_{\rm CE}=0.1$'},
    {'name': 'H', 'directory': 'alpha0_5', 'optimistic': False, 'label': r'$\alpha_{\rm CE}=0.5$'},
    {'name': 'I', 'directory': 'alpha2_0', 'optimistic': False, 'label': r'$\alpha_{\rm CE}=2$'},
    {'name': 'J', 'directory': 'alpha10', 'optimistic': False, 'label': r'$\alpha_{\rm CE}=10$'},
    {'name': 'K', 'directory': 'fiducial', 'optimistic': True, 'label': 'optimistic CE'},
    {'name': 'L', 'directory': 'rapid', 'optimistic': False, 'label': 'rapid SN'},
    {'name': 'M', 'directory': 'maxNSmass2_0', 'optimistic': False, 'label': r'max $m_{\rm NS}=2.0\,M_\odot$'},
    {'name': 'N', 'directory': 'maxNSmass3_0', 'optimistic': False, 'label': r'max $m_{\rm NS}=3.0\,M_\odot$'},
    {'name': 'O', 'directory': 'noPISN', 'optimistic': False, 'label': 'no PISN'},
    {'name': 'P', 'directory': 'ccSNkick_100km_s', 'optimistic': False, 'label': r'$\sigma_{\rm rms}^{\rm 1D}=100\,{\rm km/s}$'},
    {'name': 'Q', 'directory': 'ccSNkick_30km_s', 'optimistic': False, 'label': r'$\sigma_{\rm rms}^{\rm 1D}=30\,{\rm km/s}$'},
    {'name': 'R', 'directory': 'noBHkick', 'optimistic': False, 'label': r'$v_{\rm k,BH}=0\,{\rm km/s}$'},
    {'name': 'S', 'directory': 'wolf_rayet_multiplier_0_1', 'optimistic': False, 'label': r'$f_{\rm WR}=0.1$'},
    {'name': 'T', 'directory': 'wolf_rayet_multiplier_5', 'optimistic': False, 'label': r'$f_{\rm WR}=5$'},
]


# ---------------------------------------------------------------------------
# Kroupa IMF sampling (three-part broken power law), used to work out which
# fraction of all star-forming mass is represented by the binaries that
# COMPAS actually simulates (i.e. those with primary mass in
# [m1_min, m1_max]). This reproduces the Monte-Carlo approach used in
# legacy_compas_scripts/totalMassEvolvedPerZ.py.
# ---------------------------------------------------------------------------

def _broken_power_law_cdf(m: np.ndarray, m_breaks: Sequence[float], slopes: Sequence[float],
                           C1: float = 1.0) -> np.ndarray:
    """CDF of a three-part broken power law p(m) ~ m^slope, defined on [m_breaks[0], m_breaks[-1]]."""
    m1, m2, m3, m4 = m_breaks
    a1, a2, a3 = slopes
    C2 = C1 * m2 ** (a1 - a2)
    C3 = C2 * m3 ** (a2 - a3)

    def integral(C, a, lo, hi):
        return (C / (a + 1)) * (hi ** (a + 1) - lo ** (a + 1))

    N1 = integral(C1, a1, m1, m2)
    N2 = integral(C2, a2, m2, m3)
    N3 = integral(C3, a3, m3, m4)
    norm = N1 + N2 + N3

    cdf = np.zeros_like(m)
    mask1 = (m >= m1) & (m < m2)
    cdf[mask1] = integral(C1, a1, m1, m[mask1]) / norm
    mask2 = (m >= m2) & (m < m3)
    cdf[mask2] = (N1 + integral(C2, a2, m2, m[mask2])) / norm
    mask3 = (m >= m3) & (m <= m4)
    cdf[mask3] = (N1 + N2 + integral(C3, a3, m3, m[mask3])) / norm
    return cdf


def _invert_broken_power_law_cdf(cdf: np.ndarray, m_breaks: Sequence[float], slopes: Sequence[float],
                                  C1: float = 1.0) -> np.ndarray:
    """Inverse CDF of `_broken_power_law_cdf`, used to draw masses from the IMF."""
    m1, m2, m3, m4 = m_breaks
    a1, a2, a3 = slopes
    C2 = C1 * m2 ** (a1 - a2)
    C3 = C2 * m3 ** (a2 - a3)

    def integral(C, a, lo, hi):
        return (C / (a + 1)) * (hi ** (a + 1) - lo ** (a + 1))

    def invert(y, C, a, lo):
        return ((y + (C / (a + 1)) * lo ** (a + 1)) / (C / (a + 1))) ** (1.0 / (a + 1))

    N1 = integral(C1, a1, m1, m2)
    N2 = integral(C2, a2, m2, m3)
    N3 = integral(C3, a3, m3, m4)
    norm = N1 + N2 + N3

    cdf_at_m2 = N1 / norm
    cdf_at_m3 = (N1 + N2) / norm

    m = np.zeros_like(cdf)
    mask1 = cdf < cdf_at_m2
    m[mask1] = invert(cdf[mask1] * norm, C1, a1, m1)

    mask2 = (cdf_at_m2 <= cdf) & (cdf < cdf_at_m3)
    m[mask2] = invert(cdf[mask2] * norm - N1, C2, a2, m2)

    mask3 = cdf_at_m3 <= cdf
    m[mask3] = invert(cdf[mask3] * norm - (N1 + N2), C3, a3, m3)

    return m


def _sample_kroupa_masses(n_samples: int, binary_fraction: float = 1.0,
                           m_breaks: Sequence[float] = (0.01, 0.08, 0.5, 200.0),
                           slopes: Sequence[float] = (-0.3, -1.3, -2.3),
                           C1: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Draw primary/secondary masses (M1, M2) for a population of single stars + binaries.

    Primary masses are drawn from a Kroupa-like three-part broken power law
    IMF over [m_breaks[0], m_breaks[-1]]. A fraction `binary_fraction` of
    stars are binaries; for those, the secondary mass is M2 = q * M1 with q
    drawn uniformly in [0, 1]. Non-binaries get M2 = 0.
    """
    cdf_min = _broken_power_law_cdf(np.array([m_breaks[0]]), m_breaks, slopes, C1)[0]
    cdf_max = _broken_power_law_cdf(np.array([m_breaks[-1]]), m_breaks, slopes, C1)[0]

    draw_cdf_m1 = np.random.uniform(cdf_min, cdf_max, n_samples)
    draw_is_binary = np.random.uniform(0, 1, n_samples)
    draw_q = np.random.uniform(0, 1, n_samples)

    m1 = _invert_broken_power_law_cdf(draw_cdf_m1, m_breaks, slopes, C1)
    m2 = np.zeros(n_samples)
    is_binary = draw_is_binary < binary_fraction
    m2[is_binary] = draw_q[is_binary] * m1[is_binary]

    return m1, m2


# ---------------------------------------------------------------------------
# Total stellar mass evolved per metallicity
# ---------------------------------------------------------------------------

def get_total_mass_evolved_per_metallicity(
        path: str,
        m1_min: float,
        m1_max: float,
        binary_fraction: float = 1.0,
        m_breaks: Sequence[float] = (0.01, 0.08, 0.5, 200.0),
        slopes: Sequence[float] = (-0.3, -1.3, -2.3),
        n_imf_samples: int = 2_000_000,
        random_seed: Optional[int] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Total stellar mass formed (in Msun) at each metallicity of a COMPAS run.

    COMPAS only simulates binaries with a primary mass in
    [`m1_min`, `m1_max`]. To convert the (weighted) number of binaries
    simulated at a given metallicity into a total star-forming mass, we
    Monte-Carlo sample a full stellar population from the IMF and work out
    what fraction of the total mass formed falls in the simulated range.

    Parameters
    ----------
    path : str
        Path to a COMPAS ``COMPASOutput.h5`` file.
    m1_min, m1_max : float
        Primary-mass sampling bounds used by COMPAS (`pythonSubmit`
        ``--initial-mass-min`` / ``--initial-mass-max``).
    binary_fraction : float
        Fraction of stars that are in binaries (the rest are single stars
        with M2 = 0 and do not count towards the simulated mass).
    m_breaks, slopes : sequence of float
        Breakpoints and power-law slopes of the assumed IMF.
    n_imf_samples : int
        Number of Monte-Carlo IMF samples used to estimate the mass
        fraction. Larger values reduce the (small) Monte-Carlo noise.
    random_seed : int, optional
        Seed for the IMF sampling, for reproducibility.

    Returns
    -------
    metallicity_grid : np.ndarray
        Sorted, unique metallicities (Z) simulated by COMPAS.
    total_mass_evolved_per_Z : np.ndarray
        Total stellar mass formed (Msun) at each metallicity in
        `metallicity_grid`.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    with h5.File(path, 'r') as f:
        systems = f['systems']
        metallicities = systems['Metallicity1'][...].squeeze()
        weights = systems['weight'][...].squeeze()

        metallicity_grid, n_binaries_per_Z = np.unique(metallicities, return_counts=True)

        # A "weighted" simulation draws importance-sampled binaries (weight != 1);
        # an "unweighted" simulation draws every binary from the IMF directly.
        is_weighted = not np.all(weights[metallicities == metallicity_grid[0]] == 1)

        if is_weighted:
            mass_formed_per_Z = None  # filled in below using the IMF sample
        else:
            mass1 = systems['mass1'][...].squeeze()
            mass2 = systems['mass2'][...].squeeze()
            _, inverse = np.unique(metallicities, return_inverse=True)
            mass_formed_per_Z = (np.bincount(inverse, weights=mass1)
                                  + np.bincount(inverse, weights=mass2))

    # Monte-Carlo sample a full stellar population from the IMF, and work out
    # what fraction of its mass lies in the [m1_min, m1_max] range that COMPAS
    # actually simulates.
    m1_sample, m2_sample = _sample_kroupa_masses(
        n_imf_samples, binary_fraction=binary_fraction, m_breaks=m_breaks, slopes=slopes)

    total_mass_all = m1_sample.sum() + m2_sample.sum()

    in_compas_range = (m1_sample >= m1_min) & (m1_sample <= m1_max) & (m2_sample != 0)
    total_mass_in_compas_range = m1_sample[in_compas_range].sum() + m2_sample[in_compas_range].sum()

    compas_mass_fraction = total_mass_in_compas_range / total_mass_all

    if is_weighted:
        mean_mass_per_compas_binary = total_mass_in_compas_range / np.sum(in_compas_range)
        mass_formed_per_Z = mean_mass_per_compas_binary * n_binaries_per_Z

    total_mass_evolved_per_Z = mass_formed_per_Z / compas_mass_fraction

    return metallicity_grid, total_mass_evolved_per_Z


# ---------------------------------------------------------------------------
# DCO selection
# ---------------------------------------------------------------------------

VALID_DCO_TYPES = ('BBH', 'BHBH', 'BNS', 'NSNS', 'BHNS')


def get_dco_mask(
        h5file: h5.File,
        dco_type: str = 'BBH',
        within_hubble_time: bool = True,
        optimistic: bool = False,
        no_rlof_after_cee: bool = True,
) -> np.ndarray:
    """Boolean mask (length = full DCO table) selecting the requested DCOs.

    Parameters
    ----------
    h5file : h5py.File
        Open COMPAS HDF5 file.
    dco_type : {'BBH', 'BHBH', 'BNS', 'NSNS', 'BHNS'}
        Type of double compact object to select.
    within_hubble_time : bool
        If True, only keep DCOs that merge within a Hubble time.
    optimistic : bool
        If False (default), exclude binaries that survived a common-envelope
        phase only under the "optimistic" assumption (i.e. require the
        standard, pessimistic CE prescription used in Broekgaarden+22). If
        True, all binaries are kept regardless of their CE outcome.
    no_rlof_after_cee : bool
        If True (default), exclude binaries that immediately start RLOF right
        after a common envelope (these are not expected to survive).

    Returns
    -------
    mask : np.ndarray of bool
        Mask into the ``doubleCompactObjects`` table.
    """
    dco_type = dco_type.upper()
    fdco = h5file['doubleCompactObjects']
    type1 = fdco['stellarType1'][...].squeeze()
    type2 = fdco['stellarType2'][...].squeeze()

    if dco_type in ('BBH', 'BHBH'):
        mask = (type1 == 14) & (type2 == 14)
    elif dco_type in ('BNS', 'NSNS'):
        mask = (type1 == 13) & (type2 == 13)
    elif dco_type == 'BHNS':
        mask = ((type1 == 14) & (type2 == 13)) | ((type1 == 13) & (type2 == 14))
    else:
        raise ValueError(f"dco_type={dco_type!r} not recognised; choose from {VALID_DCO_TYPES}")

    if within_hubble_time:
        mask &= fdco['mergesInHubbleTimeFlag'][...].squeeze()

    if not optimistic:
        mask &= ~fdco['optimisticCEFlag'][...].squeeze()

    if no_rlof_after_cee:
        mask &= ~fdco['RLOFSecondaryAfterCEE'][...].squeeze()

    return mask


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def calculate_formation_efficiency_per_binary(
        path: str,
        dco_type: str = 'BBH',
        m1_min: float = 5.0,
        m1_max: float = 150.0,
        binary_fraction: float = 1.0,
        within_hubble_time: bool = True,
        optimistic: bool = False,
        no_rlof_after_cee: bool = True,
        n_imf_samples: int = 2_000_000,
        random_seed: Optional[int] = None,
        extra_dco_columns: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Formation efficiency of every individual DCO in a legacy COMPAS file.

    The formation efficiency of a binary is its sampling weight divided by
    the total stellar mass formed at its metallicity, i.e. the number of
    such binaries formed per Msun of star formation.

    Parameters
    ----------
    path : str
        Path to a COMPAS ``COMPASOutput.h5`` file.
    dco_type : {'BBH', 'BHBH', 'BNS', 'NSNS', 'BHNS'}
        Type of double compact object to select.
    m1_min, m1_max : float
        Primary-mass sampling bounds used by COMPAS (see
        `get_total_mass_evolved_per_metallicity`).
    binary_fraction : float
        Fraction of stars assumed to be in binaries.
    within_hubble_time, optimistic, no_rlof_after_cee : bool
        DCO selection options, see `get_dco_mask`.
    n_imf_samples : int
        Number of Monte-Carlo IMF samples used to normalise the weights
        (see `get_total_mass_evolved_per_metallicity`).
    random_seed : int, optional
        Seed for the IMF sampling, for reproducibility.
    extra_dco_columns : sequence of str, optional
        Extra column names to read from the ``doubleCompactObjects`` table
        and include in the returned DataFrame (e.g. ``['M1', 'M2', 'tform',
        'tc']``), useful for combining this output with other (e.g.
        cosmic-integration) pipelines.

    Returns
    -------
    binaries : pandas.DataFrame
        One row per selected DCO, with columns ``seed``, ``metallicity``,
        ``weight``, ``formation_efficiency`` (in Msun^-1), plus any columns
        requested via `extra_dco_columns`.
    metallicity_grid : np.ndarray
        Sorted, unique metallicities simulated by COMPAS.
    total_mass_evolved_per_Z : np.ndarray
        Total stellar mass formed (Msun) at each metallicity in
        `metallicity_grid`. Useful e.g. for plotting the underlying star
        formation per metallicity bin.
    """
    if dco_type.upper() not in VALID_DCO_TYPES:
        raise ValueError(f"dco_type={dco_type!r} not recognised; choose from {VALID_DCO_TYPES}")

    metallicity_grid, total_mass_evolved_per_Z = get_total_mass_evolved_per_metallicity(
        path, m1_min=m1_min, m1_max=m1_max, binary_fraction=binary_fraction,
        n_imf_samples=n_imf_samples, random_seed=random_seed)

    with h5.File(path, 'r') as f:
        mask = get_dco_mask(f, dco_type=dco_type, within_hubble_time=within_hubble_time,
                             optimistic=optimistic, no_rlof_after_cee=no_rlof_after_cee)

        fdco = f['doubleCompactObjects']
        data = {
            'seed': fdco['seed'][...].squeeze()[mask],
            'metallicity': fdco['Metallicity1'][...].squeeze()[mask],
            'weight': fdco['weight'][...].squeeze()[mask],
        }
        for column in extra_dco_columns or []:
            data[column] = fdco[column][...].squeeze()[mask]

    binaries = pd.DataFrame(data)

    z_index = np.searchsorted(metallicity_grid, binaries['metallicity'].to_numpy())
    if not np.array_equal(metallicity_grid[z_index], binaries['metallicity'].to_numpy()):
        raise ValueError(
            "Found DCOs with a metallicity that is not in this simulation's metallicity grid.")

    binaries['formation_efficiency'] = binaries['weight'] / total_mass_evolved_per_Z[z_index]

    return binaries, metallicity_grid, total_mass_evolved_per_Z


# ---------------------------------------------------------------------------
# Formation efficiency vs. metallicity, for one or more DCO types at once
# ---------------------------------------------------------------------------

def calculate_formation_efficiency_vs_metallicity(
        path: str,
        dco_types: Sequence[str] = ('BBH', 'BHNS', 'BNS'),
        m1_min: float = 5.0,
        m1_max: float = 150.0,
        binary_fraction: float = 1.0,
        within_hubble_time: bool = True,
        optimistic: bool = False,
        no_rlof_after_cee: bool = True,
        n_imf_samples: int = 2_000_000,
        random_seed: Optional[int] = None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Total formation efficiency eta(Z) for one or more DCO types.

    This sums the per-binary formation efficiencies in each metallicity bin,
    i.e. ``eta(Z) = sum(weight_i) / total_mass_evolved_at_Z`` for binaries i
    of the requested type at metallicity Z. It is a cheaper alternative to
    calling `calculate_formation_efficiency_per_binary` separately for each
    DCO type, since the (slow) IMF-normalisation and the DCO table are each
    only read once.

    Parameters
    ----------
    path : str
        Path to a COMPAS ``COMPASOutput.h5`` file.
    dco_types : sequence of {'BBH', 'BHBH', 'BNS', 'NSNS', 'BHNS'}
        DCO types to compute eta(Z) for.
    m1_min, m1_max, binary_fraction, n_imf_samples, random_seed :
        See `get_total_mass_evolved_per_metallicity`.
    within_hubble_time, optimistic, no_rlof_after_cee : bool
        DCO selection options, see `get_dco_mask`.

    Returns
    -------
    metallicity_grid : np.ndarray
        Sorted, unique metallicities simulated by COMPAS.
    formation_efficiency_per_Z : dict[str, np.ndarray]
        Maps each requested DCO type to its eta(Z) array (Msun^-1), aligned
        with `metallicity_grid`.
    """
    for dco_type in dco_types:
        if dco_type.upper() not in VALID_DCO_TYPES:
            raise ValueError(f"dco_type={dco_type!r} not recognised; choose from {VALID_DCO_TYPES}")

    metallicity_grid, total_mass_evolved_per_Z = get_total_mass_evolved_per_metallicity(
        path, m1_min=m1_min, m1_max=m1_max, binary_fraction=binary_fraction,
        n_imf_samples=n_imf_samples, random_seed=random_seed)

    formation_efficiency_per_Z = {}
    with h5.File(path, 'r') as f:
        fdco = f['doubleCompactObjects']
        metallicity = fdco['Metallicity1'][...].squeeze()
        weight = fdco['weight'][...].squeeze()

        for dco_type in dco_types:
            mask = get_dco_mask(f, dco_type=dco_type, within_hubble_time=within_hubble_time,
                                 optimistic=optimistic, no_rlof_after_cee=no_rlof_after_cee)

            z_index = np.searchsorted(metallicity_grid, metallicity[mask])
            weight_per_Z = np.bincount(z_index, weights=weight[mask], minlength=len(metallicity_grid))
            formation_efficiency_per_Z[dco_type] = weight_per_Z / total_mass_evolved_per_Z

    return metallicity_grid, formation_efficiency_per_Z
