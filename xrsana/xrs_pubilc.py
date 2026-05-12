import numpy as np
from scipy import constants
import math
from re import findall
import os
import numbers
from functools import lru_cache
from scipy import interpolate
data_installation_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"resources",  'data')

def element(z):
    """
    Convert element symbol to atomic number, or atomic number to element symbol.

    Examples:
        element("Si") -> 14
        element(14)   -> "Si"
    """

    symbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf"
    ]

    symbol_to_z = {symbol: i + 1 for i, symbol in enumerate(symbols)}

    # Optional compatibility with old notation.
    symbol_to_z["Ku"] = 104

    if isinstance(z, str):
        symbol = z.strip()

        if not symbol:
            raise ValueError("Empty element symbol.")

        symbol = symbol[0].upper() + symbol[1:].lower()

        if symbol not in symbol_to_z:
            raise ValueError(f"Unknown element symbol: {z}")

        return symbol_to_z[symbol]

    if isinstance(z, numbers.Integral):
        atomic_number = int(z)

        if atomic_number < 1 or atomic_number > len(symbols):
            raise ValueError(f"Unknown atomic number: {z}")

        return symbols[atomic_number - 1]

    raise TypeError(f"Unsupported type for element(): {type(z)}")



@lru_cache(maxsize=None)
def _load_logtable(logtablefile):
    return np.loadtxt(logtablefile)


def _logpoly(le, coeffs):
    """
    Evaluate exp(a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3).
    coeffs should contain four coefficients.
    """
    return np.exp(
        coeffs[0]
        + le * (
            coeffs[1]
            + le * (
                coeffs[2]
                + le * coeffs[3]
            )
        )
    )


def prho(
    energy,
    Z,
    logtablefile=os.path.join(data_installation_dir, "logtable.dat"),
):
    """
    Calculates photoelectric, elastic, and inelastic mass attenuation
    coefficients of an element.

    Args:
        energy:
            Energy in keV. Scalar or array-like.
        Z:
            Atomic number or element symbol.
        logtablefile:
            Path to logtable.dat.

    Returns:
        murho:
            Mass attenuation coefficients, shape = (len(energy), 3).
            Columns:
                0 = photoelectric
                1 = elastic
                2 = inelastic
            Unit is probably cm^2/g.
        rho:
            Element density, probably g/cm^3.
        m:
            Atomic mass, probably g/mol.
    """

    en = np.atleast_1d(energy).astype(float)

    if np.any(en <= 0):
        raise ValueError("energy must be strictly positive because log(energy) is used.")

    if isinstance(Z, numbers.Integral):
        atomic_number = int(Z)
    else:
        atomic_number = element(Z)

    logtable = _load_logtable(logtablefile)

    matches = np.where(logtable[:, 0] == atomic_number)[0]

    if len(matches) == 0:
        raise ValueError(f"No such element in logtable.dat: Z={Z}")

    ind = matches[0]

    if ind + 5 > logtable.shape[0]:
        raise ValueError(f"Incomplete logtable block for Z={Z}")

    c = np.array(logtable[ind:ind + 5, :], dtype=float)

    le = np.log(en)

    # Photoelectric absorption, piecewise.
    photo = _logpoly(le, c[1:5, 3])

    mask = en <= c[0, 3]
    photo[mask] = _logpoly(le[mask], c[1:5, 2])

    mask = en < c[0, 2]
    photo[mask] = _logpoly(le[mask], c[1:5, 1])

    mask = en < c[0, 1]
    photo[mask] = _logpoly(le[mask], c[1:5, 0])

    # Elastic and inelastic scattering.
    elastic = _logpoly(le, c[1:5, 4])
    inelastic = _logpoly(le, c[1:5, 5])

    mu = np.column_stack((photo, elastic, inelastic))

    atomic_mass = c[0, 4]
    density = c[0, 5]

    if atomic_mass <= 0:
        raise ValueError(f"Invalid atomic mass for Z={Z}: {atomic_mass}")

    murho = mu * 0.602252 / atomic_mass

    return murho, density, atomic_mass

def parseformula(formula):
    """Parses a chemical sum formula.

    Parses the constituing elements and stoichiometries from a given 
    chemical sum formula.

    Args:
      * formula (string): string of a chemical formula (e.g. 'SiO2', 'Ba8Si46', etc.)

    Returns:
      * elements (list): list of strings of constituting elemental symbols.
      * stoichiometries (list): list of according stoichiometries in the same order as 'elements'.
    """
    elements = []
    stoichiometries = []
    splitted = findall(r'([A-Z][a-z]*)(\d*)',formula)
    elements.extend([element[0] for element in splitted])
    stoichiometries.extend([(int(element[1]) if element[1] else 1) for element in splitted])
    return elements,stoichiometries

def mpr(energy, compound):
    """
    Calculates the mass attenuation coefficient of a chemical compound.

    Args:
        energy:
            Energy scale in keV. Can be scalar or array-like.
        compound:
            Chemical formula, for example 'SiO2'.

    Returns:
        murho_total:
            Total mass attenuation coefficient of the compound.
        murho_components:
            Separate contributions:
            column 0 = photoelectric,
            column 1 = elastic,
            column 2 = inelastic.
        element_densities:
            Element densities returned by myprho().
        element_mass_contributions:
            Stoichiometry-weighted atomic masses.
        molar_mass:
            Total formula mass of the compound.
    """

    en = np.atleast_1d(energy).astype(float)

    z, stoich = parseformula(compound)

    murho_components_weighted = np.zeros((len(en), 3), dtype=float)
    element_densities = np.zeros(len(z), dtype=float)
    element_mass_contributions = np.zeros(len(z), dtype=float)

    for i, zi in enumerate(z):
        tmp, rho, atomic_mass = prho(en, zi)

        mass_contribution = atomic_mass * stoich[i]

        element_mass_contributions[i] = mass_contribution
        element_densities[i] = rho

        murho_components_weighted += tmp * mass_contribution

    molar_mass = np.sum(element_mass_contributions)

    murho_components = murho_components_weighted / molar_mass
    murho_total = np.sum(murho_components, axis=1)

    return (
        murho_total,
        murho_components,
        element_densities,
        element_mass_contributions,
        molar_mass,
    )

def e2pz(w1, w2, th, *, eps=1e-14):
    """
    Calculates the momentum scale and the relativistic Compton cross section
    correction according to P. Holm, PRA 37, 3706 (1988).

    Args:
        w1 : float or array-like
            Incident energy in keV.
        w2 : float or array-like
            Scattered energy in keV.
        th : float or array-like
            Scattering angle 2 theta in degrees.
        eps : float
            Numerical tolerance for zero-division checks.

    Returns:
        pz : float or np.ndarray
            Momentum scale in atomic units.
        cf : float or np.ndarray
            Cross section correction factor.
    """

    scalar_input = np.isscalar(w1) and np.isscalar(w2) and np.isscalar(th)

    w1 = np.asarray(w1, dtype=float)
    w2 = np.asarray(w2, dtype=float)
    th = np.asarray(th, dtype=float)

    # Broadcast inputs to a common shape
    try:
        w1, w2, th = np.broadcast_arrays(w1, w2, th)
    except ValueError as exc:
        raise ValueError("w1, w2, and th must be broadcast-compatible.") from exc

    if np.any(w1 <= 0):
        raise ValueError("w1 must be positive, in keV.")

    if np.any(w2 <= 0):
        raise ValueError("w2 must be positive, in keV.")

    # Physical constants
    m = constants.value("electron mass energy equivalent in MeV") * 1e3  # keV
    alp = constants.value("fine-structure constant")
    r0 = constants.value("classical electron radius")  # meter

    th_rad = np.deg2rad(th)
    cos_th = np.cos(th_rad)
    one_minus_cos = 1.0 - cos_th

    if np.any(np.abs(one_minus_cos) < eps):
        raise ValueError("Scattering angle th must not be 0 degrees or too close to 0.")

    # Momentum transfer
    q2 = w1**2 + w2**2 - 2.0 * w1 * w2 * cos_th

    # Avoid tiny negative values caused by floating point roundoff
    q2 = np.where((q2 < 0) & (q2 > -eps), 0.0, q2)

    if np.any(q2 <= eps):
        raise ValueError("Momentum transfer q is zero or too close to zero.")

    q = np.sqrt(q2)

    sqrt_arg = 0.25 + m**2 / (2.0 * w1 * w2 * one_minus_cos)

    if np.any(sqrt_arg < 0):
        raise ValueError("Invalid square-root argument while computing pz.")

    # pz in natural units
    pz = q / 2.0 - (w1 - w2) * np.sqrt(sqrt_arg)

    E = np.sqrt(m**2 + pz**2)

    A = ((w1 - w2) * E - w1 * w2 * one_minus_cos) / q
    D = (w1 - w2 * cos_th) * A / q

    R = w1 * (E - D)
    R2 = R - w1 * w2 * one_minus_cos

    if np.any(np.abs(R) < eps):
        raise ValueError("R is zero or too close to zero; calculation is unstable.")

    if np.any(np.abs(R2) < eps):
        raise ValueError("R2 is zero or too close to zero; calculation is unstable.")

    chi = (
        R / R2
        + R2 / R
        + 2.0 * m**2 * (1.0 / R - 1.0 / R2)
        + m**4 * (1.0 / R - 1.0 / R2) ** 2
    )

    if np.any(np.abs(chi) < eps):
        raise ValueError("chi is zero or too close to zero; calculation is unstable.")

    cf = 2.0 * w1 * q * E / (m**2 * r0**2 * w2 * chi)

    # Unit conversions
    cf = cf * (1.0e-28 * (m * alp))
    pz = pz / (m * alp)

    if scalar_input:
        return float(pz), float(cf)

    return pz, cf










def mpr_compds(energy, formulas, concentrations, E0, rho_formu):
    """
    Calculate linear absorption coefficients for a mixture of compounds.

    Args:
        energy : float or array-like
            Incident/scanned energy scale in keV.
        formulas : list of str
            Chemical formulas, e.g. ["SiO2", "Al2O3"].
        concentrations : list of float
            Relative fractions of each compound.
        E0 : float or array-like
            Outgoing or fixed energy in keV.
        rho_formu : float or list of float
            Density for each compound in g/cm^3.

    Returns:
        mu_tot_in : np.ndarray
            Total linear absorption coefficient at `energy`.
        mu_tot_out : np.ndarray
            Total linear absorption coefficient at `E0`.
    """

    en = np.atleast_1d(np.asarray(energy, dtype=float))
    e0 = np.atleast_1d(np.asarray(E0, dtype=float))

    formulas = list(formulas)
    concentrations = np.atleast_1d(np.asarray(concentrations, dtype=float))
    rho_formu = np.atleast_1d(np.asarray(rho_formu, dtype=float))

    if not (len(formulas) == len(concentrations) == len(rho_formu)):
        raise ValueError(
            "formulas, concentrations, and rho_formu must have the same length."
        )

    if np.any(en <= 0):
        raise ValueError("All incident energies must be positive.")

    if np.any(e0 <= 0):
        raise ValueError("All outgoing energies must be positive.")

    mu_tot_in = np.zeros_like(en, dtype=float)
    mu_tot_out = np.zeros_like(e0, dtype=float)

    for formula, conc, rho in zip(formulas, concentrations, rho_formu):
        mu_tot_in += mpr(en, formula)[0] * conc * rho
        mu_tot_out += mpr(e0, formula)[0] * conc * rho

    # If outgoing energy is scalar but incident energy is an array,
    # broadcast mu_tot_out to the same shape as mu_tot_in.
    if mu_tot_out.size == 1 and mu_tot_in.size > 1:
        mu_tot_out = np.full_like(mu_tot_in, mu_tot_out.item())

    return mu_tot_in, mu_tot_out



def abscorr2(mu1, mu2, alpha, beta, samthick):
    """
    Calculates absorption correction for given mu1 and mu2.
    Multiply the measured spectrum with this correction factor.
    """

    mu1 = np.asarray(mu1, dtype=float)
    mu2 = np.asarray(mu2, dtype=float)

    cosa = math.cos(math.radians(alpha))
    cosb = math.cos(math.radians(beta))

    if beta >= 0:
        # Reflection geometry
        ac = (
            cosa * (mu1 / cosa + mu2 / cosb)
            / (
                1.0
                - np.exp(-mu1 * samthick / cosa - mu2 * samthick / cosb)
            )
        )

    else:
        # Transmission geometry
        x = mu1 / cosa
        y = mu2 / cosb

        general = (
            -cosa * (x - y)
            / (
                np.exp(-x * samthick)
                - np.exp(-y * samthick)
            )
        )

        limit = cosa / (
            samthick * np.exp(-mu1 * samthick / cosa)
        )

        close = np.isclose(x, y)

        ac = np.where(close, limit, general)

    return ac
def spline2(x,y,x2):
    """
    Extrapolates the smaller and larger valuea as a constant
    """
    xmin = np.min(x)
    xmax = np.max(x)
    imin = x == xmin
    imax = x == xmax
    f  = interpolate.interp1d(x,y, bounds_error=False, fill_value=0.0)
    y2 = f(x2)
    i     = np.where(x2<xmin)
    y2[i] = y[imin]
    i     = np.where(x2>xmax)
    y2[i] = y[imax]
    return y2

def pz2e1(w2,pz,th):
    """Calculates the incident energy for a specific scattered photon and momentum value.

    Returns the incident energy for a given photon energy and scattering angle.
    This function is translated from Keijo Hamalainen's Matlab implementation (KH 29.05.96).

    Args:
      * w2 (float): scattered photon energy in [keV]
      * pz (np.array): pz scale in [a.u.]
      * th (float): scattering angle two theta in [deg]

    Returns:
      * w1 (np.array): incident energy in [keV]
    """
    pz  = np.array(pz)
    w   = np.array(np.arange(np.array(w2)/4.0,4.0*np.array(w2),np.array(w2)/5000.0))
    p   = e2pz(w,w2,th)[0]
    if ( p[1]-p[0] <0) :
        tck = interpolate.UnivariateSpline(p[::-1],w[::-1])
    else:
        tck = interpolate.UnivariateSpline(p,w)
    w1  = tck(pz)
    return w1

def gauss(x,x0,fwhm):
    # area-normalized gaussian
    sigma = fwhm/(2*np.sqrt(2*np.log(2)));
    y = np.exp(-(x-x0)**2/2/sigma**2)/sigma/np.sqrt(2*np.pi)
    return y

def convg(x,y,fwhm):
    """
    Convolution with Gaussian
    x  = x-vector
    y  = y-vector
    fwhm = fulll width at half maximum of the gaussian with which y is convoluted
    """
    dx = np.min(np.absolute(np.diff(x)))
    x2 = np.arange(np.min(x)-1.5*fwhm, np.max(x)+1.5*fwhm, dx)
    xg = np.arange(-np.floor(2.0*fwhm/dx)*dx, np.floor(2.0*fwhm/dx)*dx, dx)
    yg = gauss(xg,0,fwhm)
    yg = yg/np.sum(yg)
    y2 = spline2(x,y,x2)
    c  = np.convolve(y2,yg, mode='full')
    n  = int( np.floor(np.max(np.shape(xg))/2))
    c  = c[n:len(c)-n+1] # not sure about the +- 1 here
    f  = interpolate.interp1d(x2,c)
    return f(x)