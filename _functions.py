from _icd import *


def modified_friedmann(time, var: list, m: float, r: float, l: float, k: float,
                       p: float):
    """
    friedmann equations from alternative model
    Parameters
    ----------
    time : float
        time
    var : list
        list of variables [a, adot]
    m : float
        matter density
    r : float
        radiation density
    l : float
        lambda density
    k : float
        constant
    p : float
        power
    Returns
    -------
    list
        list of derivatives [a, adot]
    """

    dvar = np.zeros(2)

    dvar[0] = var[1]           
    g_lcdm = (
        var[1]**2 * var[0] * m**-1 - 
        r * m**-1 * var[0]**-1 -
        l * m**-1 * var[0]**3 - 
        1
    )

    numerator = (
        np.sign(k) * np.abs(k)**p * np.abs(g_lcdm)**(1 - p) - 
        g_lcdm * var[1]**2 - var[1]**4 * var[0] * p * m**-1 - 
        r * m**-1 * var[1]**2 * var[0]**-1 * p + 
        3 * var[1]**2 * var[0]**3 * l * m**-1 * p
    )
    
    denominator = 2 * var[1]**2 * var[0]**2 * p * m**-1 - var[0] * g_lcdm
    dvar[1] = numerator / denominator

    return dvar

###############################################################################

def acceleration(
        a: list, ap: list, m: float, r: float, l: float, k: float, p: float
):
    """
    Acceleration as a function of scale factor and its derivative

    Parameters
    ----------
    a : list
        scale factor
    ap : list
        derivative of scale factor
    m : float
        matter density
    r : float
        radiation density
    l : float
        lambda density
    k : float
        constant
    p : float
        power

    Returns
    -------
    float
        acceleration
    """

    g = ap**2 * a * m**-1 - r * m**-1 * a**-1 - l * m**-1 * a**3 - 1

    numerator = (
        np.sign(k) * np.abs(k)**p * np.sign(g) * np.abs(g)**(1-p) -
        g * ap**2 - 
        ap**4 * a * p * m**-1 -
        r * m**-1 * ap**2 * a**-1 * p +
        3 * ap**2 * a**3 * l * m**-1 * p
    )
    
    denominator = (2 * ap**2 * a**2 * p * m**-1) - (a * g)

    return numerator / denominator

###############################################################################

def dm_z_o4(z: list = z_sn, q: float = -0.55, j: float = 1., s: float = 0.):
    """
    Distance modulus as a function of redshift

    :param z: redshift
    :param q: deceleration parameter
    :param j: jerk
    :param s: snap
    :return: array, distance modulus
    """
    dl =  (
        c * z / h0 * (1 + 0.5 * z * (1 - q) - 
                      (1/6) * z**2 * (1 - q - 3 * q**2 + j) + 
                      (1/24) * z**3 * (2 - 2 * q - 15 * q**2 - 15 * q**3 + 
                                       5 * j + 10*j*q + s))
    )

    return 5 * np.log10(dl) + 25

###############################################################################

def rchi2(obs: list, exp: list = sndat, method: str = 'formula'):
    """
    Reduced chi squared

    Parameters
    ----------
    obs : list
        observed data
    exp : list
        expected data
    method : str
        method to calculate chi squared

    Returns
    -------
    float
        reduced chi squared
    """

    if method == 'formula':
        chi = np.sum((obs - exp)**2)
        answer = chi / len(obs)

    elif method == 'poly':
        p, residuals, _, _, _ = np.polyfit(z_sn, exp, 2, full=True)
        answer = residuals / (len(exp) - 2)
    else:
        raise Exception('Bad "method" input for rchi2 method. \
                        Options are "formula" or "poly"')
    
    return float(answer)

###############################################################################

def read_model_data(fname: str, fdir: str = '../../Data/model_data/'):
    """
    Read model data quickly

    Parameters
    ----------
    fname : str
        name of file
    fdir : str
        directory of file

    Returns
    -------
    pandas dataframe
        data
    """

    file = pd.read_csv(fdir+fname, delim_whitespace=True, comment='#')

    return file

###############################################################################

def specific_function(array: list, number: int):
    """
    Function to do something specific and time consuming

    Parameters
    ----------
    array : list
        list of data
    number : int
        number to choose from array

    Returns
    -------
    float
        number chosen from array
    """

    chi = array['chi']
    beta = array['beta']
    kappa = array['kappa']

    chi_chosen = np.sort(chi)[number]
    beta_chosen = float(np.array(beta)[np.where(chi == chi_chosen)])
    kappa_chosen = float(np.array(kappa)[np.where(chi == chi_chosen)])

    return beta_chosen, kappa_chosen

###############################################################################

def timer(func):
    """
    Decorator to time functions

    Parameters
    ----------
    func : function
        function to time

    Returns
    -------
    function
        wrapper
    """

    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        delta = end - start
        time_elapsed = delta / 60 if delta > 60 else delta
        print('\n\t{} took {} m {} s\n'.format(
            func.__name__, int(time_elapsed) if delta > 60 else 0,
            round(time_elapsed % 1 * 60, 2) if delta > 60 else round(delta, 1)
        ))
        return result
    
    return wrapper

###############################################################################