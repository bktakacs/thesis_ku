from _icd import *


def modified_friedmann(time, var, m, r, l, k, p):
    """
    friedmann equations from alternative model
    :param time: time
    :param var: array of variables a, adot
    :param m: matter density
    :param r: radiation density
    :param l: lambda density
    :param k: constant
    :param p: power
    :return: array of derivatives [a, adot]
    """
    dvar = np.zeros(2)

    dvar[0] = var[1]           
    g_lcdm = (var[1]**2 * var[0] * m**-1 -
              r * m**-1 * var[0]**-1 -
              l * m**-1 * var[0]**3 -
              1)

    numerator = (np.sign(k) * np.abs(k)**p * np.abs(g_lcdm)**(1 - p) - 
                 g_lcdm * var[1]**2 - var[1]**4 * var[0] * p * m**-1 - 
                 r * m**-1 * var[1]**2 * var[0]**-1 * p + 
                 3 * var[1]**2 * var[0]**3 * l * m**-1 * p)

    # numerator = 10
    
    denominator = 2 * var[1]**2 * var[0]**2 * p * m**-1 - var[0] * g_lcdm
    dvar[1] = numerator / denominator

    return dvar

###############################################################################

def acceleration(a, ap, m, r, l, k, p):
    g = ap**2 * a * m**-1 - r * m**-1 * a**-1 - l * m**-1 * a**3 - 1

    numerator = (np.sign(k) * np.abs(k)**p * np.sign(g) * np.abs(g)**(1-p) -
                 g * ap**2 -
                 ap**4 * a * p * m**-1 -
                 r * m**-1 * ap**2 * a**-1 * p +
                 3 * ap**2 * a**3 * l * m**-1 * p)
    
    denominator = (2 * ap**2 * a**2 * p * m**-1) - (a * g)

    return numerator / denominator

###############################################################################

def dm_z_o4(z=z_sn, q=-0.55, j=1, s=0):
    """
    Distance modulus as a function of redshift

    :param z: array, redshift
    :param q: float, deceleration parameter
    :param j: float, jerk
    :param s: float, snap
    :return: array, distance modulus
    """
    dl =  (c * z / h0 * (1 + 0.5 * z * (1 - q) - 
                         (1/6) * z**2 * (1 - q - 3 * q**2 + j) + 
                         (1/24) * z**3 * (2 - 2 * q - 15 * q**2 - 15 * q**3 + 
                                          5 * j + 10*j*q + s)))

    return 5 * np.log10(dl) + 25

###############################################################################

def rchi2(obs, exp=sndat, method='formula'):
    """
    Reduced chi squared
    :param obs: array, observed data
    :param exp: array, expected data
    :param method: string, method to calculate chi squared
    :return: float, reduced chi squared
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

def read_model_data(fname, fdir='../../Data/model_data/'):
    """
    Read model data quickly

    :param fname: string, name of file
    :param fdir: string, directory of file
    :return: pandas dataframe, data
    """

    file = pd.read_csv(fdir+fname, delim_whitespace=True, comment='#')

    return file

###############################################################################

def specific_function(array, number):
    """
    Function to do something specific and time consuming

    :param array: array, array of data
    :param number: int, which sorted chi^2 value to choose from which beta and
    kappa are returned
    :return: float, beta, kappa
    """

    chi = array['chi']
    beta = array['beta']
    kappa = array['kappa']

    chi_chosen = np.sort(chi)[number]
    beta_chosen = float(np.array(beta)[np.where(chi == chi_chosen)])
    kappa_chosen = float(np.array(kappa)[np.where(chi == chi_chosen)])

    return beta_chosen, kappa_chosen