from _icd import *


def modified_friedmann(time, var, m, r, l, k, p):
    dvar = np.zeros(2)

    dvar[0] = var[1]
    #                   m                               r                       l               
    g_lcdm  = (var[1]**2 * var[0] * m**-1) - (r * m**-1 * var[0]**-1) - (l * m**-1 * var[0]**3) - 1

    numerator = np.sign(k) * np.abs(k)**p * np.abs(g_lcdm)**(1 - p) - g_lcdm * var[1]**2 - var[1]**4 * var[0] * p * m**-1 - r * m**-1 * var[1]**2 * var[0]**-1 * p + 3 * var[1]**2 * var[0]**3 * l * m**-1 * p
    denominator = 2 * var[1]**2 * var[0]**2 * p * m**-1 - var[0] * g_lcdm
    dvar[1] = numerator / denominator

    return dvar

##############################################################################

def acceleration(a, ap, m, r, l, k, p):
    g = ap**2 * a * m**-1 - r * m**-1 * a**-1 - l * m**-1 * a**3 - 1

    numerator = (np.sign(k) * np.abs(k)**p * np.sign(g) * np.abs(g)**(1-p)) - (g * ap**2) - (ap**4 * a * p * m**-1) - (r * m**-1 * ap**2 * a**-1 * p) + (3 * ap**2 * a**3 * l * m**-1 * p)
    denominator = (2 * ap**2 * a**2 * p * m**-1) - (a * g)

    return numerator / denominator

##############################################################################

def inverse(x, y):
    return y**(-0.5)

##############################################################################

def adot_inverse(a, ap):
    return ap**-2

##############################################################################

def integrand(z, m, r, l, f, k):
    return (m * (1 + z)**3 * (1 + k * f) + r * (1 + z)**4 + l)**(-0.5)

##############################################################################

def dm_z_o4(z=z_sn, q=-0.55, j=1, s=0):
    dl =  c * z / h0 * (1 + 0.5 * z * (1 - q) - (1/6) * z**2 * (1 - q - 3 * q**2 + j) + (1/24) * z**3 * (2 - 2 * q - 15 * q**2 - 15 * q**3 + 5 * j + 10*j*q + s)) 
    return 5 * np.log10(dl) + 25

##############################################################################

def rchi2(obs, exp=sndat, method='formula'):
    if method == 'formula':
        chi = np.sum((obs - exp)**2)
        answer = chi / len(obs)

    elif method == 'poly':
        p, residuals, _, _, _ = np.polyfit(z_sn, exp, 2, full=True)
        answer = residuals / (len(exp) - 2)
    else:
        raise Exception('Bad "method" input for rchi2 method. Options are "formula" or "poly"')
    
    return float(answer)

##############################################################################

def read_model_data(fname, fdir='../Data/model_data/'):
    """
    Read model data quickly
    """
    file = pd.read_csv(fdir+fname, delim_whitespace=True, comment='#')
    return file

#

def specific_function(array, number):
    """
    Function to do something specific and time consuming
    """
    chi = array['chi']
    beta = array['beta']
    kappa = array['kappa']

    chi_chosen = np.sort(chi)[number]
    beta_chosen = float(np.array(beta)[np.where(chi == chi_chosen)])
    kappa_chosen = float(np.array(kappa)[np.where(chi == chi_chosen)])

    return beta_chosen, kappa_chosen