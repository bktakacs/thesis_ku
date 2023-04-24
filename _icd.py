# Imports, Constants, Data

# Imports
from scipy.integrate import solve_ivp, odeint, quad, RK45, cumtrapz, simpson
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import LogNorm, NoNorm
import pandas as pd
from tqdm import tqdm
from tqdm.contrib import itertools
# from astropy.cosmology import Planck18
from astropy.cosmology import FlatLambdaCDM
from functools import lru_cache
from time import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Plot settings
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 200
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.figsize'] = (10, 5)


# Constants
h0 = 73.6
c  = 3e5        # celeritas [km/s]
dh = c / h0
mat0 = 0.334
rad0 = 1e-4
lam0 = 0.666


# Data
df      = pd.read_csv('../../Data/Pantheon+SH0ES.dat.txt', sep=' ')
sndat   = np.array(df['MU_SH0ES'])
snerr   = np.array(df['MU_SH0ES_ERR_DIAG'])
z_sn    = np.array(df['zHD'])
z_err   = np.array(df['zHDERR'])

# dat     = pd.read_csv('../Data/lcparam_full_long.txt', sep=' ')

cosmo = FlatLambdaCDM(H0=h0, Om0=mat0)
dm_astro = cosmo.distmod(z_sn).value