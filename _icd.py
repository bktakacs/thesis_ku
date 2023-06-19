# Imports, Constants, Data

# Imports
from scipy.integrate import solve_ivp, odeint, quad
from scipy.stats import chisquare
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import LogNorm, NoNorm
from mpl_toolkits.axes_grid1 import ImageGrid, make_axes_locatable
import pandas as pd
from tqdm import tqdm
from tqdm.contrib import itertools
from itertools import product
from astropy.cosmology import FlatLambdaCDM
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
c  = 3e5        # celeritas [km/s]
h0 = 73.6       # hubbles constant [km/s/Mpc]
dh = c / h0     # hubble distance [Mpc]
mat0 = 0.334    # matter density
rad0 = 1e-4     # radiation density
lam0 = 0.666    # dark energy density


# Data
df      = pd.read_csv('../../Data/Pantheon+SH0ES.dat.txt', sep=' ')
sndat   = np.array(df['MU_SH0ES'])  # supernova data
snerr   = np.array(df['MU_SH0ES_ERR_DIAG']) # supernova error
z_sn    = np.array(df['zHD'])    # redshift
z_err   = np.array(df['zHDERR'])    # redshift error

# dat     = pd.read_csv('../../Data/lcparam_full_long.txt', sep=' ')

# distance modulus from astropy for comparison
cosmo = FlatLambdaCDM(H0=h0, Om0=mat0)
dm_astro = cosmo.distmod(z_sn).value