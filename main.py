# Main file

# Import
from _icd import *
from _functions import modified_friedmann, acceleration, dm_z_o4, rchi2,\
                       read_model_data, specific_function


# Classes
class model():
    """
    class which creates a cosmological model with given parameters
    
    This class takes initial inputs and then integrates the modified Friedmann
    equations to get a, adot (a1) and t. From there we derive further time
    derivatives as well as other cosmological parameters.

    Calling this class with default arguments, model(), returns a flat LCDM
    model. A matter dominated model is returned by model(lam=0.)
    """

    def __init__(self, a_start: float = 1e-3, mat: float = mat0, 
                 rad: float = rad0, lam: float = lam0, beta: float = 3.,
                 kappa: float = 0., n: float = 1, xaxis: str = 'a', 
                 xmax: float = 1.5, xlen: int = 50000):
        """
        initialization method for model class

        Parameters
        ----------
        a_start : float, optional
            initial value of scale factor at which to begin integration,
            chosen 1e-3 to correspond to time of recombination (z ~ 1100)
            The default is 1e-3.
        mat : float, optional
            dimensionless matter density parameter (baryons + dark matter),
            defined in _icd. The default is mat0.
        rad : float, optional
            dimensionless radiation density parameter (photons + neutrinos),
            defined in _icd. The default is rad0.
        lam : float, optional
            dimensionless dark energy density parameter (cosmological
            constant), defined in _icd. The default is lam0.
        beta : float, optional
            power to which mass is raised in new DM force. The default is 3.
        kappa : float, optional
            constant in front of new DM force. The default is 0..
        n : float, optional
            not sure, but has to do with power series in Press-Schechter
            formalism, used in variable gamma. The default is 1.
        xaxis : str, optional
            axis upon which we plan to plot, can be time 't' or scale factor
            'a', scale factor and its derivatives only run until xaxis
            variable reaches unity. The default is 'a'.
        xmax : float, optional
            upper limit of integration in units of time tH0, chosen as 1.5 to
            ensure scale factor can reach unity if xaxis = 'a'. The default is
            1.5.
        xlen : int, optional
            length of time array for integration, chosen to be sufficiently
            large so that when we go to redshift there are values similar to
            those of the z_sn array (1e-3, 2.2). The default is 50000.
        """
        # properties of this model
        self.m = mat
        self.r = rad
        self.l = lam
        self.b = beta
        self.k = kappa
        self.axis = xaxis

        # gamma and p variables used in acceleration equation
        gamma  = 1 + n / 3
        self.p = gamma / (2 * beta)
        
        # set initial condtions from a_start
        a1_start = np.sqrt(self.m * a_start**-1 + self.r * a_start**-2 + \
                                                        self.l * a_start**2)
        initial_conditions = (a_start, a1_start)

        # time array upon which to integrate
        time_array = np.linspace(0, xmax, xlen)
        # time_array = np.logspace(-5, np.log10(xmax), xlen)

        # integrate modified friedmann equations
        solution = solve_ivp(modified_friedmann, y0=initial_conditions,
                             t_span=(0, xmax), t_eval=time_array,
                             args=(self.m, self.r, self.l, self.k, self.p))
        # following lines are for odeint
        # solution = odeint(modified_friedmann, y0=initial_conditions,
        #                   t=time_array,
        #                   args=(self.m, self.r, self.l, self.k, self.p),
        #                   tfirst=True)

        # extract model parameters
        self.a  = solution.y[0, :]
        self.a1 = solution.y[1, :]
        self.t  = solution.t
        # following lines are for odeint
        # self.a  = solution[:, 0]
        # self.a1 = solution[:, 1]
        # self.t  = np.copy(time_array)

        # time derivatives
        self.a2 = acceleration(self.a, self.a1, self.m, self.r, 
                               self.l, self.k, self.p)
        self.a3 = np.diff(self.a2) / np.diff(self.t)
        self.a4 = np.diff(self.a3) / np.diff(self.t[:-1])

        # find where today is in terms of a or t
        hoy = np.argmin(np.abs(1 - self.a)) if xaxis == 'a' else \
                                                  np.argmin(np.abs(1 - self.t))
        
        # truncate arrays to today
        self.a = self.a[:hoy]
        self.a1 = self.a1[:hoy]
        self.a2 = self.a2[:hoy]
        self.a3 = self.a3[:hoy]
        self.a4 = self.a4[:hoy]
        self.t = self.t[:hoy]

        # parameters from time derivatives
        self.q = - self.a2[-1] * self.a[-1] * self.a1[-1]**-2  # decel param q
        self.j = self.a3[-1] * self.a[-1]**2 * self.a1[-1]**-3 # jerk j
        self.s = self.a4[-1] * self.a[-1]**3 * self.a1[-1]**-4 # snap s


    def norm(self, matter: object):
        """
        norm() method which normalizes the scale factor acceleration
        to that of a matter only model

        Parameters
        ----------
        matter : object
            matter class object
        """

        if type(matter) != model:
            raise TypeError("matter must be a class object")

        self.a2norm = self.a2 / np.interp(self.a, matter.a, matter.a2)

    
    def distance_modulus(self, effort: bool = True):
        """
        distance_modulus() method which calculates the distance modulus of
        scale factor values calculated in the initialization in two ways.
        Returns array

        Parameters
        ----------
        effort : bool, optional
            if True, uses the "true" method of calculating the distance
            modulus, if False, uses a faster method. The default is True.
        
        Note on effort parameter:
        "True" way is calculating f and using the integrand function. "False"
        way method calculates E(z) as an array of data and simply integrates
        that. It is a bit mysterious however because it uses the inverse
        function which returns x**(-0.5) when it should just be x**(-1). For
        some reason though, -0.5 just works.
        
        The distance modulus is calculated through two methods, here called
        integration and taylor: DM through integration is calculated through a
        series of equations as follows: dc = dh * int_0^z dz/E(z) where
        E(z) = sqrt(m(1+z)(1+kf) + r(1+z)^4 + l) ; dl = (1+z) * dh ;
        dm = 5 * np.log10(dl) + 25 (for dl given in Mpc). Because z from the
        model will not exactly equal z from the SN data, we interpolate along z
        from SNe
        
        DM through taylor is calculated using a Taylor approximation of the
        scale factor which is a long equation that's written below

        At the end we add the dm_int and dm_tay attributes corresponding to the
        DM calculated from integration and taylor respectively
        """

        if type(effort) != bool:
            raise ValueError('effort must be True or False')

        z = (1 / self.a - 1)
        idx0 = np.argmin(np.abs(z - z_sn[-1]))
        idxf = np.argmin(np.abs(z - z_sn[0]))

        if effort:
            # calculate f
            integral_f = np.zeros_like(z)
            for index, scale in enumerate(self.a):
                integral_f[index] = quad(lambda x, y: y**-2, 0, scale,
                                         (self.a1[index],))[0]
            f = (self.a1 / self.a * integral_f)**(self.p)

            # calculate dl
            integrand = lambda x, m, r, l, f, k: \
                (m * (1 + x) * (1 + k * f) + r * (1 + x)**4 + l)**-0.5
            dl = np.zeros_like(z)
            for index, redshift in enumerate((z)):
                dl[index] = (1 + redshift) * dh * quad(integrand, 0, redshift,
                                                       (self.m, self.r, self.l,
                                                        f[index], self.k))[0]

        elif not effort:
            # calculate E(z)
            ez = self.a1 / self.a
            
            # calculate dl
            dl = np.zeros_like(z)
            for index, redshift in enumerate(z):
                dl[index] = (1 + redshift) * dh * quad(lambda x, y: y**-0.5, 0,
                                                       redshift,
                                                       (ez[index],))[0]
    
        # calculate dm
        dm_int = 5 * np.log10(dl) + 25
        dm_int = np.flip(dm_int[idx0:idxf])
        z = np.flip(z[idx0:idxf])

        # try to interpolate, else nan
        try:
            self.dm_int = np.interp(z_sn, z, dm_int)
        except:
            self.dm_int = np.nan

        # calculate dm_tay
        self.dm_tay = dm_z_o4(z=z_sn, q=self.q, j=self.j, s=self.s)

    
    def chi_value(self):
        """
        Method to calculate the chi squared value of the model compared to the
        SN data set. Returns chi^2 value fo both integration and Taylor methods
        in the chi_int and chi_tay attributes respectively
        """

        self.chi_int = np.nan if np.isnan(self.dm_int).all() else \
                                rchi2(obs=self.dm_int)
        self.chi_tay = np.nan if np.isnan(self.dm_tay).all() else \
                                rchi2(obs=self.dm_tay)


    def plot(self, which: str, lcdm: object, matter: object):
        """
        This method plots the acceleration of the scale factor or the distance
        modulus, depending on the input "which". When plotting acceleration it
        normalizes the acceleration of the model to the acceleration of the
        matter model and plots the acceleration of the LCDM model for
        comparison

        Parameters
        ----------
        which : str
            options are 'acc' or 'dm' corresponding to acceleration or
            distance modulus
        lcdm : object
            model object, used for comparison
        matter : object
            model object, used for normalization of acc plot
        """

        if which not in ('acc', 'dm'):
            raise Exception('Bad input for "which" in "plot" method in "model"'
                            'class. Input must be "acc" or "dm".')
        
        if type(lcdm) != model and type(matter) != model:
            raise Exception('Bad input for "lcdm" and "matter" in "plot"'
                            'method in "model" class. Input must be a model'
                            'object.')

        # normalize acceleration
        self.norm(matter)
        lcdm.norm(matter)

        x_mod = self.a if self.axis == 'a' else self.t
        x_lcdm = lcdm.a if lcdm.axis == 'a' else lcdm.t

        if which == 'acc':

            plt.figure()
            plt.plot(x_mod, self.a2norm, c='r', ls='-', 
                     label=r'Alt. Model, $\beta={:.2f}, k={:.2f}$'
                     ''.format(self.b, self.k))
            
            plt.plot(x_lcdm, lcdm.a2norm, c='k', ls='--',
                     label=r'$\Lambda$CDM Model')
            
            plt.xlabel(r'$a$')
            plt.ylabel(r'$\ddot{a}/\ddot{a}_{\mathrm{M}}$')
            plt.ylim([-5, 2])
            plt.tick_params(axis='both', which='both', direction='in',
                            bottom=True, top=True, left=True, right=True)
            plt.legend(loc='lower left')
            plt.grid()
            plt.show()

        elif which == 'dm':
            
            fig = plt.figure(constrained_layout=False)
            frame1 = fig.add_axes((.1, .3, .8, .6))
            plt.errorbar(z_sn, sndat, yerr=snerr, lw=0.5, ls='', marker='.',
                         markersize=2, label=r'Pantheon+SH0ES', zorder=0)

            plt.plot(z_sn, dm_astro, c='orange', ls='-',
                     label=r'$DM$ from flat-$\Lambda$CDM Cosmology, '\
                    r'$\chi^{{2}}_{{r}}={:.4f}$'.format(rchi2(obs=dm_astro)))
            
            plt.plot(z_sn, self.dm_int, c='k', ls='-.', label=r'$DM(z, E(z))$,'
                    r' $\chi^{{2}}_{{r}}={:.4f}$'.format(
                                                    rchi2(obs=self.dm_int)))
            
            plt.plot(z_sn, self.dm_tay, c='r', ls='--',
                     label=r'$DM(z, q, j, s)$, '
                        r'$\chi^{{2}}_{{r}}={:.4f}$'.format(
                                                    rchi2(obs=self.dm_tay)))
            
            plt.ylabel(r'$DM$ [mag]')
            plt.xscale('log')
            plt.tick_params(axis='both', which='both', direction='in',
                            bottom=True, top=True, left=True, right=True)
            plt.legend()
            plt.grid()
            frame1.set_xticklabels([])

            frame2 = fig.add_axes((.1, .1, .8, .2))
            plt.errorbar(z_sn, sndat - dm_astro, yerr=snerr, lw=0.5, ls='',
                         marker='.', markersize=2, label=r'Supernova',zorder=0)
            plt.plot(z_sn, (self.dm_int - dm_astro)/dm_astro*100, c='k',
                     ls='-.')
            plt.plot(z_sn, (self.dm_tay - dm_astro)/dm_astro*100, c='r',
                     ls='--')
            plt.xlabel(r'$z$')
            plt.xscale('log')
            plt.ylabel(r'$\%\Delta{DM}_{\mathrm{flat}\Lambda\mathrm{CDM}}$')
            plt.tick_params(axis='both', which='both', direction='in',
                            bottom=True, top=True, left=True, right=True)
            plt.tick_params(axis='x', pad=7)
            extraticks = [-1, 1]
            plt.yticks(list(plt.yticks()[0]) + extraticks)
            plt.ylim([-2, 2])
            # plt.ylim([-1, 1])
            plt.grid(which='major', axis='both')
            plt.show()


# Functions
def chi_comp(parameter: str, space: list, beta: float = 3., kappa: float = 0.,
             lam: float = lam0, dm_effort: bool = False, dm_method: str='int',
             plot: bool = True):
    """
    This function calculates the chi^2 value for a given parameter space and
    plots the chi^2 value as a function of the parameter space. It also prints
    the lowest chi^2 value and the corresponding parameter value.

    Parameters
    ----------
    parameter : str
        string, options are 'l', 'b' or 'k', corresponding to lambda, beta or
        kappa
    space : list
        parameter space to be explored
    beta : float
        beta value
    kappa : float
        kappa value
    lam : float
        lambda value
    dm_effort : bool
        if True, uses the "true" method of calculating the distance modulus,
        if False, uses a faster method
    dm_method : str
        options are 'int' or 'tay', corresponding to the integration method or
        Taylor expansion method
    plot : bool
        if True, plots the chi^2 value as a function of the parameter space
    """

    if parameter not in ('l', 'b', 'k'):
        raise Exception('Bad "parameter" input in "chi_comp" function. '
                        'Inputs are "l", "b", or "k".')

    array = np.zeros_like(space)

    if parameter == 'l':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=value, beta=beta, kappa=kappa)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi_value()
            array[index] = tmod.chi_int if dm_method == 'int' \
                                        else tmod.chi_tay
        
        model_optimized = model(lam=space[np.argmin(array)], beta=beta,
                                kappa=kappa)

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f}'
              ' and Omg_Lambda = {:.3f} for Omg_Lambda in the range '
              '[{:.1f}, {:.1f}]'.format(
                        np.min(array), beta, kappa, space[np.argmin(array)],
                        space[0], space[-1]))

    elif parameter == 'b':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=lam, beta=value, kappa=kappa)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi_value()
            array[index] = tmod.chi_int if dm_method == 'int' \
                                        else tmod.chi_tay
        
        model_optimized = model(lam=lam, beta=space[np.argmin(array)],
                                kappa=kappa)

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f}'
              ' and Omg_Lambda = {:.3f} for beta in the range [{:.1f}, {:.1f}]'
                ''.format(np.min(array), space[np.argmin(array)], kappa, lam,
                          space[0], space[-1]))
    
    elif parameter == 'k':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=lam, beta=beta, kappa=value)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi_value()
            array[index] = tmod.chi_int if dm_method == 'int' \
                                        else tmod.chi_tay
        
        model_optimized = model(lam=lam, beta=beta,
                                kappa=space[np.argmin(array)])

        print('The lowest chi^2 value is {:.5f} for beta = {:.3f}, k = {:.3f} '
              'and Omg_Lambda = {:.3f} for kappa in the range [{:.1f}, {:.1f}]'
              ''.format(np.min(array), beta, space[np.argmin(array)], lam,
                        space[0], space[-1]))

    if plot:
        xlab = r'$\beta$' if parameter == 'b'\
                else r'$k$' if parameter == 'k'\
                else r'$\Omega_{\Lambda}$'
        yscale = 'log' if np.max(array)/np.min(array) > 50 else 'linear'
        plotlab = r'$k={:.3f},\Omega_{{\Lambda}}={:.1f}$' if parameter == 'b'\
                    else r'$\beta={:.3f},\Omega_{{\Lambda}}={:.1f}$'\
                    if parameter == 'k' else r'$\beta={:.3f},k={:.3f}$'
        plotform = (kappa, lam) if parameter == 'b'\
                    else (beta, lam) if parameter == 'k'\
                        else (beta, kappa)

        plt.figure()
        plt.plot(space, array, c='k', ls='-', label=plotlab.format(*plotform))
        plt.plot([space[np.argmin(array)], space[np.argmin(array)]],
                 [0.9*np.min(array), 1.01*np.max(array)], c='r', ls='--',
                 label=xlab + r' $= {:.3f}$'.format(space[np.argmin(array)]))
        plt.xlabel(xlab)
        plt.ylabel(r'$\chi^{2}_{r}$')
        plt.yscale(yscale)
        plt.grid()
        plt.legend(loc='best', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        plt.show()

    model_optimized.distance_modulus(effort=dm_effort)
    model_optimized.chi_value()

    return model_optimized


def chi_search(fname: str, length: int = 10, blim: tuple = (2., 4.),
               klim: tuple = (1., 10.), lam: int = 0., dm_effort: bool = False,
               dm_method: str = 'int', plot: bool = True, round: int = 1,
               scale = LogNorm(), double_eval: bool = False,
               fdir: str='../../Data/model_data/'):
    """
    chi_search() function. Calculates chi^2 value for models with different
    combinations (beta, kappa). This function creates a linear range of beta
    values using length & blim, same for kappa, and then loops over
    combinations of those. For each combination a model is created and the
    distance_modulus method called using dm_effort input. Then chi_value()
    method is called. Chi^2 value is stored and shaped into (length, length)
    array for plotting. Finally stores data using inputted file name and then
    prints statement of results.

    Parameters
    ----------
    fname : str
        Name of file to save data to
    length : int, optional
        Length of array with beta or kappa values. Length^2 is number of
        iterations in loop. The default is 10.
    blim : tuple, optional
        Upper and lower bounds of beta. The default is (2., 4.).
    klim : tuple, optional
        Upper and lower bounds of kappa. The default is (1., 10.).
    lam : int, optional
        Lambda value to be used for each model. The default is 0.
    dm_effort : bool, optional
        Whether to use effort or not in calculation of distance modulus (used
        later in calculating chi^2). The default is False.
    dm_method : str, optional
        'int' or 'tay', which method to use in evaluating distance modulus.
        The default is 'int'.
    plot : bool, optional
        Plot chi^2 heat map or no. The default is True.
    round : int, optional
        Used when plot==True, what number of decimal places to round x & y
        labels to, for visual purposes. The default is 1.
    scale : matplotlib.colors.Normalize, optional
        Used when plot==True, what scale to use for colour map. The default is
        LogNorm().
    double_eval : bool, optional
        Whether to evaluate on both chi^2 values for each model.
    fdir : str, optional
        Directory to save data to. The default is '../../Data/model_data/'.
    """

    if fname == 'nosave':
        save = False
    else:
        fname = fname + '.txt'
        save = True

    brange = np.linspace(np.min(blim), np.max(blim), length)
    krange = np.linspace(np.min(klim), np.max(klim), length)

    chival_int = np.zeros(length**2)
    chival_tay = np.zeros(length**2)

    nan_count = 0

    matter = model(lam=0.)
    top_mod = model(lam=lam, beta=brange[0], kappa=krange[0])
    top_mod.distance_modulus(effort=dm_effort)
    top_mod.chi_value()
    running_chi = (top_mod.chi_int + top_mod.chi_tay) if double_eval \
            else top_mod.chi_int if dm_method == 'int' else top_mod.chi_tay

    # Iterate over all (b, k), store chi value for each
    for index, param in enumerate(itertools.product(brange, krange)):
        tmod = model(lam=lam, beta=param[0], kappa=param[1])
        tmod.norm(matter=matter)
        if (np.max(tmod.a2norm) > 3 or np.min(tmod.a2norm) < -10):
            chival_int[index] = np.nan
            chival_tay[index] = np.nan
            nan_count += 1
        else:
            # if model is physical, store chi values
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi_value()
            chival_int[index] = tmod.chi_int
            chival_tay[index] = tmod.chi_tay
            # update top model if desired chi < running_chi
            top_mod = tmod if ((tmod.chi_int + tmod.chi_tay) < running_chi and
                               double_eval) \
                            else tmod if (tmod.chi_int < running_chi and
                                          dm_method == 'int') \
                            else tmod if (tmod.chi_tay < running_chi and
                                          dm_method == 'tay') else top_mod
    
    # Raise exception if only NaN values were returned
    if nan_count == length**2:
        raise Exception('No real chi values found')
    nan_ratio = nan_count / length**2 * 100

    chival = np.copy(chival_int) if dm_method == 'int' else \
             np.copy(chival_tay) if dm_method == 'tay' else np.full_like(
                                                            chival_int, np.nan)
    
    # Convert NaNs to 1e6 because otherwise it's dumb (number not important)
    chival_nonan = np.nan_to_num(chival, nan=1e6)
    chival_int_nonan = np.nan_to_num(chival_int, nan=1e6)
    chival_tay_nonan = np.nan_to_num(chival_tay, nan=1e6)
    chival_com_nonan = chival_int_nonan + chival_tay_nonan
    lowest = np.argmin(chival_com_nonan) if double_eval else \
                                                        np.argmin(chival_nonan)

    # Find beta and kappa values for lowest chi^2 value
    chi_low_int = chival_int_nonan[lowest]
    chi_low_tay = chival_tay_nonan[lowest]
    chi_low = chival[lowest]
    beta_low  = np.repeat(brange, length)[lowest]
    kappa_low = np.tile(krange, length)[lowest]

    # Print results of function call
    print('The lowest chi^2 values are {:.5f} (int) {:.3f} (tay) for beta'
          ' = {:.3f} & mu = {:.3f} in the range\n\t'
          '{:.2f} < beta < {:.2f} and {:.2f} < mu < {:.2f}\n\t'
          '{:.1f} % of models had a chi^2 value of NaN. \n'
          ''.format(chi_low_int, chi_low_tay, beta_low, kappa_low,
                    np.min(blim), np.max(blim), np.min(klim), np.max(klim),
                    nan_ratio) if double_eval else \
          'The lowest chi^2 value is {:.5f} for beta = {:.3f} & mu = '
          '{:.3f} in the range'
          '\n{:.2f} < beta < {:.2f} and {:.2f} < mu < {:.2f}\n\t'
          '{:.1f} % of models had a chi^2 value of NaN. \n'
          ''.format(chi_low, beta_low, kappa_low, np.min(blim), np.max(blim),
                    np.min(klim), np.max(klim), nan_ratio))

    # Plot heat map of chi values
    if plot and double_eval: 
        chi_plot_z_int = np.reshape(chival_int, (length, length)).T
        chi_plot_z_tay = np.reshape(chival_tay, (length, length)).T 
        chi_plot_z = np.reshape(chival_int + chival_tay, (length, length)).T

        fig, ax = plt.subplots(2, 2)
        fig.delaxes(ax[1,1])
        f1 = ax[0,0]; f2 = ax[0,1]; f3 = ax[1,0]
        cmap = matplotlib.cm.get_cmap('viridis_r').copy()
        cmap.set_bad(color='r')

        # int
        im1 = f1.imshow(chi_plot_z_int, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f1.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='r', ha='center', va='center', fontsize=16)
        f1.set_xticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), round),
                      rotation=45, fontsize=12)
        f1.set_yticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(krange[0], krange[-1], 10), round),
                      fontsize=12)
        f1.set_xlabel(r'$\beta$')
        f1.set_ylabel(r'$k$')
        f1.set_title(r'Integration method')
        f1.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        f1.grid()
        cbar = fig.colorbar(im1, ax=f1)
        cbar.set_label(r'$\chi^2$', rotation='90')

        # tay
        im2 = f2.imshow(chi_plot_z_tay, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f2.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='r', ha='center', va='center', fontsize=16)
        f2.set_xticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), round),
                      rotation=45, fontsize=12)
        f2.set_yticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(krange[0], krange[-1], 10), round),
                      fontsize=12)
        f2.set_xlabel(r'$\beta$')
        f2.set_ylabel(r'$k$')
        f2.set_title(r'Taylor expansion')
        f2.tick_params(axis='both', which='both', direction='in',
                       bottom=True, top=True, left=True, right=True)
        f2.grid()
        fig.colorbar(im2, label=r'$\chi^{2}_{r}$', ax=f2)

        # both
        im3 = f3.imshow(chi_plot_z, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f3.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='r', ha='center', va='center', fontsize=16)
        f3.set_xticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), round),
                      rotation=45, fontsize=12)
        f3.set_yticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(krange[0], krange[-1], 10), round),
                      fontsize=12)
        f3.set_xlabel(r'$\beta$')
        f3.set_ylabel(r'$k$')
        f3.set_title(r'Combined')
        f3.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        f3.grid()      
        fig.colorbar(im3, label=r'$\chi^{2}_{r}$', ax=f3)
        plt.show()

    elif plot:
        chi_plot_z = np.reshape(chival, (length, length)).T

        fig, ax = plt.subplots()
        cmap = matplotlib.cm.get_cmap('viridis_r').copy()
        # cmap = matplotlib.cm.get_cmap('binary').copy()
        cmap.set_bad(color='r')
        im = ax.imshow(chi_plot_z, cmap=cmap, origin='lower',
                       interpolation='nearest', norm=scale)
        ax.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='r', ha='center', va='center', fontsize=20)
        ax.set_xticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), 1),
                      rotation=45, fontsize=14)
        ax.set_yticks(np.linspace(0, length-1, 10),
                      np.round(np.linspace(krange[0], krange[-1], 10), 1),
                      fontsize=14)
        ax.set_xlabel(r'$\beta$')
        ax.set_ylabel(r'$k$')
        cbar = fig.colorbar(im)
        cbar.set_label(r'$\chi^{2}_{r}$', rotation='90')
        ax.tick_params(axis='both', which='both', direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.grid()
        plt.show()

    # Save chi, beta, kappa values to file
    if save:
        f_chi = np.copy(chival) if not double_eval else np.copy(chival_int
                                                                +chival_tay)
        f_beta = np.repeat(brange, length)
        f_kappa = np.tile(krange, length)
        f_save = np.vstack((f_chi, f_beta, f_kappa)).T
        f_comment = '#Results of "chi_search" called with the following' +\
                'inputs:\n' +\
                '#length={}, blim=({}, {}), klim=({}, {}), lambda={}, '.format(
                    length, blim[0], blim[1], klim[0], klim[1], lam) +\
                'effort={}, dm_method={}, \n'.format(
                                        dm_effort, dm_method) +\
                '#Lowest chi^2 was with beta = {} & k = {}\n'.format(
                                                        beta_low, kappa_low)
        np.savetxt(fname=fdir+fname, X=f_save, header='chi beta kappa',
                   delimiter=' ', comments=f_comment)

    # Return the optimized model
    return top_mod

def chi_search_a(fname: str, length: int = 10, blim: tuple = (2., 4.),
                 klim: tuple = (1., 10.), lam: int = 0., plot: bool = True, 
                 fdir: str = '../../Data/model_data/'):
    """
    chi_search_a() function. Iterates over all combinations (beta, kappa) and
    for each creates a model and calculates the chi^2 fit of the acceleration
    to that of the LCDM model. Just a test and I'm curious to see what comes
    out of it.

    Parameters
    ----------
    fname : str
        file name to save results to
    length : int, optional
        length of beta and kappa ranges, by default 10
    blim : tuple, optional
        limit of beta values to plot, by default (2., 4.)
    klim : tuple, optional
        limit of kappa values to plot, by default (1., 10.)
    lam : int, optional
        lambda value to use, by default 0.
    plot : bool, optional
        if True, plot best fit acceleration, by default True
    fdir : str, optional
        directory to save results to, by default '../../Data/model_data/'
    """

    if len(fname) == 0:
        raise Exception('fname must be a string of at least one character')
    
    save = False if fname == 'nosave' else True

    if blim[0] > blim[1] or klim[0] > klim[1]:
        raise Exception('blim and klim must be increasing tuples')

    fname += '.txt'

    # Create beta and kappa ranges
    brange = np.linspace(blim[0], blim[1], length)
    krange = np.linspace(klim[0], klim[1], length)

    # Create empty arrays to store chi^2 values
    chival = np.zeros(length**2)

    # Initialize some things
    matter = model(lam=0.)
    lcdm = model()
    lcdm.norm(matter)
    nan_count = 0
    top_mod = model(lam=lam, beta=blim[0], kappa=klim[0])
    top_mod.norm(matter=matter)
    top_mod_a2norm_interp = np.interp(lcdm.a, top_mod.a, top_mod.a2norm)
    running_chi = rchi2(top_mod_a2norm_interp, lcdm.a2norm)

    # Iterate over all combinations of beta and kappa
    for index, param in enumerate(itertools.product(brange, krange)):
        tmod = model(lam=lam, beta=param[0], kappa=param[1])
        tmod.norm(matter=matter)
        if (np.max(tmod.a2norm) > 3 or np.min(tmod.a2norm) < -10):
            chival[index] = np.nan
            nan_count += 1
        else:
            tmod_a2norm_interp = np.interp(lcdm.a, tmod.a, tmod.a2norm)
            chival[index] = rchi2(tmod_a2norm_interp, lcdm.a2norm)
            top_mod = tmod if chival[index] < running_chi else top_mod

    # Raise exception if only Nans are returned
    if nan_count == length**2:
        raise Exception('Only NaNs returned. Check input parameters.')
    nan_ratio = nan_count / length**2 * 100

    # Find lowest chi^2 value and corresponding beta and kappa values
    chi_low = np.sort(chival)[0]
    lowest = np.argmin(np.abs(chival - chi_low))
    beta_low = np.repeat(brange, length)[lowest]
    kappa_low = np.tile(krange, length)[lowest]

    # Print results
    print('The lowest chi^2 value is {:.5f} for beta = {:.3f} & mu = {:.3f} in'
          ' the range\n{:.2f} < beta < {:.2f} and {:.2f} < mu < {:.2f}\n\t'
          '{:.1f} % of models had a chi^2 value of NaN. \n'
          ''.format(chi_low, beta_low, kappa_low, blim[0], blim[1], 
                    klim[0], klim[1], nan_ratio))
    
    # Plot results
    if plot:
        plt.figure()
        plt.plot(top_mod.a, top_mod.a2norm, c='r', ls='-', 
                    label=r'Alt. Model, $\beta={:.2f}, k={:.2f}$'.format(
                                                        top_mod.b, top_mod.k))
        
        plt.plot(lcdm.a, lcdm.a2norm, c='k', ls='--',
                 label=r'$\Lambda$CDM Model')

        plt.xlabel(r'$a$')
        plt.ylabel(r'$\ddot{a}/\ddot{a}_{\mathrm{M}}$')
        plt.ylim([-5, 2])
        plt.tick_params(axis='both', which='both', direction='in', bottom=True,
                        top=True, left=True, right=True)
        plt.legend(loc='lower left')
        plt.grid()
        plt.show()

    # Save results (maybe)
    if save:
        f_chi = chival
        f_beta = np.repeat(brange, length)
        f_kappa = np.tile(krange, length)
        f_save = np.vstack((f_chi, f_beta, f_kappa)).T
        f_comment = '#Results of "chi_search_a" called with the following' +\
                'inputs:\n' +\
                '#length={}, blim=({}, {}), klim=({}, {}), lambda={}\n'.format(
                    length, blim[0], blim[1], klim[0], klim[1], lam) +\
                '#Lowest chi^2 was with beta = {} & k = {}\n'.format(
                    beta_low, kappa_low)
        np.savetxt(fname=fdir+fname, X=f_save, header='chi beta kappa',
                   delimiter=' ', comments=f_comment)
    
    top_mod.norm(matter=matter)

    return top_mod


def q_surface(length: int = 20, blim: tuple = (2., 4.),
              klim: tuple = (1., 10.), qlim: tuple = (-1., 0.),lam: float = 0.,
              dm_method: str = 'int', dm_effort: bool = False,
              splot: bool = True, mplot: bool = True):
    """
    q_surface() function. Plots a surface of q values for a given range of
    beta and kappa values.

    Parameters
    ----------
    length : int
        number of points to plot in each direction
    blim : tuple
        range of beta values to plot
    klim : tuple
        range of kappa values to plot
    qlim : tuple
        range of q values to plot
    lam : float
        lambda value to use for model
    dm_method : str
        method to use for distance modulus calculation
    dm_effort : bool
        if True, use more accurate dist mod calculation
    splot : bool
        if True, plot surface plot of q values
    mplot : bool
        if True, plot acceleration and dist mod of opt model
    """

    if (blim[0] > blim[1] or klim[0] > klim[1] or qlim[0] > qlim[1]):
        raise Exception('blim, klim, and qlim tuples must be of form (a, b) '
                        'where a < b')

    brange = np.linspace(np.min(blim), np.max(blim), length)
    krange = np.linspace(np.min(klim), np.max(klim), length)
    q_save = np.zeros(length**2)

    qcd = []
    bcd = []
    kcd = []
    xcd = []
    current_chi = 1

    for index, param in enumerate(itertools.product(brange, krange)):
        temp_mod = model(lam=lam, beta=param[0], kappa=param[1])
        q_save[index] = temp_mod.q

        if (qlim[0] < temp_mod.q < qlim[1]):
            temp_mod.distance_modulus(effort=dm_effort)
            temp_mod.chi_value()

            qcd.append(temp_mod.q)
            bcd.append(temp_mod.b)
            kcd.append(temp_mod.k)
            xcd.append(temp_mod.chi_int if dm_method == 'int'
                                        else temp_mod.chi_tay)

            top_mod = temp_mod if xcd[-1] < current_chi else top_mod
            current_chi = xcd[-1] if xcd[-1]< current_chi else current_chi

    if len(qcd) < 1:
        raise Exception('No q values found within the range {} < q < {}'
                        ''.format(np.min(qlim), np.max(qlim)))
    
    q_save = np.reshape(q_save, (length, length))

    qround = np.round(qcd, 2)   # round down to 2 decimal places
    qred = []
    bred = []
    kred = []
    xred = []

    for i, x in enumerate(qround):
        if x not in qred:
            qred.append(x)
            bred.append(bcd[i])
            kred.append(kcd[i])
            xred.append(xcd[i])

    sortorder = np.copy(qred)
    qred,   bred = zip(*sorted(zip(sortorder, bred)))
    _,      kred = zip(*sorted(zip(sortorder, kred)))
    _,      xred = zip(*sorted(zip(sortorder, xred)))

    for i in range(len(qred)):
        if xred[i] == np.min(xred):
            print('   **** q = {:.3f}\tfor\tbeta = {:.3f}\tk = {:.3f}\twith'
                  'chi^2 = {:.4f} *****'.format(
                                        qred[i], bred[i], kred[i], xred[i]))
        else:
            print('\tq = {:.3f}\tfor\tbeta = {:.3f}\tk = {:.3f}\twith'
                  'chi^2 = {:.4f}'.format(
                                        qred[i], bred[i], kred[i], xred[i]))
            
    
    if splot:
        bplot = np.reshape((np.tile(brange, length)), (length, length)).T
        kplot = np.reshape((np.tile(krange, length)), (length, length)).T
        X, Y = np.meshgrid(bplot[:, 0], kplot[:, 0])
        fig = plt.figure(figsize=plt.figaspect(0.3), constrained_layout=False)
        ax = fig.add_subplot(1, 3, 1, projection='3d')
        ax.plot_surface(X, Y, q_save, rstride=1, cstride=1, antialiased=True)
        ax.set_xlabel(r'$\beta$', labelpad=15)
        ax.set_ylabel(r'$k$', labelpad=15)
        ax.set_zlabel(r'$q$')
        ax.set_zlim([-2, 0])
        ax.set_title(r'$\beta, k, q$')
        ax.view_init(elev=20, azim=-45)

        ax1 = fig.add_subplot(1, 3, 2, projection='3d')
        ax1.set_yticklabels('')
        ax1.plot_surface(X, Y, q_save, rstride=1, cstride=1, antialiased=True)
        ax1.set_xlabel(r'$\beta$', labelpad=15)
        ax1.set_zlabel(r'$q$')
        ax1.set_zlim([-2, 0])
        ax1.set_title(r'$\beta, q$')
        ax1.view_init(elev=0, azim=-90)#, roll=0)

        ax2 = fig.add_subplot(1, 3, 3, projection='3d')
        ax2.set_xticklabels('')
        ax2.plot_surface(X, Y, q_save, rstride=1, cstride=1, antialiased=True)
        ax2.set_ylabel(r'$k$', labelpad=15)
        ax2.set_zlabel(r'$q$')
        ax2.set_zlim([-2, 0])
        ax2.set_title(r'$k, q$')
        ax2.view_init(elev=0, azim=0)
        plt.show()


    if mplot:
        plot_mod = model(beta=bred[np.argmin(xred)],
                         kappa=kred[np.argmin(xred)],lam=lam)                 
        plot_mod.distance_modulus(effort=dm_effort)
        plot_mod.plot(which='acc', lcdm=model(), matter=model(lam=0.))
        plot_mod.plot(which='dm', lcdm=model(), matter=model(lam=0.))

    return top_mod


def auto_optimize(fname: str, it_num: int = 2,
                  search_method: str = 'acc', length: int = 20,
                  beta_lim_init: tuple = (1, 4),
                  kappa_lim_init: tuple = (1, 10), o_lambda: float = 0.,
                  dm_effort: bool = False, dm_method: str = 'int',
                  plot: int = 1, double_eval: bool = False,
                  require_decreasing_chi: bool = False,
                  fdir: str = '../../Data/model_data/', verbose=True):
    """
    Automatically optimize model parameters for a given number of iterations
    it_num >=3  using search_method method. First call search_method with
    initial inputs and then if it_num > 1, call search_method with inputs 80%
    nearer the best model from the it_num - 1 iteration.

    Parameters
    ----------
    fname : str
        Name of file to save data to. Use 'nosave' to not save data. Only saves
        data for final iteration.
    it_num : int
        Number of iterations to run.
    search_method : function
        Method to use for searching for best model. Options are chi_search or
        chi_search_a
    beta_lim_init : tuple
        Initial beta limits for search.
    kappa_lim_init : tuple
        Initial kappa limits for search.
    o_lambda : float
        Cosmological constant.
    dm_effort : bool
        Whether to use effort to calculate distance modulus.
    dm_method : str
        Method to use for calculating distance modulus.
    plot : int
        Whether to plot results. 0 = no plot, 1 = plot best, 2 = plot all.
    double_eval : bool
        Whether to evaluate on both chi^2 values for each model.
    require_decreasing_chi : bool
        Whether to require chi^2 to decrease with each iteration.
    fdir : str
        Directory to save data to.
    verbose : bool
        Whether to print progress.
    """

    # Validate arguments
    if fname == 'nosave':
        save = False
    elif len(fname) < 1:
        raise ValueError('fname must be a string')
    
    if it_num < 2:
        raise ValueError('it_num must be a positive integer greater than 1')
    
    if plot not in (0, 1, 2):
        raise ValueError('plot must be 0, 1 or 2')
    
    if search_method != 'acc' and search_method != 'dm':
        raise ValueError('search_method must be "acc" or "dm"')
    
    if search_method == 'acc':
        acc = True
    else:
        acc = False
    
    matter = model(lam=0.)
    lcdm = model()
    lcdm.norm(matter=matter)

    # Limits on initial search
    li = 0.6
    ui = 1 + (1 - li)
    # Limits on middle search  
    lm = 0.7
    um = 1 + (1 - lm)
    # Limits on final search
    lf = 0.8
    uf = 1 + (1 - lf)
    
    # Plot parameter bools
    plot_notfinal = True if plot == 2 else False
    plot_final = False if plot == 0 else True

    # Initial search
    if verbose:
        print('Initial search (1 / {})'.format(it_num))

    model_initial = chi_search_a(fname='nosave', length=length,
                                 blim=beta_lim_init, klim=kappa_lim_init,
                                 lam=o_lambda,
                                 plot=plot_notfinal) \
        if acc else chi_search(fname='nosave', length=length,
                               blim=beta_lim_init, klim=kappa_lim_init,
                               lam=o_lambda, dm_method=dm_method, round=2,
                               dm_effort=dm_effort, plot=plot_notfinal,
                               double_eval=double_eval)
    
    # Get correct chi^2 value for initial model
    if acc:
        model_init_int = np.interp(lcdm.a, model_initial.a,
                                           model_initial.a2norm)

    running_chi = rchi2(model_init_int, lcdm.a2norm) if acc \
                    else (model_initial.chi_int + model_initial.chi_tay) \
                    if double_eval else model_initial.chi_int \
                        if dm_method == 'int' else model_initial.chi_tay
    
    # Subsequent searches
    if verbose:
        print('Subsequent search (2 / {})'.format(it_num))
    model_mid = chi_search_a(fname='nosave', length=length, 
                                blim=(li*model_initial.b, ui*model_initial.b),
                                klim=(li*model_initial.k, ui*model_initial.k),
                                lam=o_lambda, plot=plot_notfinal) \
        if acc else chi_search(fname='nosave', length=length,
                               blim=(li*model_initial.b, ui*model_initial.b),
                               klim=(li*model_initial.k, ui*model_initial.k),
                               lam=o_lambda, dm_method=dm_method, round=2,
                               dm_effort=dm_effort, plot=plot_notfinal,
                               double_eval=double_eval)
    
    model_mid_int = np.interp(lcdm.a, model_mid.a, model_mid.a2norm) if acc \
                    else None
    model_mid_chi = rchi2(model_mid_int, lcdm.a2norm) if acc \
                    else (model_mid.chi_int + model_mid.chi_tay) \
                    if double_eval else model_mid.chi_int \
                        if dm_method == 'int' else model_mid.chi_tay
    
    if model_mid_chi > running_chi and require_decreasing_chi:
        raise Exception('model_mid has higher chi^2 than model_initial\n'
                        '\t{:.3f} > {:.3f}\n'
                        'Check initial search parameters'.format(model_mid_chi,
                                                                 running_chi))
    
    running_chi = model_mid_chi

    # More middle searches
    if it_num > 3:
        for i in range(it_num - 3):
            if verbose:
                print('Subsequent search ({} / {})'.format(2 + i, it_num))
            model_mid = chi_search_a(fname='nosave', length=length,
                                    blim=(lm*model_mid.b, um*model_mid.b),
                                    klim=(lm*model_mid.k, um*model_mid.k),
                                    lam=o_lambda, plot=plot_notfinal) \
                    if search_method == 'acc' \
                    else chi_search(fname='nosave', length=length,
                                    blim=(lm*model_mid.b, um*model_mid.b),
                                    klim=(lm*model_mid.k, um*model_mid.k),
                                    lam=o_lambda, dm_method=dm_method,
                                    dm_effort=dm_effort, round=2,
                                    plot=plot_notfinal,double_eval=double_eval)
            
            # compare chi^2 values
            model_mid_int = np.interp(lcdm.a,model_mid.a,model_mid.a2norm) if \
                            acc else None
            model_mid_chi = rchi2(model_mid_int, lcdm.a2norm) if acc \
                            else (model_mid.chi_int + model_mid.chi_tay) \
                            if double_eval else model_mid.chi_int \
                                if dm_method == 'int' else model_mid.chi_tay
            
            if model_mid_chi > running_chi and require_decreasing_chi:
                raise Exception('model_mid has higher chi^2 than previous '
                                'iteration\n\t{:.3f} > {:.3f}\n'
                                'Check initial search parameters'.format(
                                                model_mid_chi, running_chi))
            
            running_chi = model_mid_chi

    # Final search
    if verbose:
        print('Final search ({} / {})'.format(it_num, it_num))
    model_final = chi_search_a(fname=fname if save else 'nosave',length=length,
                                blim=(lf*model_mid.b, uf*model_mid.b),
                                klim=(lf*model_mid.k, uf*model_mid.k),
                                lam=o_lambda, plot=plot_final, fdir=fdir)\
            if acc else chi_search(fname=fname, length=length,
                                   blim=(lf*model_mid.b, uf*model_mid.b),
                                   klim=(lf*model_mid.k, uf*model_mid.k),
                                   lam=o_lambda, dm_method=dm_method,
                                   round = 3 if it_num > 3 else 2,
                                   dm_effort=dm_effort, plot=plot_final,
                                   double_eval=double_eval, fdir=fdir)
    
    if acc:
        model_fin_int = np.interp(lcdm.a, model_final.a, model_final.a2norm)

    model_fin_chi = rchi2(model_fin_int, lcdm.a2norm) if acc \
            else (model_final.chi_int + model_final.chi_tay) if double_eval \
            else model_final.chi_int if dm_method == 'int' \
            else model_final.chi_tay

    if model_fin_chi > running_chi and require_decreasing_chi:
        raise Exception('model_final has higher chi^2 than model_mid\n'
                        '\t{:.3f} > {:.3f}\nCheck initial search parameters'
                        ''.format(model_mid_chi, running_chi))

    running_chi = model_fin_chi
    
    print('After {} iterations, the best fit model has parameters '
          'beta = {:.4f} and kappa = {:.4f}\n'
          'with a chi^2 value of {:.5f}'.format(it_num, model_final.b,
                                                model_final.k, running_chi))
    
    return model_final