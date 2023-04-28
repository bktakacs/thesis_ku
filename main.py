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

    def __init__(self, a_start=1e-3, mat=mat0, rad=rad0, lam=lam0, beta=3.,
                 kappa=0., n=1, xaxis='a', xmax=1.5, xlen=50000):
        """
        initialization method for model class
        :param a_start: initial value of scale factor at which to begin
        integration, chosen 1e-3 to correspond to time of recombination
        (z ~ 1100)
        :param mat: dimensionless matter density parameter
        (baryons + dark matter), defined in _icd
        :param rad: dimensionless radiation density parameter
        (photons + neutrinos), defined in _icd
        :param lam: dimensionless dark energy density parameter
        (cosmological constant), defined in _icd
        :param beta: power to which mass is raised in new DM force
        :param kappa: constant in front of new DM force
        :param n: not sure, but has to do with power series in Press-Schechter
        formalism, used in variable gamma
        :param xaxis: axis upon which we plan to plot, can be time 't' or scale
        factor 'a', scale factor and its derivatives only run until xaxis
        variable reaches unity
        :param xamax: upper limit of integration in units of time tH0, chosen
        as 1.5 to ensure scale factor can reach unity if xaxis = 'a'
        :param xlen: length of time array for integration, chosen to be
        sufficiently large so that when we go to redshift there are values
        similar to those of the z_sn array (1e-3, 2.2)
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

    def norm(self, matter):
        """
        norm() method which normalizes the scale factor acceleration
        to that of a matter only model
        :param matter: matter class object
        :return: None
        """
        self.a2norm = self.a2 / np.interp(self.a, matter.a, matter.a2)

    
    def distance_modulus(self, effort=True):
        """
        distance_modulus() method which calculates the distance modulus of
        scale factor values calculated in the initialization in two ways

        :param effort: boolean, if True, uses the "true" method of calculating
        the distance modulus, if False, uses a faster method
        :return: distance modulus array
        
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
        
        else:
            raise Exception('Bad "effort" input for "distance_luminosity" \
                            method in "model" class. Input is boolean.')
    
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

    
    def chi_value(self, dm_method='int', chi_method='formula', eval_both=True):
        """
        Method to calculate the chi squared value of the model compared to the
        SN data set

        :param dm_method: string, either 'int' or 'tay' for distance modulus
        method
        :param chi_method: string, either 'formula' or 'poly' for chi squared
        method
        :param eval_both: boolean, if True, calculates chi2 for both distance
        modulus methods, if False, calculates for only one, defaulting to
        whichever is specified in dm_method
        :return: chi squared value of model compared to SN data set
        """

        if eval_both:
            self.chi_int = np.nan if np.isnan(self.dm_int).all() else \
                                    rchi2(obs=self.dm_int, method=chi_method)
            self.chi_tay = np.nan if np.isnan(self.dm_tay).all() else \
                                    rchi2(obs=self.dm_tay, method=chi_method)
        else:
            distmod = self.dm_int if dm_method == 'int' else self.dm_tay

            self.chi = np.nan if np.isnan(distmod).all() else \
                                        rchi2(obs=distmod, method=chi_method)


    def plot(self, which, lcdm, matter):
        """
        This method plots the acceleration of the scale factor or the distance
        modulus, depending on the input "which". When plotting acceleration it
        normalizes the acceleration of the model to the acceleration of the
        matter model and plots the acceleration of the LCDM model for
        comparison

        :param which: string, options are 'acc' or 'dm' corresponding to
        acceleration or distance modulus
        :param lcdm: model object, used for comparison
        :param matter: model object, used for normalization of acc plot
        :return: plot
        """

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
                     label=r'$DM$ from flat-$\Lambda$CDM Cosmology,\
                     $\chi^{{2}}_{{r}}={:.4f}$'.format(rchi2(obs=dm_astro)))
            
            plt.plot(z_sn, self.dm_int, c='k', ls='-.', label=r'$DM(z, E(z))$,\
                     $\chi^{{2}}_{{r}}={:.4f}$'.format(rchi2(obs=self.dm_int)))
            
            plt.plot(z_sn, self.dm_tay, c='r', ls='--',
                     label=r'$DM(z, q, j, s)$,\
                        $\chi^{{2}}_{{r}}={:.4f}$'.format(
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

        else:
            raise Exception('Bad "which" input in "plot" method in "model"\
                            class. Inputs are "acc" or "dm".')


# Functions
def chi_comp(parameter, space, beta=3., kappa=0., lam=lam0, dm_effort=False,
              dm_method='int', chi_method='formula', plot=True):
    """
    This function calculates the chi^2 value for a given parameter space and
    plots the chi^2 value as a function of the parameter space. It also prints
    the lowest chi^2 value and the corresponding parameter value.
    :param parameter: string, options are 'l', 'b' or 'k', corresponding to
    lambda, beta or kappa
    :param space: array, parameter space to be explored
    :param beta: float, beta value
    :param kappa: float, kappa value
    :param lam: float, lambda value
    :param dm_effort: bool, if True, uses the "true" method of calculating the
    distance modulus, if False, uses a faster method
    :param dm_method: string, options are 'int' or 'tay', corresponding to the
    integration method or Taylor expansion method
    :param chi_method: string, options are 'formula' or 'poly' corresponding to
    the formula method or polynomial method
    :param plot: bool, if True, plots the chi^2 value as a function of the
    parameter space
    :return: model object, model with the lowest chi^2 value   
    """

    array = np.zeros_like(space)

    if parameter == 'l':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=value, beta=beta, kappa=kappa)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi_value(dm_method=dm_method, chi_method=chi_method)
            array[index] = tmod.chi
        
        model_optimized = model(lam=space[np.argmin(array)], beta=beta,
                                kappa=kappa)

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f}'
              'and Omg_Lambda = {:.3f} for Omg_Lambda in the range'
              '[{:.1f}, {:.1f}]'.format(
                        np.min(array), beta, kappa, space[np.argmin(array)],
                        space[0], space[-1]))

    elif parameter == 'b':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=lam, beta=value, kappa=kappa)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
            array[index] = tmod.chi
        
        model_optimized = model(lam=lam, beta=space[np.argmin(array)],
                                kappa=kappa)

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f}'
              'and Omg_Lambda = {:.3f} for beta in the range [{:.1f}, {:.1f}]'
                ''.format(np.min(array), space[np.argmin(array)], kappa, lam,
                          space[0], space[-1]))
    
    elif parameter == 'k':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=lam, beta=beta, kappa=value)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
            array[index] = tmod.chi
        
        model_optimized = model(lam=lam, beta=beta,
                                kappa=space[np.argmin(array)])

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f}'
              'and Omg_Lambda = {:.3f} for kappa in the range [{:.1f}, {:.1f}]'
              ''.format(np.min(array), beta, space[np.argmin(array)], lam,
                        space[0], space[-1]))
            
    else:
        raise Exception('Bad "parameter" input in chi2_comp function')

    if plot:
        xlab = r'$\beta$' if parameter == 'b'\
                else r'$k$' if parameter == 'k'\
                else r'$\Omega_{\Lambda}$'
        yscale = 'log' if np.max(array)/np.min(array) > 50 else 'linear'
        plotlab = r'$k={:.3f},\Omega_{{\Lambda}}={:.1f}$' if parameter == 'b'\
                    else r'$\beta={:.3f},\Omega_{{\Lambda}}={:.1f}$'\
                        if parameter == 'k'\
                            else r'$\beta={:.3f},k={:.3f}$'
        plotform = (kappa, lam) if parameter == 'b'\
                    else (beta, lam) if parameter == 'k'\
                        else (beta, kappa)

        plt.figure()
        plt.plot(space, array, c='k', ls='-', label=plotlab.format(*plotform))
        plt.plot([space[np.argmin(array)], space[np.argmin(array)]],
                 [0.9*np.min(array), 1.01*np.max(array)], c='r', ls='--')
        plt.xlabel(xlab)
        plt.ylabel(r'$\chi^{2}_{r}$')
        plt.yscale(yscale)
        plt.grid()
        plt.legend(loc='best', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        plt.show()

    model_optimized.distance_modulus(effort=dm_effort)
    model_optimized.chi2value(dm_method=dm_method, chi_method=chi_method)

    return model_optimized


def chi_search(fname, length=10, blim=(2., 4.), klim=(1., 10.), l=0.,
               dm_effort=False, dm_method='int', chi_method='formula',
               plot=True, round=1, scale=LogNorm(), double_eval=False,
               fdir='../../Data/model_data/'):
    """
    chi_search() function. Calculates chi^2 value for models with different
    combinations (beta, kappa). This function creates a linear range of beta
    values using length & blim, same for kappa, and then loops over
    combinations of those. For each combination a model is created and the
    distance_modulus method called using dm_effort input. Then chi2_value()
    method is called using dm_method and chi_method inputs. Chi^2 value is
    stored and shaped into (length, length) array for plotting. Finally stores
    data using inputted file name and then prints statement of results.

    :param fname: string, name of file to save data to
    :param length: integer, length of array with beta or kappa values.Length^2
    is number of iterations in loop
    :param blim: tuple, upper and lower bounds of beta
    :param klim: tuple, upper and lower bounds of kappa
    :param l: integer, lambda value to be used for each model
    :param dm_effort: boolean, whether to use effort or not in calculation of
    distance modulus (used later in calculating chi^2)
    :param dm_method: 'int' or 'tay', which method to use in evaluating
    distance modulus
    :param chi_method: 'formula' or 'poly', which method to use in evaluating
    chi^2 value
    :param plot: boolean, plot chi^2 heat map or no
    :param round: integer, used when plot==True, what number of decimal places
    to round x & y labels to, for visual purposes
    :param scale: NoNorm() or LogNorm(), used when plot==True, scale of heat
    map. NoNorm() (linear) for fine scale grid when chi^2 values are within an
    order of magnitude, LogNorm() (log) for coarse scale grid when chi^2
    values are not within an order of magnitude
    :param fdir: string, file directory for storing data. Default stores in
    /Data/model_data/
    :param double_eval: boolean, whether to evaluate chi^2 value for each
    model twice (once for each method) or not
    :return: optimized model object
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

    # Iterate over all (b, k), store chi value for each
    for index, param in enumerate(itertools.product(brange, krange)):
        tmod = model(lam=l, beta=param[0], kappa=param[1])
        tmod.norm(matter=matter)
        if (np.max(tmod.a2norm) > 3 or np.min(tmod.a2norm) < -10):
            chival_int[index] = np.nan
            chival_tay[index] = np.nan
            nan_count += 1
        else:
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
            chival_int[index] = tmod.chi_int
            chival_tay[index] = tmod.chi_tay
    
    # Raise exception if only NaN values were returned
    if nan_count == length**2:
        raise Exception('No real chi values found')
    nan_ratio = nan_count / length**2 * 100

    chival = np.copy(chival_int) if dm_method == 'int' else \
             np.copy(chival_tay) if dm_method == 'tay' else None
    
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
    print('The lowest chi^2 values are {:.3f} (int) {:.3f} (tay) for beta'
          '= {:.3f} & mu = {:.3f} in the range\n\t'
          '{:.0f} < beta < {:.0f} and {:.0f} < mu < {:.0f}\n\t'
          '{:.1f} % of models had a chi^2 value of NaN. \n'
          ''.format(chi_low_int, chi_low_tay, beta_low, kappa_low,
                    np.min(blim), np.max(blim), np.min(klim), np.max(klim),
                    nan_ratio) if double_eval else \
          'The lowest chi^2 value is {:.3f} for beta = {:.3f} & mu = '
          '{:.3f} in the range'
          '\n{:.0f} < beta < {:.0f} and {:.0f} < mu < {:.0f}\n\t'
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

        # intf
        im1 = f1.imshow(chi_plot_z_int, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f1.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='r', ha='center', va='center', fontsize=16)
        f1.set_xticks(np.linspace(0, length-0.5, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), round),
                      rotation=45, fontsize=12)
        f1.set_yticks(np.linspace(0, length-0.5, 10),
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
        f2.set_xticks(np.linspace(0, length-0.5, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), round),
                      rotation=45, fontsize=12)
        f2.set_yticks(np.linspace(0, length-0.5, 10),
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
        f3.set_xticks(np.linspace(0, length-0.5, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), round),
                      rotation=45, fontsize=12)
        f3.set_yticks(np.linspace(0, length-0.5, 10),
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
        ax.set_xticks(np.linspace(0, length, 10),
                      np.round(np.linspace(brange[0], brange[-1], 10), 1),
                      rotation=45, fontsize=14)
        ax.set_yticks(np.linspace(0, length, 10),
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
                '#length={}, blim=({}, {}), klim=({}, {}), lambda={},'.format(
                    length, blim[0], blim[1], klim[0], klim[1], l) +\
                'effort={}, dm_method={}, chi_method={}\n'.format(
                                        dm_effort, dm_method, chi_method) +\
                '#Lowest chi^2 was with beta = {} & k = {}\n'.format(
                                                        beta_low, kappa_low)
        np.savetxt(fname=fdir+fname, X=f_save, header='chi beta kappa',
                   delimiter=' ', comments=f_comment)

    # Compute optimal model based on chi results
    model_optimized = model(lam=l, beta=beta_low, kappa=kappa_low)
    model_optimized.distance_modulus(effort=dm_effort)
    model_optimized.chi2value(dm_method=dm_method, chi_method=chi_method)

    return model_optimized

def chi_search_a(fname, length=10, blim=(2., 4.), klim=(1., 10.),
                 lam=0., plot=True, fdir='../../Data/model_data/'):
    """
    chi_search_a() function. Iterates over all combinations (beta, kappa) and
    for each creates a model and calculates the chi^2 fit of the acceleration
    to that of the LCDM model. Just a test and I'm curious to see what comes
    out of it.

    :param fname: file name to save results to
    :param length: length of beta and kappa ranges
    :param blim: limit of beta values to plot
    :param klim: limit of kappa values to plot
    :param plot: if True, plot best fit acceleration
    :param fdir: directory to save results to
    """

    if len(fname) == 0:
        raise Exception('fname must be a string of at least one character')
    
    if fname == 'nosave':
        nosave = True
    
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
    top_mod = model(lam=lam, beta=brange[0], kappa=krange[0])
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
    nan_ratio = nan_count / length**2

    # Find lowest chi^2 value and corresponding beta and kappa values
    chi_low = np.sort(chival)[0]
    lowest = np.argmin(np.abs(chival - chi_low))
    beta_low = np.repeat(brange, length)[lowest]
    kappa_low = np.tile(krange, length)[lowest]

    # Print results
    print('The lowest chi^2 value is {:.3f} for beta = {:.3f} & mu = {:.3f} in'
          'the range\n{:.0f} < beta < {:.0f} and {:.0f} < mu < {:.0f}\n\t'
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
    if not nosave:
        f_chi = chival
        f_beta = np.repeat(brange, length)
        f_kappa = np.tile(krange, length)
        f_save = np.vstack((f_chi, f_beta, f_kappa)).T
        f_comment = '#Results of "chi_search_a" called with the following' +\
                'inputs:\n' +\
                '#length={}, blim=({}, {}), klim=({}, {}), lambda={},'.format(
                    length, blim[0], blim[1], klim[0], klim[1], lam) +\
                '#Lowest chi^2 was with beta = {} & k = {}\n'.format(
                    beta_low, kappa_low)
        np.savetxt(fname=fdir+fname, X=f_save, header='chi beta kappa',
                   delimiter=' ', comments=f_comment)
    
    return top_mod


def q_surface(length=20, blim=(2., 4.), klim=(1., 10.), qlim=(-1.0, 0.0),
              lam=0., dm_method='int', dm_effort=False, chi_method='formula', 
              splot=True, mplot=True):
    """
    q_surface() function. Plots a surface of q values for a given range of
    beta and kappa values.

    :param length: int, number of points to plot in each direction
    :param blim: tuple, range of beta values to plot
    :param klim: tuple, range of kappa values to plot
    :param qlim: tuple, range of q values to plot
    :param lam: float, lambda value to use for model
    :param dm_method: string, method to use for distance modulus calculation
    :param dm_effort: bool, if True, use more accurate dist mod calculation
    :param chi_method: string, method to use for chi^2 calculation
    :param splot: bool, if True, plot surface plot of q values
    :param mplot: bool, if True, plot acceleration and dist mod of opt model
    :return: class, optimized model top_mod
    """

    if (blim[0] > blim[1] or klim[0] > klim[1] or qlim[0] > qlim[1]):
        raise Exception('blim, klim, and qlim tuples must be of form (a, b)\
                        where a < b')

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
            temp_mod.chi_value(dm_method=dm_method, chi_method=chi_method,
                               eval_both=False)

            qcd.append(temp_mod.q)
            bcd.append(temp_mod.b)
            kcd.append(temp_mod.k)
            xcd.append(temp_mod.chi)

            top_mod = temp_mod if temp_mod.chi < current_chi else top_mod
            current_chi = temp_mod.chi if temp_mod.chi < current_chi \
                                        else current_chi

    q_save = np.reshape(q_save, (length, length))

    if len(qcd) < 1:
        raise Exception('No q values found within the range {} < q < {}'
                        ''.format(np.min(qlim), np.max(qlim)))

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
        plot_mod.plot('acc', model(), model(lam=0.))
        plot_mod.plot('dm', model(), model(lam=0.))

    return top_mod


# def optimize(fname, length=20, beta_lim_init=(1, 4), kappa_lim_init=(1, 10),
#              o_lambda=0., dm_effort=False, dm_method='int', plot=2,
#              double_eval=False, fdir='../../Data/model_data/'):
#     """
#     Fill in later
#     """

#     m1 = chi_search('nosave', length=length, blim=beta_lim_init,
#                     klim=kappa_lim_init, lam=o_lambda, dm_effort=dm_effort,
#                     dm_method=dm_method, double_eval=double_eval, splot=False,
#                     mplot=False)

    


# Set up lcdm and matter models
lcdm   = model()
matter = model(lam=0.)