# Main file

# Import
from _icd import *
from _functions import modified_friedmann, acceleration, inverse, adot_inverse, integrand, dm_z_o4, rchi2, read_model_data, specific_function


# Classes
class model():
    """
    class which creates a cosmological model with given parameters
    
    This class takes initial inputs and then integrates the modified Friedmann equations to get a, adot (a1) and t. From there we derive further
    time derivatives as well as other cosmological parameters.

    Calling this class with default arguments, model(), returns a flat LCDM model. A matter dominated model is returned by model(lam=0.)
    """

    def __init__(self, a_start=1e-3, mat=mat0, rad=rad0, lam=lam0, beta=3., kappa=0., n=1, xaxis='a', xmax=1.5, xlen=50000):
        """
        initialization method for model class
        :param a_start: initial value of scale factor at which to begin integration, chosen 1e-3 to correspond to time of recombination (z ~ 1100)
        :param mat: dimensionless matter density parameter (baryons + dark matter), defined in main_icd
        :param rad: dimensionless radiation density parameter (photons + neutrinos), defined in main_icd
        :param lam: dimensionless dark energy density parameter (cosmological constant), defined in main_icd
        :param beta: power to which mass is raised in new DM force
        :param kappa: constant in front of new DM force
        :param n: not sure, but has to do with power series in Press-Schechter formalism, used in variable gamma
        :param xaxis: axis upon which we plan to plot, can be time 't' or scale factor 'a', scale factor and its derivatives only run until
        xaxis variable reaches unity
        :param xamax: upper limit of integration in units of time tH0, chosen as 1.5 to ensure scale factor can reach unity if xaxis = 'a'
        :param xlen: length of time array for integration, chosen to be sufficiently large so that when we go to redshift there are values similar
        to those of the z_sn array (1e-3, 2.2)
        """
        # properties of this model
        self.m = mat
        self.r = rad
        self.l = lam
        self.b = beta
        self.k = kappa
        self.axis = xaxis

        # gamma and p variables used in acceleration equation
        self.gamma  = 1 + n / 3
        self.p = self.gamma / (2 * beta)
        
        # set initial condtions from a_start
        a1_start = np.sqrt(self.m * a_start**-1 + self.r * a_start**-2 + self.l * a_start**2)
        initial_conditions = (a_start, a1_start)

        # time array upon which to integrate
        time_array = np.linspace(0, xmax, xlen)
        # time_array = np.logspace(-5, np.log10(xmax), xlen)

        # integrate modified friedmann equations
        # solution = odeint(modified_friedmann, y0=initial_conditions, t=time_array, args=(self.m, self.r, self.l, self.k, self.p), tfirst=True)
        solution = solve_ivp(modified_friedmann, y0=initial_conditions, t_span=(0, xmax), t_eval=time_array, args=(self.m, self.r, self.l, self.k, self.p))

        # extract model parameters
        # self.a  = solution[:, 0]
        # self.a1 = solution[:, 1]
        # self.t  = np.copy(time_array)
        self.a  = solution.y[0, :]
        self.a1 = solution.y[1, :]
        self.t  = solution.t

        # time derivatives
        self.a2 = acceleration(self.a, self.a1, self.m, self.r, self.l, self.k, self.p)
        self.a3 = np.diff(self.a2) / np.diff(self.t)
        self.a4 = np.diff(self.a3) / np.diff(self.t[:-1])

        # find where today is in terms of a or t
        hoy = np.argmin(np.abs(1 - self.a)) if xaxis == 'a' else np.argmin(np.abs(1 - self.t))
        
        # truncate arrays to today
        self.a = self.a[:hoy]
        self.a1 = self.a1[:hoy]
        self.a2 = self.a2[:hoy]
        self.a3 = self.a3[:hoy]
        self.a4 = self.a4[:hoy]
        self.t = self.t[:hoy]

        # parameters from time derivatives
        self.q = - self.a2[-1] * self.a[-1] * self.a1[-1]**-2    # deceleration parameter q
        self.j = self.a3[-1] * self.a[-1]**2 * self.a1[-1]**-3    # jerk j
        self.s = self.a4[-1] * self.a[-1]**3 * self.a1[-1]**-4    # snap s

    
    def distance_modulus(self, effort=True):
        """
        distance_modulus() method which calculates the distance modulus of scale factor values calculated in the initialization in two ways

        :param effort: boolean, if True, uses the "true" method of calculating the distance modulus, if False, uses a faster method
        :return: distance modulus array
        
        Note on effort parameter:
        "True" way is calculating f and using the integrand function. "False" way method calculates E(z) as an array of data and simply
        integrates that. It is a bit mysterious however because it uses the inverse function which returns x**(-0.5) when it should
        just be x**(-1). For some reason though, -0.5 just works.
        
        The distance modulus is calculated through two methods, here called integration and taylor:
        DM through integration is calculated through a series of equations as follows:
        dc = dh * int_0^z dz/E(z) where E(z) = sqrt(m(1+z)(1+kf) + r(1+z)^4 + l)
        dl = (1+z) * dh ;   dm = 5 * np.log10(dl) + 25 (for dl given in Mpc)
        Because z from the model will not exactly equal z from the SN data, we interpolate along z from SNe
        
        DM through taylor is calculated using a Taylor approximation of the scale factor which is a long equation that's written below

        At the end we add the dm_int and dm_tay attributes corresponding to the DM calculated from integration and taylor respectively
        """

        z = (1 / self.a - 1)
        idx0 = np.argmin(np.abs(z - z_sn[-1]))
        idxf = np.argmin(np.abs(z - z_sn[0]))

        if effort:

            integral_f = np.zeros_like(z)
            for index, scale in enumerate(self.a):
                integral_f[index] = quad(adot_inverse, 0, scale, (self.a1[index],))[0]
            f = (self.a1 / self.a * integral_f)**(self.p)

            dl = np.zeros_like(z)
            for index, redshift in enumerate((z)):
                dl[index] = (1 + redshift) * dh * quad(integrand, 0, redshift, (self.m, self.r, self.l, f[index], self.k))[0]

        elif not effort:

            ez = self.a1 / self.a
            
            dl = np.zeros_like(z)
            for index, redshift in enumerate(z):
                dl[index] = (1 + redshift) * dh * quad(inverse, 0, redshift, (ez[index],))[0]
        
        else:
            raise Exception('Bad "effort" input for "distance_luminosity" method in "model" class. Input is boolean.')
    
        dm_int = 5 * np.log10(dl) + 25
        dm_int = np.flip(dm_int[idx0:idxf])
        z = np.flip(z[idx0:idxf])

        try:
            self.dm_int = np.interp(z_sn, z, dm_int)
        except:
            self.dm_int = np.nan

        self.dm_tay = dm_z_o4(z=z_sn, q=self.q, j=self.j, s=self.s)

    
    def chi2value(self, dm_method='int', chi_method='formula', eval_both=True):
        """
        Method to calculate the chi squared value of the model compared to the SN data set

        :param dm_method: string, either 'int' or 'tay' for distance modulus method
        :param chi_method: string, either 'formula' or 'poly' for chi squared method
        :param eval_both: boolean, if True, calculates chi2 for both distance modulus methods, if False, calculates for only one,
        defaulting to whichever is specified in dm_method
        :return: chi squared value of model compared to SN data set
        """

        if eval_both:
            self.chi_int = np.nan if np.isnan(self.dm_int).all() else rchi2(obs=self.dm_int, method=chi_method)
            self.chi_tay = np.nan if np.isnan(self.dm_tay).all() else rchi2(obs=self.dm_tay, method=chi_method)
        else:
            if dm_method == 'int':
                distmod = self.dm_int
            elif dm_method == 'tay':
                distmod = self.dm_tay
            else:
                raise Exception('Bad "dm_method" input for "chi2value" method in "model" class. Inputs are "int" or "tay".')

            self.chi = np.nan if np.isnan(distmod).all() else rchi2(obs=distmod, method=chi_method)


    def plot(self, which, lcdm, matter):
        """
        This method plots the acceleration of the scale factor or the distance modulus, depending on the input "which" 
        When plotting acceleration it normalizes the acceleration of the model to the acceleration of the matter model
        and plots the acceleration of the LCDM model for comparison

        :param which: string, options are 'acc' or 'dm' corresponding to acceleration or distance modulus
        :param lcdm: model object, used for comparison
        :param matter: model object, used for normalization of acc plot
        :return: plot
        """

        self.a2_norm = self.a2 / np.interp(self.a, matter.a, matter.a2)

        if self.axis == 'a':
            x_mod = self.a
            x_lcdm = lcdm.a
        elif self.axis == 't':
            x_mod = self.t
            x_lcdm = lcdm.t
        else:
            raise Exception('Something here')

        if which == 'acc':

            plt.figure()
            plt.plot(x_mod, self.a2/np.interp(self.a, matter.a, matter.a2), c='r', ls='-',
                     label=r'Alt. Model, $\beta={:.2f}, k={:.2f}$'
                     ''.format(self.b, self.k))
            
            plt.plot(x_lcdm, lcdm.a2/np.interp(lcdm.a, matter.a, matter.a2), c='k', ls='--',
                     label=r'$\Lambda$CDM Model')
            
            plt.xlabel(r'$a$')
            plt.ylabel(r'$\ddot{a}/\ddot{a}_{\mathrm{M}}$')
            plt.ylim([-5, 2])
            plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
            plt.legend(loc='lower left')
            plt.grid()
            plt.show()

        elif which == 'dm':
            
            fig = plt.figure(constrained_layout=False)
            frame1 = fig.add_axes((.1, .3, .8, .6))
            plt.errorbar(z_sn, sndat, yerr=snerr, lw=0.5, ls='', marker='.', markersize=2, label=r'Pantheon+SH0ES', zorder=0)

            plt.plot(z_sn, dm_astro, c='orange', ls='-',
                     label=r'$DM$ from flat-$\Lambda$CDM Cosmology, $\chi^{{2}}_{{r}}={:.4f}$'.format((rchi2(obs=dm_astro))))
            
            plt.plot(z_sn, self.dm_int, c='k', ls='-.',
                     label=r'$DM(z, E(z))$, $\chi^{{2}}_{{r}}={:.4f}$'.format((rchi2(obs=self.dm_int))))
            
            plt.plot(z_sn, self.dm_tay, c='r', ls='--',
                     label=r'$DM(z, q, j, s)$, $\chi^{{2}}_{{r}}={:.4f}$'.format((rchi2(obs=self.dm_tay))))
            
            plt.ylabel(r'$DM$ [mag]')
            plt.xscale('log')
            plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
            plt.legend()
            plt.grid()
            frame1.set_xticklabels([])

            frame2 = fig.add_axes((.1, .1, .8, .2))
            plt.errorbar(z_sn, sndat - dm_astro, yerr=snerr, lw=0.5, ls='', marker='.', markersize=2, label=r'Supernova', zorder=0)
            plt.plot(z_sn, (self.dm_int - dm_astro)/dm_astro*100, c='k', ls='-.')
            plt.plot(z_sn, (self.dm_tay - dm_astro)/dm_astro*100, c='r', ls='--')
            plt.xlabel(r'$z$')
            plt.xscale('log')
            plt.ylabel(r'$\%\Delta{DM}_{\mathrm{flat}\Lambda\mathrm{CDM}}$')
            plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
            plt.tick_params(axis='x', pad=7)
            extraticks = [-1, 1]
            plt.yticks(list(plt.yticks()[0]) + extraticks)
            plt.ylim([-2, 2])
            # plt.ylim([-1, 1])
            plt.grid(which='major', axis='both')
            plt.show()

        else:
            raise Exception('Bad "which" input in "plot" method in "model" class. Inputs are "acc" or "dm".')


# Functions
def chi2_comp(parameter, space, beta=3., kappa=0., lam=lam0, dm_effort=False, dm_method='int', chi_method='formula', plot=True):
    """
    This function calculates the chi^2 value for a given parameter space and plots the chi^2 value as a function of the
    parameter space. It also prints the lowest chi^2 value and the corresponding parameter value.
    :param parameter: string, options are 'l', 'b' or 'k', corresponding to lambda, beta or kappa
    :param space: array, parameter space to be explored
    :param beta: float, beta value
    :param kappa: float, kappa value
    :param lam: float, lambda value
    :param dm_effort: bool, if True, uses the "true" method of calculating the distance modulus, if False, uses a faster method
    :param dm_method: string, options are 'int' or 'tay', corresponding to the integration method or Taylor expansion method
    :param chi_method: string, options are 'formula' or 'poly' corresponding to the formula method or polynomial method
    :param plot: bool, if True, plots the chi^2 value as a function of the parameter space
    :return: model object, model with the lowest chi^2 value   
    """

    array = np.zeros_like(space)

    if parameter == 'l':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=value, beta=beta, kappa=kappa)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
            array[index] = tmod.chi
        
        model_optimized = model(lam=space[np.argmin(array)], beta=beta, kappa=kappa)

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f} and Omg_Lambda = {:.3f} for Omg_Lambda in the range [{:.1f}, {:.1f}]'
                ''.format(np.min(array), beta, kappa, space[np.argmin(array)], space[0], space[-1]))

    elif parameter == 'b':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=lam, beta=value, kappa=kappa)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
            array[index] = tmod.chi
        
        model_optimized = model(lam=lam, beta=space[np.argmin(array)], kappa=kappa)

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f} and Omg_Lambda = {:.3f} for beta in the range [{:.1f}, {:.1f}]'
                ''.format(np.min(array), space[np.argmin(array)], kappa, lam, space[0], space[-1]))
    
    elif parameter == 'k':
        for index, value in enumerate(tqdm(space)):
            tmod = model(lam=lam, beta=beta, kappa=value)
            tmod.distance_modulus(effort=dm_effort)
            tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
            array[index] = tmod.chi
        
        model_optimized = model(lam=lam, beta=beta, kappa=space[np.argmin(array)])

        print('The lowest chi^2 value is {:.3f} for beta = {:.3f}, k = {:.3f} and Omg_Lambda = {:.3f} for kappa in the range [{:.1f}, {:.1f}]'
                ''.format(np.min(array), beta, space[np.argmin(array)], lam, space[0], space[-1]))
            
    else:
        raise Exception('Bad "parameter" input in chi2_comp function')

    if plot:
        xlab = r'$\beta$' if parameter == 'b' else r'$k$' if parameter == 'k' else r'$\Omega_{\Lambda}$'
        yscale = 'log' if np.max(array)/np.min(array) > 50 else 'linear'
        plotlab = r'$k={:.3f},\Omega_{{\Lambda}}={:.1f}$' if parameter == 'b' else r'$\beta={:.3f},\Omega_{{\Lambda}}={:.1f}$' if parameter == 'k' else r'$\beta={:.3f},k={:.3f}$'
        plotform = (kappa, lam) if parameter == 'b' else (beta, lam) if parameter == 'k' else (beta, kappa)

        plt.figure()
        plt.plot(space, array, c='k', ls='-', label=plotlab.format(*plotform))
        plt.plot([space[np.argmin(array)], space[np.argmin(array)]], [0.9*np.min(array), 1.01*np.max(array)], c='r', ls='--')
        plt.xlabel(xlab)
        plt.ylabel(r'$\chi^{2}_{r}$')
        plt.yscale(yscale)
        plt.grid()
        plt.legend(loc='best', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
        plt.show()

    else:
        pass

    model_optimized.distance_modulus(effort=dm_effort)
    model_optimized.chi2value(dm_method=dm_method, chi_method=chi_method)

    return model_optimized


def chi_search(fname, length=10, blim=(2., 4.), klim=(1., 10.), l=0., dm_effort=False, dm_method='int', chi_method='formula', 
               plot=True, round=1, scale=LogNorm(), fdir='../../Data/model_data/', double_eval=True):
    """
    chi_search() function. Calculates chi^2 value for models with different combinations (beta, kappa). This function creates a linear range of
    beta values using length & blim, same for kappa, and then loops over combinations of those. For each combination a model is created and the
    distance_modulus method called using dm_effort input. Then chi2_value() method is calledusing dm_method and chi_method inputs. Chi^2 value
    is stored and shaped into (length, length) array for plotting. Finally stores data using inputted file name and then prints statement of results.

    :param fname: string, name of file to save data to
    :param length: integer, length of array with beta or kappa values. Length^2 is number of iterations in loop
    :param blim: tuple, upper and lower bounds of beta
    :param klim: tuple, upper and lower bounds of kappa
    :param l: integer, lambda value to be used for each model
    :param dm_effort: boolean, whether to use effort or not in calculation of distance modulus (used later in calculating chi^2)
    :param dm_method: 'int' or 'tay', which method to use in evaluating distance modulus
    :param chi_method: 'formula' or 'poly', which method to use in evaluating chi^2 value
    :param plot: boolean, plot chi^2 heat map or no
    :param round: integer, used when plot==True, what number of decimal places to round x & y labels to, for visual purposes
    :param scale: NoNorm() or LogNorm(), used when plot==True, scale of heat map. NoNorm() (linear) for fine scale grid when chi^2 values are within an order of magnitude,
                    LogNorm() (log) for coarse scale grid when chi^2 values are not within an order of magnitude
    :param fdir: string, file directory for storing data. Default stores in /Data/model_data/
    :param double_eval: boolean, whether to evaluate chi^2 value for each model twice (once for each method) or not
    :return: optimized model object
    """

    if fname == 'nosave':
        save = False
    else:
        fname = fname + '.txt'
        save = True

    brange = np.linspace(np.min(blim), np.max(blim), length)
    krange = np.linspace(np.min(klim), np.max(klim), length)

    # if double_eval:
    chival_int = np.zeros(length**2)
    chival_tay = np.zeros(length**2)

    # Iterate over all (b, k), store chi value for each
    for index, param in enumerate(itertools.product(brange, krange)):
        tmod = model(lam=l, beta=param[0], kappa=param[1])
        tmod.distance_modulus(effort=dm_effort)
        tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
        if np.max(tmod.a2/np.interp(tmod.a, matter.a, matter.a2)) > 3 or np.min(tmod.a2/np.interp(tmod.a, matter.a, matter.a2)) < -10:
            chival_int[index] = np.nan
            chival_tay[index] = np.nan
        else:
            chival_int[index] = tmod.chi_int
            chival_tay[index] = tmod.chi_tay
    
    # Convert NaNs to 1e6 because otherwise it's dumb (number not important)
    chival_int_nonan = np.nan_to_num(chival_int, nan=1e6)
    chival_tay_nonan = np.nan_to_num(chival_tay, nan=1e6)
    chival_com_nonan = chival_int_nonan + chival_tay_nonan
    lowest = np.argmin(chival_com_nonan) if double_eval else np.argmin(chival_int_nonan) if dm_method == 'int' else np.argmin(chival_tay_nonan)

    if double_eval:
        nan_ratio = len(chival_com_nonan[chival_com_nonan >= 1e6]) / len(chival_com_nonan) * 100
    else:
        nan_ratio = len(chival_int_nonan[chival_int_nonan == 1e6]) / len(chival_int_nonan) * 100 if dm_method == 'int' else len(chival_tay_nonan[chival_tay_nonan == 1e6]) / len(chival_tay_nonan) * 100

    # Raise exception if only NaN values were returned
    if np.all(chival_com_nonan >= 1e6) or np.all(chival_int_nonan == 1e6) or np.all(chival_tay_nonan == 1e6):
        raise Exception('No real chi values found')

    # Find beta and kappa values for lowest chi^2 value
    chi_low_int = chival_int_nonan[lowest]
    chi_low_tay = chival_tay_nonan[lowest]
    beta_low  = np.repeat(brange, length)[lowest]
    kappa_low = np.tile(krange, length)[lowest]

    # Print results of function call
    if double_eval:
        print('The lowest chi^2 values are {:.3f} (int) {:.3f} (tay) for beta = {:.3f} & mu = {:.3f} in the range\n{:.0f} < beta < {:.0f} and {:.0f} < mu < {:.0f}\n\t'
            '{:.1f} % of models had a chi^2 value of NaN. \n'
            ''.format(chi_low_int, chi_low_tay, beta_low, kappa_low, np.min(blim), np.max(blim), np.min(klim), np.max(klim), nan_ratio))
    else:
        print('The lowest chi^2 value is {:.3f} for beta = {:.3f} & mu = {:.3f} in the range\n{:.0f} < beta < {:.0f} and {:.0f} < mu < {:.0f}\n\t'
            '{:.1f} % of models had a chi^2 value of NaN. \n'
            ''.format(chi_low_int if dm_method == 'int' else chi_low_tay, 
                      beta_low, kappa_low, np.min(blim), np.max(blim), np.min(klim), np.max(klim), nan_ratio))

    # Plot heat map of chi values
    if plot and double_eval: 
        chi_plot_z_int = np.reshape(chival_int, (length, length)).T
        chi_plot_z_tay = np.reshape(chival_tay, (length, length)).T 
        chi_plot_z = np.reshape(chival_int + chival_tay, (length, length)).T

        fig, ax = plt.subplots(2, 2, figsize=(10, 8))
        fig.delaxes(ax[1,1])
        f1 = ax[0,0]; f2 = ax[0,1]; f3 = ax[1,0]
        cmap = matplotlib.cm.get_cmap('viridis_r').copy()
        cmap.set_bad(color='r')

        # intf
        im1 = f1.imshow(chi_plot_z_int, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f1.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='k', ha='center', va='center', fontsize=16)
        f1.set_xticks(np.linspace(0, length-0.5, 10),
                        np.round(np.linspace(brange[0], brange[-1], 10), 1),
                        rotation=45, fontsize=12)
        f1.set_yticks(np.linspace(0, length-0.5, 10),
                        np.round(np.linspace(krange[0], krange[-1], 10), 1),
                        fontsize=12)
        f1.set_xlabel(r'$\beta$')
        f1.set_ylabel(r'$k$')
        f1.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        f1.grid()
        fig.colorbar(im1, label=r'$\chi^{2}_{r}$', ax=f1)

        # tay
        im2 = f2.imshow(chi_plot_z_tay, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f2.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='w', ha='center', va='center', fontsize=16)
        f2.set_xticks(np.linspace(0, length-0.5, 10),
                        np.round(np.linspace(brange[0], brange[-1], 10), 1),
                        rotation=45, fontsize=12)
        f2.set_yticks(np.linspace(0, length-0.5, 10),
                        np.round(np.linspace(krange[0], krange[-1], 10), 1),
                        fontsize=12)
        # f2.text(x=np.array(range(length))[np.where(brange == beta_low)],
        #         y=np.array(range(length))[np.where(krange == kappa_low)],
        #         s=r'$\ast$', color='k', ha='center', va='center', fontsize=14)
        # f2.set_xticks(range(11), np.sort(np.append(np.linspace(blim[0],
        #                                                        blim[-1], 10),
        #                                                        beta_low)),
        #               rotation=45, fontsize=14)
        # f2.set_yticks(range(11), np.sort(np.append(np.linspace(klim[0],
        #                                                         klim[-1], 10),
        #                                                         kappa_low)),
        #                  fontsize=14)
        f2.set_xlabel(r'$\beta$')
        f2.set_ylabel(r'$k$')
        f2.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        f2.grid()
        fig.colorbar(im2, label=r'$\chi^{2}_{r}$', ax=f2)

        # both
        im3 = f3.imshow(chi_plot_z, cmap=cmap, origin='lower',
                        interpolation='nearest', norm=scale)
        f3.text(x=np.array(range(length))[np.where(brange == beta_low)],
                y=np.array(range(length))[np.where(krange == kappa_low)],
                s=r'$\ast$', color='w', ha='center', va='center', fontsize=16)
        f3.set_xticks(np.linspace(0, length-0.5, 10),
                        np.round(np.linspace(brange[0], brange[-1], 10), 1),
                        rotation=45, fontsize=12)
        f3.set_yticks(np.linspace(0, length-0.5, 10),
                        np.round(np.linspace(krange[0], krange[-1], 10), 1),
                        fontsize=12)


        # f3.text(x=np.array(range(length))[np.where(brange == beta_low)],
        #         y=np.array(range(length))[np.where(krange == kappa_low)],
        #         s=r'$\ast$', color='k', ha='center', va='center', fontsize=14)
        # f3.set_xticks(range(11), np.sort(np.append(np.linspace(blim[0],
        #                                                        blim[-1], 10),
        #                                                        beta_low)),
        #               rotation=45, fontsize=14)
        # f3.set_yticks(range(11), np.sort(np.append(np.linspace(klim[0],
        #                                                         klim[-1], 10),
        #                                                         kappa_low)),
        #                  fontsize=14)
        f3.set_xlabel(r'$\beta$')
        f3.set_ylabel(r'$k$')
        f3.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        f3.grid()      
        fig.colorbar(im3, label=r'$\chi^{2}_{r}$', ax=f3)
        plt.show()
    else:
        pass

    # Save chi, beta, kappa values to file
    if save:
        f_chi = np.copy(chival_int + chival_tay)
        f_beta = np.repeat(brange, length)
        f_kappa = np.tile(krange, length)
        f_save = np.vstack((f_chi, f_beta, f_kappa)).T
        f_comment = '#Results of "chi_search" called with the following inputs:\n' +\
                    '#length={}, blim=({}, {}), klim=({}, {}), lambda={}, effort={}, dm_method={}, chi_method={}\n'.format(
                        length, np.min(blim), np.max(blim), np.min(klim), np.max(klim), l, dm_effort, dm_method, chi_method) +\
                    '#Lowest chi^2 was with beta = {} & k = {}\n'.format(beta_low, kappa_low)
        np.savetxt(fname=fdir+fname, X=f_save, header='chi  beta    kappa', delimiter='   ', comments=f_comment)

    # Compute optimal model based on chi results
    model_optimized = model(lam=l, beta=beta_low, kappa=kappa_low)
    model_optimized.distance_modulus(effort=dm_effort)
    model_optimized.chi2value(dm_method=dm_method, chi_method=chi_method)

    return model_optimized
    
    # else:

    #     fname = fname + '.txt'

    #     brange = np.linspace(np.min(blim), np.max(blim), length)
    #     krange = np.linspace(np.min(klim), np.max(klim), length)
    #     chival = np.zeros(length**2)

    #     # Iterate over all (b, k), store chi value for each
    #     for index, param in enumerate(itertools.product(brange, krange)):
    #         beta, kappa = param
    #         tmod = model(lam=l, beta=beta, kappa=kappa)
    #         tmod.distance_modulus(effort=dm_effort)
    #         tmod.chi2value(dm_method=dm_method, chi_method=chi_method)
    #         if np.max(tmod.a2/np.interp(tmod.a, matter.a, matter.a2)) > 3 or np.min(tmod.a2/np.interp(tmod.a, matter.a, matter.a2)) < -10:
    #             chival[index] = np.nan
    #         else:            
    #             chival[index] = tmod.chi
        
    #     # Raise exception if only NaN values were returned
    #     if np.isnan(chival).all():
    #         raise Exception('No real chi values found')
    #     else:
    #         pass
        
    #     # Save chi, beta, kappa values to file
    #     f_chi = np.copy(chival)
    #     f_beta = np.repeat(brange, length)
    #     f_kappa = np.tile(krange, length)
    #     f_save = np.vstack((f_chi, f_beta, f_kappa)).T
    #     fcomment = '#Results of "chi_search" called with the following inputs:\n' +\
    #                 '#length={}, blim=({}, {}), klim=({}, {}), lambda={}, effort={}, dm_method={}, chi_method={}\n'.format(
    #                     length, np.min(blim), np.max(blim), np.min(klim), np.max(klim), l, dm_effort, dm_method, chi_method) +\
    #                 '#Lowest chi^2 was with beta = {} & k = {}\n'.format(beta_low, kappa_low) +\
    #                 '#NaN ratio = {:.1f} %'.format(nan_ratio)
    #     np.savetxt(fname=fdir+fname, X=f_save, header='chi  beta    kappa', delimiter='   ', comments=fcomment)
    #     print('Data saved in {}'.format(fdir+fname))

    #     # Convert NaNs to 1e6 because otherwise it's dumb (number not important)
    #     chival_nonan = np.nan_to_num(chival, nan=1e6)
    #     lowest = np.argmin(chival_nonan)
    #     nan_ratio = len(chival_nonan[chival_nonan == 1e6]) / len(chival_nonan) * 100

    #     # Find lowest chi^2 value and corresponding beta, kappa
    #     chi_low   = chival_nonan[lowest]
    #     beta_low  = np.repeat(brange, length)[lowest]
    #     kappa_low = np.tile(krange, length)[lowest]
        
    #     # Print results of function call
    #     print('The lowest chi^2 value is {:.3f} for beta = {:.3f} & mu = {:.3f} in the range {:.0f} < beta < {:.0f} and {:.0f} < mu < {:.0f}\n\t'
    #         '{:.1f} % of models had a chi^2 value of NaN. \n'
    #         ''.format(chi_low, beta_low, kappa_low, np.min(blim), np.max(blim), np.min(klim), np.max(klim), nan_ratio))


    #     # Plot heat map of chi values
    #     if plot:
    #         chi_plot_z = np.reshape(chival, (length, length)).T

    #         fig, ax = plt.subplots(figsize=(7, 5))
    #         cmap = matplotlib.cm.get_cmap('viridis_r').copy()
    #         cmap.set_bad(color='r')
    #         im = ax.imshow(chi_plot_z, cmap=cmap, origin='lower', interpolation='nearest', norm=scale)
    #         ax.text(x=np.array(range(length))[np.where(brange == beta_low)],
    #                 y=np.array(range(length))[np.where(krange == kappa_low)],
    #                 s=r'$\ast$', color='k', ha='center', va='center', fontsize=20)
            
    #         ax.set_xticks(range(length), np.round(brange, round), rotation=45, fontsize=14)
    #         ax.set_yticks(range(length), np.round(krange, round), fontsize=14)
    #         ax.set_xlabel(r'$\beta$')
    #         ax.set_ylabel(r'$k$')
    #         fig.colorbar(im, label=r'$\chi^{2}_{r}$')
    #         ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
    #         ax.grid()
    #         plt.show()
    #     else:
    #         pass


    #     # Compute optimal model based on chi results
    #     model_optimized = model(lam=l, beta=beta_low, kappa=kappa_low)
    #     model_optimized.distance_modulus(effort=dm_effort)
    #     model_optimized.chi2value(dm_method=dm_method, chi_method=chi_method)

    #     return model_optimized


def q_surface(length=20, blim=(2, 4), klim=(1, 10), qlim=(-1.0, 0.0), lam=0., dm_method='int', dm_effort=False, chi_method='formula',
              splot=True, mplot=True):
    """
    q_surface() function. Plots a surface of q values for a given range of beta and kappa values.

    :param length: number of points to plot in each direction
    :param blim: range of beta values to plot
    :param klim: range of kappa values to plot
    :param qlim: range of q values to plot
    :param lam: lambda value to use for model
    :param dm_method: method to use for distance modulus calculation
    :param dm_effort: if True, use more accurate distance modulus calculation
    :param chi_method: method to use for chi^2 calculation
    :param splot: if True, plot surface plot of q values
    :param mplot: if True, plot acceleration and distance modulus of optimized model
    :return: optimized model top_mod
    """

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

        if (temp_mod.q > np.min(qlim)) and (temp_mod.q < np.max(qlim)):
            temp_mod.distance_modulus(effort=dm_effort)
            temp_mod.chi2value(dm_method=dm_method, chi_method=chi_method, eval_both=False)

            qcd.append(temp_mod.q)
            bcd.append(temp_mod.b)
            kcd.append(temp_mod.k)
            xcd.append(temp_mod.chi)

            if temp_mod.chi < current_chi:
                current_chi = temp_mod.chi
                top_mod = temp_mod
            else:
                pass

        else:
            pass

    q_save = np.reshape(q_save, (length, length))

    if len(qcd) < 1:
        raise Exception('No q values found within the range {} < q < {}'.format(np.min(qlim), np.max(qlim)))
    else:
        pass

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
        else:
            continue

    sortorder = np.copy(qred)
    qred,   bred = zip(*sorted(zip(sortorder, bred)))
    _,      kred = zip(*sorted(zip(sortorder, kred)))
    _,      xred = zip(*sorted(zip(sortorder, xred)))

    for i in range(len(qred)):
        if xred[i] == np.min(xred):
            print('   **** q = {:.3f}\tfor\tbeta = {:.3f}\tk = {:.3f}\twith chi^2 = {:.4f} *****'
                  ''.format(qred[i], bred[i], kred[i], xred[i]))
        else:
            print('\tq = {:.3f}\tfor\tbeta = {:.3f}\tk = {:.3f}\twith chi^2 = {:.4f}'
                  ''.format(qred[i], bred[i], kred[i], xred[i]))
            
    
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
        plot_mod = model(lam=lam, beta=bred[np.argmin(xred)], kappa=kred[np.argmin(xred)])
        plot_mod.distance_modulus(effort=dm_effort)
        plot_mod.plot('acc', model(), model(lam=0.))
        plot_mod.plot('dm', model(), model(lam=0.))

    return top_mod



# Call
lcdm   = model()
matter = model(lam=0.)

# space = np.linspace(1, 10, 100)
# chi2_comp(parameter='k', space=space, beta=1.0, lam=0.)

# model_opt = chi_search(length=20, blim=[0.8, 1.5], klim=[4.3, 5.3], l=0.)
# model_opt.distance_modulus()
# model_opt.plot('acc', lcdm, matter)
# model_opt.plot('dm',  lcdm, matter)



# k_test = 0.6
# b_test = 3

# plt.figure(figsize=fs)
# # for i in [1, 1.5, 2, 2.5, 3]:
# for i in [0.1, 0.3, 0.5, 0.7, 1, 2, 5, 10]:
#     tmod = model(lam=0., beta=b_test, kappa=i)
#     accmat = np.interp(tmod.a, matter.a, matter.a2)
#     plt.plot(tmod.a, tmod.a2/accmat, label=r'$k={}$'.format(i))

# plt.plot(lcdm.a, lcdm.a2/np.interp(lcdm.a, matter.a, matter.a2), label=r'$\Lambda$CDM')
# plt.text(x=0.67, y=1.3, s=r'$k={}, \Omega_{{\Lambda}}=0$'.format(b_test))
# plt.grid()
# plt.legend()
# plt.ylim([-5, 2])
# plt.xlim([-0.1, 1.1])
# plt.xlabel(r'$a$')
# plt.ylabel(r'$\ddot{a}/\ddot{a}_{\mathrm{M}}$')
# plt.show()