# # Run a bunch of definitive stuff. I'm keeping this in a separate file
from main import *


# @timer
def main():
    """
    Main function
    """

    if __name__ == '__main__':

        # lcdm = model(beta=2.770, kappa=2.704, lam=0.)
        # lcdm = model(beta=2.530, kappa=1.0, lam=0.)
        # lcdm = model(beta=2.833, kappa=0.567, lam=0.)
        # lcdm = model(beta=2.685, kappa=1.939, lam=0.)
        # lcdm = model(beta=2.238, kappa=15.687, lam=0.)
        # lcdm = model()
        lcdm.distance_modulus(effort=True)

        fig = plt.figure(constrained_layout=False)
        frame1 = fig.add_axes((.1, .3, .8, .6))
        plt.errorbar(z_sn, df['m_b_corr'], yerr=snerr, lw=0.5, ls='', marker='.',
                        markersize=2, label=r'Pantheon+SH0ES', zorder=0)

        mb_astro, Mb_astro = mbcorr_fit(dm_astro)
        plt.plot(z_sn, mb_astro, c='orange', ls='-',
                    label=r'$m_{{B, \Lambda\mathrm{{CDM}}}}, M_{{B}} = {:.1f}$, '\
                r'$\chi^{{2}}_{{r}}={:.4f}$'.format(
                        Mb_astro, rchi2(obs=mb_astro, exp=df['m_b_corr'])))
        
        mb_int, Mb_int = mbcorr_fit(lcdm.dm_int)
        plt.plot(z_sn, mb_int, c='k', ls='-.', label=r'$m_{{B}}(z, E(z)), M_{{B}}={:.1f}$,'
                r' $\chi^{{2}}_{{r}}={:.4f}$'.format(
                                Mb_int, rchi2(obs=mb_int, exp=df['m_b_corr'])))
        
        mb_tay, Mb_tay = mbcorr_fit(lcdm.dm_tay)
        plt.plot(z_sn, mb_tay, c='r', ls='--',
                    label=r'$m_{{B}}(z, q, j, s), M_{{B}} = {:.1f}$, '
                    r'$\chi^{{2}}_{{r}}={:.4f}$'.format(
            Mb_tay, rchi2(obs=mb_tay, exp=df['m_b_corr'])))
        
        plt.ylabel(r'$m_{B}$ [mag]')
        plt.xscale('log')
        plt.tick_params(axis='both', which='both', direction='in',
                        bottom=True, top=True, left=True, right=True)
        plt.legend()
        plt.grid()
        plt.ylim([7, 28])
        frame1.set_xticklabels([])

        frame2 = fig.add_axes((.1, .1, .8, .2))
        plt.errorbar(
            z_sn, df['m_b_corr'] - mb_astro, yerr=snerr, lw=0.5, ls='', marker='.',
            markersize=2, label=r'Supernova',zorder=0
        )
        plt.plot(
            z_sn, (mb_int - mb_astro), c='k', ls='-.'
        )
        plt.plot(
            z_sn, (mb_tay - mb_astro), c='r', ls='--'
        )
        plt.xlabel(r'$z$')
        plt.xscale('log')
        plt.ylabel(r'$\Delta{m_{B}}_{\Lambda\mathrm{CDM}}$')
        plt.tick_params(
            axis='both', which='both', direction='in',
            bottom=True, top=True, left=True, right=True
        )
        plt.tick_params(axis='x', pad=7)
        extraticks = [-1, 1]
        plt.yticks(list(plt.yticks()[0]) + extraticks)
        plt.ylim([-2, 2])
        # plt.ylim([-1, 1])
        plt.grid(which='major', axis='both')
        plt.show()

main()