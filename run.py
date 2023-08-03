# # Run a bunch of definitive stuff. I'm keeping this in a separate file
from main import *


# @timer
def main():
    """
    Main function
    """

    if __name__ == '__main__':

        m = chi_search(
            'new-friedmann-eq30-bdf-extended.txt',
            length=50, blim=(1, 4), klim=(0.01, 1000), 
        )
        m.plot('acc')
        m.distance_modulus()
        m.plot('dm')

        # b = np.linspace(1, 4, 50)
        # # k = np.linspace(0.1, 1000, 51)
        # # k = np.logspace(-1, 3, 50)
        # k = np.geomspace(0.1, 1000, 50)

        # data = read_model_data('new-friedmann-eq30-bdf-1.txt')
        # b, k = specific_function(data, 0)
        # m = model(beta=b, kappa=k, lam=0.)
        # m.distance_modulus()
        # # print(b, k)
        # # m.plot('acc')
        # # m.distance_modulus()
        # # m.plot('dm')


        # m1 = model(beta=1.979591836734693855, kappa=15.99858719606057278, lam=0.)
        # m2 = model(beta=1.612244897959183687, kappa=86.85113737513520960, lam=0.)
        # m3 = model(beta=2.224489795918367374, kappa=10.98541141987558412, lam=0.)

        # m1.distance_modulus()
        # m2.distance_modulus()
        # m3.distance_modulus()

        # fig = plt.figure(constrained_layout=False)
        # frame1 = fig.add_axes((.1, .3, .8, .6))
        # plt.errorbar(
        #     z_sn, df['m_b_corr'], yerr=snerr, lw=0.5, ls='', marker='.',
        #     markersize=2, label=r'Pantheon+SH0ES', zorder=0
        # )

        # mb_astro, _ = mbcorr_fit(dm_astro)
        # plt.plot(
        #     z_sn, mb_astro, c='orange', ls='-',
        #     label=r'$\Lambda$CDM',
        # )
        
        # plt.plot(
        #     z_sn, m.mb_int, c='k', ls='-.',
        #     label=r'Alt. Model, $\beta = {:.2f}$'.format(m.b),
        # )

        # # plt.plot(
        # #     z_sn, m2.mb_int, c='r', ls='--',
        # #     label=r'Model 2',
        # # )

        # # plt.plot(
        # #     z_sn, m3.mb_int, c='b', ls=':',
        # #     label=r'Model 3',
        # # )


        # # plt.plot(
        # #     z_sn, m2.mb_int, c='r', ls='--',
        # #     label=r'$\beta=1.98, \mu = $',
        # # )
        
        # plt.ylabel(r'$m_{B}$ [mag]')
        # plt.xscale('log')
        # plt.tick_params(axis='both', which='both', direction='in',
        #                 bottom=True, top=True, left=True, right=True)
        # plt.legend()
        # plt.grid()
        # plt.ylim([7, 28])
        # frame1.set_xticklabels([])

        # frame2 = fig.add_axes((.1, .1, .8, .2))
        # plt.errorbar(
        #     z_sn, df['m_b_corr'] - mb_astro, yerr=snerr, lw=0.5, ls='',
        #     marker='.', markersize=2, label=r'Supernova',zorder=0
        # )
        # plt.plot(z_sn, (m.mb_int - mb_astro), c='k', ls='-.')
        # # plt.plot(z_sn, (m2.mb_int - mb_astro), c='r', ls='--')
        # # plt.plot(z_sn, (m3.mb_int - mb_astro), c='b', ls=':')

        # plt.xlabel(r'$z$')
        # plt.xscale('log')
        # plt.ylabel(r'$\Delta{m_{B}}_{\Lambda\mathrm{CDM}}$')
        # plt.tick_params(
        #     axis='both', which='both', direction='in',
        #     bottom=True, top=True, left=True, right=True
        # )
        # plt.tick_params(axis='x', pad=7)
        # # extraticks = [-1, 1]
        # # plt.yticks(list(plt.yticks()[0]) + extraticks)
        # # plt.ylim([-2, 2])
        # plt.ylim([-1, 1])
        # plt.grid(which='major', axis='both')
        # plt.show()

main()