# # Run a bunch of definitive stuff. I'm keeping this in a separate file
from main import *


# @timer
def main():
    """
    Main function
    """

    if __name__ == '__main__':

        # b = np.linspace(1, 4, 50)
        # # k = np.linspace(0.1, 1000, 51)
        # # k = np.logspace(-1, 3, 50)
        # k = np.geomspace(0.1, 1000, 50)

        # data = read_model_data('new-friedmann-1.txt')
        # chival = np.array(data['chi'])
        # brange = np.unique(np.array(data['beta']))
        # krange = np.unique(np.array(data['kappa']))
        # length = int(np.sqrt(len(chival)))
        # beta_low, kappa_low = specific_function(data, 0)
        # acc = False
        # round = 1

        # chi_plot_z = np.reshape(chival, (length, length)).T

        # fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        # ax.set_yscale('log')
        # cmap = matplotlib.cm.get_cmap('viridis_r').copy()
        # cmap.set_bad(color='w')
        # # im = ax.imshow(chi_plot_z, cmap=cmap, origin='lower',
        # #                interpolation='nearest',
        # #                norm=scale,
        # #                )
        # im = ax.pcolormesh(
        #     b, k, chi_plot_z, norm=None, edgecolors='w',
        #     cmap='viridis_r'
        # )
        # index = np.argmin(np.abs(k - kappa_low))
        # highlight_cell(
        #     beta_low-0.5*np.diff(b)[0],
        #     kappa_low-0.5*(k[index] - k[index-1]),
        #     np.diff(b)[0],
        #     (k[index] - k[index-1]),
        #     ax=ax, color='k', linewidth=1
        # )
        # ax.set_xlabel(r'$\beta$')
        # ax.set_ylabel(r'$\mu$')
        # fig.colorbar(im, ax=ax, label=r'$\chi^{{2}}_{{r, {}}}$'.format(
        #     '\ddot{a}/\ddot{a}_{\mathrm{M}}' if acc
        #     else 'E(z)'
        # ))
        # ax.tick_params(
        #     axis='both', which='both', direction='in',
        #     bottom=True, top=True, left=True, right=True
        # )
        # plt.show()

        m1 = chi_search(
            # 'nosave', length=50, blim=(1, 4), klim=(0.1, 10), solver='BDF',
            'new-friedmann-before30', length=50, blim=(1, 4), klim=(0.1, 1000),
            solver='BDF', acc=False,
        )

        # m = model(beta=0.6, kappa=1, lam=0.)
        # m.plot('acc')

main()