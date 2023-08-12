# # Run a bunch of definitive stuff. I'm keeping this in a separate file
from main import *


# @timer
def main():
    """
    Main function
    """

    if __name__ == '__main__':

        m = chi_search(
            'new-friedmann-eq30-bdf-extended-super.txt',
            length=50, blim=(1, 3), klim=(0.01, 1000), solver='BDF'
        )
        # m.plot('acc')
        # m.distance_modulus()
        # m.plot('dm')

        # data = read_model_data('new-friedmann-eq30-bdf-extended.txt')
        # chival = np.array(data['chi'])
        # length = int(np.sqrt(len(chival)))
        # beta_low, kappa_low = specific_function(data, 0)
        # brange = np.linspace(1, 4, 50)
        # krange = np.logspace(-2, 3, 50)
        # # print((chival), length)
        # chi_plot_z = np.reshape(chival, (length, length)).T

        # fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        # ax.set_yscale('log')
        # cmap = matplotlib.cm.get_cmap('viridis_r').copy()
        # cmap.set_bad(color='w')
        # im = ax.pcolormesh(
        #     brange, krange, chi_plot_z, norm=None, edgecolors='w',
        #     cmap='viridis_r'
        # )
        # index = np.argmin(np.abs(krange - kappa_low))
        # highlight_cell(
        #     beta_low-0.5*np.diff(brange)[0],
        #     kappa_low-0.5*(krange[index] - krange[index-1]),
        #     np.diff(brange)[0],
        #     (krange[index+1] - krange[index]),
        #     ax=ax, color='k', linewidth=1
        # )
        # for index in np.argwhere(mask_lower == 1):
        #     plus_cell(
        #         np.repeat(brange, length)[index],
        #         np.tile(krange, length)[index],
        #         0.8 * np.diff(brange)[0],
        #         0.8 * (krange[index+1] - krange[index]),
        #         ax=ax, color='k', linewidth=1
        #     )
        # ax.set_xlabel(r'$\beta$')
        # ax.set_ylabel(r'$\mu$')
        # fig.colorbar(im, ax=ax, label=r'$\chi^{{2}}_{{r, {}}}$'.format(
        #     'E(z)'
        # ))
        # ax.tick_params(
        #     axis='both', which='both', direction='in',
        #     bottom=True, top=True, left=True, right=True
        # )
        # plt.show()

        # chi_search('nosave', length=10)



main()