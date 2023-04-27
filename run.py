# # Run a bunch of definitive stuff. I'm keeping this in a separate file
from main import *

# # Run chi_search 4 times: 1: int-noeff, 2:int-eff, 3:tay-noeff, 4:tay-eff

# m1 = chi_search(fname='optimization-int-noeff-1', length=50, blim=(1, 10), klim=(0.1, 100), dm_effort=False, dm_method='int', plot=False)

# m2 = chi_search(fname='optimization-int-eff-1', length=50, blim=(1, 10), klim=(0.1, 100), dm_effort=True, dm_method='int', plot=False)

# m3 = chi_search(fname='optimization-tay-noeff-1', length=50, blim=(1, 10), klim=(0.1, 100), dm_effort=False, dm_method='tay', plot=False)

# m4 = chi_search(fname='optimization-tay-eff-1', length=50, blim=(1, 10), klim=(0.1, 100), dm_effort=True, dm_method='tay', plot=False)

# m1.plot('acc', lcdm, matter)
# m1.distance_modulus(effort=False)
# m1.plot('dm', lcdm, matter)

# m2.plot('acc', lcdm, matter)
# m2.distance_modulus(effort=False)
# m2.plot('dm', lcdm, matter)

# m3.plot('acc', lcdm, matter)
# m3.distance_modulus(effort=True)
# m3.plot('dm', lcdm, matter)

# m4.plot('acc', lcdm, matter)
# m4.distance_modulus(effort=False)
# m4.plot('dm', lcdm, matter)

# Testing double optimization

# m1 = chi_search(fname='double_optimization_test', length=50, blim=(1, 10), klim=(0.1, 100), dm_effort=False, dm_method='int', plot=True)
# m1 = chi_search(fname='double_optimization_int_noeff_fine', length=50, blim=(1.4, 1.7), klim=(0.1, 4.2),
#                 dm_effort=False, dm_method='int', plot=True)

# m1.distance_modulus(effort=False)

# m1.plot('acc', lcdm, matter)
# m1.plot('dm', lcdm, matter)


# Testing q surface function
# m1 = q_surface(length=30, blim=(1, 2), klim=(0.1, 7))




# Steen has asked about fixing one value (beta, for example) and then plotting
# a bunch of values of kappa and seeing what happens or trying to get close
# to a decent looking plot. Let's try that.

# Start by fixing beta and looking at different kappa
beta = 3
# klook = (0.1, 0.5, 1, 2, 3, 4, 5)
klook = (0.1, 0.2, 0.3, 0.4, 0.5)

for i in tqdm(klook):
    temp_mod = model(beta=beta, kappa=i, lam=0.)
    temp_mod.norm(matter=matter)
    plt.plot(temp_mod.a, temp_mod.a2norm, label='kappa = {}'.format(i))

plt.xlabel(r'$a$')
plt.ylabel(r'$\ddot{a}$')
plt.ylim([-5, 2])
plt.legend()
plt.grid()
plt.tick_params(axis='both', which='both', direction='in', bottom=True,
                top=True, left=True, right=True)
plt.show()