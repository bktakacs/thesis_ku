# Run a bunch of definitive stuff. I'm keeping this in a separate file
from main import *

# Run chi_search 4 times: 1: int-noeff, 2:int-eff, 3:tay-noeff, 4:tay-eff

m1 = chi_search(fname='optimization-int-noeff', lenght=50, blim=(1, 10), klim=(0.1, 100), dm_effort=False, dm_method='int')

m2 = chi_search(fname='optimization-int-eff', lenght=50, blim=(1, 10), klim=(0.1, 100), dm_effort=True, dm_method='int')

m3 = chi_search(fname='optimization-tay-noeff', lenght=50, blim=(1, 10), klim=(0.1, 100), dm_effort=False, dm_method='tay')

m4 = chi_search(fname='optimization-tay-eff', lenght=50, blim=(1, 10), klim=(0.1, 100), dm_effort=True, dm_method='tay')

m1.plot('acc', lcdm, matter)
m1.distance_modulus(effort=False)
m1.plot('dm', lcdm, matter)

m2.plot('acc', lcdm, matter)
m2.distance_modulus(effort=False)
m2.plot('dm', lcdm, matter)

m3.plot('acc', lcdm, matter)
m3.distance_modulus(effort=True)
m3.plot('dm', lcdm, matter)

m4.plot('acc', lcdm, matter)
m4.distance_modulus(effort=False)
m4.plot('dm', lcdm, matter)

