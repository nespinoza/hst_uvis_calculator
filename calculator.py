import astropy.units as q

import matplotlib.pyplot as plt
import numpy as np

import utils

# Numbers for HP41:
texp_HP41 = 53.332#190. # Seconds
SNR_HP41 = 419.# 710.

# Teff, V and J magnitudes of the star you want to use for the calculation:
Teff = 6776. 
J = 9.625

# Transit parameters of your target:
transit_duration_correction = True
tdur = 2.9 # Hours

# Texp to reach SNR of SNR_HP41 (results from WFC3/UVIS ETC):
texp = 30 # seconds

# What follows, is based on the following idea: the error on the transit depth is given by 2 x (Photometric Precision) / sqrt(N), 
# where the photometric precision is in ppm, and N is the number of total datapoints in a time-series. The (Photometric Precision) 
# term is wavelength-dependant, and depends on the SED of the star as compared to HAT-P-41, after going through the instrument.
# N, on the other hand, is wavelength-independant and depends on comparing both the transit duration comparison between our planet 
# and HAT-P-41b *and* the exposure time of the observations.
#
# What we want to do, thus, is to scale the errorbars of HAT-P-41's transmission spectrum by the ratio of those two numbers. Let's 
# First gather the data for HAT-P-41b:
#
# First, V-magnitude of HAT-P-41b:
# J-magnitude for HAT-P-41b:
J_HP41 = 10.006  # From 2MASS.

# Effective temperature of HAT-P-41b:
Teff_HP41 = 6390.

# And transit parameters for HAT-P-41b:
tdur_HP41 = 4.09 # Hours

# Let's now compute the N on sqrt(N):
if transit_duration_correction:

    N = (texp_HP41 / texp) * (tdur / tdur_HP41)

else:

    N = (texp_HP41 / texp)

# And let's do the photometric precision term. For this, we are assuming both reach SNR_HP41 at 4500 angstroms, but this might 
# actually be better/worse at other wavelengths. So what we do is as follows. First, get the ratio between the spectra of HP41 
# and the target star, assuming both are solar-like:

w, spectrum_HP41 = utils.get_stellar_model(teff = Teff_HP41, jmag = J_HP41)
w, spectrum = utils.get_stellar_model(teff = Teff, jmag = J)

ratio = spectrum / spectrum_HP41

# We assume the conversion between flux ratio and electron ratio is more or less the same between the two, and thus, assume any 
# improvement in photometric precision enters as a sqrt(ratio) improvement on the precision. So, first, we normalize the sqrt of the 
# ratio to 1 at 0.45 microns (4500 angstroms), which is where we now at this exposure time both reach the same SNR:
idx = np.where(np.abs(w.value - 0.45) == np.min(np.abs(w.value - 0.45)))[0][0]

sqrt_ratio = np.sqrt(ratio)

norm_sqrt_ratio = sqrt_ratio / sqrt_ratio[idx] 

# And now the final improvement on the transit depth precision then is:
improvement = np.sqrt(N) * norm_sqrt_ratio

# Now, apply this to the transit depths obtained by Wakeford et al. (2018). First, extract those:
wavelength, depth, error = np.loadtxt('data/hat-p-41-depths.dat', unpack = True, usecols = (0,1,2))

# Convert depths and errors from % to ppm:
depth = (depth / 100. ) * 1e6
error = (error / 100. ) * 1e6

# Multiply errors by sqrt(2) because two visits were used originally:
error = error * np.sqrt(2)

# Now save in machine readable format:
fout = open('simulated_errors.dat', 'w')

fout.write('# Simulated precisions for star of '+str(Teff)+' K, J = '+str(J)+'\n')
fout.write('# \n')
fout.write('# Column 1 : Central wavelength (nm)\n')
fout.write('# Column 2 : Wavelength bin half-width (nm)\n')
fout.write('# Column 3 : Error on bin (parts-per-million)\n')

for i in range(len(wavelength)):

    if i == 0:

        delta_w = ( wavelength[1] - wavelength[0] ) * 0.5

    else:

        delta_w = np.abs( wavelength[i - 1] - wavelength[i] ) * 0.5

    w_nm = w.to(q.nm).value


    idx = np.where((w_nm > wavelength[i] - delta_w)&(w_nm < wavelength[i] + delta_w))[0]

    print(wavelength[i], delta_w)
    print(w_nm[idx])    

    average_improvement = np.mean(improvement[idx]) 

    fout.write('{0:.2f} {1:.2f} {2:.2f}\n'.format(wavelength[i], delta_w, error[i] / average_improvement))

fout.close()
