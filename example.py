import numpy as np
import ratios


# Load the data
names = ['l', 'n', 'freq', 'err']
fmts = [int, int, float, float]
obsFreq = np.genfromtxt('16CygA.txt', dtype={'names': names, 'formats': fmts})

# Compute ratios
r02, r01, r10 = ratios.ratios(obsFreq)

# Compute ratios as well as corresponding covariance matrix
obsr02, covr02 = ratios.ratio_and_cov(obsFreq, rtype='R02', nrealizations=10000)
