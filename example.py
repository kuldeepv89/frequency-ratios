import numpy as np
import ratios


# Load the data
names = ['l', 'n', 'freq', 'err']
fmts = [int, int, float, float]

obsFreq = np.genfromtxt('16CygA.txt', dtype={'names': names, 'formats': fmts})


# Compute r01 as well as the corresponding covariance matrix
print ()
print ('Calculate r01...')

obsr01, covr01 = ratios.ratio_and_cov(obsFreq, rtype='R01', nrealizations=10000)

for i in range(obsr01.shape[0]):
    print (4 * '%12.5f' %(obsr01[i, 0], obsr01[i, 1], obsr01[i, 2], obsr01[i, 3]))


# Compute r02 as well as the corresponding covariance matrix
print ()
print ('Calculate r02...')

obsr02, covr02 = ratios.ratio_and_cov(obsFreq, rtype='R02', nrealizations=10000)

for i in range(obsr02.shape[0]):
    print (4 * '%12.5f' %(obsr02[i, 0], obsr02[i, 1], obsr02[i, 2], obsr02[i, 3]))
