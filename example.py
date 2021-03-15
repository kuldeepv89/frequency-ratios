import numpy as np
import ratios


# whether it is observation, futh path to the data 
#-------------------------------------------------
isObs, fname = True, '16CygA.txt'


if isObs:

    # Load the observed data
    names = ['l', 'n', 'freq', 'err']
    fmts = [int, int, float, float]
    
    freq = np.genfromtxt(fname, dtype={'names': names, 'formats': fmts})
    
    
    # Compute r02 as well as the corresponding covariance matrix
    r02, covr02 = ratios.ratio_and_cov(freq, rtype='R02', nrealizations=10000)
    
    print ()
    print ('Printing observed r02...')
    for i in range(r02.shape[0]):
        print (4 * '%12.5f' %(r02[i, 0], r02[i, 1], r02[i, 2], r02[i, 3]))
    
    
    # Compute r01 as well as the corresponding covariance matrix
    r01, covr01 = ratios.ratio_and_cov(freq, rtype='R01', nrealizations=10000)
    
    print ()
    print ('Printing observed r01...')
    for i in range(r01.shape[0]):
        print (4 * '%12.5f' %(r01[i, 0], r01[i, 1], r01[i, 2], r01[i, 3]))

else:

    # Load the model data
    names = ['l', 'n', 'freq', 'iner']
    fmts = [int, int, float, float]
    
    freq = np.genfromtxt(fname, dtype={'names': names, 'formats': fmts})
    

    # Compute r02, r01 and r10
    r02, r01, r10 = ratios.ratios(freq)

    print ()
    print ('Printing model r02...')
    for i in range(r02.shape[0]):
        print (4 * '%12.5f' %(r02[i, 0], r02[i, 1], r02[i, 2], r02[i, 3]))

    print ()
    print ('Printing model r01...')
    for i in range(r01.shape[0]):
        print (4 * '%12.5f' %(r01[i, 0], r01[i, 1], r01[i, 2], r01[i, 3]))
