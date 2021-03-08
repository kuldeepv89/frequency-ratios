import numpy as np
import ratios



fname = '/Users/kuldeep/Dropbox/seismology/obs_data/kic012069424/freq/freq.dat'
names = ['l', 'n', 'freq', 'err']
fmts = [int, int, float, float]
obsFreq = np.genfromtxt(fname, dtype={'names': names, 'formats': fmts})

r02, r01, r10, r010, r012, r102 = ratios.ratios(obsFreq)
print(r02)

#obsr012, covr012, icovr012 = ratios.ratio_and_cov(obsFreq, rtype='R02', nrealizations=10000)
#print (obsr012)
