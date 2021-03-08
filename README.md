# The frequency ratios
Certain combinations of the stellar oscillation frequencies are particularly interesting because they are independent of uncertain aspects of the models. For the definition of such combinations of frequencies, known as frequency ratios, see Roxburgh & Vorontsov (2003). 

This repository can be used to cumpute various types of frequency ratios, for example r01, r10, r02, r012 and r102. The code automatically computes **all possible ratios** of a given type. The code also allows the computation of their **covariance matrices** assuming uncertainties on the oscillation frequencies are normally distributed.

The script *example.py* explains the functioning of the code through examples.

# Troubleshooting 
The current implimentation fails if the input frequencies have missing radial order.
