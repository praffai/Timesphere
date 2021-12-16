'''
This is the code you should run for fitting the cosmic chronometer dataset to the reference LambdaCDM and Timesphere cosmologies. Chrono_fit.py calculates the measure of goodness of fit (in our case, the weighted sum of squared deviations between data and the reference model, a.k.a. the chi-squared statistic) for different input values of cosmological parameters H0 (Hubble constant, for both the Timesphere and the LambdaCDM models) and Omega_m (dimensionless matter density parameter, for the LambdaCDM model only) by fitting the corresponding model curves to z (redshift), H(z) (Hubble parameter), and sigma_H (errors on H(z)) data arrays obtained from the cosmic chronometer data table in Sharov, G. S., & Vasiliev, V. O., MMG 6, 1, 1 (2018), eprint arXiv:1807.07323 (see data under title 'DA method' in Table II). We obtain the chi-square values for the LambdaCDM and Timesphere models of the universe using two different functions: chisqfunc_lcdm and chisqfunc_timesphere. We then search for the Hubble parameter and matter density parameter values for which chi-square is the smallest in the two cosmologies.


Inputs for Chrono_fit:

Hz_1807.07323.dat - input ascii file with three columns containing the z (redshift), H(z) (Hubble parameter), and sigma_H (errors on H(z)) values, respectively, obtained from the cosmic chronometer data table in Sharov, G. S., & Vasiliev, V. O., MMG 6, 1, 1 (2018), eprint arXiv:1807.07323 (see data under title 'DA method' in Table II). The Hubble parameter and sigma_H values are given in (km/s)/Mpc units. The file should only contain a single-line description of columns in the header, otherwise input data arrays are not read in properly.

c - a hardcoded parameter, it is the speed of light in vacuum expressed in km/s units (thus, c=299792.458).

Harray, H2array - arrays of potential Hubble parameter values given in (km/s)/Mpc units, which we use as input parameters in the chi-square computing functions. We use the H2array only for the Timesphere fit.

OMarray - an array of potential Omega_m values we use as an input parameter array for the LambdaCDM chi-square computing function.


Outputs of Chrono_fit:

chisq - the measure of goodness of fit (in our case, the weighted sum of squared deviations between data and the reference model, a.k.a. the chi-squared statistic)

E_z - the cosmological expansion function for the LambdaCDM model. It is a defined function of redshift, matter density parameter, cosmological constant density parameter and curvature density parameter (note that we always choose Omega_k=0 for the curvature density parameter). We use this function for the computation of chi-square values for the LambdaCDM model.


Credits: Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A. 
All rights reserved. (2021) 
Contact: peter.raffai@ttk.elte.hu
'''

# We import the necessary Python packages for running the code.
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# The best-fit Planck values from ”Planck 2018 results. VI. Cosmological parameters” (A&A 641, A6 (2020)):
# Hubble constant in (km/s)/Mpc: H0 = 67.66   
# dimensionless matter density parameter today: OM = 0.3111
# dimensionless cosmological constant density parameter today: OL = 0.6889
# dimensionless curvature density parameter today: OK = 0
c = 299792.458 # speed of light in km/s

# Importing input z, H(z), and sigma_H data from input file Hz_1807.07323.dat

f = open("Hz_1807.07323.dat", "r") # loading the given values
l = f.readlines()[1:]

z = []
Hz = []
sigma = []

for x in l:
	z.append(float(x.split(' ')[0]))
	Hz.append(float(x.split(' ')[1]))
	sigma.append(float(x.split(' ')[2]))

z = np.asarray(z)
Hz = np.asarray(Hz)
sigma = np.asarray(sigma)

# Defining the function for determining chisq values in Timesphere cosmology

def chisqfunc_timesphere(arr):
	a = arr[0]
	b = arr[1]
	model = a*z+b
	chisq = np.sum( ((Hz - model)/sigma)**2 )
	return chisq

# Defining E(z) for use in further calculations
	
def E_z(z, OM, OK, OL):
    return np.sqrt((1 + z)**3 * OM + (1 + z)**2 * OK + OL)
	
# Defining the function for determining chisq values in LambdaCDM cosmology
# the total density parameter equals 1
# the curvature density parameter is 0
# the radiation density parameter is negligible and ignored

def chisqfunc_lcdm(arr):
	H0 = arr[0]
	Omega_m0 = arr[1]
	Omega_l0 = arr[2]
	Omega_k0 = 0
	model = H0*E_z(z, Omega_m0, Omega_k0, Omega_l0)
	chisq = np.sum( ((Hz - model)/sigma)**2 )
	return chisq
	
# Creating arrays (2 with potential Hubble parameter and 1 with potential matter density parameter values) from which the ones resulting with the best fits are chosen

Harray = np.linspace(68, 69, num=1001)
H2array = np.linspace(62, 63, num=1001)
OMarray = np.linspace(0.31, 0.33, num=1001)

# These three parameters can be equal to 0 or 1, depending on which cosmological model we want to fit and with how many parameters.

tsfit = 0  # 1 if you want to fit Timesphere
twoparam = 0 # 1 if you want to fit Timesphere with 2 parameters
lcdmfit = 1 # 1 if you want to fit LambdaCDM


# Creating matrices with appropriate size

if (twoparam == 1):
	res_ts = np.zeros((len(Harray), len(H2array)))
else:
	res_ts = []

res_lcdm = np.zeros((len(Harray), len(OMarray)))

# Determining the values of matrix elements

for i in range(len(Harray)):

	if (tsfit == 1):

		if (twoparam == 1):
			for j in range(len(H2array)):
				res_ts[i][j] = chisqfunc_timesphere((Harray[i], H2array[j]))
		else:
			res_ts.append(chisqfunc_timesphere((Harray[i], Harray[i])))
			
	if (lcdmfit == 1):	
	
		for j in range(len(OMarray)):
			res_lcdm[i][j] = chisqfunc_lcdm((Harray[i], OMarray[j], 1-OMarray[j]))

res_ts = np.asarray(res_ts)

# Finding the parameter value(s) for which chisq is the smallest in Timesphere cosmology

if (tsfit == 1):
	if (twoparam == 1):
		(inda, indb) = np.unravel_index(res_ts.argmin(), res_ts.shape)
		print(Harray[inda], H2array[indb], res_ts[inda][indb])
	else:
		print(Harray[np.argmin(res_ts)], min(res_ts))

# Finding the Hubble parameter and matter density parameter values for which chisq is the smallest in LambdaCDM cosmology

if (lcdmfit == 1):
	(inda, indb) = np.unravel_index(res_lcdm.argmin(), res_lcdm.shape)
	print(Harray[inda], OMarray[indb], res_lcdm[inda][indb])
		
