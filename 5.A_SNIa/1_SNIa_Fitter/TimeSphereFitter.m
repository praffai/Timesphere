function Chi=TimeSphereFitter(Params,Data,r_vect,sigma_sample_data)
% This function calculates the measure of goodness of fit (in our case, the
% weighted sum of squared deviations between data and the reference model, 
% a.k.a. the chi-squared statistic) for parameter values input in vector
% 'Params'. The reference cosmological model is a Timesphere model of the
% universe, with the only parameter, the h=H0/(100 km/(s*Mpc)) reduced 
% Hubble constant set to h=0.62339 based on its best-fit value we
% obtained from the cosmic chronometer fit (see. Sharov, G. S., & Vasiliev,
% V. O., MMG 6, 1, 1 (2018), eprint arXiv:1807.07323). We also set the
% value of sigma_ext to sigma_ext=0.05, and the values of sigma_sample to
% the ones given in input matrix sigma_sample_data (for a detailed 
% description of sigma_ext and sigma_sample, please check paper Suzuki et 
% al., ApJ 746, 1, 85, 24 (2012), eprint arXiv:1105.3470). 
%
%
% Inputs for TimeSphereFitter:
%
% Params - a vector that contains the values of Alpha, Beta, Delta, and MB,
% (for a detailed description of these parameters, see the header of 
% 'parameter_limits.dat' and paper Suzuki et al. 2012) for which the 
% measure of goodness of fit is to be calculated by TimeSphereFitter
%
% Data - the matrix of SNe data loaded from 'SNIa_data.dat', where each row
% corresponds to a single SNIa event. For a description of the matrix 
% columns, see the header of the ascii text file 'SNIa_data.dat'.
%
% r_vect - a vector of Pearson correlation coefficients between various 
% SNIa parameters. The coefficients are between mB and x1, mB and color,
% and x1 and color, respectively. They are used in calculating the 
% covariant terms of the uncertainties of measured SNIa distance moduli.
%
% sigma_sample_data - the matrix of sigma_sample values for the various 
% SNIa survey samples, loaded from 'sigma_sample_data.dat'. The first
% column contains integer numbers from 1 to 19 corresponding to the various
% SNIa surveys (the names of the surveys are listed in a reversed order in 
% Figure 3 of Suzuki et al. 2012). The second column contains the
% sigma_sample values for the surveys. These values start from 
% sigma_sample=0.15 (the median of sigma_sample values in Suzuki et al. 
% 2012), and then are recalculated and changed in every iterations of the 
% fitting process until they and the fitted parameters converge to a stable
% set of values.
%
% sigma_ext, h - hardcoded parameters of the SNIa surveys and the reference
% cosmological model. 'sigma_ext' is described in Suzuki et al. 2012, and
% should be kept to be sigma_ext=0.05 (this is an approximate value we
% chose based on the various terms that contribute to sigma_ext). 
% h=H0/(100 km/(s*Mpc)) is the reduced Hubble constant that we set to
% h=0.62339 based on its best-fit value we obtained from the cosmic 
% chronometer fit (see. Sharov, G. S., & Vasiliev, V. O., MMG 6, 1, 1 
% (2018), eprint arXiv:1807.07323).
%
% c - a hardcoded parameter, it is the speed of light in vacuum expressed 
% in m/s units (thus, c=299792458).
%
%
% Output of TimeSphereFitter:
% Chi - the measure of goodness of fit (in our case, the weighted sum of 
% squared deviations between data and the reference model, a.k.a. the 
% chi-squared statistic)
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the values of sigma_ext and h. These should be left unchanged
% throughout the iterative fitting process.
sigma_ext=0.05;
h=0.62339;

% Setting the value of the speed of light
c=299792458; % [in m/s]

%  Extracting the values of SNIa parameters from input vector 'Params'
Alpha=Params(1);
Beta=Params(2);
Delta=Params(3);
MB=Params(4);

% Defining the different SNIa parameters from the different columns of 
% input data matrix 'Data'. For a detailed description of these parameters,
% see the header of data file 'SNIa_data.dat'.
z=Data(:,1);
mB=Data(:,2);
sigma_mB=Data(:,3);
x1=Data(:,4);
sigma_x1=Data(:,5);
color=Data(:,6);
sigma_color=Data(:,7);
P=Data(:,8);
sample=Data(:,9);

% Creating a vector of sigma_sample values, where each vector element 
% corresponds to a single SNIa within a survey. Parameter 'sample' contains
% the survey identification number (1-19) for each SNIa.
sigma_sample=sigma_sample_data(sample,2);

% Extracting the Pearson correlation coefficients from the input vector
r_mB_x1=r_vect(1);
r_mB_color=r_vect(2);
r_x1_color=r_vect(3);

% Calculating the distance moduli for all SNe. For details about this, see
% Suzuki et al. 2012.
mu_meas=mB+Alpha*x1-Beta*color+Delta*P-MB;

% Calculating the distance moduli for the reference cosmological model at
% all redshifts of the full SNIa sample.
mu_th=5*log10(c*(1+z).*sin(log(1+z)))-5*log10(h);

% Calculating the uncertainties of all mu_meas values. For details about
% this, see Suzuki et al. 2012 and references therein.
sigma_lc_sq=sigma_mB.^2+(Alpha^2)*sigma_x1.^2+(Beta^2)*sigma_color.^2+2*Alpha*r_mB_x1*sigma_mB.*sigma_x1-2*Beta*r_mB_color*sigma_mB.*sigma_color-2*Alpha*Beta*r_x1_color*sigma_x1.*sigma_color;

% Calcuting the measure of goodness of fit (in our case, the weighted sum 
% of squared deviations between data and the reference model, a.k.a. the 
% chi-squared statistic)
Chi=sum(((mu_meas-mu_th).^2)./(sigma_lc_sq+sigma_ext^2+sigma_sample.^2));
