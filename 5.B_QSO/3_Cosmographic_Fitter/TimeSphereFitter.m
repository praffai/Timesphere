function Chi=TimeSphereFitter(Params,Data)
% This function calculates the measure of goodness of fit (in our case, the
% weighted sum of squared deviations between data and the reference model, 
% a.k.a. the chi-squared statistic) for parameter values input in vector
% 'Params'. The reference model is a cosmographic approximation of the 
% distance modulus vs. redshift curve, and the QSO distance modulus data
% has previously been fitted to the baseline Timesphere model of the
% universe with the only parameter, the h=H0/(100 km/(s*Mpc)) reduced 
% Hubble constant set to h=0.62339 based on its best-fit value we
% obtained from the cosmic chronometer fit (see. Sharov, G. S., & Vasiliev,
% V. O., MMG 6, 1, 1 (2018), eprint arXiv:1807.07323). The cosmographic
% parameters fitted in this process are 'a2' and 'a3, described in Risaliti
% G., & Lusso, E., Nature Astronomy 3, p. 272-277 (2019), eprint 
% arXiv:1811.02590.
%
%
% Inputs for TimeSphereFitter:
%
% Params - a vector that contains the value of parameters for which the 
% measure of goodness of fit is to be calculated by TimeSphereFitter.
%
% Data - the matrix of QSO data loaded from 'QSO_TimeSphere_Data.dat', 
% where each row corresponds to a single QSO event. For a description of 
% the matrix columns, see the header of the ascii text file 
% 'QSO_TimeSphere_Data.dat'.
%
% h - hardcoded parameter of the reference cosmological model. 
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

% Setting the values of fitted parameters. 
a2=Params(1);
a3=Params(2);

% Setting the value of h. This should be left unchanged throughout the 
% iterative fitting process.
h=0.62339;

% Setting the value of the speed of light.
c=299792458; % [in m/s]

% Defining the different QSO parameters from the different columns of 
% input data matrix 'Data'. For a detailed description of these parameters,
% see the header of data file 'QSO_TimeSphere_Data.dat'.
z=Data(:,1);
mu_meas=Data(:,2);
sigma_mu=Data(:,3);

% Calculating the distance moduli for the reference cosmological model at
% all redshifts of the QSO sample.
x=log10(1+z);
mu_th=5*log10(x+a2*x.^2+a3*x.^3)+5*log10(log(10)*c)-5*log10(h);

% Calcuting the test statistic
Chi=sum(((mu_meas-mu_th).^2)./(sigma_mu.^2));

