function Chi=TimeSphereFitter(Params,Data)
% This function calculates the measure of goodness of fit (in our case, the
% weighted sum of squared deviations between data and the reference model, 
% a.k.a. the chi-squared statistic) for parameter values input in vector
% 'Params'. The reference cosmological model is a Timesphere model of the
% universe, with the only parameter, the h=H0/(100 km/(s*Mpc)) reduced 
% Hubble constant set to h=0.62339 based on its best-fit value we
% obtained from the cosmic chronometer fit (see. Sharov, G. S., & Vasiliev,
% V. O., MMG 6, 1, 1 (2018), eprint arXiv:1807.07323). The only parameter
% fitted in this process is 'mu_shift', which is a vertical shift applied 
% on all distance moduli of QSOs.
%
%
% Inputs for TimeSphereFitter:
%
% Params - a vector that contains the value of mu_shift for which the 
% measure of goodness of fit is to be calculated by TimeSphereFitter.
%
% Data - the matrix of QSO data loaded from 'QSO_data.dat', where each row
% corresponds to a single SNIa event. For a description of the matrix 
% columns, see the header of the ascii text file 'QSO_data.dat'.
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

% Setting the value of fitted parameter mu_shift. 
mu_shift=Params(1);

% Setting the value of h. This should be left unchanged throughout the 
% iterative fitting process.
h=0.62339;

% Setting the value of the speed of light
c=299792458; % [in m/s]

% Defining the different QSO parameters from the different columns of 
% input data matrix 'Data'. For a detailed description of these parameters,
% see the header of data file 'QSO_data.dat'.
z=Data(:,1);
mu=Data(:,2);
sigma_mu=Data(:,3);

% Calculating the distance moduli for all QSOs.
mu_meas=mu+mu_shift;

% Calculating the distance moduli for the reference cosmological model at
% all redshifts of the QSO sample.
mu_th=5*log10(c*(1+z).*sin(log(1+z)))-5*log10(h);

% Calculating the uncertainties of all mu_meas values.
sigma_lc_sq=sigma_mu.^2;

% Calcuting the test statistic
Chi=sum(((mu_meas-mu_th).^2)./(sigma_lc_sq));

