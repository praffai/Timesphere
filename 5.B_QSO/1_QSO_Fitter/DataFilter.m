function [Inds,Data]=DataFilter(Data,r_vect)
% This function filters out QSOs (rows) from 'Data' that deviate from a
% reference model by more than 3 sigma. This process is usually referred to
% as 'sigma clipping', and is commonly applied when fitting a cosmological
% model to various databases of standardized candles (e.g. Ia SNe, GRB,
% QSOs, etc.). Here our reference model is a Timesphere model of the
% universe, with the only parameter, the h=H0/(100 km/(s*Mpc)) reduced 
% Hubble constant set to h=0.62339 based on its best-fit value we obtained 
% from the cosmic chronometer fit (see. Sharov, G. S., & Vasiliev, V. O., 
% MMG 6, 1, 1 (2018), eprint arXiv:1807.07323). We also set the values of 
% gamma, beta, and their uncertainties sigma_gamma and sigma_beta. Our test 
% statistic on which the 3 sigma cut is applied is the weighted 
% root-sum-square deviation between the distance modulus calculated from 
% QSOs data and from the reference cosmological model.
%
%
% Inputs for DataFilter:
%
% Data - the matrix of QSO data loaded from 'QSO_data.dat', where each row
% corresponds to a single QSO event. For a description of the matrix 
% columns, see the header of the ascii text file 'QSO_data.dat'.
%
% r_vect - the Pearson correlation coefficient between QSO parameters 
% logFUV and logFX. It is used in calculating the covariant term of the 
% uncertainties of measured QSO distance moduli.
%
% gamma, beta, sigma_gamma, sigma_beta - hardcoded parameters obtained from
% the fitting process. Their values should be updated after every 
% iterations of the process.
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
% Outputs of DataFilter:
% 
% Inds - a vector that contains the row indices of QSOs that DataFilter 
% filtered out. If no QSO was filtered out, then Inds is an empty vector.
%
% Data - data matrix that contains QSO data that was left untouched by 
% DataFilter's filtering process
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the values of fitted parameters gamma, beta, and their 
% uncertainties. These are obtained from the previous iteration of the 
% fitting process (thus the data filtering is not applied when the fitting 
% process is run first). The values of these parameters should be changed 
% after each iteration of the fitting process.
gamma=0.707;%0.663;%0.634;
beta=57.4;%-7.06;%-7.33;
sigma_gamma=0.001;
sigma_beta=0.2;

% Setting the value of h. This should be left unchanged throughout the 
% iterative fitting process.
h=0.62339;

% Setting the value of the speed of light
c=299792458; % [in m/s]

% Defining the different QSO parameters from the different columns of 
% input data matrix 'Data'. For a detailed description of these parameters,
% see the header of data file 'QSO_data.dat'.
z=Data(:,1);
logFUV=Data(:,2);
sigma_logFUV=Data(:,3);
logFX=Data(:,4);
sigma_logFX=Data(:,5);

% Extracting the Pearson correlation coefficient from input vector
r_logFUV_logFX=r_vect(1);

% Calculating the distance moduli for all QSOs.
mu_meas=(5/(2*(gamma-1)))*(logFX-gamma*logFUV)-beta;

% Calculating the distance moduli for the reference cosmological model at
% all redshifts of the QSO sample.
mu_th=5*log10(c*(1+z).*sin(log(1+z)))-5*log10(h);

% Calculating the uncertainties of all mu_meas values.
sigma_lc_sq=((5/(2*(gamma-1)))^2)*sigma_logFX.^2+((5*gamma/(2*(gamma-1)))^2)*sigma_logFUV.^2+sigma_beta^2+((5/(2*(gamma-1)^2))^2)*((logFX-logFUV).^2)*sigma_gamma^2-(5*gamma/(gamma-1))*r_logFUV_logFX*sigma_logFX.*sigma_logFUV;

% Calcuting the test statistic
Chi=sqrt(((mu_meas-mu_th).^2)./(sigma_lc_sq));

% Filtering out all QSOs from the Data matrix that deviates from the
% reference model by more than 3 sigma.
Inds=find(Chi>3);
Data(Inds,:)=[];

% Making sure that 'Inds' is returned as a column vector
Inds=Inds(:);

