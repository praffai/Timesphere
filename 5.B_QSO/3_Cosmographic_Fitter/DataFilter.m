function [Inds,Data]=DataFilter(Data)
% This function filters out QSOs (rows) from 'Data' that deviate from a
% reference model by more than 3 sigma. This process is usually referred to
% as 'sigma clipping', and is commonly applied when fitting a cosmological
% model to various databases of standardized candles (e.g. Ia SNe, GRB,
% QSOs, etc.). Here The reference model is a cosmographic approximation of 
% the distance modulus vs. redshift curve, and the QSO distance modulus 
% data has previously been fitted to the baseline Timesphere model of the
% universe with the only parameter, the h=H0/(100 km/(s*Mpc)) reduced 
% Hubble constant set to h=0.62339 based on its best-fit value we
% obtained from the cosmic chronometer fit (see. Sharov, G. S., & Vasiliev,
% V. O., MMG 6, 1, 1 (2018), eprint arXiv:1807.07323). We also set the 
% values of a2, a3, and their uncertainties sigma_a2 and sigma_a3. Our test 
% statistic on which the 3 sigma cut is applied is the weighted 
% root-sum-square deviation between the distance modulus calculated from 
% QSOs data and from the reference cosmological model.
%
%
% Inputs for DataFilter:
%
% Data - the matrix of QSO data loaded from 'QSO_TimeSphere_Data.dat', 
% where each row corresponds to a single QSO event. For a description of 
% the matrix columns, see the header of the ascii text file 
% 'QSO_TimeSphere_Data.dat'.
%
% a2, a3, sigma_a2, sigma_a3 - hardcoded parameters obtained from the 
% fitting process. Their values should be updated after every iterations of
% the process.
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

% Setting the values of fitted parameters a2, a3, and their uncertainties. 
% These are obtained from the previous iteration of the fitting process 
% (thus the data filtering is not applied when the fitting process is run 
% first). The values of these parameters should be changed after each 
% iteration of the fitting process.
a2=2.20;% this is the final value of 'a2' we obtained from the iterative fitting process. a2=ln(10)=2.30 for the Timesphere model.
a3=2.31;% this is the final value of 'a3' we obtained from the iterative fitting process. a3=(ln(10)^2)/3=1.78 for the Timesphere model.
sigma_a2=0.11;
sigma_a3=0.31;

% Setting the value of h. This should be left unchanged throughout the 
% iterative fitting process.
h=0.62339;

% Setting the value of the speed of light
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

% Calculating the combined sigma of parameters.
y=1+a2*x+a3*x.^2;
sigma_p_sq=((25*(x.^2)*log(10)^2)./y.^2).*(sigma_a2^2+(sigma_a3^2)*x.^2);

% Calcuting the test statistic.
Chi=sqrt(((mu_meas-mu_th).^2)./(sigma_mu.^2+sigma_p_sq));

% Filtering out all QSOs from the Data matrix that deviates from the
% reference model by more than 3 sigma.
Inds=find(Chi>3);
Data(Inds,:)=[];

% Making sure that 'Inds' is returned as a column vector
Inds=Inds(:);

