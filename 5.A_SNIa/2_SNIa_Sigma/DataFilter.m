function [Inds,Data]=DataFilter(Data,r_vect,FixParams)
% This function filters out SNe (rows) from 'Data' that deviate from a
% reference model by more than 3 sigma. This process is usually referred to
% as 'sigma clipping', and is commonly applied when fitting a cosmological
% model to various databases of standardized candles (e.g. Ia SNe, GRB,
% QSOs, etc.). Here our reference model is a Timesphere model of the
% universe, with the only parameter, the h=H0/(100 km/(s*Mpc)) reduced 
% Hubble constant set to h=0.62339 based on its best-fit value we obtained 
% from the cosmic chronometer fit (see. Sharov, G. S., & Vasiliev, V. O., 
% MMG 6, 1, 1 (2018), eprint arXiv:1807.07323). We also set the values of 
% Alpha, Beta, Delta, MB, their uncertainties (sigma_[parameter]) and 
% sigma_ext (for a detailed description of these, please check the comments 
% included in file 'parameter_limits.dat' and paper Suzuki et al., ApJ 746, 
% 1, 85, 24 (2012), eprint arXiv:1105.3470). Values of sigma_sample for the 
% various SNIa surveys are also set from input parameter 
% 'sigma_sample_data'. Our test statistic on which the 3 sigma cut is 
% applied is the weighted root-sum-square deviation between the distance 
% modulus calculated from SNe data and from the reference cosmological 
% model. Weights applied in the test statistic are the sum of squared 
% sigma_lc, sigma_ext, and sigma_sample, where these terms are described in 
% details in Suzuki et al. 2012.
%
%
% Inputs for DataFilter:
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
% FixParams - SNIa parameters Alpha, Beta, Delta, MB, sigma_Alpha, 
% sigma_Beta, sigma_Delta, sigma_MB, and sigma_sample, organized in a row 
% vector, and taken from the fitting process preceding this run. For a
% description of these parameters, see the header of 
% 'parameter_limits.dat', paper Suzuki et al. (2012), and the header of
% 'sigma_sample_data.dat' in folder '1_SNIa_Fitter'.
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
% Outputs of DataFilter:
% 
% Inds - a vector that contains the row indices of SNe that DataFilter 
% filtered out. If no SNIa was filtered out, then Inds is an empty vector.
%
% Data - data matrix that contains SNIa data that was left untouched by 
% DataFilter's filtering process
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the values of SNIa parameters and their uncertainties. These are
% obtained from the previous iteration of the fitting process. The values 
% of these parameters should be changed after each iteration of the 
% fitting process.
Alpha=FixParams(1);
Beta=FixParams(2);
Delta=FixParams(3);
MB=FixParams(4);
sigma_Alpha=FixParams(5);
sigma_Beta=FixParams(6);
sigma_Delta=FixParams(7);
sigma_MB=FixParams(8);
sigma_sample=FixParams(9);

% Setting the values of sigma_ext and h. These should be left unchanged
% throughout the iterative fitting process.
sigma_ext=0.05;
h=0.62339;

% Setting the value of the speed of light
c=299792458; % [in m/s]

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
sigma_lc_sq=sigma_MB^2+(P.^2)*sigma_Delta^2+(color.^2)*sigma_Beta^2+(x1.^2)*sigma_Alpha^2+sigma_mB.^2+(Alpha^2)*sigma_x1.^2+(Beta^2)*sigma_color.^2+2*Alpha*r_mB_x1*sigma_mB.*sigma_x1-2*Beta*r_mB_color*sigma_mB.*sigma_color-2*Alpha*Beta*r_x1_color*sigma_x1.*sigma_color;

% Calcuting the test statistic
Chi=sqrt(((mu_meas-mu_th).^2)./(sigma_sample.^2+sigma_lc_sq+sigma_ext^2));

% Filtering out all SNe from the Data matrix that deviates from the
% reference model by more than 3 sigma.
Inds=find(Chi>3);
Data(Inds,:)=[];

% Making sure that 'Inds' is returned as a column vector
Inds=Inds(:);
