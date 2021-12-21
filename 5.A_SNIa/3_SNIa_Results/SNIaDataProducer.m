function SNIaDataProducer
% This code produces ascii data file 'SNIa_TimeSphere_Data.dat' containing 
% SNe redshifts (Column #1), distance moduli according to the best-fit
% Timesphere model (Column #2), and distance moduli errors (Column #3). The
% output file is an analog of 'SCPUnion2.1_mu_vs_z.txt' (available for
% download on the Union2.1 Supernova Cosmology Project website: 
% http://supernova.lbl.gov/Union/; see under title 'Union2.1 Compilation 
% Magnitude vs. Redshift Table (for your own cosmology fitter)') in which
% the distance moduli are calculated for the best-fit LCDM cosmological
% model.
%
%
% Inputs for SNIaDataProducer:
%
% SNIa_data.dat - ascii data file produced by combining data in 
% 'SCPUnion2.1_AllSNe.tex' and 'SCPUnion2.1_mu_vs_z.txt', both of them 
% downloaded from the Union2.1 Supernova Cosmology Project website: 
% http://supernova.lbl.gov/Union/ (see under titles 'Full Table of All SNe'
% and 'Union2.1 Compilation Magnitude vs. Redshift Table (for your own 
% cosmology fitter)'). The data was originally used in Suzuki et al.
% (2012).
%
% final_sigma_sample_data.dat - an external ascii file defining the 
% sigma_sample values for the various SNIa survey samples (for details, see
% Suzuki et al., ApJ 746, 1, 85, 24 (2012); e-print: 
% https://arxiv.org/abs/1105.3470). The file contains the following 
% columns:
% Column #1 - an integer number from 1-19 indicating which survey sample 
% the SNIa data was taken from. The names of the surveys are listed in a 
% reversed order in Figure 3 of Suzuki et al. 2012.
% Column #2 - values of sigma_sample for the different surveys. 
% The file may potentially contain a third column, which gives the number 
% of SNe filtered out from the different survey samples by DataFilter.
% Column #3 (not used) - the number of SNe filtered out from the different 
% survey samples by DataFilter. A data point is filtered out if it deviates
% from a reference model by more than 3 sigma (this process is commonly 
% used under the name 'sigma clipping'). The reference model is a 
% Timesphere model of the universe, with the only parameter, the 
% h=H0/(100 km/(s*Mpc)) reduced Hubble constant set to h=0.62339 based on 
% its best-fit value we obtained from the cosmic chronometer fit (see. 
% Sharov, G. S., & Vasiliev, V. O., MMG 6, 1, 1 (2018), eprint 
% arXiv:1807.07323).
%
% Alpha,Beta,Delta,MB,sigma_Alpha,sigma_Beta,sigma_Delta,sigma_MB -
% hardcoded parameters of SNe. Their values are taken from the fitting
% process preceding this run. For a detailed description of them, see the 
% header of file 'sigma_sample_data.dat' in folder '1_SNIa_Fitter'.
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
% Outputs of SNIaDataProducer:
% SNIa_TimeSphere_Data.dat - ascii data file 'SNIa_TimeSphere_Data.dat' 
% containing SNe redshifts (Column #1), distance moduli according to the 
% best-fit Timesphere model (Column #2), and distance moduli errors (Column
% #3). The output file is an analog of 'SCPUnion2.1_mu_vs_z.txt' (available
% for download on the Union2.1 Supernova Cosmology Project website: 
% http://supernova.lbl.gov/Union/; see under title 'Union2.1 Compilation 
% Magnitude vs. Redshift Table (for your own cosmology fitter)') in which
% the distance moduli are calculated for the best-fit LCDM cosmological
% model.
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
% 

% Loading SNIa and sigma_sample data
Data=load('SNIa_data.dat');
sigma_sample_data=load('final_sigma_sample_data.dat');

% Setting the values of SNIa parameters and their uncertainties. These are
% obtained from the previous iterations of the fitting process.
Alpha=0.122;
Beta=2.417;
Delta=-0.181;
MB=-19.497;
sigma_Alpha=0.008;
sigma_Beta=0.066;
sigma_Delta=0.024;
sigma_MB=0.015;

% Setting the values of sigma_ext and h. These should be left unchanged
% throughout the iterative fitting process.
sigma_ext=0.05;
h=0.62339;

% Setting the value of the speed of light
c=299792458; % [in m/s]

% Calculating the Pearson correlation coefficients between various SNIa
% parameters (indicated in the parameter names). They are used in 
% calculating the covariant terms of the uncertainties of distance moduli. 
r_mB_x1=corrcoef(Data(:,2),Data(:,4));
r_mB_color=corrcoef(Data(:,2),Data(:,6));
r_x1_color=corrcoef(Data(:,4),Data(:,6));
r_vect=[r_mB_x1(1,2),r_mB_color(1,2),r_x1_color(1,2)];

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

% Calculating the distance moduli for all SNe. For details about this, see
% Suzuki et al. 2012.
mu_meas=mB+Alpha*x1-Beta*color+Delta*P-MB;

% Calculating the distance moduli for the reference cosmological model at
% all redshifts of the full SNIa sample.
mu_th=5*log10(c*(1+z).*sin(log(1+z)))-5*log10(h);

% Extracting the Pearson correlation coefficients from the input vector
r_mB_x1=r_vect(1);
r_mB_color=r_vect(2);
r_x1_color=r_vect(3);

% Calculating the uncertainties of all mu_meas values. For details about
% this, see Suzuki et al. 2012 and references therein.
sigma_lc_sq_Ref=sigma_ext^2+sigma_sample.^2+sigma_MB^2+(P.^2)*sigma_Delta^2+(color.^2)*sigma_Beta^2+(x1.^2)*sigma_Alpha^2+sigma_mB.^2+(Alpha^2)*sigma_x1.^2+(Beta^2)*sigma_color.^2+2*Alpha*r_mB_x1*sigma_mB.*sigma_x1-2*Beta*r_mB_color*sigma_mB.*sigma_color-2*Alpha*Beta*r_x1_color*sigma_x1.*sigma_color;

% Creating matrix of output data. All SNe are sorted by their redshift in
% an ascending order.
OutData(:,1)=z;
OutData(:,2)=mu_meas;
OutData(:,3)=sqrt(sigma_lc_sq_Ref);
OutData=sortrows(OutData,1);

% Saving output file
save('SNIa_TimeSphere_Data.dat','OutData','-ascii','-double','-tabs');

