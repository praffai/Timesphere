function QSODataProducer
% This code produces ascii data file 'QSO_TimeSphere_Data.dat' containing 
% QSO redshifts (Column #1), distance moduli according to the best-fit
% Timesphere model (Column #2), and distance moduli errors (Column #3). The
% output file is an analog of 'QSO_data.dat' in which the distance moduli 
% and their errors are calculated for the best-fit LCDM cosmological model.
% Function QSODataProducer also provide the same QSO data in a binned form
% using logarithmic bins in ascii file 'QSO_Binned_TimeSphere_Data.dat',
% containing central redshifts of bins (Column #1), mean distance moduli of
% QSOs in bins (Column #2), errors on mean distance moduli (Column #3). The
% binned data is analog with the data visualized with red points in 
% Figure 7 of Bargiacchi, G., et al., A&A 649, A65, 10 (2021); e-print: 
% https://arxiv.org/abs/2101.08278.
%
%
% Inputs for QSODataProducer:
%
% QSO_data.dat - ascii data file produced by extracting Columns #4
% (redshifts), #12 (distance moduli) and #13 (errors on distance moduli)
% from 'table3.dat' made available under title 'Quasars as standard 
% candles. III.: J/A+A/642/A150' in VizieR 
% (URL: http://cdsarc.u-strasbg.fr/viz-bin/cat/J/A+A/642/A150), which is
% the data used by Lusso, E., et al., A&A 642, A150, 24 (2020); e-print: 
% https://arxiv.org/abs/2008.08586.
%
% mu_shift, sigma_shift - hardcoded parameters obtained from the fitting 
% process. Their values should be refreshed after every iterations of the 
% process.
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
% N - the number of bins to be used in binning output QSO data. The default
% value we set is N=18.
%
%
% Outputs of QSODataProducer:
%
% QSO_TimeSphere_Data.dat - ascii data file containing QSO redshifts 
% (Column #1), distance moduli according to the best-fit Timesphere model 
% (Column #2), and distance moduli errors (Column #3). The output file is 
% an analog of 'QSO_data.dat' in which the distance moduli and their errors
% are calculated for the best-fit LCDM cosmological model.
%
% QSO_Binned_TimeSphere_Data.dat - ascii data file containing binned QSO
% data as central redshifts of bins (Column #1), mean distance moduli of
% QSOs in bins (Column #2), errors on mean distance moduli (Column #3). The
% binned data is analog with the data visualized with red points in 
% Figure 7 of Bargiacchi, G., et al., A&A 649, A65, 10 (2021); e-print: 
% https://arxiv.org/abs/2101.08278.
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
% 

% Loading QSO data
Data=load('QSO_data.dat');

% Setting the best-fit value of mu_shift and its uncertainty. These are 
% obtained from the fitting process.
mu_shift=-0.149;
sigma_shift=0.008;

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
sigma_lc_sq=sigma_mu.^2+sigma_shift^2;

% Creating matrix of output data. All QSOs are sorted by their redshift in
% an ascending order.
OutData(:,1)=z;
OutData(:,2)=mu_meas;
OutData(:,3)=sqrt(sigma_lc_sq);
OutData=sortrows(OutData,1);

% Saving output QSO data file
save('QSO_TimeSphere_Data.dat','OutData','-ascii','-double','-tabs');

% Setting up logarithmic binning for the binned QSO data
N=18;
zmin=log(min(z));
zmax=log(max(z));
step=(zmax-zmin)/N;
logscale=zmin:step:zmax;

% Binning QSO data
for i=1:length(logscale)-1
    
    % Calculating central redshift value associated to a bin
    bin_z(i)=exp(zmin+(i-1)*step+step/2);
    
    % Selecting QSOs falling into the bin
    MinLimZ=exp(zmin+(i-1)*step);
    MaxLimZ=exp(zmin+i*step);
    Inds=find((z>=MinLimZ)&(z<=MaxLimZ));
    
    % Calculating mean distance modulus value in the bin
    bin_mu(i)=mean(mu(Inds));
    
    % Calculating error on the mean distance modulus value
    bin_sigma(i)=sqrt(sum(sigma_mu(Inds).^2))/sqrt(length(Inds));
    
end

% Deleting empty bins from output binned QSO data
DelInds=find(isnan(bin_mu));
bin_z(DelInds)=[];
bin_mu(DelInds)=[];
bin_sigma(DelInds)=[];

% Creating matrix of output binned data.
OutBinnedData(:,1)=bin_z;
OutBinnedData(:,2)=bin_mu;
OutBinnedData(:,3)=bin_sigma;

% Saving output binned QSO data file
save('QSO_Binned_TimeSphere_Data.dat','OutBinnedData','-ascii','-double','-tabs');

