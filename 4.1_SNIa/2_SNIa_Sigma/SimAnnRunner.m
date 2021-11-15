function OutFiltered=SimAnnRunner(TrialNum,Sample)
% This function prepares running the simulated annealing based fitting
% process 'TrialNum' number of times on data taken from survey sample 
% number 'Sample'. SimAnnRunner loads all data stored in external files and
% needed for the process, and provides them as inputs for the SimAnn 
% function. This way data files are loaded only once instead of 'TrialNum' 
% number of times, which significantly reduces the running time. 
% SimAnnRunner also runs 'DataFilter', a data filtering process applied 
% before the simulated annealing, which filters out all SNIa from the 
% survey sample that deviate from a reference model by more than 3 sigma (a
% process named as 'sigma clipping'). The number of SNe filtered out are
% output in integer number 'OutFiltered'.
%
%
% Inputs for SimAnnRunner:
%
% TrialNum - the number of times the fitting process is run. The value of
% it, in the current format of the code, should be an integer times 100.
%
% Sample - an integer number between 1-19, which is the identification 
% number of the survey sample processed. The names of the surveys are 
% listed in a reversed order in Figure 3 of in Suzuki et al., ApJ 746, 1, 
% 85, 24 (2012); e-print: https://arxiv.org/abs/1105.3470.
%
% Alpha,Beta,Delta,MB,sigma_Alpha,sigma_Beta,sigma_Delta,sigma_MB -
% hardcoded parameters of SNe. Their values are taken from the fitting
% process preceding this run that aims to calculate the proper sigma_sample
% values. For a detailed description of them, see the header of file
% 'sigma_sample_data.dat' in folder '1_SNIa_Fitter'.
%
% parameter_limits.dat - an external ascii file defining the boundaries of
% the parameter space to be explored in the fitting process, as well as the
% step length applied by simulated annealing along the different dimensions
% of the parameter space. Each row in the file corresponds to one parameter
% to be fitted; the number of rows define the number of parameters to be 
% fitted (or equivalently, the number of dimensions of the parameter 
% space). The columns in the file correspond to the lower (Column #1) and 
% the upper boundary (Column #2) of a parameter interval to be explored, 
% and the step length (Column #3) the simulated annealing applies within 
% the interval. As in this process, only sigma_sample is fit,
% parameter_limits.dat should only contain one single row.
%
% sigma_sample_data.txt - an external ascii file defining the sigma_sample 
% values for the various SNIa survey samples, applied both in the data
% filtering and in the fitting process (for details, see Suzuki et al., ApJ 
% 746, 1, 85, 24 (2012); e-print: https://arxiv.org/abs/1105.3470). The 
% file contains the following columns:
% Column #1 - an integer number from 1-19 indicating which survey sample 
% the SNIa data was taken from. The names of the surveys are listed in a 
% reversed order in Figure 3 of Suzuki et al. 2012.
% Column #2 - values of sigma_sample for the different surveys. 
% The file may potentially contain a third column, which gives the number 
% of SNe filtered out from the different survey samples by DataFilter.
%
% SNIa_data.dat - ascii data file produced by combining data in 
% 'SCPUnion2.1_AllSNe.tex' and 'SCPUnion2.1_mu_vs_z.txt', both of them 
% downloaded from the Union2.1 Supernova Cosmology Project website: 
% http://supernova.lbl.gov/Union/ (see under titles 'Full Table of All SNe'
% and 'Union2.1 Compilation Magnitude vs. Redshift Table (for your own 
% cosmology fitter)'). The data was originally used in Suzuki et al.
% (2012).
% 
%
% Outputs of SimAnnRunner:
% - The main output is produced in the format of an ascii text file 
% ('Optimal_Params_[TrialNum/100].txt'). Each row of the output file 
% contains the values of the parameters fitted in a single trial. The last 
% column gives the measure of goodness of fit corresponding to the 
% parameter values in the previous columns. The measure of goodness of fit 
% in this case is the absolute deviation of the reduced chi statistic (the 
% weighted sum of absolute deviations between data and the reference model,
% normalized by the degrees of freedom of the fit, which in our case is 
% equal to the total number of data points in the sample) from 1. The 
% output file is saved after every 100th trial under the name 
% 'Optimal_Params_[i/100].txt', and the previous output file is deleted.
% - SimAnnRunner also outputs the number of SNe outfiltered by DataFilter
% from the survey sample. This number is stored under the name 
% 'OutFiltered'.
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the values of SNe parameters using the results of the preceding
% fitting process run.
Alpha=0.123;
Beta=2.433;
Delta=-0.172;
MB=-19.502;
sigma_Alpha=0.008;
sigma_Beta=0.057;
sigma_Delta=0.026;
sigma_MB=0.014;

% Loading input files. See the detailed description of them in the header 
% of this file.
InputParams=load('parameter_limits.dat');
sigma_sample_data=load('sigma_sample_data.txt');
RawData=load('SNIa_data.dat');

% We make sure that the SNIa data table is sorted by redshifts in an 
% ascending order.
RawData=sortrows(RawData,1);

% Extracting the sigma_sample value of the processed survey sample (having
% identification number 'Sample') from input data matrix
% 'sigma_sample_data'.
sigma_sample=sigma_sample_data(Sample,2);

% Organizing all parameters needed by DataFilter and SimAnn into one row
% vector.
FixParams=[Alpha,Beta,Delta,MB,sigma_Alpha,sigma_Beta,sigma_Delta,sigma_MB,sigma_sample];

% Calculating the Pearson correlation coefficients between various SNIa
% parameters (indicated in the parameter names). The coefficients are
% passed on as input parameters for DataFilter.m and SimAnn.m. They are
% used in calculating the covariant terms of the uncertainties of distance
% moduli in the fitting process. Note that these coefficients are
% calculated using the full set of SNIa data, and not for individual survey
% samples.
r_mB_x1=corrcoef(RawData(:,2),RawData(:,4));
r_mB_color=corrcoef(RawData(:,2),RawData(:,6));
r_x1_color=corrcoef(RawData(:,4),RawData(:,6));
r_vect=[r_mB_x1(1,2),r_mB_color(1,2),r_x1_color(1,2)];

% Extracting only the SNIa survey sample data that we process. All other
% SNIa data are discarded.
SampleData=RawData(find(RawData(:,9)==Sample),:);
clear RawData;

% Filtering out (a.k.a. applying sigma clipping on) SNe for which the 
% distance moduli deviate from the reference model by more than 3 sigma.
% Output parameter 'Data' contains SNIa data that was left untouched by the
% filtering process. Output parameter 'Inds' contains the row indices of
% SNe that have been filtered out. If no SNIa is filtered out, then
% Data=RawData, and Inds is an empty vector with a length of zero. The
% length of Inds, i.e. the number of SNe filtered out, is stored and ouput
% in integer variable 'OutFiltered'.
[Inds,Data]=DataFilter(SampleData,r_vect,FixParams);
OutFiltered=length(Inds);

% This is the main cycle of the program, that we run 'TrialNum' number of
% times.
for i=1:TrialNum
    
    % We run the simulated annealing defined in function SimAnn.m. Each
    % trial results with a row of optimal parameters, and a corresponding
    % measure of goodness of fit, which in this case, is the absolute 
    % deviation of the reduced chi statistic (the weighted sum of absolute 
    % deviations between data and the reference model, normalized by the 
    % degrees of freedom of the fit, which in our case is equal to the 
    % total number of data points in the sample) from 1.
    OutMatrix(i,:)=SimAnn(InputParams,Data,r_vect,FixParams);

    % We save the output matrix 'OutMatrix' after every 100th trial in file 
    % 'Optimal_Params_[i/100].txt', and delete the previous output file 
    % (if there is one).    
    if((i/100)==floor(i/100))
        FileName=sprintf('Optimal_Params_%i.dat',i/100);
        save(FileName,'OutMatrix','-ascii');
        if((i/100)>1)
            DelFileName=sprintf('Optimal_Params_%i.dat',(i/100)-1);
            delete(DelFileName);
        end
    end
    
end

