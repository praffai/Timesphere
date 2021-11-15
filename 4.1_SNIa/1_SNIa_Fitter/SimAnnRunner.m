function SimAnnRunner(TrialNum)
% This function prepares running the simulated annealing based fitting
% process 'TrialNum' number of times. SimAnnRunner loads all data stored in
% external files and needed for the process, and provides them as inputs 
% for the SimAnn function. This way data files are loaded only once instead 
% of 'TrialNum' number of times, which significantly reduces the running 
% time. SimAnnRunner also runs 'DataFilter', a data filtering process
% applied before the simulated annealing, which filters out all SNIe that
% deviate from a reference model by more than 3 sigma (a process named as
% 'sigma clipping'). 
%
%
% Inputs for SimAnnRunner:
%
% TrialNum - the number of times the fitting process is run. The value of
% it, in the current format of the code, should be an integer times 100.
%
% FilterSwitch - a hardcoded parameter which, if set to 1, terminates the 
% run if DataFilter filters out the exact same SNe that were filtered out 
% in the last fitting. In such a case, the row indices of these SNe should
% be available from the last run, saved in ascii file 'outfiltered.txt'. 
% FilterSwitch should be set to any value other than 1 when the fitting 
% process is run for the first time.
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
% the interval.
%
% sigma_sample_data.dat - an external ascii file defining the sigma_sample 
% values for the various SNIa survey samples, applied both in the data
% filtering and in the fitting process (for details, see Suzuki et al., ApJ 
% 746, 1, 85, 24 (2012); e-print: https://arxiv.org/abs/1105.3470). The 
% file contains the following columns:
% Column #1 - an integer number from 1-19 indicating which survey sample 
% the SNIa data was taken from. The names of the surveys are listed in a 
% reversed order in Figure 3 of Suzuki et al. 2012.
% Column #2 - values of sigma_sample for the different surveys. The
% starting value is sigma_sample=0.15 (the median of sigma_sample values in 
% Suzuki et al. 2012) for all surveys. The values are then recalculated
% after each fit in an iterative process until the fitting process 
% converges to a stable set of parameters and sigma_sample values.
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
% column gives the measure of goodness of fit (in our case, the weighted 
% sum of squared deviations between data and the reference model, a.k.a. 
% the chi-squared statistic) corresponding to the parameter values in the 
% previous columns. The output file is saved after every 100th trial under 
% the name 'Optimal_Params_[i/100].txt', and the previous output file is 
% deleted.
% - SimAnnRunner.m also saves the row indices of SNe filtered out by 
% DataFilter.m in ascii file 'outfiltered.txt'. If no SNIa is filtered out,
% then 'outfiltered.txt' will be an empty file.
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the value of FilterSwitch. See the detailed description of the
% parameter in the header of this file.
FilterSwitch=0;

% Loading input files. See the detailed description of them in the header 
% of this file.
InputParams=load('parameter_limits.dat');
sigma_sample_data=load('sigma_sample_data.dat');
RawData=load('SNIa_data.dat');

% We make sure that the SNIa data table is sorted by redshifts in an 
% ascending order.
RawData=sortrows(RawData,1);

% Calculating the Pearson correlation coefficients between various SNIa
% parameters (indicated in the parameter names). The coefficients are
% passed on as input parameters for DataFilter.m and SimAnn.m. They are
% used in calculating the covariant terms of the uncertainties of distance
% moduli in the fitting process.
r_mB_x1=corrcoef(RawData(:,2),RawData(:,4));
r_mB_color=corrcoef(RawData(:,2),RawData(:,6));
r_x1_color=corrcoef(RawData(:,4),RawData(:,6));
r_vect=[r_mB_x1(1,2),r_mB_color(1,2),r_x1_color(1,2)];

% Filtering out (a.k.a. applying sigma clipping on) SNe for which the 
% distance moduli deviate from the reference model by more than 3 sigma.
% Output parameter 'Data' contains SNIa data that was left untouched by the
% filtering process. Output parameter 'Inds' contains the row indices of
% SNe that have been filtered out. If no SNIa is filtered out, then
% Data=RawData, and Inds is an empty vector.
[Inds,Data]=DataFilter(RawData,r_vect,sigma_sample_data);

% If FilterSwitch is set to 1, we load 'outfiltered.txt' from the previous
% run, and check if the indices of the outfiltered SNe from the previous
% run match the ones DataFilter.m filtered out in this run. The two being
% identical means that the iterative fit has converged to stable
% parameters, and thus the run is terminated.
if(FilterSwitch==1)
    Inds_Ref=load('outfiltered.txt');
    Inds_Ref=Inds_Ref(:);
    Inds=Inds(:);
    if(isequal(Inds,Inds_Ref))
        error('Iteration converged.');
    end
end

% If the run is not terminated in the previous step, we save the row 
% indices of the SNe filtered out by DataFilter.m in ascii file 
% 'outfiltered.txt'.
save('outfiltered.txt','Inds','-ascii');

% We start measuring the running time here. Whenever command 'toc' is used, 
% the running time until that point is output to the screen.
tic 

% This is the main cycle of the program, that we run 'TrialNum' number of
% times.
for i=1:TrialNum
    
    % We run the simulated annealing defined in function SimAnn.m. Each
    % trial results with a row of optimal parameters, and a corresponding
    % measure of goodness of fit (in our case, the weighted sum of squared 
    % deviations between data and the reference model, a.k.a. the 
    % chi-squared statistic)
    OutMatrix(i,:)=SimAnn(InputParams,Data,r_vect,sigma_sample_data);
    
    % We save the output matrix 'OutMatrix' after every 100th trial in file 
    % 'Optimal_Params_[i/100].txt', and delete the previous output file 
    % (if there is one).
    if((i/100)==floor(i/100))
        [i,100*i/TrialNum,toc] % checking the status of the run on the screen
        FileName=sprintf('Optimal_Params_%i.txt',i/100);
        save(FileName,'OutMatrix','-ascii');
        if((i/100)>1)
            DelFileName=sprintf('Optimal_Params_%i.txt',(i/100)-1);
            delete(DelFileName);
        end
    end
    
end

