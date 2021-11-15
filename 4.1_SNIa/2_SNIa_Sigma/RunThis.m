function RunThis
% This is the main code for the modified version of the code package
% provided in folder '1_SNIa_Fitter'. The purpose of this version is to
% calculate the sigma_sample values after each fit that properly normalize
% the measure of goodness of fit in the fitting process. A description of
% how sigma_sample values for the different SNIa surveys should be
% calculated can be found in Suzuki et al., ApJ 746, 1, 85, 24 (2012); 
% e-print: https://arxiv.org/abs/1105.3470. We apply the same simulated
% annealing based fitting process to find the values of sigma_sample-s, but
% here the fitting process is run 19 times, where the data is restricted to
% data within an individual SNIa survey (all other data are filtered out in
% SimAnnRunner.m), and the measure of goodnees of fit is now the absolute
% deviation of the reduced chi statistic (the weighted sum of absolute 
% deviations between data and the reference model, normalized by the 
% degrees of freedom of the fit, which in our case is equal to the total 
% number of data points in the sample) from 1. By minimizing this measure
% of goodness of fit we make sure that the sigma_sample of a survey
% normalizes the reduced chi statistic to 1 (a requirement set by Suzuki et
% al. 2012).
%
% The hierarchy of functions in the code package is the following:
% (1) RunThis.m calls SimAnnRunner.m
% (2) SimAnnRunner.m needs 
% - parameter_limits.dat, 
% - sigma_sample_data.txt, and 
% - SNIa_data.dat 
% as input data files in the same directory. 
% (3) SimAnnRunner.m calls DataFilter.m and SimAnn.m
% (4) SimAnn.m calls TimeSphereFitter.m
%
% Inputs for RunThis:
% TrialNum - a hardcoded parameter defining the number of times the fitting
% process is run. The value of it, in the current format of the code, 
% should be an integer times 100. Due to sigma_sample being the only
% parameter fitted in the process, 'TrialNum' can be much less than what we
% use in '1_SNIa_Fitter' (e.g. 300 versus 10000).
%
% SampleNum - a hardcoded integer parameter defining the number of SNIa
% survey samples. For data taken from Suzuki et al. 2012, SampleNum=19; the
% names of the surveys are listed in a reversed order in Figure 3 of Suzuki
% et al. 2012.
%
%
% Output of RunThis:
% sigma_sample_data.txt - an ascii file that stores the calculated 
% sigma_sample values for the various SNIa survey samples (for details 
% about what sigma_sample values are, see Suzuki et al. 2012). The file 
% contains the following columns:
% Column #1 - an integer number from 1-19 indicating which survey sample 
% the SNIa data was taken from. The names of the surveys are listed in a 
% reversed order in Figure 3 of Suzuki et al. 2012.
% Column #2 - calculated values of sigma_sample for the different surveys. 
% Column #3 - the number of SNe filtered out from the different survey 
% samples by DataFilter. A data point is filtered out if it deviates from a
% reference model by more than 3 sigma (this process is commonly used under
% the name 'sigma clipping'). Note that 'sigma_sample_data.txt' from a
% previous run is also used in the process as an input (for DataFilter.m).
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the number of trials in the fitting process
TrialNum=300;

% Setting the number of SNIa survey samples
SampleNum=19;

% Setting the name of the file output by the simulated annealing based
% fitting process
OutFileName=sprintf('Optimal_Params_%i.dat',TrialNum/100);

% This is the main cycle of the code
for i=1:SampleNum
    
    % Setting the identification number of the survey sample processed in a
    % cycle. The number runs from 1 to SampleNum=19, and are stored as
    % elements of vector 'Sample'.
    Sample(i)=i;
    
    % Running the simulated annealing-based fitting process 'TrialNum' 
    % number of times on a SNIa survey sample identified by number 
    % Sample(i). Output vector OutFiltered stores the number of SNIa
    % datapoints filtered out by DataFilter from the different survey
    % samples.
    OutFiltered(i)=SimAnnRunner(TrialNum,Sample(i));
    
    % Loading the SimAnnRunner output file, and extracting the best-fit
    % sigma_sample value from it. The sigma_sample value is then stored as
    % the element of vector 'Sigma_Sample'.
    FitResults=load(OutFileName);
    [MinVal,MinInd]=min(FitResults(:,2));
    Sigma_Sample(i)=FitResults(MinInd,1);
    
    % Deleting the output file of SimAnnRunner, as it is no longer needed.
    delete(OutFileName);
    
    % Displaying the survey number and the corresponding sigma_sample value
    % on the screen, so that the user can follow where the process is at.
    [i,Sigma_Sample(i)]
    
end

% Creating the output matrix to be stored in output file 
% 'sigma_sample_data.txt'. For the detailed description of the matrix
% columns, see the header of this file.
OutMatrix(:,1)=Sample;
OutMatrix(:,2)=Sigma_Sample;
OutMatrix(:,3)=OutFiltered;

% Saving the output file in an ascii format. Elements in the ascii file are
% saved as tab-separated double-precision floating-point numbers.
save('sigma_sample_data.txt','OutMatrix','-ascii','-double','-tabs');
