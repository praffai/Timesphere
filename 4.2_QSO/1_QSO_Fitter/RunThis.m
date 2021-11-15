function RunThis
% This is the main code you should run for fitting the QSO dataset from
% Lusso, E., et al. A&A 642, A150, 24 (2020); e-print: 
% https://arxiv.org/abs/2008.08586, to the reference Timesphere cosmology. 
% RunThis.m runs a simulated annealing based fitting process 'TrialNum' 
% number of times, each resulting with a set of fitted parameters given in 
% a row of output file 'Optimal_Params_[TrialNum/100].txt'. The last column
% in the output file gives the measure of goodness of fit (in our case, the 
% weighted sum of squared deviations between data and the reference model, 
% a.k.a. the chi-squared statistic) corresponding to the fitted parameter 
% values in the first columns.
%
% The hierarchy of functions in the code package is the following:
% (1) RunThis.m calls SimAnnRunner.m
% (2) SimAnnRunner.m needs 
% - parameter_limits.dat, and
% - QSO_data.dat (and potentially outfiltered.txt) 
% as input data files in the same directory. 
% (3) SimAnnRunner.m calls DataFilter.m and SimAnn.m
% (4) SimAnn.m calls TimeSphereFitter.m
%
% Inputs for RunThis:
% TrialNum - a hardcoded parameter defining the number of times the fitting
% process is run. The value of it, in the current format of the code, 
% should be an integer times 100.
%
% Outputs of RunThis:
% None. The output is produced by SimAnnRunner.m in the format of an ascii
% file ('Optimal_Params_[TrialNum/100].txt').
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the number of trials in the fitting process
TrialNum=200;

% Running the simulated annealing based fitting process 'TrialNum' number
% of times
SimAnnRunner(TrialNum); 
