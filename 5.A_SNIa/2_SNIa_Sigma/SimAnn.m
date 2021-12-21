function Params=SimAnn(InputParams,Data,r_vect,FixParams)
% This function runs the simulated annealing based fitting process on input
% SNIa data stored in matrix 'Data'. The number of dimensions of parameter
% space is set by the number of rows in input matrix 'InputParams'. The
% boundaries of the parameter space along the different dimensions are 
% set in the first (lower boundary) and second (upper boundary) columns of 
% 'InputParams'. In the fitting process, we try to find the minimum of the 
% measure of goodness of fit within the boundaries of the parameter space. 
% The measure of goodness of fit in this case is the absolute deviation of 
% the reduced chi statistic (the weighted sum of absolute deviations 
% between data and the reference model, normalized by the degrees of 
% freedom of the fit, which in our case is equal to the total number of 
% data points in the sample) from 1. First, we start a wandering process at
% a randomly chosen allowed point of parameter space, whose coordinates are
% randomized from uniform distributions defined between the boundaries of 
% individual parameters. We then lower the value of a temperature parameter
% 'T' with each step in the wandering, from 'T_start' (which is set to be 
% the measure of goodness of fit at the starting point) to 
% T_crit=0.01*T_start by a factor of (1-T_step), where T_step=0.001. The 
% wandering finishes when T<=T_crit, and then the final set of parameters 
% (coordinates of the point where the wandering terminates) as well as 
% corresponding measure of goodness of fit is output in vector 'Params'. 
% Each step of the wandering starts with randomly choosing a parameter to 
% change (i.e. the axis which we step along), and randomizing a step length
% from a Gaussian distribution with zero mean and sigma equal to the value 
% set in the third column of 'InputParams'. If the step moves us out from 
% the interval of allowed parameters, we step toward the opposite direction
% along the axis of the same parameter (in the exceptional case when this 
% also moves us out from the boundaries - which means that the step length 
% in InputParams was set to be too large -, then we terminate the whole 
% fitting process, and display an error message on the screen). The step is
% accepted with a probability of 1 if the measure of goodness of fit is 
% smaller than it was before, and it is accepted with a probability of 
% p=exp(-Diff/T) if the measure is larger, where 'Diff' is the absolute 
% difference between the old and the new measure (one can see that as T 
% decreases with every step, the same Diff results with a decreasing 
% probability of the new position being accepted).
%
%
% Inputs for SimAnn:
% 
% InputParams - a matrix that defines the boundaries of the parameter space
% to be explored in the fitting process, as well as the step length applied
% by simulated annealing along the different dimensions of parameter space. 
% Each row in the matrix corresponds to one parameter to be fitted; the 
% number of rows define the number of parameters to be fitted (or 
% equivalently, the number of dimensions of the parameter space). The 
% columns in the matrix correspond to the lower (Column #1) and the upper 
% boundary (Column #2) of a parameter interval to be explored, and the step
% length (Column #3) the simulated annealing applies within the interval.
%
% Data - the matrix of SNe data loaded from 'SNIa_data.dat' and then 
% restricted only to data from a single SNIa survey, where each row
% corresponds to a single SNIa event. For a description of the matrix 
% columns, see the header of the ascii text file 'SNIa_data.dat'.
%
% r_vect - a vector of Pearson correlation coefficients between various 
% SNIa parameters. The coefficients are between mB and x1, mB and color,
% and x1 and color, respectively. They are used in calculating the 
% covariant terms of the uncertainties of measured SNIa distance moduli. 
% For a description of these quantities, see the header of the ascii text 
% file 'SNIa_data.dat'.
%
% FixParams - SNIa parameters Alpha, Beta, Delta, MB, sigma_Alpha, 
% sigma_Beta, sigma_Delta, sigma_MB, and sigma_sample, organized in a row 
% vector, and taken from the fitting process preceding this run. For a
% description of these parameters, see the header of 
% 'parameter_limits.dat', paper Suzuki et al. (2012), and the header of
% 'sigma_sample_data.dat' in folder '1_SNIa_Fitter'.
%
% T_step - a hardcoded parameter that sets how fast the temperature
% parameter of the simulated annealing is decreased with every step. By
% default it is set to T_step=0.001.
% 
% T_crit - a hardcoded parameter that sets the temperature at which the
% simulated annealing is stopped. By default it is set to 
% T_crit=0.01*T_start, where T_start is the starting temperature set to be
% equal to the measure of goodness of fit at the starting point of the
% wandering.
%
%
% Outputs of SimAnn:
% Params - a vector containing the best fit parameters obtained by the
% simulated annealing process. Th parameters in Params follow the same
% order as the one defined by the order of rows in InputParams. The last
% element of Params gives the measure of goodness of fit (in our case, the 
% weighted sum of squared deviations between data and the reference model, 
% a.k.a. the chi-squared statistic) corresponding to the best fit parameter
% values stored in the previous vector elements.
%
%
% Credits: 
% Peter Raffai, Gergely Dalya, Alexandra Karsai; Institute of Physics, 
% Eotvos Lorand University, H-1117 Budapest, Pazmany P. s. 1/A.
% All rights reserved. (2021)
% Contact: peter.raffai@ttk.elte.hu
%

% Setting the value of T_step
T_step=0.001;

% Checking the number of parameters to be fitted in the process (the number
% of rows in matrix InputParams), and storing the number in integer 
% NumParams
NumParams=size(InputParams,1);

% Setting seed of randomization to computer clock (uncomment and use one of
% the following three methods, depending on the version number of your
% Matlab)
% s = RandStream('mt19937ar','Seed','shuffle');
% RandStream.setGlobalStream(s);
% % (For earlier Matlab versions:) RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% % (For even more earlier Matlab versions:) 
rand('state',sum(100*clock));

% Randomizing the coordinates for the starting point in parameter space
for i=1:NumParams        
    Params(i)=InputParams(i,1)+(InputParams(i,2)-InputParams(i,1))*rand(1,1);        
end

% Calculating the measure of goodness of fit for the starting point of the
% wandering, and setting the starting temperature 'T_start' to be equal 
% with it.
OutStat=TimeSphereFitter(Params,Data,r_vect,FixParams);
T_start=OutStat;

% Initializing the starting temperature
T=T_start;

% Defining the critical temperature at or below which the simulated 
% annealing in stopped
T_crit=0.01*T_start;

% This is the main cycle of the code
while(T>T_crit)
      
    % Randomly choosing a parameter to change in the step
    ParamInd=randi(NumParams,1);        
    
    % Creating a temporary vector of new parameters
    Test_Params=Params;
    
    % Calculating the step length along the chosen parameter axis
    RandStepParam=InputParams(ParamInd,3)*randn(1,1);
    
    % Stepping along the chosen parameter axis 
    Test_Params(ParamInd)=Test_Params(ParamInd)+RandStepParam;   
    
    % Checking if we stayed within the parameter boundaries with the step.
    % If not, then we reverse the step to the opposite direction along the
    % same parameter axis.
    if((Test_Params(ParamInd)<InputParams(ParamInd,1))|(Test_Params(ParamInd)>InputParams(ParamInd,2)))            
        Test_Params=Params;            
        Test_Params(ParamInd)=Test_Params(ParamInd)-RandStepParam;        
    end 
    
    % If we stayed within the allowed parameter range with the step, we
    % calculate the measure of goodness of fit at the new position of
    % parameter space. Otherwise we terminate the process and display an
    % error message on the screen (in this case, the user should decrease
    % the step length along this parameter axis in file
    % 'parameter_limits.dat').
    if((Test_Params(ParamInd)>=InputParams(ParamInd,1))&(Test_Params(ParamInd)<=InputParams(ParamInd,2)))            
        Test_OutStat=TimeSphereFitter(Test_Params,Data,r_vect,FixParams);        
    else        
        error('Step exceeded parameter boundaries.');        
    end
        
    % We accept the new position in parameter space with a probability of 1 
    % if the measure of goodness of fit at the new position is smaller than 
    % it was before, and it is accepted with a probability of 
    % p=exp(-Diff/T) if the measure is larger, where 'Diff' is the absolute 
    % difference between the old and the new measure.
    if(Test_OutStat<OutStat)  
        OutStat=Test_OutStat;   
        Params=Test_Params;
    else        
        Diff=abs(Test_OutStat-OutStat);            
        p=exp(-Diff/T);                                   
        Test_P=rand(1,1);               
        if(Test_P<p)               
            OutStat=Test_OutStat;                
            Params=Test_Params;            
        end        
    end
       
    % Decreasing the temperature parameter after each step
    T=T*(1-T_step);            
    
end

% Adding the final measure of goodness of fit to the last element of the 
% output vector
Params(length(Params)+1)=OutStat;

