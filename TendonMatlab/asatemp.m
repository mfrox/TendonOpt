function temperature = asatemp(optimValues,options)
%SATEMPERATUREFCNTEMPLATE Template to write custom temperature function
%   TEMPERATURE = SATEMPERATUREFCNTEMPLATE(optimValues,options) updates the
%   current temperature
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter 'k'
%
%   OPTIONS: options object created by using OPTIMOPTIONS.

%   Copyright 2006-2015 The MathWorks, Inc.

% The body of this function should simply return a new temperature vector 
% The new temperature will generally depend on either
% options.InitialTemperature or optimValues.temperature and on
% optimValues.k 
% temperature = options.InitialTemperature./log(optimValues.k);
tempRatioScale = 1.0e-5;
tempAnnealScale = 100.0;
currentx = optimValues.x;
nvar = numel(currentx);
m = -log(tempRatioScale);
n = log(tempAnnealScale);
c =  m.*exp(-n ./ nvar);


% Make sure temperature is of same length as optimValues.k (nvars)
temperature = options.InitialTemperature.*exp(-c.*optimValues.k.^(1/nvar));