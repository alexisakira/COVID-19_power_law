%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discreteLogistic_ML
% (c) 2020 Alexis Akira Toda
% 
% Purpose: 
%       Estimate parameters of truncated discrete logistic distribution by
%       maximum likelihood
%
% Usage:
%       [param,PMF] = discreteLogistic_ML(counts,trunc)
%
% Inputs:
% counts    - counts data
%
% Optional:
% trunc     - 0 (truncate at t=0, default) or 1 (truncate at t=1)
%
% Output:
% param     - parmeters of discrete logistic distribution
% PMF       - probability mass function
%
% Version 1.0: April 25, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [param,PMF] = discreteLogistic_ML(counts,trunc)

if nargin < 2
    trunc = 0; % default is no truncation
end

if (any(counts < 0))||(any(rem(counts,1) ~= 0))
    error('counts must be a vector of nonnegative integers')
end

if size(counts,1) > size(counts,2)
    counts = counts'; % convert to row vector
end

% construct initial guesses
phi0 = 1/2;
p0 = sum(counts)/sum([1:length(counts)].*counts); % implied value from geometric distribution
param0 = [phi0 p0];

func = @(theta)discreteLogistic_obj(theta,counts,trunc);
param = fmincon(func,param0,[],[],[],[],zeros(1,2),ones(1,2)); % ML estimator

phi = param(1);
p = param(2);

if trunc == 0
    Tmax = length(counts) - 1; % maximum age
    F = (1-phi)./(1-phi + phi*(1-p).^[0:Tmax]); % CDF
    PMF = [1-phi diff(F)]; % fitted probability mass function
else % trunc = 1;
    Tmax = length(counts);
    F = (1-phi)./(1-phi + phi*(1-p).^[0:Tmax]); % CDF
    PMF = diff(F)/phi; % fitted probability mass function
end

end

