%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discreteLogistic_obj
% (c) 2020 Alexis Akira Toda
% 
% Purpose: 
%       Objective function of maximum likelihood estimation of truncated
%       discrete logistic distribution
%
% Usage:
%       Out = discreteLogistic_obj(param,counts,trunc)
%
% Inputs:
% param     - parameter
% counts    - counts data
%
% Optional:
% trunc     - 0 (truncate at t=0, default) or 1 (truncate at t=1)
%
% Output:
% Out       - negative of log likelihood
%
% Version 1.0: April 25, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function Out = discreteLogistic_obj(param,counts,trunc)

if nargin < 3
    trunc = 0; % default is no truncation
end

phi = param(1);
p = param(2);

if size(counts,1) > size(counts,2)
    counts = counts'; % convert to row vector
end

if trunc == 0
    Tmax = length(counts) - 1; % maximum age
    F = (1-phi)./(1-phi + phi*(1-p).^[0:Tmax]); % CDF
    temp = [1-phi diff(F)]; % bin probabilities
    Out = -sum(counts.*log(temp)); % negative of log likelihood
else % trunc = 1
    Tmax = length(counts); % maximum age
    F = (1-phi)./(1-phi + phi*(1-p).^[0:Tmax]); % CDF
    temp = diff(F)/phi; % bin probabilities conditional on T>0
    Out = -sum(counts.*log(temp)); % negative of log likelihood
end

end

