%% Implementation of the fast ridge regression, specifically for the case when p > n
%
% function [beta, b0, tau2, DOF, lambda, score] = fastridge(X, y, 'param1', value1, ...)
%
% Parameters:
%   X           = regressor matrix [n x p]
%   y           = response vector  [n x 1]
%
% Optional arguments:
%
% If no path, lambda or dof is specified, fastridge will search for the
% optimal ridge model using the specified criterion (default is 'mmlu')
%
%   'path'      = generate a path of beta's, with the number of lambda's given by the argument  [ integer scalar ]
%   'criterion' = specify which model selection criterion to use [ 'mmlu', 'aicc', 'kicc' ]
%   'lambda'    = fit a ridge model with a specified lambda [ positive scalar ]
%   'dof'       = fit a ridge model with specified degrees-of-freedom [ postive scalar ]
%
% Returns:
%   beta        = regression parameters  [p x nlambda]
%   b0          = regression param. for constant [1 x 1]
%   tau2        = noise variance sigma^2 [1 x nlambda]
%   DOF         = degrees-of-freedom of fitted models [1 x nlambda]
%   score       = model selection criterion score of fitted models [1 x nlambda]
%  
%
% Example:
%  % Generate a regularisation path with 100 lambda's
%  [beta, b0, tau2, DOF, lambda, score] = fastridge(X, y, 'path', 1e2);
%
%  % Fit with a specified DOF
%  [beta, b0, tau2, DOF, lambda, score] = fastridge(X, y, 'dof', 88);
%
% See fastridge_example.m for a more extensive example with timings
%
% References:
%
% Efficient quadratic regularization for expression arrays
% Trevor Hastie and Robert Tibshirani
% Biostatistics, Vol. 5, No. 3, pp. 329--340, 2004
%
% MML Invariant Linear Regression
% Daniel F. Schmidt and Enes Makalic
% Lecture Notes in Computer Science, Vol. 5866, pp. 312--321, 2009
%
% Regression and Time Series Model Selection in Small Samples
% Clifford M. Hurvich and Chih-Ling Tsai
% Biometrika, Vol. 76, No. 2, pp. 297--307, 1989
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2015
%
function [beta, b0, tau2, DOF, lambda, score] = fastridge(X, y, varargin)

%% Standardise the input data
[n, p] = size(X);
[X, muX, normX, y, muY] = standardise(X, y);

% Determine the maximum DOF
maxDOF   = min(n-3/2, p);

%% Parse options
inParser = inputParser;  

%% Default parameter values
defaultLambda     = NaN;
defaultDOF        = NaN;
defaultCriterion  = 'mmlu';
defaultPath       = NaN;

expectedCriterion = {'mmlu', 'aicc'};

%% Define parameters
addParameter(inParser,'lambda',defaultLambda,@(x) isnumeric(x) && (x >= 0) && isscalar(x));
addParameter(inParser,'dof',defaultDOF,@(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= maxDOF);
addParameter(inParser,'criterion',defaultCriterion,@(x) any(validatestring(x,expectedCriterion)));
addParameter(inParser,'path',defaultPath,@(x) isnumeric(x) && isscalar(x) && x > 0);

%% Parse options
parse(inParser, varargin{:});  

% General options for all samplers
DOF       = inParser.Results.dof; 
lambda    = inParser.Results.lambda;
path      = inParser.Results.path;
criterion = inParser.Results.criterion;

%% Compute the SVD and sufficient statistics
[u,s,v] = svd(X, 'econ');

ds  = diag(s);
r   = bsxfun(@times, u, ds'); % r = u*s;
rtr = ds.^2;                  % r'*r
ry  = r'*y;
yty = y'*y;

%% Error checking
if (~isnan(DOF) && ~isnan(lambda))
    error('Cannot specify both DOF and lambda');
end
if ((~isnan(DOF) || ~isnan(lambda)) && ~isnan(path))
    error('Cannot specify a lambda/DOF and request a path');
end

%% Initialise
score = NaN;

%% If a DOF was specified, convert it to a lambda
if (~isnan(DOF))
    lambda = DOF2lambda(DOF, rtr);
end

%% If no DOF/lambda, or path option was specified, search for the best model using the specified criterion
if (isnan(DOF) && isnan(lambda) && isnan(path))
    minLambda    = DOF2lambda(maxDOF, rtr, 0);
    maxLambda    = DOF2lambda(1e-04, rtr, 0);

    % Search for the smallest model by criterion score
    [lambda, score] = fmincon( @(l) fitfastridge_wrapper(exp(l), X, y, yty, v, rtr, ry, criterion), minLambda*1.05, [], [], [], [], log(minLambda), log(maxLambda), [], optimoptions('fmincon','Display','off'));
    lambda = exp(lambda);
end

%% If a lambda was specified, fit the specified model
if (~isnan(lambda))
    [beta, tau2, DOF] = fitfastridge(X, y, yty, v, rtr, ry, lambda);
end

%% If a path was requested
if (~isnan(path))
    % Use a linear grid over the DOF's
    DOF   = linspace(0, maxDOF, path);
    score = zeros(1, path);
    beta     = zeros(p, path);
    tau2  = zeros(1, path);

    % Try each model in the lambda grid
    for i = 1:path,
        L = DOF2lambda(DOF(i), rtr, 0);
        lambda(i) = L;

        % Fit the model
        [beta(:,i), tau2(i), DOF(i), mmlu, aicc] = fitfastridge(X, y, yty, v, rtr, ry, lambda(i));
        switch (criterion)
            case 'mmlu'
                score(i) = mmlu;
            case 'aicc'
                score(i) = aicc;
        end
    end
end

%% Re-scale coefficients
beta = bsxfun(@rdivide, beta, normX');
b0 = muY-muX*beta;

end


%% Convert degrees-of-freedom to a lambda regularisation parameter
function [lambda, f] = DOF2lambda(dof, rtr, lambda)

if (~exist('lambda','var'))
    lambda = 0;
end
f = inf;

% Newton-Raphson to find lambda
while (abs(f-dof)>1e-2)
    lambda = lambda + (sum( rtr ./ (rtr + lambda) ) - dof) / sum(rtr./(rtr + lambda).^2);
    f = sum(rtr ./ (rtr + lambda));
end

end


%% Wrapper function for fast ridge fitting
function f = fitfastridge_wrapper(lambda, X, y, yty, v, rtr, ry, criterion)

[~, ~, ~, mmlu, aicc] = fitfastridge(X, y, yty, v, rtr, ry, lambda);
switch lower(criterion)
    case 'mmlu'
        f = mmlu;    
    case 'aicc'
        f = aicc;
end
    
end


%% Fit a ridge model in a fast manner
function [b, tau2, dof, mmlu, aicc] = fitfastridge(X, y, yty, v, rtr, ry, lambda)

n = length(y);
b = bsxfun(@times, v, 1./(rtr + lambda)') * ry;
dof = sum(rtr ./ (rtr + lambda));

e = (y - X*b);
tau2 = e'*e/n;

%% Model selection criteria
mmlu = (n-dof)/2*log(2*pi) + (n-dof)/2*log(tau2) + (n-dof)/2 + dof/2*log(pi*yty) - gammaln(dof/2+1) + 1/2*log(dof+1);
aicc = n/2*log(tau2) + n/2 + dof*n/(n-dof-1);

end


%% Standardise the covariates to have zero mean and x_i'x_i = 1
% function [X,meanX,stdX,y,meany]=standardise(X,y)
function [X,meanX,stdX,y,meany]=standardise(X,y)

%% params
n=size(X,1);
meanX=mean(X);
stdX=std(X,0,1);

%% Standardise Xs
X=bsxfun(@minus,X,meanX);
X=bsxfun(@rdivide,X,stdX);

%% Standardise ys (if neccessary)
if(nargin == 2)
    meany=mean(y);
    y=y-meany;
end;

end