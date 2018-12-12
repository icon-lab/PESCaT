function parseInputs(obj,varargin)
% PARSEINPUTS parses the required and optional inputs to the PESCaT
% and replaces the remaining parameters with the default values. 
%
%  (c) ICON Lab 2018

validType = {'PES','SURE','flatPES','STD'};
checkType = @(x) any(validatestring(x,validType));
validDataset = {'bSSFP','T1','TOF','phantom'};
checkDataset = @(x) any(validatestring(x,validDataset));
p = inputParser;
p.KeepUnmatched = true;

addRequired(p,'data');                              % Input k-space data
addRequired(p,'sampling',@isstruct);                % Structure containg sampling mask
addParameter(p,'dataset','bSSFP',checkDataset);     % Dataset type, default: bSSFP
addParameter(p,'TVType','PES',checkType);           % Method used for TV projection
addParameter(p,'sparsityType','PES',checkType);     % Method used for sparsity projection
addParameter(p,'TykWeight',0.05,@isnumeric);        % Tykhonov regularization parameter used for kernel estimation
addParameter(p,'slcs',[]);                          % Slices to reconstruct - selects midslice if empty
addParameter(p,'sparsityLambda',[]);                % Would use as sparsity regularization parameter if not empty
addParameter(p,'TVLambda',[]);                      % Would use as TV regularization parameter if not empty
addParameter(p,'verbose',0,@isnumeric);             % Print out reconstruction info
addParameter(p,'den',5,@isnumeric);                 % Used for cropping calibration area
addParameter(p,'TVNit',5,@isnumeric);               % Number of iterations in TV projection using SURE
addParameter(p,'iterLim',50,@isnumeric);            % Maximum number of iterations
addParameter(p,'L1beta',0.2,@isnumeric);            % Scaling parameter of L1 epigraph
addParameter(p,'TVbeta',0.3,@isnumeric);            % Scaling parameter of TV epigraph
addParameter(p,'lambda',0.0001,@isnumeric);         % L2 regularization parameter for LSQR 
addParameter(p,'timeLim',[],@isnumeric);            % Limit of the reconstruction time
addParameter(p,'thresh',20,@isnumeric);             % The threshold used in stopping criterion
addParameter(p,'measurments',0,@isnumeric);         % Measure the separate cost terms and the intermediate images
parse(p,varargin{:})
pResults = p.Results;
pResults.TVType = validatestring(pResults.TVType,validType);
pResults.sparsityType = validatestring(pResults.sparsityType,validType);
obj.data = pResults.data;
obj.sampling = pResults.sampling;
pResults = rmfield(pResults,{'data','sampling'});
obj.simParams = pResults;

if isempty(obj.simParams.slcs)
    slcs = floor(length([1:size(obj.data,3)])/2)+1;
else
    slcs = obj.simParams.slcs;
end

if length(pResults.iterLim) == 1
    obj.simParams.iterLim = repmat(pResults.iterLim,[1,length(slcs)]);
end
obj.simParams.slcs = slcs;

[sx, sy, ~, ~, ~] = size(obj.data);

ssx = 2^ceil(log2(sx));
ssy = 2^ceil(log2(sy));
obj.simParams.ss = max(ssx, ssy);