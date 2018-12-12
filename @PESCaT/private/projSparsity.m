function projSparsity(obj,varargin)
% PROJSPARSITY performes the sparsity projection while solving the PESCAT
% reconstruction problem. The method used for determining the
% regularization parameter is read from PESCAT object properties.
%
% see also: PESCAT
%
%  (c) ICON Lab 2018

X = ifft2c(obj.trecon);
W = Wavelet('Daubechies',4,4);
switch obj.simParams.sparsityType
    case 'PES' %joint sparsity by projecting in separate subbands from each level
        X = reshape(X,[size(X,1),size(X,2),size(X,3)*size(X,4)]);
        obj.PES_L1(X,W);

    case 'SURE' %joint sparsity using SURE criterion
        X = reshape(X,[size(X,1),size(X,2),size(X,3)*size(X,4)]);
        obj.SURE_L1(X,W);
end
end
