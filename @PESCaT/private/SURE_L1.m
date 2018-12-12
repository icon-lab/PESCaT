function SURE_L1(obj,X,W,varargin)
% SURE_L1 performes the auto-tuning sparsity projection using SURE criterion
%
% see also: PESCAT, PROJSPARSITY
%
%  (c) ICON Lab 2018

[sx,sy,npc] = size(X);
X = zpad(X,obj.simParams.ss,obj.simParams.ss,size(X,3),size(X,4));

if obj.iter == 1
    for pidx=1:size(X,3)
        obj.SURE.noiseSigma(pidx) = find_sigma(X(:,:,pidx));
    end
    obj.SURE.f0 = zeros(size(X));
    initialRange = linspace(2e-5,2e-2,100);
    obj.SURE.nextRange = permute(repmat(initialRange,[size(X,3),1]),[1 3 2]);
    obj.SURE.msk = (obj.data==0);
    obj.SURE.cdata = obj.data.*(~obj.SURE.msk);
end
Xi = X;
X = W*(X); 
obj.SURE.fn = Xi;
obj.SURE.fne1 = obj.SURE.f0;
for pidx=1:size(obj.SURE.fn,3)
    [best_lambda(pidx,obj.iter),~,obj.SURE.nextRange(pidx,obj.iter+1,:)] = autoSparsity(obj.SURE.fn(:,:,pidx),...
        obj.SURE.fne1(:,:,pidx),...
        obj.iter,...
        obj.SURE.noiseSigma(pidx),...
        obj.SURE.nextRange(pidx,obj.iter,:),...
        obj.SURE.cdata(:,:,pidx),...
        ~obj.SURE.msk(:,:,pidx));
    if isempty(obj.simParams.sparsityLambda)
        [X(:,:,pidx)] = softThresh(X(:,:,pidx),best_lambda(pidx,obj.iter));
    else
        [X(:,:,pidx)] = softThresh(X(:,:,pidx),obj.simParams.sparsityLambda);
    end
    obj.SURE.f0(:,:,pidx) = obj.SURE.fn(:,:,pidx);
end
X = W'*X;
obj.trecon = fft2c(crop(X,sx,sy,npc));
end

