function optimize(obj,varargin)
% OPTIMIZE solves the optimization problem of PESCaT using a POCS approach
% 
% see also: PESCAT/RECONPESCAT
%  (c) ICON Lab 2018


[sx, sy, sl, N, D] = size(obj.data);

data = obj.data;

p_acq = 4; % p value for p-norm combination over phase-cycles
p_coils = 2; % p value for p-norm combination over coils

for sidx = 1:sl
    fprintf('reconstructing slice %d...',sidx)
    
    datas = reshape(data(:,:,sidx,:,:),[sx,sy,N,D]); %slice data
    kernel = obj.kernel(end:-1:1,end:-1:1,:,:,:,:,sidx);
    
    imkernel = zeros(sx, sy, N, D, N, D);
    for n=1:N
        for d=1:D
            imkernel(:,:,:,:,n,d) = ifft2c(zpad(kernel(:,:,:,:,n,d)*sqrt(sx*sy), sx, sy, N, D));
        end
    end
    
    idx_acq = find(abs(datas)>0); %index of acquired data
    idx_nacq = find(abs(datas)==0); %index of non-acquired data
    empNum = length(idx_nacq(:));
    
    target = interpolate(datas,imkernel,0);
    target = [-target(:); zeros(empNum,1)];
    xn = zeros(empNum,1);
    niter = 10; %initial number of lsqr iterations
    obj.trecon = datas; 
    obj.optimParams.relMSEdiff = [];
    obj.optimParams.relMSE(1) = realmin;
    obj.optimParams.relMSE(2) = realmax;
    obj.optimParams.relMSEdiff(1) = 0;
    obj.iter = 1;
    elapsed = 0;
    while 1
        tstart = tic;
        XPre = reshape(ifft2c(obj.trecon),[sx,sy,N,D]);
        [xn, ~] = lsqr(@interp,target,1e-6,niter,speye(empNum,empNum),...
            speye(empNum,empNum),xn,imkernel,sx,sy,N,D,idx_nacq,...
            obj.simParams.lambda); %interpolation projection
        niter = 5;
        obj.trecon(idx_nacq) = xn;  
        
        obj.projSparsity(); %sparsity projection
        obj.trecon(idx_acq) = datas(idx_acq); 

        obj.projTV();   %TV projection
        obj.trecon(idx_acq) = datas(idx_acq); 
        X = reshape(ifft2c(obj.trecon),[sx,sy,N,D]);
        
        obj.optimParams.relMSE(obj.iter+1) = mean(mean(mean(mean(abs(abs(X)-abs(XPre)).^2)))); %MSE w.r.t. the previous iteration
        obj.optimParams.relMSEdiff(obj.iter+1) = 100*(obj.optimParams.relMSE(obj.iter) - ...
            obj.optimParams.relMSE(obj.iter+1))/obj.optimParams.relMSE(obj.iter); %relative change in inter-iteration MSE
        
            if obj.iter==obj.simParams.iterLim(sidx) % stop if maximum iterations reached
                break;
            end
            obj.optimParams.totalCost = abs(obj.optimParams.relMSEdiff); %relative MSE
            if obj.optimParams.totalCost(end) < obj.simParams.thresh % stop if realative MSE is below threshold
                break
            end
        
        if obj.simParams.verbose
            obj.outputInfo();
        end
        
        obj.optimParams.elapsed(sidx,obj.iter) = toc(tstart); 
        elapsed = elapsed + obj.optimParams.elapsed(sidx,obj.iter); %update elapsed time
        
        if ~isempty(obj.simParams.timeLim)
            if elapsed > obj.simParams.timeLim(sidx) %stop if elapsed time reaches the limit
                break;
            end
        end
        xn = obj.trecon(idx_nacq);
        obj.iter = obj.iter+1;
    end
    obj.recon(:,:,sidx) = sos(sos(X,4,p_coils),3,p_acq);
    obj.optimParams.niter(sidx) = obj.iter;
    fprintf('done\n')
end
end