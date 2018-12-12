function PES_L1(obj,X,W)
% PES_L1(X,J) performs joint PES-L1 projection algorithm on
% high-frequency subbands of the image X in the wavelet domain.
% obj: PESCaT object
% X: input image
% W: Wavelet object
%
% see also: PESCAT, PROJSPARSITY
%
%  (c) ICON Lab 2018

params = obj.simParams;
[sx,sy,npc] = size(X);
X = zpad(X,params.ss,params.ss,size(X,3),size(X,4));
maskbase = [size(X,1),size(X,2)]/16;
lambda = zeros(4,3);
xwi = W*X;
cost = [];
wden = zeros(size(xwi));
for level = 1:4
    wdenl = zeros([maskbase*2^level,size(X,3)]);
    for subband = 1:3
        mask = zeros(size(wdenl));
        switch subband
            case 1
                mask(1:size(mask,1)/2,1+size(mask,2)/2:end,:) = 1;
            case 2
                mask(1+size(mask,1)/2:end,1:size(mask,2)/2,:) = 1;
            case 3
                mask(1+size(mask,1)/2:end,1+size(mask,2)/2:end,:) = 1;
        end
        Xsubt = xwi(1:size(mask,1),1:size(mask,2),:).*mask;
        cost = [cost;Xsubt(:)];
        if isempty(obj.simParams.sparsityLambda)
            oPhase = angle(Xsubt); % store original phase
            Xabs = reshape(abs(Xsubt),[],1);% perform projection on abs
            Xsub = padarray(Xabs, 1, 'symmetric');
            wn = [Xsub; 0];
            wsgn = [sign(wn(1:end-1)); -1];
            wsgn(wsgn==0) = 1;
            ww = params.L1beta*(wn'*wsgn);
            wp = ProjectOntoL1Ball(wn,ww);
            lambda(level,subband) = (ww/(wp'*sign(wp)))^2;
            wp = wp(2:end-2);
            wd = reshape(wp,size(mask));
            wd = wd.*exp(1i*oPhase); %restore the phase
        else
            wd = softThresh(Xsubt,obj.simParams.sparsityLambda(level,subband));
        end
        switch subband
            case 1
                wdenl(1:size(mask,1)/2,1+size(mask,2)/2:end,:) = wd(1:size(mask,1)/2,1+size(mask,2)/2:end,:);
            case 2
                wdenl(1+size(mask,1)/2:end,1:size(mask,2)/2,:) = wd(1+size(mask,1)/2:end,1:size(mask,2)/2,:);
            case 3
                wdenl(1+size(mask,1)/2:end,1+size(mask,2)/2:end,:) = wd(1+size(mask,1)/2:end,1+size(mask,2)/2:end,:);
        end
    end
    wden(1:size(wdenl,1),1:size(wdenl,2),:) = wdenl + wden(1:size(wdenl,1),1:size(wdenl,2),:);
end
wden(1:maskbase(1),1:maskbase(2),:) = xwi(1:maskbase(1),1:maskbase(2),:);
wden = W'*wden;
Wden(:,:,:, 1) = wden;
Xp = squeeze(mean(Wden, 4));
Xp = crop(Xp,sx,sy,npc);
obj.trecon = fft2c(Xp);
end