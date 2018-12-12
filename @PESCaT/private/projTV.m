function projTV(obj)
% PROJTV performes the TV projection while solving the PESCAT
% reconstruction problem. The method used for determining the
% regularization parameter is read from PESCAT object properties.
%
% see also: PESCAT
%
%  (c) ICON Lab 2018

X = ifft2c(obj.trecon);
switch obj.simParams.TVType
    case 'PES' %TV projection using PES
        X = reshape(X,[size(X,1),size(X,2),size(X,3)*size(X,4)]);
        for pidx=1:size(X,3)
            img = X(:,:,pidx);
            [resTV,l,itv] = PES_TV(img,obj.simParams.TVbeta);
            X(:,:,pidx) = resTV;
        end
    case 'STD' %TV projection using local STD
        X = reshape(X,[size(X,1),size(X,2),size(X,3)*size(X,4)]);
        for pidx=1:size(X,3)
            img = X(:,:,pidx);
            TVWeight = autoTV(img);
            if isempty(obj.simParams.TVLambda)
                [resTV]=denoiseTV(img,TVWeight,obj.simParams.TVNit,0);
            else
                [resTV]=denoiseTV(img,obj.simParams.TVLambda,obj.simParams.TVNit,0);
            end
            X(:,:,pidx) = resTV;
        end
end
obj.trecon = fft2c(X);
end
