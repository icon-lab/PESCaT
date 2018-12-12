function train (obj)
% TRAIN trains the tensor interpolation kernel used in PESCaT
%
% see also: PESCAT/RECONPESCAT
%  (c) ICON Lab 2018

fprintf('training the kernel...')
data = obj.data(:,:,obj.simParams.slcs,:,:);
[~,~,sl,N,D] = size(data);
[calibsize,~] = getCalibSize(obj.sampling.mask(:,:,1,1,1));
for sidx = 1:sl
        fdata = data(:,:,sidx,:,:); %data in fourier domain
        ksize = round(calibsize/obj.simParams.den)+1-mod(round(calibsize/obj.simParams.den),2); %calibration area size
        kcalib = crop(fdata,[calibsize,1,N,D]); % get data from calibration area
        [YtY,Y] = dat2AtA(kcalib,ksize); %prepare the data to be used in kernel calculation
        for pidx = 1:N
            for cidx = 1:D
                if obj.simParams.verbose
                    fprintf('calibrating for slice: %d, acq. %d, coil %d\n',sidx,pidx,cidx);
                end
                kernel(:,:,:,:,pidx,cidx,sidx) = calibrate(YtY,Y,ksize,N,D,pidx,cidx,obj.simParams.TykWeight);
            end
        end
end
obj.data = data;
obj.simParams.isTrained = 1;
obj.kernel = kernel;
fprintf('done\n')
end