function [res,tflag] = interp(x,imkernel,sx,sy,N,D,idx_nacq,lambda,tflag)
% Used by the LSQR to perform the calibration consistency projection during
% PESCAT.
%
%  (c) ICON Lab 2018
if ~strcmp(tflag,'transp')
    tmpx = zeros(sx,sy,N,D);
    tmpx(idx_nacq) = x;
    res = interpolate(tmpx,imkernel,0);
    res = [res(:) ; sqrt(lambda)*x(:)];
else
    tmpy = reshape(x(1:sx*sy*N*D),sx,sy,N,D);
    res = interpolate(tmpy,imkernel,1);
    res = res(idx_nacq) + x(sx*sy*N*D+1:end)*sqrt(lambda);
end
end