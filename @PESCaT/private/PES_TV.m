function [out,lambda,itv] = PES_TV(imIN,TValpha)
% PES_TV performes the PES-TV algorithm on the input image.
% Code is borrowed from: http://signal.ee.bilkent.edu.tr/PES-TV.html
%
%  (c) ICON Lab 2018
im1 = imIN;

iter = 3000;
[ny nx] = size(im1);
y0 = im1(:);

x = [im1(:); 0 ];
imn = im1(:);
snr_oi(1) = -inf;
for i = 2:iter
    
    cost0 = costNorm(y0, x(1:end-1));
    
    x1 = [x(1:end-1); 0];
    x2 = [x(1:end-1); TValpha*tv(x(1:end-1), size(im1))];
    x3 = proj(x1, x2, ny);
    z = x3(end);
    
    cost1 = costNorm(y0, x3(1:end-1));
    
    if (cost1 > cost0)
        xLeft = [x3(1:end-1); TValpha*tv(x(1:end-1), size(im1))];
        xRight = x2;
        x3 = (xLeft + (1*i)*xRight)/(1*i+1);
        z = x3(end);
        x3(end) = 0;
        cost1 = costNorm(y0, x3(1:end-1));
    else
        xLeft = [x3(1:end-1); tv(x(1:end-1), size(im1))];
        xRight = x2;
        x3 = ((1*i)*xLeft + xRight)/(1*i+1);
        z = x3(end);
        x3(end) = 0;
        cost1 = costNorm(y0, x3(1:end-1));
    end
    
    x = x3;
    snr_oi(i) = snr(imn(:),x(1:end-1));
    if (snr_oi(i)<snr_oi(i-1))
        break
    end
    imn = x(1:end-1);
end

xx = x(1:end-1);
lambda = (z/tv(xx, size(im1)))^2;
xx = reshape(xx,size(im1));
out = xx;
itv = lambda*tv(xx, size(im1)).^2;
