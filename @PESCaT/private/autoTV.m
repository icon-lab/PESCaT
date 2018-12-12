function [TVWeight] = autoTV(img)
%
% AUTO TV Calculates the regularization term for total variation penalty.
%
% Input:
% img:      Image to calculate the TV penalty from.
%
% Output:
% TVWeight: Weight for total variaton regularization term.

% Filter Image
HSIZE = round((size(img))/25);
SIGMA = 3;
h = fspecial('gaussian', HSIZE, SIGMA);
tmp = imfilter(img,h);

% Find and normalize local variance of filtered image
nsize(1) = round(size(img,1)/100)+1-mod(round(size(img,1)/100),2);
nsize(2) = round(size(img,2)/100)+1-mod(round(size(img,2)/100),2);
NHOOD = ones(nsize);
% LV = stdfilt(abs(tmp),NHOOD);
% LVn = 0.5*norm2d(LV);
tmp = norm2d(abs(tmp));
LV = stdfilt(tmp,NHOOD);
LVn = sqrt(LV);

% Obtain best TV regularization term as the 1/3 of median of local variance
TVWeight = (median(LVn(:))/13);