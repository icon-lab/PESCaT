function [ tv, df ] = tv( im,rs )
% Computes the total-variation norm of the input image.
% Code is borrowed from: http://signal.ee.bilkent.edu.tr/PES-TV.html
%
%  (c) ICON Lab 2018

% Total variation
im = reshape(im, rs);

zVec1 = zeros(size(im,1),1);
im1 = [zVec1  im  zVec1];

zVec2 = zeros(1,size(im1,2));
im = [zVec2;  im1;  zVec2];

df1 = conv2(im,[1 -1],'valid');
df2 = conv2(im,[1; -1],'valid');
tv = sum(abs(df1(:))) + sum(abs(df2(:)));

end