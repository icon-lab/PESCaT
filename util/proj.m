function [proj, v] = proj(x1, x2, ny)
% PROJ is used by the projTV method to project onto the supporting
% hyperplane during PES-TV operation
% 
% see also: PROJTV
%
%  (c) ICON Lab 2018
ox = x1;

ox2 = x2;
x2 = x2(1:end-1);

im = reshape(x2, ny, []);
[ny nx] = size(im);

%% Gradient 
zVec1 = zeros(size(im,1),1);
im1 = [zVec1  im  zVec1];

zVec2 = zeros(1,size(im1,2));
im = [zVec2;  im1;  zVec2];

df1 = sign(conv2(im,[1 -1],'valid'));
df2 = sign(conv2(im,[1; -1],'valid'));

zVec1 = zeros(size(df1,1),1);
zVec2 = zeros(1,size(df2,2));

df1 = [zVec1  df1  zVec1];
df2 = [zVec2;  df2;  zVec2];

grad_tv = conv2(df1,[-1 +1],'valid') + conv2(df2,[-1; 1],'valid');
grad_tv = grad_tv(2:ny+1,2:nx+1);

%% Projection 
normal = grad_tv(:);

v = [normal; -1];
v = v./norm(v);

AA = (v*(v'*(ox-ox2)));

proj = ox - (AA);
end