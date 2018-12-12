function res = norm2d(fin)
% Normalize absolute value of data between 0 and 1
res = (fin-min(abs(fin(:))))/(max(abs(fin(:)))-min(abs(fin(:))));