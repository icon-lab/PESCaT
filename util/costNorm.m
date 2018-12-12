function [costNorm] = costNorm(y0, x)

costNorm = sum((x(:) - y0(:)).^2);

end