function [coeff,fitY,peakX,peakY] = quadFit(x,y)
%quadFit Quadratic fit
%   x = row vec
%   y = row vec
%   Outputs the following:
%       Polynomial coefficient
%       Y value from fitting
%       X value of the extrema
%       Y value of the extrema

coeff = quadFitM(x,y); % fit parabola
fitY = polyval(coeff, x);
peakX = -0.5*coeff(2)/coeff(1);
peakY = coeff(1)*peakX.^2 + coeff(2)*peakX + coeff(3);
end

function c = quadFitM(x, y)
V = [x.^2 x ones(numel(x),1)];
c = V \ y;
end