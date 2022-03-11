function y = hrfGamma(t,T,W,A)
%hrfGamma Hemodynamic response gamma fit
%   t = time points
%   T = time to peak
%   W = full width half max
%   A = amplitude


alpha = (T/W)^2*8*log(2);
beta = W^2/(T*8*log(2));
y = A*(t/T).^alpha.*exp((t-T)/(-beta));
end