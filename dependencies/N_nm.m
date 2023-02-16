function [out] = N_nm(n, m, beta)
% Eq. (3) from 
%
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%
% Written by Jens Ahrens, 2022

Lnm = asslegendre(n, abs(m), cos(beta));

factor_1 = (2*n + 1) / (4*pi);
factor_2 = factorial(n - abs(m)) ./ factorial(n + abs(m));

out = (-1).^m .* sqrt(factor_1 .* factor_2) .* Lnm;

end

