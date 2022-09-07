function [out] = N_nm(n, m, beta)
% Eq. (3) from 
%
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%
% Written by Jens Ahrens, 2022

% Don't use 'real' here!
out = sphharm(n, m, beta, 0, 'complex'); 

% For the beauty: Remove stray imaginary part
out = real(out);

end

