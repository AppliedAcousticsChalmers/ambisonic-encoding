function [coefficients] = least_squares_sh_fit(N, data, azimuth, colatitude, sphharm_type, reg_type, reg_parameter)
% Performs an Nth-order least squares fit of spherical harmonics on data 
% that is available in the directions (azimuth, colatitude) for each 
% frequency bin separately using the pseudo inverse
%
%   data: frequency goes downwards, position to the right
%   azimuth, colatitude: row vectors in radians
%   sphharm_type: see file sphharm.m for the options
%   reg_type (optional): '' (default), 'constant' or 'quadratic'
%                        weighting of reg_parameter w.r.t. the order n
%   reg_parameter (optional): regularization parameter; defaults to 0
%
%   coefficients: frequency goes downwards, order to the right (sorted
%                 according to (n^2+n+m+1)
%
% Author: Jens Ahrens, March 2020

if nargin < 6
    reg_type = '';
end

if nargin < 7
    reg_parameter = 0;
end

azimuth    = azimuth.';
colatitude = colatitude.';
data       = data.';

% generate a series of n: (0, 1, 1, 1, 2, 2, 2, 2, 2, 3, ...)
ns = zeros(0, 1);

for n = 0 : N
    for m = -n : n
        ns = [ns; n];
    end
end

% -------- compute basis functions ----------
Y_nm = zeros(size(azimuth, 1), (N+1)^2);

for n = 0 : N
    for m = -n : n
        Y_nm(:, n^2+n+m+1) = sphharm(n, m, colatitude, azimuth, sphharm_type);
    end
end

% --------- perform the ls fit w/ or w/o Tikhonov regularization ----------
fprintf('Performing LS fit... ');

if isempty(reg_type)
    coefficients = pinv(Y_nm) * data;

elseif strcmp(reg_type, 'constant')
    coefficients = inv(Y_nm' * Y_nm + reg_parameter * eye((N+1)^2)) * Y_nm' * data;

elseif strcmp(reg_type, 'quadratic')
    % as proposed in Duraiswami et al., IEEE ICASSP, 2004
    coefficients = inv(Y_nm' * Y_nm + reg_parameter * diag(1 + ns .* (ns + 1))) * Y_nm' * data;

else
    error('Unknown regularization type.');

end
    
fprintf('done.\n');

coefficients = coefficients.';

end

