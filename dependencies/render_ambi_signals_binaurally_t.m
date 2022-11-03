function [out_lr] = render_ambi_signals_binaurally_t(ambi_signals, head_orientation_azimuth, N, fit_type)
% fit_type: 'transform_integral' or 'ls'
%
% So far, this function works only for sphharm_type = 'real'. It is not
% conventient to extend it for complex SHs because this would produce
% time-domain output signals that are complex.
%
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00583
%
% and/or
%
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00584
%
% Written by Jens Ahrens, 2022

% Get Lebedev grid weights for the quadrature of the transformation 
% integral; not need if least squares SH fit is performed. See lines 68-69.
tmp = load('resources/2702_lebedev_grid.mat');
grid_weights = tmp.grid_weights;
clear tmp;

% download HRIRs (skipped automatically if the HRIR dataset already exists)
hrir_path = 'resources/HRIR_L2702.sofa';

if ~isfile(hrir_path)
    fprintf('Downloading HRTFs from https://zenodo.org/record/3928297 ... ');
    websave(hrir_path, 'https://zenodo.org/record/3928297/files/HRIR_L2702.sofa?download=1');
    fprintf('done.\n');
end

SOFAstart;

hrirs_sofa = SOFAload('resources/HRIR_L2702.sofa');

% extract irs
hrirs_l = squeeze(hrirs_sofa.Data.IR(:, 1, :)).';
hrirs_r = squeeze(hrirs_sofa.Data.IR(:, 2, :)).';

% extract incidence angles
azi_hrtfs =        hrirs_sofa.SourcePosition(:, 1).' / 180 * pi;
col_hrtfs = pi/2 - hrirs_sofa.SourcePosition(:, 2).' / 180 * pi; % ele to col

% SH decomposition in time domain
% h_l_ring and h_r_ring are defined via Eq. (13)/(17)

% ---------------------- SH decomposition of HRTFs ------------------------

if strcmp(fit_type, 'transform_integral')
    % evalute the discretized transformation integral (requires grid_weights)
    
    h_l_ring = zeros(size(hrirs_l, 1), (N+1)^2);
    h_r_ring = zeros(size(hrirs_r, 1), (N+1)^2);

    for n = 0 : N 
        for m = -n : n        

            h_l_ring(:, n^2+n+m+1) = sum(hrirs_l .* repmat(grid_weights .* sphharm(n, m, col_hrtfs, azi_hrtfs, 'real'), [size(hrirs_l, 1) 1]), 2) * 4*pi;
            h_r_ring(:, n^2+n+m+1) = sum(hrirs_r .* repmat(grid_weights .* sphharm(n, m, col_hrtfs, azi_hrtfs, 'real'), [size(hrirs_r, 1) 1]), 2) * 4*pi;

        end
    end
    
elseif strcmp(fit_type, 'ls')
    
    % perform a (potentially regularized) least-squares fit; this works
    % also with kind-of irregular grids
    h_l_ring = least_squares_sh_fit(N, hrirs_l, azi_hrtfs, col_hrtfs, 'real', 'constant', 1e-3); % do some mild regularization, just to be sure
    h_r_ring = least_squares_sh_fit(N, hrirs_r, azi_hrtfs, col_hrtfs, 'real', 'constant', 1e-3); % do some mild regularization, just to be sure
    
else 
    error('Unknown fit_type.');
end
    
out_lr = zeros(size(ambi_signals, 1), 2);

% ---------------------- Rendering, Eq. (12)/(16) -------------------------
for n = 0 : N
    for m = -n : n
        
        out_lr(:, 1) = out_lr(:, 1) + fftfilt(h_l_ring(:, n^2+n+m+1) .* cos(m*head_orientation_azimuth) - h_l_ring(:, n^2+n-m+1) .* sin(m*head_orientation_azimuth), ambi_signals(:, n^2+n+m+1));
        out_lr(:, 2) = out_lr(:, 2) + fftfilt(h_r_ring(:, n^2+n+m+1) .* cos(m*head_orientation_azimuth) - h_r_ring(:, n^2+n-m+1) .* sin(m*head_orientation_azimuth), ambi_signals(:, n^2+n+m+1));

    end
end

end
