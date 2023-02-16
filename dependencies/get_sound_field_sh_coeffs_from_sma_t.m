function [s_breve] = get_sound_field_sh_coeffs_from_sma_t(array_signals, sma_inv_rf_t, N, beta, alpha, grid_weights)
%
% The function assumes that the grid_weights are normalized such that they 
% sum up to 1.
%
% It makes no real sense to use it with anything other than sphharm_type =
% 'real'.
%
% The equation numbers refer to 
%
%   Jens Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00583
%
% Written by Jens Ahrens, 2022

% --------------------------- Evaluate Eq. (10) ---------------------------

s_breve = zeros(size(array_signals, 1), (N+1)^2);

for n = 0 : N 
    for m = -n : n        
        
        % evaluate the transformation integral discretely
        s_breve(:, n^2+n+m+1) = sum(array_signals .* repmat(grid_weights .* sphharm(n, m, beta, alpha, 'real'), [size(array_signals, 1) 1]), 2) * 4*pi;
        
        % If the microphone placement grid is irregular, you might want to try a
        % least-squares fit as performed in the script
        % render_ambi_signals_binaurally_t.m.

        s_breve(:, n^2+n+m+1) = fftfilt(sma_inv_rf_t(:, n+1), s_breve(:, n^2+n+m+1));
        
    end
end

end




