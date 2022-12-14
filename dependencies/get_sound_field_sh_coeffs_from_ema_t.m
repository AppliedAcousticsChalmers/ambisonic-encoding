function [s_breve] = get_sound_field_sh_coeffs_from_ema_t(array_signals, ema_inv_rf_t, N, alpha_ema, sphharm_type)
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00584
%
% Written by Jens Ahrens, 2022

% -------- PWD of the sound pressure along the equator, Eq. (8) -----------

s_ring_m_surf = zeros(size(array_signals, 1), 2*N+1);

if strcmp(sphharm_type, 'complex')
    
    for m = -N : N
        % discretized transformation integral
        s_ring_m_surf(:, m+N+1) = sum(array_signals .* exp(-1i .* m .* repmat(alpha_ema, size(array_signals, 1), 1)), 2) / size(alpha_ema, 2);
    end

elseif strcmp(sphharm_type, 'real')
    
    for m = -N : -1
        % discretized transformation integral
        s_ring_m_surf(:, m+N+1) = sum(array_signals .* sqrt(2) .* sin(abs(m) .* alpha_ema), 2) / size(alpha_ema, 2);
    end
    
    % m = 0
    s_ring_m_surf(:, 0+N+1) = sum(array_signals, 2) / size(alpha_ema, 2);
    
    for m = 1 : N
        % discretized transformation integral
        s_ring_m_surf(:, m+N+1) = sum(array_signals .* sqrt(2) .* cos(m .* alpha_ema), 2) / size(alpha_ema, 2);
    end

else
    error('Unknown sphharm_type.');
end


% --------------------------- Evaluate Eq. (13) ---------------------------

s_ring_bar_m = fftfilt(ema_inv_rf_t, s_ring_m_surf);

% --------------------------- Evaluate Eq. (15) ---------------------------

s_breve = zeros(size(s_ring_bar_m, 1), (N+1)^2);

for n = 0 : N
    for m = -n : n % for m = -n : 2 : n
        
        % The cases are redundant ... anyway.
        if strcmp(sphharm_type, 'complex')
            s_breve(:, n^2+n+m+1) = s_ring_bar_m(:, m+N+1) .* sphharm(n, m, pi/2, 0);
        elseif strcmp(sphharm_type, 'real')
        	s_breve(:, n^2+n+m+1) = s_ring_bar_m(:, m+N+1) .* N_nm(n, m, pi/2);
        else
            error('Unknown sphharm_type.');
        end

    end
end

end




