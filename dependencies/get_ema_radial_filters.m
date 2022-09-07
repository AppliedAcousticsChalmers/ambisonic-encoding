function [ema_inv_rf, ema_inv_rf_t] = get_ema_radial_filters(k, R, N, limit_db, reg_type, hankel_type, sphharm_type)
% ema_inv_rf and ema_inv_rf_t are the 1st factor on the right hand side of
% Eq. (13). They are causal (they comprise comprises a modeling delay).
%
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%
% Written by Jens Ahrens, 2022


% ----------------------- Compute the terms b_n, Eq. (7) ------------------
b_n = zeros(size(k, 1), N+1);

kR = k.*R;

for n = 0 : N
       
    hankel_prime = 1/(2*n+1) * (n * sphbesselh(n-1, hankel_type, kR) - (n+1) * sphbesselh(n+1, hankel_type, kR));
    
    b_n(:, n+1) = 1i./kR.^2 .* 1./hankel_prime;
    
    if hankel_type == 1
        b_n(:, n+1) =   4*pi * 1i^(-n) .* b_n(:, n+1);
    elseif hankel_type == 2
        b_n(:, n+1) = - 4*pi * 1i^n    .* b_n(:, n+1);
    else 
        error('Unknown hankel_type.');
    end
    
end

% ------------ Compute the inverse EMA radial filters, Eq. (13) -----------

ema_rf = zeros(size(b_n, 1), 2*N+1);

for m = -N : N

    for n_prime = abs(m) : N
                
        if strcmp(sphharm_type, 'complex')
            ema_rf(:, m+N+1) = ema_rf(:, m+N+1) + b_n(:, n_prime+1) .* sphharm(n_prime, m, pi/2, 0, 'complex').^2;
        elseif strcmp(sphharm_type, 'real')
            ema_rf(:, m+N+1) = ema_rf(:, m+N+1) + b_n(:, n_prime+1) .* N_nm(n_prime, m, pi/2).^2;
        else
            error('Unknown sphharm_type.');
        end
        
    end
    
end


ema_inv_rf = 1./ema_rf;

% ----------------------------- limiting ----------------------------------
if limit_db < inf
    
    % --- soft clipping ---
    if strcmp(reg_type, 'soft')
    
        ema_inv_rf = (2*10^(limit_db/20)) / pi * ema_inv_rf ./ abs(ema_inv_rf) .* atan(pi/(2*10^(limit_db/20)) .* abs(ema_inv_rf));
        
    % --- hard clipping ---
    elseif strcmp(reg_type, 'hard')
        
        ema_inv_rf(abs(ema_inv_rf) > 10^(limit_db/20)) = ema_inv_rf(abs(ema_inv_rf) > 10^(limit_db/20)) ./ abs(ema_inv_rf(abs(ema_inv_rf) > 10^(limit_db/20))) * 10^(limit_db/20); 
        
    % --- Moreau regularization (Moreau, Daniel, Bertet, AES 2006) ---
    elseif strcmp(reg_type, 'moreau')
        
       	limit          = 10^(limit_db/20);
        lambda_squared = (1 - sqrt(1 - 1/limit^2)) ./ (1 + sqrt(1 - 1/limit^2));
       
        ema_inv_rf     = conj(ema_rf) ./ (abs(ema_rf).^2 + lambda_squared);        
    
    else
        
        error('Unknown reg_type.');
        
    end

end

% --------------------------- make filters causal -------------------------

spec       = [ema_inv_rf; conj(flipud(ema_inv_rf(2:end-1, :)))];
spec(1, :) = abs(spec(2, :)); % fix DC

ema_inv_rf_t = ifft(spec, [], 1, 'symmetric'); 

filter_length = size(ema_inv_rf_t, 1);

% make causal (we assume that the filter_length is 2048 or longer)
ema_inv_rf_t = circshift(ema_inv_rf_t, [filter_length/2 0]);

% fade in and out
window = hann(400);
fade_in = window(1:end/2);
fade_out = window(end/2+1:end);

% apply fade window
ema_inv_rf_t(1:200, :)         = ema_inv_rf_t(1:200, :)         .* repmat(fade_in,  [1 size(ema_inv_rf_t, 2)]);
ema_inv_rf_t(end-200+1:end, :) = ema_inv_rf_t(end-200+1:end, :) .* repmat(fade_out, [1 size(ema_inv_rf_t, 2)]);

spec = fft(ema_inv_rf_t, [], 1);
ema_inv_rf = spec(1:end/2+1, :);

end


