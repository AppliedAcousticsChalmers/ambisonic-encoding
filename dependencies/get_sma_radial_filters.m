function [sma_inv_rf, sma_inv_rf_t] = get_sma_radial_filters(k, R, N, limit_db, reg_type, hankel_type)
% ema_inv_rf and ema_inv_rf_t are the 1st factor on the right hand side of
% Eq. (15). They are causal (they comprise comprises a modeling delay).
%
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%
% Written by Jens Ahrens, 2022

% ----------------------- compute the terms b_n, Eq. (6) ------------------
b_n = zeros(size(k, 1), N+1);

kR = k.*R;

for n = 0 : N
       
    hankel_prime = 1/(2*n+1) * (n * sphbesselh(n-1, hankel_type, kR) - (n+1) * sphbesselh(n+1, hankel_type, kR));
         
    b_n(:, n+1) = 1i./kR.^2 .* 1./hankel_prime;
    
    if hankel_type == 1
        b_n(:, n+1) =   4*pi * 1i^(-n) .* b_n(:, n+1);
    elseif hankel_type == 2
        b_n(:, n+1) = - 4*pi * 1i^n    .* b_n(:, n+1);
    end
    
end


sma_inv_rf = 1./b_n;

% ----------------------------- limiting ----------------------------------
if limit_db < inf
    
    % --- soft clipping ---
    if strcmp(reg_type, 'soft')
    
        sma_inv_rf = (2*10^(limit_db/20)) / pi * sma_inv_rf ./ abs(sma_inv_rf) .* atan(pi/(2*10^(limit_db/20)) .* abs(sma_inv_rf));
        
    % --- hard clipping ---
    elseif strcmp(reg_type, 'hard')
        
        sma_inv_rf(abs(sma_inv_rf) > 10^(limit_db/20)) = sma_inv_rf(abs(sma_inv_rf) > 10^(limit_db/20)) ./ abs(sma_inv_rf(abs(sma_inv_rf) > 10^(limit_db/20))) * 10^(limit_db/20); 
        
    % --- Moreau regularization (Moreau, Daniel, Bertet, AES 2006) ---
    elseif strcmp(reg_type, 'moreau')
        
       limit          = 10^(limit_db/20);
       lambda_squared = (1 - sqrt(1 - 1/limit^2)) ./ (1 + sqrt(1 - 1/limit^2));
       
       sma_inv_rf     = conj(b_n) ./ (abs(b_n).^2 + lambda_squared);      
    
    else
        
        error('Unknown reg_type.');
        
    end

end

% --------------------------- make filters causal -------------------------

spec       = [sma_inv_rf; conj(flipud(sma_inv_rf(2:end-1, :)))];
spec(1, :) = abs(spec(2, :)); % fix DC

sma_inv_rf_t = ifft(spec, 'symmetric'); 

filter_length = size(sma_inv_rf_t, 1);

% make causal (we assume that the filter_length is 2048 or longer)
sma_inv_rf_t = circshift(sma_inv_rf_t, [filter_length/2 0]);

% fade in and out
window = hann(400);
fade_in = window(1:end/2);
fade_out = window(end/2+1:end);

% apply fade window
sma_inv_rf_t(1:200, :)         = sma_inv_rf_t(1:200, :)         .* repmat(fade_in,  [1 size(sma_inv_rf_t, 2)]);
sma_inv_rf_t(end-200+1:end, :) = sma_inv_rf_t(end-200+1:end, :) .* repmat(fade_out, [1 size(sma_inv_rf_t, 2)]);

spec = fft(sma_inv_rf_t);
sma_inv_rf = spec(1:end/2+1, :);

end


