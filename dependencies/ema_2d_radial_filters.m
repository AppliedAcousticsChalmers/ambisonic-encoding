function [b_nm_inv, b_mn_inv_t] = ema_2d_radial_filters(N, kr, radial_filter_limit_dB, array_type)
% array_type: 'cardioid', 'omni'
% Includes SMA regularization mask
% f_limit_per_order: for omni arrays, length N+1 vector

% assure correct dimensions
kr = kr(:);

% regularize to avoid nuerical issues
limit_tmp = 10^(200/20);

% get the radial filter regularization (convert limit_db to SMA convention)
% See AES AVARIG, 2026 paper
% Give it 15 dB extra room
[~, ~, ~, radial_filter_mask] = get_sma_radial_filters(kr, N, radial_filter_limit_dB-22+15, 'tikhonov', 2);

b_nm_inv = zeros(length(kr), (N+1)^2);

for n = 0 : N
    for m = -n : n
        
        % compute radial filter
        if strcmp(array_type, 'omni')
            b_m = 1i^m .* besselj(m, kr);
        elseif strcmp(array_type, 'cardioid')
            % compute radial filter
            b_m = 1i^m .* (besselj(m, kr) - 1i .* besselj_prime(m, kr));
        else
            error('Unknown array type.');
        end

        % tikhonov regularizarion and inversion
        lambda_squared = (1 - sqrt(1 - 1/limit_tmp^2)) ./ (1 + sqrt(1 - 1/limit_tmp^2));
        b_nm_inv(:, n^2+n+m+1) = conj(b_m) ./ (abs(b_m).^2 + lambda_squared);


        % apply regularization mask
        b_nm_inv(:, n^2+n+m+1) = b_nm_inv(:, n^2+n+m+1) .* radial_filter_mask(:, n+1);

    end
end

% extra regularisation for omni-filter to obtain a useful time-damain
% representation
% (set all values at f > 1000 to f = 1000)
if strcmp(array_type, 'omni')
    f = kr/(2*pi*0.06)*343;
    index = sum(f < 2000);
    b_nm_inv(index+1:end, :) = repmat(b_nm_inv(index, :), size(b_nm_inv, 1)-index, 1);
end

% fix DC
b_nm_inv(1, :) = real(b_nm_inv(2, :));

% fix Nyquist
b_nm_inv(end, :) = real(b_nm_inv(end, :));

b_mn_inv_t = ifft([b_nm_inv; conj(flipud(b_nm_inv(2:end-1, :)))], 'symmetric');

% make causal
b_mn_inv_t = circshift(b_mn_inv_t, size(b_mn_inv_t, 1)/2);

end



