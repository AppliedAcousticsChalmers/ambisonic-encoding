function [out_binaural] = render_ambi_signals_binaurally_custom_hrtf_t(ambi_signals, hrirs_sh, head_orientation_azimuth_rad, N)
% So far, this function works only for sphharm_type = 'real'. 
% hrirs_sh: SH coefficients of the HRIRs

% ------------------------- Rendering, Eq. (12) ---------------------------

% zero pad for the convolution
ambi_signals = [ambi_signals; zeros(size(hrirs_sh.h_l_ring, 1), size(ambi_signals, 2))];

out_binaural = zeros(size(ambi_signals, 1), 2);

for n = 0 : N
    for m = -n : n

        out_binaural(:, 1) = out_binaural(:, 1) + fftfilt(hrirs_sh.h_l_ring(:, n^2+n+m+1) .* cos(m*head_orientation_azimuth_rad) - hrirs_sh.h_l_ring(:, n^2+n-m+1) .* sin(m*head_orientation_azimuth_rad), ambi_signals(:, n^2+n+m+1));
        out_binaural(:, 2) = out_binaural(:, 2) + fftfilt(hrirs_sh.h_r_ring(:, n^2+n+m+1) .* cos(m*head_orientation_azimuth_rad) - hrirs_sh.h_r_ring(:, n^2+n-m+1) .* sin(m*head_orientation_azimuth_rad), ambi_signals(:, n^2+n+m+1));

    end
end


end
