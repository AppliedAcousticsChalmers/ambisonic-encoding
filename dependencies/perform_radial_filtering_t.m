function [s_breve] = perform_radial_filtering_t(s_breve, radial_filter_ir)

N = sqrt(size(s_breve, 2)) - 1;

% ----------------------- check if EMA or SMA -----------------------------
if size(radial_filter_ir, 2)-1 == N
    array_type = 'SMA';
elseif size(radial_filter_ir, 2) == size(s_breve, 2)
    array_type = 'EMA';
else
    error('Radial filters have wrong number of channels.');
end

% --------------------- perform radial filtering --------------------------
for n = 0 : N 
    for m = -n : n        
    
        if strcmp(array_type, 'SMA')
            s_breve(:, n^2+n+m+1) = fftfilt(radial_filter_ir(:, n+1      ), s_breve(:, n^2+n+m+1));
        elseif strcmp(array_type, 'EMA')
            s_breve(:, n^2+n+m+1) = fftfilt(radial_filter_ir(:, n^2+n+m+1), s_breve(:, n^2+n+m+1));
        end

    end
end

end




