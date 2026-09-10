function [y_lp, y_hp] = apply_crossover(sig_lp, sig_hp, fs, fc)

% Linkwitz Riley
crossFilt = crossoverFilter( ...
                    NumCrossovers=1, ...
                    CrossoverFrequencies=fc, ...
                    CrossoverSlopes=48, ...
                    SampleRate=fs);

% use loooong impulse responses for a smooth result
imp = zeros(8*2048, 1);
imp(1) = 1;

[ir_lp, ir_hp] = crossFilt(imp);

y_lp = conv_multichannel(ir_lp, sig_lp);
y_hp = conv_multichannel(ir_hp, sig_hp);

end

