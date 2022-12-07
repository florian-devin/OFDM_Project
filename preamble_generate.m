function [preamble] = preamble_generate(length)
% preamble_generate() 
% input : length: a scaler value, desired length of preamble.
% output: preamble: preamble bits
preamble = zeros(length, 1);
LSFR_state = ones(8, 1);
for i = 1:length
    new = mod(sum(LSFR_state([4 5 6 8])), 2);

    preamble(i) = LSFR_state(end);
    LSFR_state(2:end) = LSFR_state(1:end-1);
    LSFR_state(1) = new;
end
end