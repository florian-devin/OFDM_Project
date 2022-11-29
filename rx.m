function [rxbits conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

% dummy 

time = 1:3/length(rxsignal):4-3/length(rxsignal);
rxsignal_shifted = rxsignal .* exp(1i*2*pi*400*time)'; 
%rxsignal_filtered = lowpass(rxsignal,conf);
rxsignal_filtered = rxsignal_shifted;
downsampled_rxsignal = rxsignal_filtered(5:5:end);
figure(2);
plot(time(5:5:end), downsampled_rxsignal);
rolloff = 0.22;
pulse = rrc(conf.os_factor,rolloff, 20);

filtered_rx_signal = conv(downsampled_rxsignal, conj(pulse))

preamble = zeros(conf.npreamble, 1);
LSFR_state = ones(8, 1);
for i = 1:conf.npreamble
    new = mod(sum(LSFR_state([4 5 6 8])), 2);

    preamble(i) = LSFR_state(end);
    LSFR_state(2:end) = LSFR_state(1:end-1);
    LSFR_state(1) = new;
end

preamble_bpsk = -2*(preamble) + 1;

start = detector(preamble_bpsk, filtered_rx_signal, 0.5)-1;

filtered_rx_signal = filtered_rx_signal(start:end);

rxbits = demapper(filtered_rx_signal);

%rxbits = zeros(conf.nbits,1);
