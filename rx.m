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
close all;
t = 0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s;
rxsignal= rxsignal.* exp(-1i*2*pi*(conf.f_c*t')); 
%rxsignal_filtered = lowpass(rxsignal,conf);
f_corner = 1.05*conf.bandwidth; %TODO from where come this 1.05
rxsignal_filtered = 2 * ofdmlowpass(rxsignal, conf, f_corner);
start_idx = frame_sync(rxsignal_filtered, conf);
%rx_size = conf.os_factor_ofdm * ((conf.nbdatapertrainning*conf.cp_length)+((conf.nbdatapertrainning+1)*conf.nbcarriers));
rx_size = conf.os_factor_ofdm*(conf.nbcarriers + conf.cp_length)*ceil(((conf.nbits/conf.modulation_order)/conf.nbcarriers)); 
rxsignal_filtered = rxsignal_filtered(start_idx:start_idx+rx_size-1);

rxsignal_filtered = reshape(rxsignal_filtered, [conf.os_factor_ofdm*(conf.nbcarriers+conf.cp_length), conf.nbdatapertrainning]);

% remove CP
rxsignal_filtered(1:conf.os_factor_ofdm*conf.cp_length,:) = [];

%OFDM demod
%TODO debug this part
for data_idx = 1:conf.nbdatapertrainning
    rxsignal_filtered(:,data_idx) = osfft(rxsignal_filtered(:,data_idx), conf.os_factor_ofdm);
end
%Channel estimation

% rxsignal_filtered = rxsignal_shifted;
% downsampled_rxsignal = rxsignal_filtered(5:5:end);
% figure(2);
% plot(time(5:5:end), downsampled_rxsignal);
% rolloff = 0.22;
% pulse = rrc(conf.os_factor,rolloff, 20);
% 
% filtered_rx_signal = conv(downsampled_rxsignal, conj(pulse))
% 
% preamble = zeros(conf.npreamble, 1);
% LSFR_state = ones(8, 1);
% for i = 1:conf.npreamble
%     new = mod(sum(LSFR_state([4 5 6 8])), 2);
% 
%     preamble(i) = LSFR_state(end);
%     LSFR_state(2:end) = LSFR_state(1:end-1);
%     LSFR_state(1) = new;
% end
% 
% preamble_bpsk = -2*(preamble) + 1;
% 
% start = detector(preamble_bpsk, filtered_rx_signal, 0.5)-1;
% 
% filtered_rx_signal = filtered_rx_signal(start:end);
% 
% rxbits = demapper(filtered_rx_signal);

%rxbits = zeros(conf.nbits,1);
