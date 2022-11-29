function [txsignal conf] = tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BPSK LSFR PREAMBLE GEN
preamble = zeros(conf.npreamble, 1);
LSFR_state = ones(8, 1);
for i = 1:conf.npreamble
    new = mod(sum(LSFR_state([4 5 6 8])), 2);

    preamble(i) = LSFR_state(end);
    LSFR_state(2:end) = LSFR_state(1:end-1);
    LSFR_state(1) = new;
end
preamble_bpsk = -2*(preamble) + 1;

%%%%%%%%%%% Map 
if (conf.modulation_order == 2)
    tx_symbols = zeros(conf.nbits:2,1);
    tx_symbols(txbits(1:2:end) == 0 & txbits(2:2:end) == 0) =  1/sqrt(2) + 1i/sqrt(2);
    tx_symbols(txbits(1:2:end) == 1 & txbits(2:2:end) == 0) =  1/sqrt(2) - 1i/sqrt(2);
    tx_symbols(txbits(1:2:end) == 1 & txbits(2:2:end) == 1) = -1/sqrt(2) - 1i/sqrt(2);
    tx_symbols(txbits(1:2:end) == 0 & txbits(2:2:end) == 1) = -1/sqrt(2) + 1i/sqrt(2);
end

%%%%% Add preambule
%tx_symbols = cat(1, preamble_bpsk, tx_symbols);
%tx_symbols_upsampled = zeros(length(tx_symbols)*5, 1);
%for ii = 1:length(tx_symbols)
%    tx_symbols_upsampled(5*ii) = tx_symbols(ii);
%end


%%% OFDM
% IFFT
ofdm_ifft_training = osifft(preamble_bpsk, conf.os_factor);

% Add CP to OFDM training
ofdm_ifft_training_cp = cat(1,ofdm_ifft_training(end-(length(ofdm_ifft_training)/2):end), ofdm_ifft_training);

% Add CP
ofdm_tx_signal = zeros(conf.nbcarriers, round(length(tx_symbols_upsampled)/conf.nbcarriers));

% translate data into ofdm

for symbol_idx = 1:length(tx_symbols)/conf.nbcarriers
    ofdm_symbol = osifft(tx_symbols(symbol_idx*conf.nbcarriers-conf.nbcarriers:symbol_idx*conf.nbcarriers), conf.os_factor);
    ofdm_symbol_cp = cat(1,ofdm_tx_signal(end-(length(ofdm_tx_signal)/2):end), ofdm_tx_signal);
    



for ii = 1:length(tx_symbols_upsampled)/conf.nbcarriers 
    ofdm_tx_signal(:,ii) = tx_symbols_upsampled(conf.nbcarriers*(ii-1)+1:conf.nbcarriers*(ii));

end

cp_length = round(size(ofdm_tx_signal,2)/2);
cp = ofdm_tx_signal(:,end-cp_length+1:end);
ofdm_tx_signal_cp = zeros(conf.nbcarriers, size(ofdm_tx_signal,2)+cp_length);


ofdm_tx_signal_cp = cat(2, cp, ofdm_tx_signal);


 


%%% filter
%rolloff = 0.22;
%pulse = rrc( conf.os_factor, rolloff, 20);

%filtered_tx_signal = conv(tx_symbols_upsampled, pulse, 'full');

%  modulation
time = 1:1/length(filtered_tx_signal):4;
tx_signal_modulated = real(filtered_tx_signal * exp(1i*2*pi*conf.f_c*time));

txsignal = tx_signal_modulated;
