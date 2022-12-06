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


%%% OFDM
% Training symbols insertion
tx_symbols =  reshape(tx_symbols,[conf.nbcarriers 4]); %TODO change 4 by a variable
tx_symbols =  [preamble_bpsk, tx_symbols]; %TODO Maybe change this preamble_bpsk

ofdm_symbol = zeros((conf.nbcarriers + conf.cp_length)*conf.os_factor_ofdm,size(tx_symbols,2));
for symbol_idx = 1 : size(tx_symbols, 2)
    ofdm_symbol(conf.os_factor_ofdm*conf.cp_length+1:end,symbol_idx) = osifft(tx_symbols(:,symbol_idx),conf.os_factor_ofdm);
    ofdm_symbol(1:conf.os_factor_ofdm*conf.cp_length,symbol_idx) = ofdm_symbol(end-conf.os_factor_ofdm*conf.cp_length+1:end,symbol_idx);
end

ofdm_symbol = ofdm_symbol(:);


pulse = rrc(conf.os_factor_preambul, conf.rolloff,conf.filterlength);
preamble_filtered = conv(upsample(preamble_bpsk, conf.os_factor_preambul), pulse, 'same');
%Add preambule BPSK to ofdm signal
txsignal = [preamble_filtered/sqrt(mean(abs(preamble_filtered).^2)); ofdm_symbol/sqrt(mean(abs(ofdm_symbol).^2))];

%Modulation
t = 0:1/conf.f_s: (length(txsignal) - 1)/conf.f_s;
txsignal = real(txsignal .* exp(1j*2*pi*conf.f_c*t'));


% plots
figure('name','txsignal');
f = - conf.f_s/2 : conf.f_s/length(txsignal) : conf.f_s/2 - conf.f_s/length(txsignal);
plot(f,abs(fftshift(fft(txsignal))))
grid on
title('Transmitted signal - FFT','interpreter','latex','FontSize',16);
xlabel('frequency/Hz','interpreter','latex','FontSize',16);
ylabel('amplitude','interpreter','latex','FontSize',16);
xline(conf.f_c,'--r'); xline(- conf.f_c,'--r');


