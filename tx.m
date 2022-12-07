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


% Map data to QPSK symbols
tx_symbols = mapper(txbits,conf);

%%% OFDM Part
% Add training blocks
tx_symbols =  reshape(tx_symbols,[conf.nbcarriers conf.nbdatapertrainning*conf.nbtraining]);
for block = 0:conf.nbtraining-1
    if block == 0
        tx_symbols =  [conf.trainingseq, tx_symbols];
    else
        tx_symbols = [tx_symbols(:,1:block*conf.nbdatapertrainning+block-1) conf.trainingseq tx_symbols(:,block*conf.nbdatapertrainning+block:end)];
    end
end
%tx_symbols =  [conf.trainingseq, tx_symbols];

% Map into OFDM and add CP
ofdm_symbol = zeros((conf.nbcarriers + conf.cp_length)*conf.os_factor_ofdm,size(tx_symbols,2));
for symbol_idx = 1 : size(tx_symbols, 2)
    % Map data symbol to OFDM symbol
    ofdm_symbol(conf.os_factor_ofdm*conf.cp_length+1:end,symbol_idx) = osifft(tx_symbols(:,symbol_idx),conf.os_factor_ofdm);
    % Add CP
    ofdm_symbol(1:conf.os_factor_ofdm*conf.cp_length,symbol_idx) = ofdm_symbol(end-conf.os_factor_ofdm*conf.cp_length+1:end,symbol_idx);
end

ofdm_symbol = ofdm_symbol(:);

% Mach filter for single carrier preamble
preamble_filtered = matched_filter(upsample(conf.preamble, conf.os_factor_preambul), conf);

%Add BPSK single carrier preambule to ofdm signal
txsignal = [preamble_filtered/sqrt(mean(abs(preamble_filtered).^2)); ofdm_symbol/sqrt(mean(abs(ofdm_symbol).^2))];

%Modulation
tx_non_modulated = txsignal;
txsignal = modulate(txsignal, conf);



% plots
if (conf.plotfig == 1)
    figure;
    f = - conf.f_s/2 : conf.f_s/length(txsignal) : conf.f_s/2 - conf.f_s/length(txsignal);
    plot(f,abs(fftshift(fft(txsignal))))
    grid on
    title('Spectre of transmitted signal');
    xlabel('frequency/Hz');
    ylabel('Amplitude');
end


