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


%Demodulate
rxsignal = demodulate(rxsignal, conf);

% OFDM lowpass filter % TODO understand why we use a lowpass filter
f_corner = 1.05*conf.bandwidth; % The size of the filter is a little bit biger than the signal windows (5% here)
rxsignal_filtered = 2 * ofdmlowpass(rxsignal, conf, f_corner);

% Frame syncronisation
%TODO remplace by this line : start_idx = frame_sync(rxsignal_filtered, conf);
start_idx = conf.f_s+1 + conf.npreamble;
rx_size = (conf.nbdatapertrainning+1) * conf.os_factor_ofdm * conf.nbcarriers * (conf.cp_length / conf.nbcarriers + 1);
rxsignal_filtered = rxsignal_filtered(start_idx:start_idx+rx_size-1);

rxsignal_filtered = reshape(rxsignal_filtered, [conf.os_factor_ofdm*(conf.nbcarriers+conf.cp_length), conf.nbdatapertrainning+1]);

% Remove all CPs
rxsignal_filtered(1:conf.os_factor_ofdm*conf.cp_length,:) = [];

% OFDM demodulation
for data_idx = 1:conf.nbdatapertrainning+1 % +1 because of the training 
    rxsymbol(:,data_idx) = osfft(rxsignal_filtered(:,data_idx), conf.os_factor_ofdm);
end
disp('Rx first QPSK symbol')
rxsymbol(1,2)

if (conf.plotfig == 'yes')
    figure;
    plot(real(rxsymbol(:,1)),imag(rxsymbol(:,1)), 'o');
    title('Received training Symbols (BPSK)');
end

% Channel estimation
switch conf.estimationtype
    case 'none'
        rxsymbol(:,1) = []; %remove training seq
        rxdata = rxsymbol;
    case 'block'
        H_hat = rxsymbol(:,1)./conf.trainingseq;
        theta_hat = mod(angle(H_hat), 2*pi);
        rxsymbol(:,1) = []; %remove training seq
        for data_idx = 1:conf.nbdatapertrainning
            %amplitude correction
            rxdata(:,data_idx) = rxsymbol(:,data_idx)./abs(H_hat);
            %phase correction
            rxdata(:,data_idx) = rxdata(:,data_idx) .* exp(-1i*theta_hat);
        end
end

rxdata = rxdata(:);

if (conf.plotfig == 'yes')
    figure;
    plot(real(rxdata),imag(rxdata), 'o');
    title('Received Symbols');
end

rxbits = demapper(rxdata);


