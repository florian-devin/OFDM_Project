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

% plots rx spectre
% plotdata.figrx.spt.x     = - conf.f_s/2 : conf.f_s/length(rxsignal) : conf.f_s/2 - conf.f_s/length(rxsignal);
% plotdata.figrx.spt.y     = abs(fftshift(fft(rxsignal)));
% plotdata.figrx.spt.title = 'Spectre of received signal';
% plotdata.figrx.spt.xlabel= 'frequency/Hz';
% plotdata.figrx.spt.ylabel= 'Amplitude';


%Demodulate
rxsignal = demodulate(rxsignal, conf);

% OFDM lowpass filter
f_corner = 1.05*conf.bandwidth; % The size of the filter is a little bit biger than the signal windows (5% here)
rxsignal_filtered = 2 * ofdmlowpass(rxsignal, conf, f_corner);

% Frame syncronisation
%TODO remplace by this line : 
start_idx = frame_sync(rxsignal_filtered, conf);
%start_idx = detector(conf.preamble, rxsignal_filtered, 15);
%start_idx = 49065;
%rx_size = (conf.nbdatapertrainning+1) * conf.os_factor_ofdm * conf.nbcarriers * (conf.cp_length / conf.nbcarriers + 1);
rx_size = (conf.nbcarriers*conf.os_factor_ofdm*conf.nbtraining*(conf.nbdatapertrainning+1)*(1+conf.cp_length/conf.nbcarriers));
rxsignal_filtered = rxsignal_filtered(start_idx:start_idx+rx_size-1);

rxsignal_filtered = reshape(rxsignal_filtered, [conf.os_factor_ofdm*(conf.nbcarriers+conf.cp_length), (conf.nbdatapertrainning+1)*conf.nbtraining]);

% Remove all CPs
rxsignal_filtered(1:conf.os_factor_ofdm*conf.cp_length,:) = [];

% OFDM demodulation
for data_idx = 1:(conf.nbdatapertrainning+1)*conf.nbtraining % +1 because of the training 
    rxsymbol(:,data_idx) = osfft(rxsignal_filtered(:,data_idx), conf.os_factor_ofdm);
end


if (conf.plotfig == 1)
    rxsymbolplot = rxsymbol;
    for block = 0:conf.nbtraining-1
            block_idx = block*conf.nbdatapertrainning+1;
            rxsymbolplot(:,block_idx) = []; %remove training seq
    end
    figure;
    plot(real(rxsymbolplot),imag(rxsymbolplot), 'o');
    title('Received Symbols before correction');
end

% Channel estimation
switch conf.estimationtype
    case 'none'
        for block = 0:conf.nbtraining-1
            block_idx = block*conf.nbdatapertrainning+1;
            rxsymbol(:,block_idx) = []; %remove training seq
        end
        rxdata = rxsymbol;
    case 'block'
        for block = 0:conf.nbtraining-1
            block_idx = block*conf.nbdatapertrainning+1;
            H_hat        = rxsymbol(:,block_idx)./conf.trainingseq;
            theta_hat    = mod(angle(H_hat), 2*pi);
            rxsymbol(:,block_idx) = []; %remove training seq
            for data_idx = block*conf.nbdatapertrainning+1:(block+1)*conf.nbdatapertrainning
                %amplitude correction
                rxdata(:,data_idx) = rxsymbol(:,data_idx)./abs(H_hat);
                %phase correction
                rxdata(:,data_idx) = rxdata(:,data_idx) .* exp(-1i*theta_hat);
            end
        end
    case 'viterbi'
        rxdata = [];
        for block = 0:conf.nbtraining-1

            block_idx = block*(conf.nbdatapertrainning+1)+1;

            H_hat                       = rxsymbol(:,block_idx)./conf.trainingseq;
            theta_hat(:,block_idx)      = mod(angle(H_hat), 2*pi);

            for data_idx = block*(conf.nbdatapertrainning+1)+2:(block+1)*(conf.nbdatapertrainning+1)

                    deltaTheta = (1/4*angle(-rxsymbol(:, data_idx).^4) + pi/2*(-1:4));
                    [~, ind] = min(abs(deltaTheta - theta_hat(:, data_idx-1)),[],2);
                    indvec = (0:conf.nbcarriers-1).*6 + ind'; 
                    deltaTheta = deltaTheta';
                    theta = deltaTheta(indvec);
                    theta_hat(:, data_idx) = mod(0.01*theta' + 0.99*theta_hat(:, data_idx-1), 2*pi);

                    rxsymbol(:,data_idx) = rxsymbol(:,data_idx)./abs(H_hat);
                    rxsymbol(:,data_idx) = rxsymbol(:,data_idx) .* exp(-1i*theta_hat(:,data_idx));


            end
            rxdata = [rxdata rxsymbol(:, block*(conf.nbdatapertrainning+1)+2:(block+1)*(conf.nbdatapertrainning+1))];
           %rxdata(:,block_idx) = []; %remove training seq
            
        end
end

rxdata = rxdata(:);

rxdata = rxdata / sqrt(mean(abs(rxdata).^2));

if (conf.plotfig == 1)
    figure;
    plot(real(rxdata),imag(rxdata), 'o');
    title('Received Symbols after correction');
end

rxbits = demapper(rxdata);


