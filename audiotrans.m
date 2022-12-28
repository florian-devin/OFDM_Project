% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   4 operating modes:
%   - 'matlab'  : generic MATLAB audio routines (unreliable under Linux)
%   - 'native'  : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass'  : no audio transmission, takes txsignal as received signal
%   - 'bypass2' : white noize and phase shift is added to the signal by matlab

%TODO : make it work for multiple training sequence
%TODO : add better estimation of the chanel
%TODO : add random bits to complete all cariers


clear variables;
close all;
clc;

% Configuration Values
conf.audiosystem = 'bypass2'; % Values: 'matlab','native','bypass','bypass2'
conf.estimationtype = 'viterbi'; % For chanel estimation and correction : 'none', 'block', 'viterbi'
conf.datatype = 'image'; %type of transmited data : 'image' 'random'
conf.plotfig = 1;


% OFDM 
conf.nbcarriers = 512;
conf.carriersSpacing = 2; % Hz
conf.cp_length = 128;
conf.bandwidth = ceil((conf.nbcarriers + 1)/ 2)*conf.carriersSpacing;
conf.nbdatapertrainning = 64;

conf.f_s     = 48000;   % sampling rate  
conf.f_sym   = 100;     % symbol rate
% for rx filter (preamble filter)
conf.rolloff = 0.22;
conf.filterlength = 20;

conf.nframes = 1;       % number of frames to transmit
conf.nbits   = conf.nbdatapertrainning*conf.nbcarriers*2      * 5;    % number of bits
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 6000;

conf.npreamble  = 256;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;



% Init Section
conf.os_factor_ofdm  = conf.f_s/(conf.carriersSpacing*conf.nbcarriers); % given by >> help osifft
conf.os_factor_preambul = 4;
% Preamble generation
conf.preamble =  -2*(preamble_generate(conf.npreamble)) + 1; % BPSK (-1 or 1)
% Training generation
conf.trainingseq = -2*(randi([0 1],conf.nbcarriers,1)) + 1; % BPSK (-1 or 1)



if mod(conf.os_factor_preambul,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end


% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);



for k=1:conf.nframes

    if strcmp(conf.datatype,'random')
        % Generate random data
        txbits = randi([0 1],conf.nbits,1);
    elseif strcmp(conf.datatype,'image')
        % Load image
        load("Matterhorn.mat");
        txbits = Matterhorn_bin;

        % Add random bits at the end
        [txbits conf] = add_random_bit(txbits,conf);
        nb_rdm_bits = txbits(1:32);
    end
    conf.nbtraining = conf.nbits/ (conf.nbcarriers * conf.modulation_order * conf.nbdatapertrainning) ; % Dont touch this variable
    conf.nsyms      = ceil(conf.nbits/conf.modulation_order);
    

    % Transmit Function
    [txsignal conf] = tx(txbits,conf,k);


    % % % % % % % % % % % %
    % Begin
    % Audio Transmission
    %
    
    % normalize values
    peakvalue       = max(abs(txsignal));
    normtxsignal    = txsignal / (peakvalue + 0.3);

    
    % create vector for transmission
    rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
    rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
    txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
    
if (conf.plotfig == 1)
    figure(4);
    subplot(244);
    plot(rawtxsignal(:,1));
    title('Transmited signal');
    xlabel('Sample');
    ylabel('Amplitude');
end   
    %wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
    audiowrite('out.wav',rawtxsignal,conf.f_s)  
    
    % Platform native audio mode 
    if strcmp(conf.audiosystem,'native')
        
        % Windows WAV mode 
        if ispc()
            disp('Windows WAV');
            wavplay(rawtxsignal,conf.f_s,'async');
            disp('Recording in Progress');
            rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
            disp('Recording complete')
            rxsignal = rawrxsignal(1:end,1);

        % ALSA WAV mode 
        elseif isunix()
            disp('Linux ALSA');
            cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
            system(cmd); 
            disp('Recording in Progress');
            system('aplay  out.wav')
            pause(2);
            disp('Recording complete')
            rawrxsignal = audioread('in.wav');
            rxsignal    = rawrxsignal(1:end,1);
        end
        
    % MATLAB audio mode
    elseif strcmp(conf.audiosystem,'matlab')
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
        recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
        record(recobj);
        disp('Recording in Progress');
        playblocking(playobj)
        pause(0.5);
        stop(recobj);
        disp('Recording complete')
        rawrxsignal  = getaudiodata(recobj,'int16');
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
        
    elseif strcmp(conf.audiosystem,'bypass')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal;
    elseif strcmp(conf.audiosystem, 'bypass2')
        SNR = 30;
        SNRlin = 10^(SNR/10);
        rawrxsignal = rawtxsignal(:,1);
        sigmaDeltaTheta = 0.0005;
        theta_n = generate_phase_noise(length(rawrxsignal), sigmaDeltaTheta);
        %apply phase-shift and white noize
        rawrxsymbol =  demodulate(rawrxsignal, conf);
        rawrxsymbol = rawrxsymbol.*exp(1j*theta_n);
        rawrxsymbol = rawrxsymbol + sqrt(1/(2*SNRlin)) * (randn(size(rawrxsymbol)) + 1i*randn(size(rawrxsymbol))); 
        rawrxsignal =  modulate(rawrxsymbol, conf);
        rxsignal    = rawrxsignal;
    end
    
    % Plot received signal for debugging
    if (conf.plotfig == 1)
        figure(4);
        subplot(248);
        plot(rxsignal);
        title('Received Signal');
        xlabel('Sample');
        ylabel('Amplitude');
    end   
    %
    % End
    % Audio Transmission   
    % % % % % % % % % % % %
    
    [rxbits conf]       = rx(rxsignal,conf);
    
    res.rxnbits(k)      = length(rxbits);
    res.biterrors(k)    = sum(rxbits ~= txbits);

    if strcmp(conf.datatype,'image')
        % decode image
        % remove random bits
        %nb_random_rx_bits = bi2de(rxbits(1:32)','left-msb');
        % Use transmitted variable to avoid error during demonstration
        nb_random_rx_bits = bi2de(nb_rdm_bits','left-msb');
        %nb_random_rx_bits = bi2de(rxbits(1:32)','left-msb');
        payload_data = rxbits(33:end - nb_random_rx_bits);
        image_size = [128 128];
        figure(5);
        image = image_decoder(payload_data, image_size);
        imshow(image/255);
        title("Recovered Image")
    end




end

per = sum(res.biterrors > 0)/conf.nframes
ber = sum(res.biterrors)/sum(res.rxnbits)
