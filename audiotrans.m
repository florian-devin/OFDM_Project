% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal

%TODO : make it work for multiple training sequence
%TODO : add better estimation of the chanel
%TODO : add random bits to complete all cariers


clear variables;
close all;
clc;
plotdata.plotfig = 1;
plotdata.figrx = [];
plotdata.figtx = [];

% Configuration Values
conf.audiosystem = 'bypass'; % Values: 'matlab','native','bypass','bypass2'
conf.estimationtype = 'viterbi'; % For chanel estimation and correction : 'none', 'block', 'viterbi'
conf.plotfig = 1;


% OFDM 
conf.nbcarriers = 256;
conf.carriersSpacing = 5; % Hz
conf.cp_length = conf.nbcarriers / 2;
conf.bandwidth = ceil((conf.nbcarriers + 1)/ 2)*conf.carriersSpacing;
conf.nbdatapertrainning = 10;

conf.f_s     = 48000;   % sampling rate  
conf.f_sym   = 100;     % symbol rate
% for rx filter (preamble filter)
conf.rolloff = 0.22;
conf.filterlength = 20;

conf.nframes = 1;       % number of frames to transmit
conf.nbits   = conf.nbdatapertrainning*conf.nbcarriers*2      * 3;    % number of bits thes last 2 is for 2 training block insertion
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 10000;

conf.npreamble  = 256;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;



% Init Section
% all calculations that you only have to do once
conf.os_factor_ofdm  = conf.f_s/(conf.carriersSpacing*conf.nbcarriers); % given by >> help osifft
%conf.os_factor_preambul = conf.f_s/conf.f_sym; % os_factor for BPSK preambule (single carrier)
conf.os_factor_preambul = 4;
% Preamble generation
conf.preamble =  -2*(preamble_generate(conf.npreamble)) + 1; % BPSK (-1 or 1)
% Training generation
%preamble_generate generate a random sequence perfect for trainingseq
%conf.trainingseq = -2*(preamble_generate(conf.nbcarriers)) + 1; % BPSK (-1 or 1)
conf.trainingseq = -2*(randi([0 1],conf.nbcarriers,1)) + 1;

conf.nbtraining = conf.nbits/ (conf.nbcarriers * conf.modulation_order * conf.nbdatapertrainning) ; % Dont touch this variable

if mod(conf.os_factor_preambul,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end
conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% TODO: To speed up your simulation pregenerate data you can reuse
% beforehand.


% Results


for k=1:conf.nframes
    
    % Generate random data
    txbits = randi([0 1],conf.nbits,1);
    
    % TODO: Implement tx() Transmit Function
    [txsignal conf plotdata] = tx(txbits,conf,k,plotdata);

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
    
%     wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
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
        SNR = 1;
        SNRlin = 10^(SNR/10);
        rawrxsignal = rawtxsignal(:,1);
        %rawrxsignal = rawrxsignal + sqrt(1/(2*SNRlin)) * (randn(size(rawrxsignal)) + 1i*randn(size(rawrxsignal))); 
        sigmaDeltaTheta = 0.004;
        theta_n = generate_phase_noise(length(rawrxsignal), sigmaDeltaTheta);
        %apply phase noize
        rawrxsymbol =  demodulate(rawrxsignal, conf);
        rawrxsymbol = rawrxsymbol.*exp(1j*theta_n);
        rawrxsymbol = rawrxsymbol + sqrt(1/(2*SNRlin)) * (randn(size(rawrxsymbol)) + 1i*randn(size(rawrxsymbol))); 
        rawrxsignal =  modulate(rawrxsymbol, conf);
        rxsignal    = rawrxsignal;
    end
    
    % Plot received signal for debugging
    if (conf.plotfig == 1)
        figure;
        plot(rxsignal);
        title('Received Signal')
    end
    %
    % End
    % Audio Transmission   
    % % % % % % % % % % % %
    
    [rxbits conf]       = rx(rxsignal,conf);
    
    res.rxnbits(k)      = length(rxbits);
    res.biterrors(k)    = sum(rxbits ~= txbits);

    %%% Plot section
    figure(4);
    subplot(231);
    grid on;
    plot(plotdata.figtx.spt.x,plotdata.figtx.spt.y);
    title(plotdata.figtx.spt.title);
    xlabel(plotdata.figtx.spt.xlabel);
    ylabel(plotdata.figtx.spt.ylabel);

    subplot(232);
    hold on;
    plot(plotdata.figtx.symb.xpreamble,plotdata.figtx.symb.ypreamble,'yo');
    plot(plotdata.figtx.symb.xtraining,plotdata.figtx.symb.ytraining,'g+');
    plot(plotdata.figtx.symb.xdata    ,plotdata.figtx.symb.ydata    ,'bo');
    title(plotdata.figtx.symb.title);
    xlabel(plotdata.figtx.symb.xlabel);
    ylabel(plotdata.figtx.symb.ylabel);

    subplot(233);
    plot(rawtxsignal);
    title('Transmited signal');
    xlabel('Sample');
    ylabel('Amplitude');

    subplot(234);
    plot(rxsignal);
    title('Received signal');
    xlabel('Sample');
    ylabel('Amplitude');

    subplot(235);
%     grid on;
%     plot(plotdata.figrx.spt.x,plotdata.figrx.spt.y);
%     title(plotdata.figrx.spt.title);
%     xlabel(plotdata.figrx.spt.xlabel);
%     ylabel(plotdata.figrx.spt.ylabel);


end

per = sum(res.biterrors > 0)/conf.nframes
ber = sum(res.biterrors)/sum(res.rxnbits)
