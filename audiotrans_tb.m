
function [ber,conf] = audiotrans_tb(conf)






conf.nframes = 1;       % number of frames to transmit
conf.nbits   = conf.nbdatapertrainning*conf.nbcarriers*2      * 4;    % number of bits
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 8000;

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
conf.trainingseq = -2*(randi([0 1],conf.nbcarriers,1)) + 1;



if mod(conf.os_factor_preambul,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end


% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% TODO: To speed up your simulation pregenerate data you can reuse
% beforehand.


% Results


for k=1:conf.nframes
    

        % Generate random data 
        txbits = randi([0 1],conf.nbits,1);


    conf.nbtraining = conf.nbits/ (conf.nbcarriers * conf.modulation_order * conf.nbdatapertrainning) ; % Dont touch this variable
        

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
    
%     wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
 %   audiowrite('out.wav',rawtxsignal,conf.f_s)  
    

       
        SNR = conf.SNR;
        SNRlin = 10^(SNR/10);
        rawrxsignal = rawtxsignal(:,1);
        %sigmaDeltaTheta = 0.0007;
        %theta_n = generate_phase_noise(length(rawrxsignal), sigmaDeltaTheta);
        %apply phase noize
        rawrxsymbol =  demodulate(rawrxsignal, conf);
       % rawrxsymbol = rawrxsymbol.*exp(1j*theta_n);
        rawrxsymbol = rawrxsymbol + sqrt(1/(2*SNRlin)) * (randn(size(rawrxsymbol)) + 1i*randn(size(rawrxsymbol))); 
        rawrxsignal =  modulate(rawrxsymbol, conf);
        rxsignal    = rawrxsignal;
 
 
    %
    % End
    % Audio Transmission   
    % % % % % % % % % % % %
    
    [rxbits conf]       = rx(rxsignal,conf);
    
    res.rxnbits(k)      = length(rxbits);
    res.biterrors(k)    = sum(rxbits ~= txbits);

   
end

per = sum(res.biterrors > 0)/conf.nframes
ber = sum(res.biterrors)/sum(res.rxnbits)