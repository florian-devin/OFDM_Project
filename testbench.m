

%csvwrite(csv_filename, header, A);



% Configuration Values
conf.audiosystem = 'bypass2'; % Values: 'matlab','native','bypass','bypass2'
conf.estimationtype = 'block'; % For chanel estimation and correction : 'none', 'block', 'viterbi'
conf.datatype = 'random'; % Send an image or random bits : 'image', 'random' 
conf.plotfig = 0;


% OFDM 
conf.nbcarriers = 256;
conf.carriersSpacing = 5; % Hz
conf.cp_length = conf.nbcarriers / 2;
conf.bandwidth = ceil((conf.nbcarriers + 1)/ 2)*conf.carriersSpacing;
conf.nbdatapertrainning = 30;


conf.f_s     = 48000;   % sampling rate  
conf.f_sym   = 100;     % symbol rate
% for rx filter (preamble filter)
conf.rolloff = 0.22;
conf.filterlength = 20;

nb_run = 10;

csv_filename = "result_ber_snr.csv";
header = {'run' ,'SNR' 'delta_theta' ,'nbcarriers' ,'carriersSpacing' ,'cp_length' ,'correction methode' ,'BER'};


A = zeros(nb_run,8);
SNR_max = 2;
SNR_min = -2;
SNR_range = logspace(SNR_min,SNR_max,nb_run);
for run = 1:nb_run
    run
    conf.SNR = SNR_range(run);
    S_ber = 0;
    for i = 1:5
        [ber, conf] = audiotrans_tb(conf);
        S_ber = S_ber+ber/5;
    end
    A(run,:) = [run SNR_range(run) 0 256 5 128 1 S_ber];
end


T = array2table(A);
T.Properties.VariableNames(1:8) ={'run','SNR','delta_theta','nbcarriers','carriersSpacing','cp_length','correction_methode','BER'};
writetable(T, csv_filename);