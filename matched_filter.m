function filtered_signal = matched_filter(signal, conf)
% Create a root-raised cosine filter and filter the signal with it.

rolloff_factor = 0.22;

h = rrc(conf.os_factor_ofdm, conf.rolloff, conf.filterlength);
filtered_signal = conv(h, signal);