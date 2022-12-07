function demodulated_signal = demodulate(modulated_signal, conf)
%Demodulation
    t = 0:1/conf.f_s:(length(modulated_signal)-1)/conf.f_s;
    demodulated_signal= modulated_signal.* exp(-1i*2*pi*(conf.f_c*t'));
end