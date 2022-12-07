function modulated_signal = modulate(signal, conf)
%Modulation
    t = 0:1/conf.f_s: (length(signal) - 1)/conf.f_s;
    modulated_signal = real(signal .* exp(1j*2*pi*conf.f_c*t'));
end