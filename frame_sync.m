function [beginning_of_data, phase_of_peak] = frame_sync(rx_signal, conf)

% Frame synchronizer.
% rx_signal is the noisy received signal, and L is the oversampling factor (L=1 in chapter 2, L=4 in all later chapters).
% The returned value is the index of the first data symbol in rx_signal.

if (rx_signal(1) == 0)
    warning('Signal seems to be noise-free. The frame synchronizer will not work in this case.');
    
end

detection_threshold = 15;

% Calculate the frame synchronization sequence and map it to BPSK: 0 -> +1, 1 -> -1
%frame_sync_sequence = 1 - 2*lfsr_framesync(frame_sync_length);
frame_sync_length = conf.npreamble;
frame_sync_sequence = conf.preamble;

% When processing an oversampled signal (L>1), the following is important:
% Do not simply return the index where T exceeds the threshold for the first time. Since the signal is oversampled, so will be the
% peak in the correlator output. So once we have detected a peak, we keep on processing the next L samples and return the index
% where the test statistic takes on the maximum value.
% The following two variables exist for exactly this purpose.
current_peak_value = 0;
L = conf.os_factor_preambul;
samples_after_threshold = L;
for i = L * frame_sync_length + 1 : length(rx_signal)
    r = rx_signal(i - L * frame_sync_length : L : i - L); % The part of the received signal that is currently inside the correlator.
    c = frame_sync_sequence' * r;
    T = abs(c)^2 / abs(r' * r);
    
    if (T > detection_threshold || samples_after_threshold < L)
        samples_after_threshold = samples_after_threshold - 1;
        if (T > current_peak_value)
            beginning_of_data = i;
            % TODO
            
            phase_of_peak = angle(c);
            current_peak_value = T;
        end
        if (samples_after_threshold == 0)
            return;
        end
    end
end


warning('No synchronization sequence found.');
beginning_of_data = 1;
phase_of_peak = 0;
return