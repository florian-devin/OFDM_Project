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
index = 1;
for i = L * frame_sync_length + 1 : length(rx_signal)
    r(:,index) = rx_signal(i - L * frame_sync_length : L : i - L); % The part of the received signal that is currently inside the correlator.
    c(:,index) = frame_sync_sequence' * r(:,index);
    T(:,index) = abs(c(:,index))^2 / abs(r(:,index)' * r(:,index));
    
    if (T(:,index) > detection_threshold || samples_after_threshold < L)
        samples_after_threshold = samples_after_threshold - 1;
        if (T(:,index) > current_peak_value)
            beginning_of_data = i;
            % TODO
            
            phase_of_peak = angle(c(:,index));
            current_peak_value = T(:,index);
        end
        if (samples_after_threshold == 0)
            return;
        end
    end
    index = index +1;
end

if (conf.plotfig == 1)
    figure;
    plot(L * frame_sync_length + 1 : length(rx_signal),c(:,index));
    title('sync');
end

error('No synchronization sequence found.');
return