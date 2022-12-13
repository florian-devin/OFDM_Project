function theta_n = generate_phase_noise(length_of_noise, sigmaDeltaTheta)
    % Create phase noise
    theta_n = zeros(length_of_noise,1);
    %% TODO
    theta_n(1) = rand(1);
    for n=2:length_of_noise
        theta_n(n) = theta_n(n-1) + randn(1,1)*sigmaDeltaTheta;
    end
end