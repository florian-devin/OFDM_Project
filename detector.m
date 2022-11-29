function [start] = detector(p,r, thr)
% Input :   p: preamble (shape = (100, 1)), r: received signal (shape = (Nr, 1)), thr: scalar threshold 
% output:   start: signal start index
Np = size(p,1);
Nr = size(r,1);
c = zeros(Nr-Np, 1);
c_norm = zeros(Nr-Np, 1);
    for n = 1:Nr-Np
        % c: complex value
        c(n, 1) = sum(p.*r(n : n+ (Np-1), 1));
        demominator = sum(abs(r(n  : n + (Np-1), 1)).^2);
        c_norm(n, 1) = abs(c(n, 1))^2/demominator;
        if c_norm(n, 1) > thr
            start = n+100;
            return
        end
    end
disp('Frame start not found.')
start = -1;
end

