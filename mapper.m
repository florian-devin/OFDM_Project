function symbols = mapper(bits, conf)
% Convert bit vector into QPSK symbols.
    if (conf.modulation_order == 2)
        symbols = zeros(conf.nbits:2,1);
        symbols(bits(1:2:end) == 0 & bits(2:2:end) == 0) =  1/sqrt(2) + 1i/sqrt(2);
        symbols(bits(1:2:end) == 1 & bits(2:2:end) == 0) =  1/sqrt(2) - 1i/sqrt(2);
        symbols(bits(1:2:end) == 1 & bits(2:2:end) == 1) = -1/sqrt(2) - 1i/sqrt(2);
        symbols(bits(1:2:end) == 0 & bits(2:2:end) == 1) = -1/sqrt(2) + 1i/sqrt(2);
    end
end