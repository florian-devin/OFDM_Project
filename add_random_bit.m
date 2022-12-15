function [txsignal conf] = add_random_bit(txbits,conf)
    bit_per_tram = conf.nbdatapertrainning*conf.nbcarriers*2;
    conf.nbits = conf.nbdatapertrainning*conf.nbcarriers*2  * (floor(length(txbits)/bit_per_tram)+1);
    nb_rdm_bit_char = dec2bin(conf.nbits-length(txbits)-32,32);
    nb_rdm_bit = str2num(nb_rdm_bit_char(1));
    for bit = 2:length(nb_rdm_bit_char)
        nb_rdm_bit = [nb_rdm_bit str2num(nb_rdm_bit_char(bit))];
    end
    txbits = [nb_rdm_bit txbits];
    txbits = txbits(:);
    txbits = [txbits ; randi([0 1],conf.nbits-length(txbits),1)];
    txsignal = txbits;
end