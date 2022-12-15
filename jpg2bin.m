close all;
clear variables;
clc;

image = dec2bin(imread('Matterhorn.jpg'),8);
txbits_char = image(1,:);
for pix = 2:size(image,1)
    txbits_char = [txbits_char image(pix,:)];
end
txbits = str2num(txbits_char(1));
for bit = 2:length(txbits_char)
    txbits = [txbits str2num(txbits_char(bit))];
end

Matterhorn_bin = txbits;

save("Matterhorn.mat","Matterhorn_bin");