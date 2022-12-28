import matplotlib.pyplot as plt
import math

x_cp = [128, 64, 32, 16, 8]
x_SNR = [ 20, 15, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3]

y_SNR = [ 0.00000000000000000000000001, 0.000000001, 0.00004, 0.0005, 0.002, 0.004, 0.01, 0.02, 0.035, 0.05, 0.08, 0.11, 0.15, 0.2, 0.22, 0.25, 0.33, 0.5]
y_SNR_log =  [math.log(snr, 10) for snr in y_SNR]

x_phi = [0.00001, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004]
y_phi = [0.00000000000000000001,0.00000000000000000001,0.00000000000000000001,0.00000000000000000001,0.00000000000000000001,0.00000000000000000001,0.00000000000000000001,0.00000000000000000001,0.02,0.1,0.12,0.17,0.31,0.48]
y_phi_log =  [math.log(phi, 10) for phi in y_phi]
x_phi_mrad = [phi*1000 for phi in x_phi]
y_cp = [0.025, 0.05, 0.08, 0.1, 0.11]


x_nbsub = [5, 2.5, 1.25, 0.75, 0.375, 0.15625]
y_nbsub = [0, 0, 0.0253, 0.15, 0.20945, 0.5]

x_nb_data = [8, 16, 32, 64, 128, 256, 512, 1024]
y_nb_data_block = [0.04, 0.1, 0.16, 0.2, 0.25, 0.32, 0.43, 0.5]
y_nb_data_viterbi = [0.06, 0.0865, 0.119, 0.13, 0.25, 0.3, 0.33, 0.45]

"""
plt.plot(x_nb_data,y_nb_data_block, label='block')
plt.plot(x_nb_data,y_nb_data_viterbi, label='viterbi')
plt.legend(loc='upper left')
# Ajout d'un titre et des étiquettes aux axes
plt.title('Training spacing variation')
plt.xlabel('Number of data block per training')
plt.ylabel('BER')

plt.plot(x_nbsub,y_nbsub)
# Ajout d'un titre et des étiquettes aux axes
plt.title('Cariers spacing variation')
plt.xlabel('Cariers spacing [Hz]')
plt.ylabel('BER')
"""
plt.plot(x_phi_mrad,y_phi_log)
# Ajout d'un titre et des étiquettes aux axes
plt.title('Phase-shift variation')
plt.xlabel('Phase-shift [milirad]')
plt.ylabel('BER')

"""
plt.plot(x_SNR,y_SNR_log)
# Ajout d'un titre et des étiquettes aux axes
plt.title('SNR variation')
plt.xlabel('SNR [dB]')
plt.ylabel('log(BER)')



plt.plot(x_cp, y_cp)
# Ajout d'un titre et des étiquettes aux axes
plt.title('CP Length variation')
plt.xlabel('CP Length')
plt.ylabel('BER')
"""
# Affichage du graphique
plt.show()