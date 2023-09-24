"""
PHY407 Lab 05
Question Q2 b) Filtering an Audio Recording

Low pass filters sound from a .wav file and writes it to an output file.
Plots the original and fourier transformed signals, before and after filtering.

Author: Lisa Nasu-Yu, Oct 2021
"""

from scipy.io.wavfile import read, write
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# plot settings
S = 15
L = 15
T = 15

plt.rc('font', size=S)
plt.rc('axes', titlesize=T)
plt.rc('axes', labelsize=S)
plt.rc('xtick', labelsize=S)
plt.rc('ytick', labelsize=S)
plt.rc('legend', fontsize=L)
plt.rc('figure', titlesize=S)

# set plot fonts. runs slower

# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['figure.autolayout'] = True
#
# plt.rcParams.update({"text.usetex":True, "font.family": 'sans-serif', 'font.sans-serif': ['Helvetica']})
# plt.rcParams.update({"text.usetex":True, "font.family": 'serif', 'font.serif': ['Palatino']})


def plot(ft_channel, filtered_ft, channel, filtered, save):
    """
    Plots original and fourier transforms of a wav file, before and after filtering.
    :param ft_channel: [array] fourier transformed data
    :param filtered_ft: [array] filtered fourier transformed data
    :param channel: [array] original data
    :param filtered: [array] filtered data in original space
    :param save: [str] name of file to be saved
    :return: figure with 4 subplots
    """
    fig, axes = plt.subplots(2, 2, figsize=(11, 4))
    axes[0, 0].set_title('Fourier Transform')
    axes[0, 0].set_xlabel(r'f [Hz]')
    axes[0, 0].set_ylabel(r'Amplitude $|c_k|$')
    axes[0, 0].plot(freq, np.abs(ft_channel))

    axes[0, 1].set_title('Filtered Fourier Transform')
    axes[0, 1].set_xlabel(r'f [Hz]')
    axes[0, 1].set_ylabel(r'Amplitude $|c_k|$')
    axes[0, 1].plot(freq, np.abs(filtered_ft))

    axes[1, 0].set_title('Original')
    axes[1, 0].set_ylabel('Amplitude')
    axes[1, 0].set_xlabel('Time [s]')
    axes[1, 0].set_xlim(0, 30e-3)
    axes[1, 0].plot(np.arange(N_Points) / sample, channel)

    axes[1, 1].set_title('Filtered')
    axes[1, 1].set_ylabel('Amplitude')
    axes[1, 1].set_xlabel('Time [s]')
    axes[1, 1].set_xlim(0, 30e-3)
    axes[1, 1].plot(np.arange(N_Points) / sample, filtered)

    plt.tight_layout()
    plt.savefig(save)
    plt.show()


# read in wav file
sample, data = read('graviteaTime.wav')
# separate data into individual arrays
channel_0 = data[:, 0]
channel_1 = data[:, 1]
N_Points = len(channel_0)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
# plot Channel 0
axes[0].set_title('Channel 0')
axes[0].set_ylabel('Amplitude')
axes[0].set_xlabel('Time [s]')
axes[0].plot(np.arange(N_Points) / sample, channel_0)

# plot Channel 1
axes[1].set_title('Channel 1')
axes[1].set_ylabel('Amplitude')
axes[1].set_xlabel('Time [s]')
axes[1].plot(np.arange(N_Points) / sample, channel_1)

plt.tight_layout()
# plt.savefig('original_wav.pdf')
plt.show()

# Compute frequencies of fft
freq = np.fft.rfftfreq(len(channel_0), 1/sample)

# compute fast fourier transforms
ft_channel0 = np.fft.rfft(channel_0)
filtered_ft_0 = np.copy(ft_channel0)
filtered_ft_0[freq > 880] = 0
filtered_0 = np.fft.irfft(filtered_ft_0)

ft_channel1 = np.fft.rfft(channel_1)
filtered_ft_1 = np.copy(ft_channel1)
filtered_ft_1[freq > 880] = 0
filtered_1 = np.fft.irfft(filtered_ft_1)

# make plots
plot(ft_channel0, filtered_ft_0, channel_0, filtered_0, 'filtered_ch0.pdf')
plot(ft_channel1, filtered_ft_1, channel_1, filtered_1, 'filtered_ch1.pdf')

# create array for output data and write to file
data_out = np.empty(data.shape, dtype=data.dtype)
data_out[:, 0] = filtered_0
data_out[:, 1] = filtered_1
write('output_file.wav', sample, data_out)

