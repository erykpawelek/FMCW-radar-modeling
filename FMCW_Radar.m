clc;
clear;
close all;

% #### PART 1 MANUAL DEFINITION OF SIGNALS ####
               
c = 299792458;  % Light speed const in m/s
fc = 77e9;      % Carrier frequency (77 GHz - standard for automotive radars)

% ---- Project requirements -----

max_range = 200;                        % Maximal range in meters
range_resolution = 0.5;                 % Minimal distinguishable difference in distance in meters
max_velocity = 150;                     % Maximal velocity of measured target in km/h
bandwidth_ = c / (2*range_resolution);  % Bandwidth of our signal 
t_of_flight = 2*max_range/c;            % Time required to signal travel both ways 

% Value below is increased by some factor in order to obtain sufficient
% amount of time to compare it to beat signal
t_sweep = 10 * t_of_flight;             % Duration of 1 chirp 
slope = bandwidth_/t_sweep;             % Frequency slope of chirp
fs = 10*bandwidth_;                     % Sampling frequency

% Ploting parameters
disp(['System parameters        ', 'value']);
disp('------------------------------------');
disp(['Carrier frequency :      ', num2str(fc/1e9), ' GHz']);
disp(['Maximum target range:    ', num2str(max_range), ' m']);
disp(['Maximum target velocity: ', num2str(max_velocity), ' km/h']);
disp(['Chirp time :             ', num2str(t_sweep*1e6), ' us']);
disp(['Bandwidth :              ', num2str(bandwidth_/1e9), ' GHz']);
disp(['Chirp time (t_sweep):    ', num2str(t_sweep*1e6), ' us']);
disp(['Slope:                   ', num2str(slope/1e12), ' THz/s']);
disp(['Sampling frequency:      ', num2str(fs/1e9), ' GHz']);
disp(' ');

% ---- Target definition ---- 
target_range_0 = 100;   % In meters
target_vel_0 = 0;       % In km/h

target_range_1 = target_range_0;   % In meters
target_vel_1 = 100;     % In km/h

% ---- Signal definition ----

t = 0: (1/fs) : t_sweep;             % Time vector for one chirp
phase_tx = (2*pi*slope*t.^2)/2;      % Phase of one chirp
signal_tx = exp(1j * phase_tx);      % Definition of transitted signal

% % Diagnostic print
% figure;
% plot(t(1:500), real(signal_tx(1:500)));
% title('Begining of transmition signal (Re)');
% xlabel('Time (s)');
% ylabel('Amplittude');
freq_function = (fc + slope * t)/1e9;
figure;
plot(t, freq_function);
title('Frequency function of one chirp');
xlabel('Time (s)');
ylabel('Freq (GHz)');
grid on;
% ---- Signal propagation and echo ----

%Calculation for our first target (target_0)
tau = (2 * target_range_0) / c;                                         % This is our signal travel time (round-trip)
lambda = c / fc;                                                        % wavelenght
f_d_0 = (2 * target_vel_0 * 1000 / 3600) / lambda;                      % Doppler frequency 
phase_rx = (2 * pi * slope * (t-tau).^2) / 2;                           % Reciving signal phase
doppler_shift_signal_0 = exp(1j * 2 * pi * f_d_0 * t);                  % Generationg doppler shift signal
noise = 0.5 * (randn(size(t)) + 1j * randn(size(t)));                   % Generating complex gausian noise 
signal_rx_0 = (exp(1j * phase_rx) .* doppler_shift_signal_0) + noise;   % Return signal with noise and phase shift

%Calculation for our first target (target_1)
f_d_1 = (2 * target_vel_1 * 1000 / 3600) / lambda;                      % Doppler frequency 
doppler_shift_signal_1 = exp(1j * 2 * pi * f_d_1 * t);                  % Generationg doppler shift signal
signal_rx_1 = (exp(1j * phase_rx) .* doppler_shift_signal_1) + noise;   % Return signal with noise and phase shift


% We multiply signals to make use of exponent adding properiety of a
% signals
disp('Comparison od stationary and non stationary target');
disp('------------------------------------');
disp(['Travel time for target_0:              ', num2str(tau), ' s']);
disp(['Doppler freqency for target_0 (vel=0): ', num2str(f_d_0/1e6), ' MHz']);
disp(['Travel time for target_1:              ', num2str(tau), ' s']);
disp(['Doppler freqency for target_1 (vel>0): ', num2str(f_d_1/1e6), ' MHz']);
% ---- Beat signal ----

beat_signal_0 = signal_tx .* conj(signal_rx_0);                         % Generating beat signal for target_0
beat_signal_1 = signal_tx .* conj(signal_rx_1);                         % Generating beat signal for target_1
% Beat signal is generated using Complex Conjugate function, by that we
% ramowe our tx component of recived signal (we subtract it)

% ---- printing results ----
figure;
plot(t, real(beat_signal_0));
title('Beat Signal');
xlabel('Time (s)');
ylabel('Amplitude (Re)');
grid on;
figure;
subplot(3,1,1);
plot(t(1:4000), real(signal_tx(1:4000)),'-g');
xlabel('Time (s)');
ylabel('Amplitude (Re)');
title('Transmitted signal');
subplot(3,1,2);
t_start = round(tau*fs)+1;
plot(t(t_start:4000), real(signal_rx_0(t_start:4000)),'-r');
xlabel('Time (s)');
ylabel('Amplitude (Re)');
title('Received signal');
subplot(3,1,3);
plot(t(t_start:4000), real(beat_signal_0(t_start:4000)),'-b');
title('Beat signal');
xlabel('Time (s)');
ylabel('Amplitude (Re)');
grid on;

% ---- Fast Fourier Transform

% For target_0
beat_signal_cleaned_0 = beat_signal_0(t_start:end);        % Cutting of first samples (before received signal arivall)
N_fft_0 = 2^nextpow2(length(beat_signal_cleaned_0));       % Optymalizing fft transform function
fft_0 = fft(beat_signal_cleaned_0, N_fft_0);               % Fast Fourier transform
fft_shifted_0 = fftshift(fft_0);                           % Shifting order

% For target_1
beat_signal_cleaned_1 = beat_signal_1(t_start:end);        % Cutting of first samples (before received signal arivall)
N_fft_1 = 2^nextpow2(length(beat_signal_cleaned_1));       % Optymalizing fft transform function
fft_1 = fft(beat_signal_cleaned_1, N_fft_1);               % Fast Fourier transform
fft_shifted_1 = fftshift(fft_1);                           % Shifting order

f_axis = fs * (-N_fft_0/2 : N_fft_0/2 - 1) / N_fft_0;      % Axis for fft ploting

% FFT ploting
figure;
plot(f_axis/1e6, abs(fft_shifted_0), 'b-', 'LineWidth', 1.5); 
hold on;
plot(f_axis/1e6, abs(fft_shifted_1), 'r--', 'LineWidth', 1.5); 
hold off;
grid on;
title('FFT Spectrum Comparison: Doppler Shift Effect');
xlabel('Frequency (MHz)');
ylabel('Amplitude |Y(f)|');
legend('Stationary Target (v = 0 km/h)', 'Moving Target (v = 100 km/h)');
% Zoom settings
[~, idx_peak] = max(abs(fft_shifted_0));
f_peak_MHz = f_axis(idx_peak) / 1e6;
zoom = 0.5; 
xlim([f_peak_MHz - zoom, f_peak_MHz + zoom]);

% ---- Spectrogram of chirp ----

%Spectrogram parametrs
window_size = 1024;
overlap = 1000;
nfft = 2048;
% Spectrogram plotting 
figure;
[s, f, t_spec] = spectrogram(signal_tx, window_size, overlap, nfft, fs);    % Generating data for spectrogram drawin
imagesc(t_spec * 1e6, f / 1e6, abs(s));                                     % Plotting spectrogram
ylim([0, bandwidth_ / 1e6]);    
set(gca, 'YDir', 'normal'); 
colormap('jet');
colorbar;
title('Spetrogram of transmitted signal (TX Chirp)');
xlabel('Time (us)');
ylabel('Frequency (MHz)');