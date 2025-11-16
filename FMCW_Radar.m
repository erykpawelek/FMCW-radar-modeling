clc;
clear;
close all;

% Light speed const in m/s
c = 299792458;
%Carrier frequency (77 GHz - standard for automotive radars)
fc = 77e9; 

% ---- Project requirements -----

max_range = 200; % Maximal range in meters
range_resolution = 0.1; % Minimal distinguishable difference in distance in meters
max_velocity = 150; % Maximal velocity of measured target in km/h

bandwidth_ = c / (2*range_resolution); % Bandwidth of our signal 

T_of_flight = 2*max_range/c; % Time required to signal travel both ways 

% Value below is increased by some factor in order to obtain sufficient
% amount of time to compare it to beat signal
T_sweep = 100 * T_of_flight; % Duration of 1 chirp 

slope = bandwidth_/T_sweep;

% Ploting parameters
disp(['Bandwidth : ', num2str(bandwidth_/1e9), ' GHz']);
disp(['Chirp time (T_sweep): ', num2str(T_sweep*1e6), ' us']);
disp(['Slope: ', num2str(slope/1e12), ' THz/s']);

% ---- Target definition ---- 
target_range = 100; % In meters
target_vel = 80; % In km/h

if target_vel > max_velocity
    error('Target velocity exceeds maximum allowable velocity.');

else
    % ---- Signal definition ----

    fs = 10*bandwidth_; % Sampling frequency
    t = 0: (1/fs) : T_sweep; % Time vector for one chirp
    phase_tx = (2*pi*slope*t.^2)/2; % Phase of one chirp
    signal_tx = exp(1j * phase_tx);
    disp(['Sampling frequency: ', num2str(fs/1e9), ' GHz']);
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

    tau = (2 * target_range) / c; % This is our signal travel time (round-trip)
    lambda = c / fc; % wavelenght
    f_d = (2 * target_vel * 1000 / 3600) / lambda; % Doppler frequency 
    phase_rx = (2 * pi * slope * (t-tau).^2) / 2; % Reciving signal phase
    doppler_shift_signal = exp(1j * 2 * pi * f_d * t); % Generationg doppler shift signal
    noise = 0.5 * (randn(size(t)) + 1j * randn(size(t))); % Generating complex gausian noise 
    signal_rx = (exp(1j * phase_rx) .* doppler_shift_signal) + noise; %Return signal with noise and phase shift
    % We multiply signals to make use of exponent adding properiety of a
    % signals
    disp(['Travel time: ', num2str(tau), ' s']);
    disp(['Doppler freqency: ', num2str(f_d/1e6), ' MHz']);
    
   % ---- Beat signal ----
   
    beat_signal = signal_rx .* conj(signal_tx); % Generating beat signal
    % Beat signal is generated using Complex Conjugate function, by that we
    % ramowe our tx component of recived signal (we subtract it)

    % ---- printing results ----
    figure;
    plot(t, real(beat_signal));
    title('Beat Signal');
    xlabel('Time (s)');
    ylabel('Amplitude (Re)');
    grid on;
    figure;
    subplot(3,1,1);
    plot(t(1:50000), real(signal_tx(1:50000)),'-g');
    xlabel('Time (s)');
    ylabel('Amplitude (Re)');
    title('Transmitted signal');
    subplot(3,1,2);
    t_start = round(tau*fs)+1;
    plot(t(t_start:50000), real(signal_rx(t_start:50000)),'-r');
    xlabel('Time (s)');
    ylabel('Amplitude (Re)');
    title('Received signal');
    subplot(3,1,3);
    plot(t(t_start:50000), real(beat_signal(t_start:50000)),'-b');
    title('Beat signal');
    xlabel('Time (s)');
    ylabel('Amplitude (Re)');
    grid on;

    % ---- Fast Fourier Transform

    beat_signal_cleaned = beat_signal(t_start:end); % Cutting of first samples (before received signal arivall)
    N_fft = 2^nextpow2(length(beat_signal_cleaned)); % Optymalizing fft transform function
    fft_ = fft(beat_signal_cleaned, N_fft); % Fast Fourier transform
    fft_shifted = fftshift(fft_); % Shifting order
    f_axis = fs * (-N_fft/2 : N_fft/2 - 1) / N_fft; % Creating x axis to print fft
    

    figure;
    subplot(2,1,1)
    plot(f_axis / 1e6, abs(fft_shifted)); 
    title('FFT Spectrum of the Beat Signal');
    xlabel('Frequency (MHz)');
    ylabel('|fft(beat(signal)|');
    grid on;
    subplot(2,1,2)
    plot(f_axis / 1e6, abs(fft_shifted)); 
    title('FFT Spectrum of the Beat Signal');
    xlabel('Frequency (MHz)');
    ylabel('|fft(beat(signal)|');

    [max_amplitude, idx_peak] = max(abs(fft_shifted));
    f_peak_hz = f_axis(idx_peak); 
    
    % Peak zooming
    zoom_range_hz = max(10e6, abs(f_peak_hz * 0.5)); 
    xlim_min = (f_peak_hz - zoom_range_hz) / 1e6;
    xlim_max = (f_peak_hz + zoom_range_hz) / 1e6;
    xlim([max(xlim_min, f_axis(1)/1e6), min(xlim_max, f_axis(end)/1e6)]);
end
