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
    
    % ---- Signal propagation and echo ----

    tau = (2 * target_range) / c; % This is our signal travel time (round-trip)
    lambda = c / fc; % wavelenght
    f_d = (2 * target_vel * 1000 / 3600) / lambda; % Doppler frequency 
    phase_rx = (2 * pi * slope * (t-tau).^2) / 2; % Reciving signal phase
    doppler_shift_signal = exp(1j * 2 * pi * f_d * t); % Generationg doppler shift signal
    noise = 0.5 * (randn(size(t)) + 1j * randn(size(t))); % Generating complex gausian noise 
    signal_rx = (exp(1j * phase_rx) .* doppler_shift_signal) + noise; %Return signal with noise and phase shift
    disp(['Travel time: ', num2str(tau), ' s']);
    disp(['Doppler freqency: ', num2str(f_d/1e6), ' MHz']);
   % ---- Beat signal ----
   
    beat_signal = signal_rx .* conj(signal_tx);
    figure;
    plot(t, real(beat_signal));
    title('Beat Signal');
    xlabel('Time (s)');
    ylabel('Amplitude (Re)');
    grid on;
end
