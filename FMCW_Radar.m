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
disp(['Measurement resolution:  ', num2str(range_resolution), ' m']);
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
[s, f, t_spec] = spectrogram(signal_tx, window_size, overlap, nfft, fs);    % Generating data for spectrogram plot
s_db = 10 * log10((abs(s)).^2 +eps);                                        % Converting to db Power scale , we add eps, which is relly inrelevant
                                                                            % positive number that allows to deal with log10(0) = -inf rule and by that we can scale our spectrogram properetly 
imagesc(t_spec * 1e6, f / 1e6, s_db);                                       % Plotting spectrogram
ylim([0, bandwidth_ / 1e6]);
max_val = max(s_db(:));                                                     % Regulate renges of colorbar to get clear colour image
clim([max_val - 60, max_val]); 
set(gca, 'YDir', 'normal');                                                 % We rotate y axis to obtain growing in up direction
colormap('jet');
c = colorbar;
c.Label.String = 'Power (dB)';
title('Spetrogram of transmitted signal (TX Chirp)');
xlabel('Time (us)');
ylabel('Frequency (MHz)');

%% 
%   #### PART 2 REALIZATION USING PHASED ARRAY SYSTEM TOOLBOX ####

clc;
clear all;

c = 299792458;  % Light speed const in m/s
fc = 77e9;      % Carrier frequency (77 GHz - standard for automotive radars)
lambda = c/fc;  % Wavelength 

% ---- Project requirements -----

max_range = 200;                                % Maximal range in meters
t_sweep = 5.5 * range2time(max_range, c);       % Duration of 1 chirp 
range_resolution = 0.5;                         % Minimal distinguishable difference in distance in meters
bandwidth_ = rangeres2bw(range_resolution, c);  % Bandwidth of our signal 
slope = bandwidth_/t_sweep;                     % Frequency slope of chirp
fr_max = range2beat(max_range,slope,c);   % Obtaining max freq component of beat signal colerated to distance
max_velocity = 150*1000/3600;                   % Max velocity in meters per second
fd_max = speed2dop(2*max_velocity,lambda);      % Obtaining max freq component of beat signal colerated to speed
fb_max = fr_max+fd_max;                         % Maximum beat signal freq containing distance component as well as velocity component


% ---- HINT ---- 
%{ 
Because an FMCW signal often occupies a large bandwidth, setting the sample rate blindly to twice the bandwidth often 
stresses the capability of A/D converter hardware. To address this issue, you can often choose a lower sample rate.
Consider two things here:

* For a complex sampled signal, the sample rate can be set to the same as the bandwidth.
* FMCW radars estimate the target range using the beat frequency embedded in the dechirped signal. 
   The maximum beat frequency the radar needs to detect is the sum of the beat frequency corresponding to the
   maximum range and the maximum Doppler frequency. Hence, the sample rate only needs to be twice the maximum beat frequency.
%}

fs = max(2*fb_max, bandwidth_); % Securing sample time to meet sapling requirement

%{
 For our simulation purpuses we use bandwidth sampling frequency becausewe want to 
 generate our tx_signal. In real radar decoding rx is proceeded in analog
 mixer so we anly need sufficient sampling freq to decode beat signal.
%}

waveform = phased.FMCWWaveform('SweepTime',t_sweep, ...
    'SweepBandwidth',bandwidth_, ...
    'SampleRate',fs, ...
    'SweepInterval', 'Positive');

signal_tx = waveform();
subplot(211); plot(0:1/fs:t_sweep-1/fs,real(signal_tx));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal (aliasing)'); axis tight;
subplot(212); spectrogram(signal_tx,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');

% ---- Target definition ---- 

target_range = 100;
target_speed = -96*1000/3600;
car_rcs = db2pow(min(10*log10(target_range)+5,20));      % Simulation of changing size of our target based on a distance  

target = phased.RadarTarget( ...                           % Definition of our target system object from electromagnetic side 
    'MeanRCS',car_rcs, ...
    'PropagationSpeed',c,...
    'OperatingFrequency',fc);

target_motion = phased.Platform( ...                       % Definition another system object which describes movement of a target
    'InitialPosition',[target_range;0;0.5],...
    'Velocity',[target_speed;0;0]);

channel = phased.FreeSpace( ...                            % Definition of another system object which describes wave propagation model
    'PropagationSpeed',c,...
    'OperatingFrequency',fc, ...
    'SampleRate',fs, ...
    'TwoWayPropagation',true);

% Hardware config
ant_aperture = 6.06e-4;                         % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(5)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB

% Definition of transmiter and receiver system objects
transmitter = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,'SampleRate',fs);

% Definition of transmitter movement system object
radar_speed = 10*1000/3600;
radar_motion = phased.Platform('InitialPosition',[0;0;0.5],...
    'Velocity',[radar_speed;0;0]);
% Configuration of spectrum analizer, we chec if our signal didn't got lost
% in noise
specanalyzer = spectrumAnalyzer('SampleRate',fs, ...
    'Method','welch','AveragingMethod','running', ...
    'PlotAsTwoSidedSpectrum',true, 'FrequencyResolutionMethod','rbw', ...
    'Title','Spectrum for received and dechirped signal', ...
    'ShowLegend',true, ...
    'ChannelNames', {'Transmitted Signal', 'Dechirped (Beat) Signal'});


rng(2012);                  % Defining randomnes to be the same in every program run
Nsweep = 64;                % Numer of chirps for measurment
xr = complex(zeros(waveform.SampleRate*waveform.SweepTime,Nsweep)); % Data storage 2D matrix, 

for m = 1:Nsweep
    % Update radar and target positions
    [radar_pos,radar_vel] = radar_motion(waveform.SweepTime);
    [tgt_pos,tgt_vel] = target_motion(waveform.SweepTime);

    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);

    % Propagate the signal and reflect off the target
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    txsig = target(txsig);

    % Dechirp the received radar return
    txsig = receiver(txsig);
    dechirpsig = dechirp(txsig,sig);

    % Visualize the spectrum
    specanalyzer([txsig dechirpsig]);
    % Data collection
    xr(:,m) = dechirpsig;
end

rngdopresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
    'RangeMethod','FFT','SweepSlope',slope,...
    'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',256);

clf;
plotResponse(rngdopresp,xr);                     % Plot range Doppler map
axis([-max_velocity max_velocity 0 max_range]);
title('Range-Doppler Map');