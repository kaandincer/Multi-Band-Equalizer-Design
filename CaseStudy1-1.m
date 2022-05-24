close all; clear; clc;
%% load noisy violin/Ode data, then play recording
%[xv,xvfs] = audioread('violindirty.wav');
[xv,xvfs] = audioread('OdeSnippetNoisy.wav');
fs = xvfs; 
% sound(xv,fs)
xv_ct = xv;
figure, spectrogram(xv,1024,200,1024,fs)
title('Original Spectogram');

XV = fft(xv); 
XV = XV(1:length(XV)/2);
f = [0:length(XV)-1].*fs./length(XV);
figure, plot(f,abs(XV));
title('Original Violin FFT');
xlabel('Frequency: Hz');
ylabel('FFT Magnitude');
xlim([1,10000]);

% Bandpass Filters
cutoff1 = 800;
cutoff2 = 1800;
cutoff3 = 3000;
cutoff4 = 6000;

factor = 250;

% Change gain accoring to the signal being filtered:
% - for OdeSnippet: gain = [.001 0.2 1 1 1]
gain = [.001 0.2 1 1 1];

% - for Violin: gain = [.001 .2 1 .2 .001];
%gain = [.001 .2 1 .2 .001];

% Filter 1 
[xv1, d1] = lowpass(xv,cutoff1-factor,fs,'ImpulseResponse','iir','Steepness',0.95);

% Filter 2
[xv2, d2] = bandpass(xv,[cutoff1+factor cutoff2-factor],fs);

% Filter 3
[xv3, d3] = bandpass(xv,[cutoff2+factor cutoff3-factor],fs);

% Filter 4
[xv4, d4] = bandpass(xv,[cutoff3+factor cutoff4-factor],fs);

% Filter 5
[xv5, d5] = highpass(xv,cutoff4+factor,fs);

% Combined Filtered Signal
xv = gain(1)*xv1+gain(2)*xv2+gain(3)*xv3+gain(4)*xv4+gain(5)*xv5;

% Post-Filter Spectrogram
figure, spectrogram(xv,256,200,256,fs);
title('Post-Filter Spectogram');

sound(xv,fs)
%% FFT and Frequency Response
XV = fft(xv); 
XV = XV(1:length(XV)/2);
f = [0:length(XV)-1].*fs./length(XV);

figure();
plot(f,abs(XV));
title('Filtered Violin FFT');
xlabel('Frequency: Hz');
ylabel('FFT Magnitude');
xlim([1,10000]);


[h1, freq1] = freqz(d1,1024,fs);
[h2, freq2] = freqz(d2,1024,fs);
[h3, freq3] = freqz(d3,1024,fs);
[h4, freq4] = freqz(d4,1024,fs);
[h5, freq5] = freqz(d5,1024,fs);

% Frequency response of each band
figure();
subplot(2,1,1);
plot(freq1, mag2db(abs(h1)));
title('Band 1 Magnitude')
xlabel('Freq.');
ylabel('dB');
subplot(2,1,2);
plot(freq1, angle(h1)/pi);
title('Band 1 Phase');
xlabel('Freq.');
ylabel('Degrees');

figure();
subplot(2,1,1);
plot(freq2, mag2db(abs(h2)));
title('Band 2 Magnitude')
xlabel('Freq.');
ylabel('dB');
subplot(2,1,2);
plot(freq2, angle(h2)/pi);
title('Band 2 Phase');
xlabel('Freq.');
ylabel('Degrees');

figure();
subplot(2,1,1);
plot(freq3, mag2db(abs(h3)));
title('Band 3 Magnitude')
xlabel('Freq.');
ylabel('dB');
subplot(2,1,2);
plot(freq2, angle(h3)/pi);
title('Band 3 Phase');
xlabel('Freq.');
ylabel('Degrees');

figure();
subplot(2,1,1);
plot(freq4, mag2db(abs(h4)));
title('Band 4 Magnitude')
xlabel('Freq.');
ylabel('dB');
subplot(2,1,2);
plot(freq2, angle(h4)/pi);
title('Band 4 Phase');
xlabel('Freq.');
ylabel('Degrees');

figure();
subplot(2,1,1);
plot(freq5, mag2db(abs(h5)));
title('Band 5 Magnitude')
xlabel('Freq.');
ylabel('dB');
subplot(2,1,2);
plot(freq2, angle(h5)/pi);
title('Band 5 Phase');
xlabel('Freq.');
ylabel('Degrees');

% Combined plot
figure();
hold on;
plot(freq1, mag2db(abs(h1)));
plot(freq2, mag2db(abs(h2)));
plot(freq3, mag2db(abs(h3)));
plot(freq4, mag2db(abs(h4)));
plot(freq5, mag2db(abs(h5)));
title(['Crosstalk Comparison: Factor = ',num2str(factor)]);
legend('Lowpass','Band1','Band2','Band3','Highpass');
xlabel('Hz');
ylabel('dB');
%% Impulse Response
impulse = zeros(1,fs);
impulse(5) = fs;

% Filters
xv1 = lowpass(impulse,cutoff1-factor,fs,'ImpulseResponse','iir','Steepness',0.95);
XV1 = fft(xv1); 
XV1 = XV1(1:length(XV1)/2);
f = [0:length(XV1)-1].*fs./length(XV1);

xv2 = bandpass(impulse,[cutoff1+factor cutoff2-factor],fs);
XV2 = fft(xv2); 
XV2 = XV2(1:length(XV2)/2);


xv3 = bandpass(impulse,[cutoff2+factor cutoff3-factor],fs);
XV3 = fft(xv3); 
XV3 = XV3(1:length(XV3)/2);

xv4 = bandpass(impulse,[cutoff3+factor cutoff4-factor],fs);
XV4 = fft(xv4); 
XV4 = XV4(1:length(XV4)/2);

xv5 = highpass(impulse,cutoff4+factor,fs);
XV5 = fft(xv5); 
XV5 = XV5(1:length(XV5)/2);

xv_imp = xv1+xv2+xv3+xv4+xv5;
XV_imp = fft(xv_imp); 
XV_imp = XV_imp(1:length(XV_imp)/2);

figure();
hold on;
plot(f,abs(XV1));
plot(f,abs(XV2));
plot(f,abs(XV3));
plot(f,abs(XV4));
plot(f,abs(XV5));
title('Impulse Responses');
xlabel('Frequency: hZ');
ylabel('Amplitude');
legend('Lowpass','Band1','Band2','Band3','Highpass');

figure();
plot(f,abs(XV_imp));
title('Total Impulse Response');
xlabel('Frequency: hZ');
ylabel('Amplitude');

%% Signal of choice: rising synth
[xv_synth,xvfs] = audioread('RisingSynth.wav');
fs = xvfs; 

factor = 300;
gain = [1 1 1 .01 1];
xv1 = lowpass(xv_synth,cutoff1-factor,fs,'ImpulseResponse','iir','Steepness',0.95);
xv2 = bandpass(xv_synth,[cutoff1+factor cutoff2-factor],fs);
xv3 = bandpass(xv_synth,[cutoff2+factor cutoff3-factor],fs);
xv4 = bandpass(xv_synth,[cutoff3+factor cutoff4-factor],fs);
xv5 = highpass(xv_synth,cutoff4+factor,fs);
xv_synth_post1 = gain(1)*xv1+gain(2)*xv2+gain(3)*xv3+gain(4)*xv4+gain(5)*xv5;

% sound(xv_synth,fs)

figure();
hold on;
plot(xv_synth(:,1));
plot(xv_synth_post1(:,1));
title('Rising Synth');
legend('Input','Output');
xlabel('Time');
ylabel('Waveform');