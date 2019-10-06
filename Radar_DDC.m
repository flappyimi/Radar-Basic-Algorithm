%% DDC Chirp
clear ; close all; clc;

% parameter
N  = 4096;
f0 = 60e+6;      % 60MHz中频 
B  = 8e+6;       % 8MHz带宽
Tp = 10e-6;      %10us时宽
fs = 80e+6;      % 80MHz采样频率
fnco = 20e+6;    % 20MHz NCO频率

t  = (0:N-1)/fs;
fd = (-N/2:N/2-1)*fs/N;
kr = B/Tp;
% Generate LFM @f0
s_in  = rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2).*exp(-1i*2*pi*f0*t);
s_in_fft = fft(s_in,N);
figure(1);
subplot(221);plot(t/1e-6,real(s_in));
xlabel('t/us');ylabel('Amplitude');title('Input sampled signal（time Re）');
subplot(222);plot(fd/1e6,abs(fftshift(s_in_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('Input sampled signal(spectrum)');

s_fix = s_in.*exp(-1i*2*pi*fnco*t);
s_fix_fft = fft(s_fix,N);
subplot(223);
plot(t/1e-6,real(s_fix));hold on;
plot(t/1e-6,imag(s_fix));xlim([0, 10]);
xlabel('t/us');ylabel('Amplitude');
title('NCO signal (time)');legend('Re','Im');
subplot(224);plot(fd/1e6,abs(fftshift(s_fix_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('NCO signal (spectrum)');

coeff = fir1(127,B/(fs/2),hamming(128)); % 0.2 = B/(fs/2) 此处必须B大于信号带宽，否则信号失真
figure(2);freqz(coeff);
title('FIR filter speciation');

% fir filter
s_ddc = conv(s_fix,coeff);
scale = 4;
s_ddc = downsample(s_ddc,scale);
s_ddc_fft = fft(s_ddc,N);

figure(3);
subplot(211);
plot(real(s_ddc));hold on;
plot(imag(s_ddc));xlim([0, 1000/scale]);
xlabel('t/us');ylabel('Amplitude');
title('Baseband signal (time)');legend('Re','Im');
subplot(212);plot(fd/1e6/scale,abs(fftshift(s_ddc_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('Baseband signal (spectrum)');

%% DDC sine
clear ; close all; clc;

% parameter
N  = 4000;
f0 = 61e+6;      % 60MHz中频 
fs = 80e+6;      % 80MHz采样频率
fnco = 20e+6;    % 20MHz NCO采样频率

t  = (0:N-1)/fs;
fd = (-N/2:N/2-1)*fs/N;

s_in  = sin(2*pi*t*f0);
s_in_fft = fft(s_in,N);
figure(1);
subplot(221);plot(t/1e-6,s_in);xlim([0, 2]);
xlabel('t/us');ylabel('Amplitude');title('Input sampled signal（time）');
subplot(222);plot(fd/1e6,abs(fftshift(s_in_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('Input sampled signal(spectrum)');

s_fix = s_in.*exp(-1i*2*pi*fnco*t);
s_fix_fft = fft(s_fix,N);

subplot(223);
plot(t/1e-6,real(s_fix));hold on;
plot(t/1e-6,imag(s_fix));xlim([0, 2]);
xlabel('t/us');ylabel('Amplitude');
title('NCO output signal(time)');legend('Re','Im');
subplot(224);plot(fd/1e6,abs(fftshift(s_fix_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('NCO output signal(spectrum)');

coeff = fir1(127,2e6/(fs/2),hamming(128));
figure(2);freqz(coeff);
title('FIR filter speciation');

% fir filter
s_ddc = conv(s_fix,coeff);
scale = 10;
s_ddc = downsample(s_ddc,scale);
s_ddc_fft = fft(s_ddc,N);


figure(3);
subplot(211);
plot(real(s_ddc));hold on;
plot(imag(s_ddc));xlim([0, 1000/scale]);
xlabel('t/us');ylabel('Amplitude');
title('Baseband signal (time)');legend('Re','Im');
subplot(212);plot(fd/1e6/scale,abs(fftshift(s_ddc_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('Baseband signal (spectrum)');