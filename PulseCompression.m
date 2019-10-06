%% Pulse Compression
clear,clc,close all;
%Parameters
N  = 4000;
Tp = 10e-6;
B  = 10e6;
f0 = 10e9;
fs = 100e6;
R1 = 3000;
R2 = 4000;
c  = 3e8;

t  = (0:N-1)/fs;
rd = c/2 .* t;
fd = (-N/2:N/2-1)*fs/N;
kr = B/Tp;
tao1 = 2*R1/c;
tao2 = 2*R2/c;
%Receive Signal & Reference with Hamming
sr  = rectpuls(t-Tp/2-tao1,Tp).*exp(1i*pi*kr*(t-Tp/2-tao1).^2).*exp(-1i*2*pi*f0*tao1);
sf  = rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2);
sf_h =rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2).*[hamming(1001)',zeros(1,(N-1001))];

%FFT
sr_fft = fft(sr,N);
sf_fft = fft(sf,N);

figure(1);
subplot(321);plot(rd/1e3,real(sr));
xlabel('Range/km');ylabel('Amplitude');title('Signal£¨Re£©');
subplot(322);plot(fd/1e6,abs(fftshift(sr_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('Spectrum of signal');
subplot(323);plot(rd/1e3,real(sf));
xlabel('Range/km');ylabel('Amplitude');title('Reference signal£¨Re£©');
subplot(324);plot(fd/1e6,abs(fftshift(sf_fft)));
xlabel('Frequency/MHz');ylabel('Amplitude');title('Spectrum of reference signal');
%ylim([80,120]);
%pulse Compression
sf_hfft = fft(sf_h,N);
s_out = ifft(sr_fft.*conj(sf_hfft));

subplot(3,2,[5 6]);
plot(rd/1e3,20*log10(abs(s_out)/max(abs(s_out))));
xlabel('Range/km');ylabel('Amplitude/dB');title('Result of pulse compression');
axis ([0 6 -60 1]);

%Use pulse compression change the range resolution
sr2  = rectpuls(t-Tp/2-tao2,Tp).*exp(1i*pi*kr*(t-Tp/2-tao2).^2).*exp(-1i*2*pi*f0*tao2);
sr_in2 = sr+sr2;

figure(2);
subplot(211);plot(rd/1e3,real(sr+sr2));
xlabel('Range/km');ylabel('Amplitude/dB');title('Signal before pluse comprssion');
s_out2 = ifft(fft(sr_in2).*conj(sf_hfft));
subplot(212);plot(rd/1e3,20*log10(abs(s_out2)/max(abs(s_out2))));
xlabel('Range/km');ylabel('Amplitude/dB');title('Signal after pluse comprssion');
axis ([0 6 -60 1]);