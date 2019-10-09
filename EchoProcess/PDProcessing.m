%% Pulse Doppler Processing
clear,clc,close all;
N = 4096;
N_PRT = 64;
Tp = 10e-6;
T_PRT = 100e-6;
f0 = 10e9;
fs = 100e6;
B  = 10e6;
kr  = B/Tp;
v0  = 60;
R0 = 3000;
c  = 3e8;
t = (0:N-1)/fs;
r = c*t/2;
m = 0:N_PRT-1;
v = (m-N_PRT/2)/N_PRT /2/T_PRT *c/f0;
s = zeros(N_PRT,N);
s_pd = zeros(N_PRT,N);
noise = zeros(N_PRT,N);
np = zeros(N_PRT,N);
n_pd = zeros(N_PRT,N);
SNR = 10;                        %信噪比 dB
sigma = 1;
A = 10^(SNR/20)*sigma;  %信号幅度

%Pulse compression along range direction
s_ref = rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2).*[hamming(1001)',zeros(1,N-1001)];

for j = 1:N_PRT
    tao = 2*(R0-(j-1)*T_PRT*v0)/c;
    noise(j,:) = sigma/sqrt(2).*(randn(1,N)+1i*randn(1,N));
    s_in = A.*rectpuls(t-Tp/2-tao,Tp).*exp(1i*pi*kr*(t-Tp/2-tao).^2).*exp(-1i*2*pi*f0*tao) + noise(j,:);
    s(j,:) = ifft(fft(s_in,N).*conj(fft(s_ref,N)));
    np(j,:)= ifft(fft(noise(j,:),N).*conj(fft(s_ref,N)));
end
subplot(211);imagesc(r,m,abs(s));%?dB
xlabel('Range(m)');ylabel('PRT No');title('Pulse Compression');
for n = 1:N
    s_pd(:,n) = fftshift(fft(s(:,n),N_PRT));
    n_pd(:,n) = fftshift(fft(np(:,n),N_PRT));
end
subplot(212);imagesc(r,v,abs(s_pd));
xlabel('Range(m)');ylabel('Velocity(m/s)');title('Pulse Doppler');

[v_point,r_point] = find(abs(s_pd) == max(max(abs(s_pd))));
v_est = (v_point-N_PRT/2-1)/N_PRT./T_PRT.*c/f0/2
r_est = (r_point-1)*c/2/fs
v_uamb = 1/T_PRT/2.*c/f0/2 

SNRi  = 20*log10(abs(s(v_point,r_point)-np(v_point,r_point))/(sqrt(N_PRT)*sigma));
SNRo1 = 20*log10(abs(s_pd(v_point,r_point)-n_pd(v_point,r_point))/(sqrt(N_PRT)*sigma));
gain1 = 10^((SNRo1 - SNRi)/20)
SNRo2 = 10*log10(sum(abs(s(:,r_point)-np(:,r_point)).^2)/(N_PRT*sigma^2));
gain2 = 10^((SNRo2 - SNRi)/20)

%% PD Processing Gain



