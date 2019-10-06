%% Target Detection Pd Pfa,CA-CFAR
close all;clear all;clc;
N = 1e4;        %Points
n = 1:N;
SNR = 10;       %SNR in dB
Pfa = 1e-3;     
sigma_dB = 0;   %Noise power in dB
sigma = db2pow(sigma_dB);           %Noise Power=1w

noise = sigma/sqrt(2)*(randn(1,N)+1i*randn(1,N));
A = sqrt(sigma^2 * 10^(SNR/10));    %Target signal
T_fix =  sqrt(-sigma^2*log(Pfa));   %fix thresholding
c = 500;                            %one side,ref cell = c-guardc
guardc = 10;   
CUTIdx = N/2;
%Simulation of target detection%
sample = rand(1,N) < 0.001;
A0 = A .* sample;
%A1 = abs(A0 + noise.*sample);
signal1 = A0 + noise;              %Receive signal

figure(1);
idx_target = find(A0 > 0);
idx_detect = find(abs(signal1) > T_fix);
plot(n,abs(signal1),'-o','MarkerIndices',idx_detect,'MarkerSize',10);
hold on; plot(n(idx_target),abs(signal1(1,idx_target)),'p','MarkerEdgeColor','r','MarkerFaceColor','r');
hold on; plot(n,T_fix*ones(1,N));
xlim([c+2,N-c-1]);

%Evaluation of Pd & Pfa %
signal2 = A + noise;
disp('Detection Rate£º')
Pd_est  = sum(abs(signal2) > T_fix) / N
disp('False alarm Rate£º')
Pfa_est = sum(abs(noise(1,c+2:N-c-1)) > T_fix) / N

%Simulation of CFAR detection
%Ncells = 23;%Ntrails =1e4;


signal = noise;
signal(1,CUTIdx) = signal(1,CUTIdx) + A;
k = Pfa^(-1/(2*(c-guardc)))-1;

T_cfar = T_fix .* ones(1,N);
for j=c+1:N-c
    pn = sum(abs(signal(1,(j-c:j-guardc-1))).^2) + sum(abs(signal(1,(j+guardc+1:j+c))).^2);
    T_cfar(j) =k.*abs(pn);
end
figure(2);
idxC = find(abs(signal).^2 > T_cfar);
plot(n,abs(signal).^2,'-o','MarkerIndices',idxC,'MarkerSize',10);
hold on;stem(N/2,abs(signal(1,N/2)).^2,'-p','MarkerEdgeColor','r','MarkerFaceColor','r');
hold on;plot(n,T_cfar);
xlim([c+2,N-c-1]);
Pfa_cfar = sum(abs(noise(1,c+2:N-c-1)).^2 > T_cfar(1,c+2:N-c-1))./(N-2*c)
% %Pfa_np = sum(x(CUTIdx,:)>T_ideal)/Ntrials

%% Pd
N_trial = 10000;
sum_fix = 0;
sum_cfar = 0;
T_fix =  sqrt(-sigma^2*log(Pfa)); 
for n = 1:N_trial
    noise = sigma/sqrt(2)*(randn(1,N)+1i*randn(1,N));
    signal = noise;
    signal(1,CUTIdx) = signal(1,CUTIdx) + A;
    sum_fix = sum_fix + (abs(signal(1,CUTIdx)) > T_fix);
    
    pn = sum(abs(signal(1,(CUTIdx-c:CUTIdx-guardc-1))).^2) + sum(abs(signal(1,(CUTIdx+guardc+1:CUTIdx+c))).^2);
    T_cfar0 =k.*abs(pn);
    sum_cfar = sum_cfar + (abs(signal(1,CUTIdx)).^2 > T_cfar0);
end
pd_fix = sum_fix / N_trial
pd_cfar = sum_cfar / N_trial

%% CFAR COMPARE
close all;clear all;clc;
N = 6e3;        
n = 1:N;
SNR = 20;       
Pfa = 1e-6;     
sigma_dB = 0;   
sigma = db2pow(sigma_dB);           %Noise Power=1w
ndeta = 100;
scale = 5;
noise = sigma*[randn(1,N/2-ndeta), scale*randn(1,N/2+ndeta)];
A = sqrt((scale*sigma)^2 * 10^(SNR/10));    %Target signal amplitude
T_fix =  sqrt(-sigma^2*log(Pfa));   %fix thresholding
c = N/40;                             %one side,ref cell = c-guardc
guardc = 10;                         %guard cell number
CUTIdx = N/2;                       %position of target

%CFAR mode choose%%%%%,0 CA;1 GO;2 SO;3:OS
cfar = 3 ;

%%Simulation of CFAR detection
%signal model
signal = noise;
signal(1,CUTIdx) = signal(1,CUTIdx) + A;
%k = Pfa^(-1/(2*(c-guardc)))-1;
T_cfar = ones(1,N);
ncell = 0;
pn = 0;

for j=c+1:N-c
    a = sum(abs(signal(1,(j-c:j-guardc-1))).^2);
    b = sum(abs(signal(1,(j+guardc+1:j+c))).^2);
    if cfar==0
        pn = a + b;
        ncell = 2*(c-guardc);
    end
    if cfar==1
        pn = max(a,b);
        ncell = c-guardc;
    end
    if cfar==2
        pn = min(a,b);
        ncell = c-guardc;   
    end 
    k = Pfa^(-1/ncell)-1;
    T_cfar(j) =k.*abs(pn);
end
n_target = 1; %target numbers
if cfar==3
    s = sort(signal);
    T_cfar = ones(1,N) .* s(1,N-n_target).^2;
end

%traget detected
idxC = find(abs(signal).^2 > T_cfar);

figure;
plot(n,10*log(abs(signal).^2),'-o','MarkerIndices',idxC,'MarkerSize',10);
hold on;stem(N/2,10*log(abs(signal(1,N/2)).^2),'-p','MarkerEdgeColor','r','MarkerFaceColor','r');
hold on;plot(n,10*log(T_cfar));
xlim([c+2,N-c-1]);ylim([0,100]);
Pfa_cfar = sum(abs(noise(1,c+2:N-c-1)).^2 > T_cfar(1,c+2:N-c-1))./(N-2*c)