%% OSINR
clear,clc,close all;
L=8;
K=50;
d=1;
M=3;
sigma=1;
snr=-10:30;
N = 500;
%初始化输出信干噪比矩阵
osinr_mpdr=zeros(length(snr),N);
osinr_smi=zeros(length(snr),N);
osinr_rmi=zeros(length(snr),N);
osinr_dl=zeros(length(snr),N);
osinr_sip=zeros(length(snr),N);

%初始化失配广义夹角矩阵
angle_mismatched_mpdr=zeros(1,length(snr));
angle_mismatched_smi=zeros(1,length(snr));
angle_mismatched_rmi=zeros(1,length(snr));
angle_mismatched_dl=zeros(1,length(snr));
angle_mismatched_sisp=zeros(1,length(snr));


theta0=0;
theta1=60;
theta2=-60;
a0=exp(1i*2*pi*d*sin(theta0*pi/180)*(0:L-1)');
a1=exp(1i*2*pi*d*sin(theta1*pi/180)*(0:L-1)');
a2=exp(1i*2*pi*d*sin(theta2*pi/180)*(0:L-1)');
b0=exp(1i*2*pi*d*sin((theta0+2)*pi/180)*(0:L-1)');
for m=1:length(snr)
    sigma0=sigma*10^(snr(m)/20);
    %两个干扰源，信干噪比为-30dB
    inr1=snr(m)+30;
    inr2=snr(m)+30;
    sigma1=sigma*10^(inr1/20);
    sigma2=sigma*10^(inr2/20);
        
    Rxx = sigma0^2*(a0*a0')+sigma1^2*(a1*a1')+sigma2^2*(a2*a2')+sigma^2*eye(L);
    Rin = sigma1^2*(a1*a1')+sigma2^2*(a2*a2')+sigma^2*eye(L);
    
    for i = 1:N
        st0=(randn(1,K)+1i*randn(1,K))/sqrt(2)*sigma0;
        st1=(randn(1,K)+1i*randn(1,K))/sqrt(2)*sigma1;
        st2=(randn(1,K)+1i*randn(1,K))/sqrt(2)*sigma2;
        nt =(randn(L,K)+1i*randn(L,K))/sqrt(2)*sigma;
        
        xt = a0*st0+a1*st1+a2*st2+nt;
        sRxx = xt*xt'/ K;
        
        [eive,eiva] = eig(sRxx);         %eive是特征向量，eiva是特征值
        [temp,order] = sort(diag(eiva)); %order排序后对应原向量元素位置；temp 排序后的元素；从小到大排序
        sigma_e = sqrt(sum(temp(1:L-M)));
        Us = eive(:,order(L-M+1:L));
        Rin_e = a1*a1'/(a1'*inv(sRxx)*a1) + a2*a2'/(a2'*inv(sRxx)*a2) + sigma_e^2*eye(L);%?小特征值对应噪声功率

        w_mpdr = inv(Rxx)*a0 /(a0'*inv(Rxx)*a0);
        w_smi  = inv(sRxx)*b0/(b0'*inv(sRxx)*b0);
        w_rmi  = inv(Rin_e)*b0/(b0'*inv(Rin_e)*b0);
        w_dl   = inv(sRxx+1*sigma_e^2*eye(L))*b0/(b0'*inv(sRxx+1*sigma_e^2*eye(L))*b0); %rou=10
        w_sip  = inv(sRxx)*Us*Us'*b0;
        
        osinr_mpdr(m,i)=sigma0^2*(abs(w_mpdr'*a0))^2/abs(w_mpdr'*Rin*w_mpdr);
        osinr_smi(m,i) =sigma0^2*(abs(w_smi'*a0))^2/abs(w_smi'*Rin*w_smi);
        osinr_rmi(m,i) =sigma0^2*(abs(w_rmi'*a0))^2/abs(w_rmi'*Rin*w_rmi);
        osinr_dl(m,i)  =sigma0^2*(abs(w_dl'*a0))^2/abs(w_dl'*Rin*w_dl);
        osinr_sip(m,i) =sigma0^2*(abs(w_sip'*a0))^2/abs(w_sip'*Rin*w_sip);
        
        a0_mpdr=Rxx*inv(Rxx)*a0;
        a0_smi=Rxx*inv(sRxx)*b0;
        a0_rmi=Rxx*inv(Rin_e)*b0;
        a0_dl=Rxx*inv(sRxx+1*sigma_e^2*eye(L))*b0;
        a0_sisp=Rxx*inv(sRxx)*Us*Us'*b0;

        angle_mismatched_mpdr(m)=angle_mismatched_mpdr(m)+(abs(a0_mpdr'*inv(Rin)*a0))^2/((a0_mpdr'*inv(Rin)*a0_mpdr)*(a0'*inv(Rin)*a0));
        angle_mismatched_smi(m)=angle_mismatched_smi(m)+(abs(a0_smi'*inv(Rin)*a0))^2/((a0_smi'*inv(Rin)*a0_smi)*(a0'*inv(Rin)*a0));
        angle_mismatched_rmi(m)=angle_mismatched_rmi(m)+(abs(a0_rmi'*inv(Rin)*a0))^2/((a0_rmi'*inv(Rin)*a0_rmi)*(a0'*inv(Rin)*a0));
        angle_mismatched_dl(m)=angle_mismatched_dl(m)+(abs(a0_dl'*inv(Rin)*a0))^2/((a0_dl'*inv(Rin)*a0_dl)*(a0'*inv(Rin)*a0));
        angle_mismatched_sisp(m)=angle_mismatched_sisp(m)+(abs(a0_sisp'*inv(Rin)*a0))^2/((a0_sisp'*inv(Rin)*a0_sisp)*(a0'*inv(Rin)*a0));
    end
    osinr_mpdr(m,1) = sum(osinr_mpdr(m,:))/N ;
    osinr_smi(m,1) = sum(osinr_smi(m,:))/N ;
    osinr_rmi(m,1) = sum(osinr_rmi(m,:))/N ;
    osinr_dl(m,1) = sum(osinr_dl(m,:))/N ;
    osinr_sip(m,1) = sum(osinr_sip(m,:))/N ;
    
    angle_mismatched_mpdr(m)= abs(angle_mismatched_mpdr(m)/N);
    angle_mismatched_mpdr(m)= acos(sqrt(angle_mismatched_mpdr(m)))/pi*180;

    angle_mismatched_smi(m)=angle_mismatched_smi(m)/N;
    angle_mismatched_smi(m)= acos(sqrt(angle_mismatched_smi(m)))/pi*180;      

    angle_mismatched_rmi(m)= angle_mismatched_rmi(m)/N;
    angle_mismatched_rmi(m)= acos(sqrt(angle_mismatched_rmi(m)))/pi*180;  

    angle_mismatched_dl(m)=angle_mismatched_dl(m)/N;
    angle_mismatched_dl(m)= acos(sqrt(angle_mismatched_dl(m)))/pi*180;  

    angle_mismatched_sisp(m)=angle_mismatched_sisp(m)/N;
    angle_mismatched_sisp(m)= acos(sqrt(angle_mismatched_sisp(m)))/pi*180;
    
end

%绘制输出信干噪比随输入信噪比变化
figure(1);hold on;
plot(snr,10*log(osinr_mpdr(:,1)),'--','linewidth',2);
plot(snr,10*log(osinr_smi(:,1)),'-.','linewidth',2);
plot(snr,10*log(osinr_rmi(:,1)),'-*','linewidth',2);
plot(snr,10*log(osinr_dl(:,1)),'-o','linewidth',2);
plot(snr,10*log(osinr_sip(:,1)),'-+','linewidth',2);
xlabel('SNR/dB');ylabel('OSINR/dB');
title('输出信干噪比随输入信噪比变化趋势');
legend('MPDR','SMI','RMI','DL','SIP');

%输出失配广义夹角随输入信噪比变化
figure(2);hold on;
plot(snr,angle_mismatched_mpdr,'--','linewidth',2);
plot(snr,abs(angle_mismatched_smi),'-.','linewidth',2);
plot(snr,abs(angle_mismatched_rmi),'-*','linewidth',2);
plot(snr,abs(angle_mismatched_dl),'-o','linewidth',2);
plot(snr,abs(angle_mismatched_sisp),'-+','linewidth',2);
xlabel('SNR/dB');ylabel('失配广义夹角/(^o)');
title('失配广义夹角随输入信噪比变化趋势');
legend('MPDR','SMI','RMI','DL','SIP');

%% DL
clear,clc,close all;
L=16;
K=50;
d=1;
M=3;
sigma=1;
snr=-10:30;
N = 500;
%初始化输出信干噪比矩阵
osinr_mpdr=zeros(length(snr),N);
osinr_smi=zeros(length(snr),N);
osinr_rmi=zeros(length(snr),N);
osinr_dl=zeros(length(snr),N);
osinr_sip=zeros(length(snr),N);

%初始化失配广义夹角矩阵
angle_mismatched_mpdr=zeros(1,length(snr));
angle_mismatched_smi=zeros(1,length(snr));
angle_mismatched_rmi=zeros(1,length(snr));
angle_mismatched_dl=zeros(1,length(snr));
angle_mismatched_sisp=zeros(1,length(snr));


theta0=0;
theta1=60;
theta2=-60;
a0=exp(1i*2*pi*d*sin(theta0*pi/180)*(0:L-1)');
a1=exp(1i*2*pi*d*sin(theta1*pi/180)*(0:L-1)');
a2=exp(1i*2*pi*d*sin(theta2*pi/180)*(0:L-1)');
b0=exp(1i*2*pi*d*sin((theta0+2)*pi/180)*(0:L-1)');
for m=1:length(snr)
    sigma0=sigma*10^(snr(m)/20);
    %两个干扰源，信干噪比为-30dB
    inr1=snr(m)+30;
    inr2=snr(m)+30;
    sigma1=sigma*10^(inr1/20);
    sigma2=sigma*10^(inr2/20);
        
    Rxx = sigma0^2*(a0*a0')+sigma1^2*(a1*a1')+sigma2^2*(a2*a2')+sigma^2*eye(L);
    Rin = sigma1^2*(a1*a1')+sigma2^2*(a2*a2')+sigma^2*eye(L);
    
    for i = 1:N
        st0=(randn(1,K)+1i*randn(1,K))/sqrt(2)*sigma0;
        st1=(randn(1,K)+1i*randn(1,K))/sqrt(2)*sigma1;
        st2=(randn(1,K)+1i*randn(1,K))/sqrt(2)*sigma2;
        nt =(randn(L,K)+1i*randn(L,K))/sqrt(2)*sigma;
        
        xt = a0*st0+a1*st1+a2*st2+nt;
        sRxx = xt*xt'/ K;
        
        [eive,eiva] = eig(sRxx);         %eive是特征向量，eiva是特征值
        [temp,order] = sort(diag(eiva)); %order排序后对应原向量元素位置；temp 排序后的元素；从小到大排序
        sigma_e = sqrt(sum(temp(1:L-M)));

        w_mpdr = inv(sRxx+1*sigma_e^2*eye(L))*b0/(b0'*inv(sRxx+1*sigma_e^2*eye(L))*b0);
        w_smi  = inv(sRxx+5*sigma_e^2*eye(L))*b0/(b0'*inv(sRxx+5*sigma_e^2*eye(L))*b0);
        w_rmi  = inv(sRxx+10*sigma_e^2*eye(L))*b0/(b0'*inv(sRxx+10*sigma_e^2*eye(L))*b0);
        w_dl   = inv(sRxx+15*sigma_e^2*eye(L))*b0/(b0'*inv(sRxx+15*sigma_e^2*eye(L))*b0);
        w_sip  = inv(sRxx+20*sigma_e^2*eye(L))*b0/(b0'*inv(sRxx+20*sigma_e^2*eye(L))*b0);
        
        osinr_mpdr(m,i)=sigma0^2*(abs(w_mpdr'*a0))^2/abs(w_mpdr'*Rin*w_mpdr);
        osinr_smi(m,i) =sigma0^2*(abs(w_smi'*a0))^2/abs(w_smi'*Rin*w_smi);
        osinr_rmi(m,i) =sigma0^2*(abs(w_rmi'*a0))^2/abs(w_rmi'*Rin*w_rmi);
        osinr_dl(m,i)  =sigma0^2*(abs(w_dl'*a0))^2/abs(w_dl'*Rin*w_dl);
        osinr_sip(m,i) =sigma0^2*(abs(w_sip'*a0))^2/abs(w_sip'*Rin*w_sip);
        
        a0_mpdr=Rxx*inv(sRxx+1*sigma_e^2*eye(L))*b0;
        a0_smi=Rxx*inv(sRxx+5*sigma_e^2*eye(L))*b0;
        a0_rmi=Rxx*inv(sRxx+10*sigma_e^2*eye(L))*b0;
        a0_dl=Rxx*inv(sRxx+15*sigma_e^2*eye(L))*b0;
        a0_sisp=Rxx*inv(sRxx+20*sigma_e^2*eye(L))*b0;

        angle_mismatched_mpdr(m)=angle_mismatched_mpdr(m)+(abs(a0_mpdr'*inv(Rin)*a0))^2/((a0_mpdr'*inv(Rin)*a0_mpdr)*(a0'*inv(Rin)*a0));
        angle_mismatched_smi(m)=angle_mismatched_smi(m)+(abs(a0_smi'*inv(Rin)*a0))^2/((a0_smi'*inv(Rin)*a0_smi)*(a0'*inv(Rin)*a0));
        angle_mismatched_rmi(m)=angle_mismatched_rmi(m)+(abs(a0_rmi'*inv(Rin)*a0))^2/((a0_rmi'*inv(Rin)*a0_rmi)*(a0'*inv(Rin)*a0));
        angle_mismatched_dl(m)=angle_mismatched_dl(m)+(abs(a0_dl'*inv(Rin)*a0))^2/((a0_dl'*inv(Rin)*a0_dl)*(a0'*inv(Rin)*a0));
        angle_mismatched_sisp(m)=angle_mismatched_sisp(m)+(abs(a0_sisp'*inv(Rin)*a0))^2/((a0_sisp'*inv(Rin)*a0_sisp)*(a0'*inv(Rin)*a0));
    end
    osinr_mpdr(m,1) = sum(osinr_mpdr(m,:))/N ;
    osinr_smi(m,1) = sum(osinr_smi(m,:))/N ;
    osinr_rmi(m,1) = sum(osinr_rmi(m,:))/N ;
    osinr_dl(m,1) = sum(osinr_dl(m,:))/N ;
    osinr_sip(m,1) = sum(osinr_sip(m,:))/N ;
    
    angle_mismatched_mpdr(m)= abs(angle_mismatched_mpdr(m)/N);
    angle_mismatched_mpdr(m)= acos(sqrt(angle_mismatched_mpdr(m)))/pi*180;

    angle_mismatched_smi(m)=angle_mismatched_smi(m)/N;
    angle_mismatched_smi(m)= acos(sqrt(angle_mismatched_smi(m)))/pi*180;      

    angle_mismatched_rmi(m)= angle_mismatched_rmi(m)/N;
    angle_mismatched_rmi(m)= acos(sqrt(angle_mismatched_rmi(m)))/pi*180;  

    angle_mismatched_dl(m)=angle_mismatched_dl(m)/N;
    angle_mismatched_dl(m)= acos(sqrt(angle_mismatched_dl(m)))/pi*180;  

    angle_mismatched_sisp(m)=angle_mismatched_sisp(m)/N;
    angle_mismatched_sisp(m)= acos(sqrt(angle_mismatched_sisp(m)))/pi*180;
    
end

%绘制输出信干噪比随输入信噪比变化
figure(1);hold on;
plot(snr,10*log(osinr_mpdr(:,1)),'--','linewidth',2);
plot(snr,10*log(osinr_smi(:,1)),'-.','linewidth',2);
plot(snr,10*log(osinr_rmi(:,1)),'-*','linewidth',2);
plot(snr,10*log(osinr_dl(:,1)),'-o','linewidth',2);
plot(snr,10*log(osinr_sip(:,1)),'-+','linewidth',2);
xlabel('SNR/dB');ylabel('OSINR/dB');
title('输出信干噪比随输入信噪比变化趋势');
legend('DL ρ=1','DL ρ=5','DL ρ=10','DL ρ=15','DL ρ=20');

%输出失配广义夹角随输入信噪比变化
figure(2);hold on;
plot(snr,angle_mismatched_mpdr,'--','linewidth',2);
plot(snr,abs(angle_mismatched_smi),'-.','linewidth',2);
plot(snr,abs(angle_mismatched_rmi),'-*','linewidth',2);
plot(snr,abs(angle_mismatched_dl),'-o','linewidth',2);
plot(snr,abs(angle_mismatched_sisp),'-+','linewidth',2);
xlabel('SNR/dB');ylabel('失配广义夹角/(^o)');
title('失配广义夹角随输入信噪比变化趋势');
legend('DL ρ=1','DL ρ=5','DL ρ=10','DL ρ=15','DL ρ=20');