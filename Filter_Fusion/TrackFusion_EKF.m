

clear,clc;
close all;
%% EKF Track Fusion 2D
T_total = 200;       %Observation time s
T= 1;             
N = T_total/T;
t = T:T:T_total;
M = 1;              %Monto-carlo time
%Motion parameters
Rx = 50;
Ry = -100;
vx = -1;
vy = 2;
%观测站的位置
x1=-10;
y1=0;
x2=10;
y2=0;

r1 = sqrt((Rx-x1)^2+(Ry-y1)^2);
beta1 = atan2((Ry-y1),(Rx-x1));
r2 = sqrt((Rx-x2)^2+(Ry-y2)^2);
beta2 = atan2((Ry-y2),(Rx-x2));

%noise
sigma_u = sqrt(0.0001);     %过程噪声
sigma_R = sqrt(5);        %距离量测噪声
sigma_beta = sqrt(0.0001);    %角度量测噪声

%% Kalman filter CV 2D

%-------Kalman Parameters-------%

A1 = [cos(beta1) -r1*sin(beta1); sin(beta1) r1*cos(beta1)] ;
R1 = A1*[sigma_R^2 0;0 sigma_beta^2]*A1' ;
P1 = [R1(1,1)   R1(1,1)/T     R1(1,2)   R1(1,2)/T
    R1(1,1)/T 2*R1(1,1)/T^2 R1(1,2)/T 2*R1(1,2)/T^2
    R1(1,2)   R1(1,2)/T     R1(2,2)   R1(2,2)/T
    R1(1,2)/T 2*R1(1,2)/T^2 R1(2,2)/T 2*R1(2,2)/T^2 ];
% P1 = 100*eye(4);
A2 = [cos(beta2) -r2*sin(beta2); sin(beta2) r2*cos(beta2)] ;
R2 = A2*[sigma_R^2 0;0 sigma_beta^2]*A2' ;
P2 = [R2(1,1)   R2(1,1)/T     R2(1,2)   R2(1,2)/T
    R2(1,1)/T 2*R2(1,1)/T^2 R2(1,2)/T 2*R2(1,2)/T^2
    R2(1,2)   R2(1,2)/T     R2(2,2)   R2(2,2)/T
    R2(1,2)/T 2*R2(1,2)/T^2 R2(2,2)/T 2*R2(2,2)/T^2 ];
% P2 = 100*eye(4);
%状态转移矩阵
F = [1 T 0 0 
     0 1 0 0
     0 0 1 T
     0 0 0 1];

%过程噪声
B = [T^2/2, 0; T, 0;
     0, T^2/2; 0, T]; %过程噪声分布矩阵
v = sigma_u^2;   %x方向的过程噪声向量//相当于Q
V = B * v * B';
% %观测噪声??
% W = B * noise_x;

%------Data initial-------%
X_real = zeros(4,N);
X = zeros(4,N);

Z1 = zeros(2,N);
Z_polar1 = zeros(2,N);
X_EKF1 = zeros(4,N);
Z2 = zeros(2,N);
Z_polar2 = zeros(2,N);
X_EKF2 = zeros(4,N);

X_CC = zeros(4,N);
X_BC = zeros(4,N);
bias = zeros(8,N,M);

%-------Track Initial-------%
X_real(:,1) = [Rx, vx, Ry, vy]'; %x: km,km/s

X_EKF1(:,1) = [Rx, 0, Ry, 0];
X_EKF2(:,1) = [Rx, 0, Ry, 0];
X_CC(:,1) = [Rx, 0, Ry, 0];
X_BC(:,1) = [Rx, 0, Ry, 0];

%Monto-carlo
for m=1:M
    
    noise_u = randn(2,N).*sigma_u;
    noise_R = randn(1,N).*sigma_R; %过程噪声
    noise_beta = randn(1,N).*sigma_beta; %观测噪声
    
    %构造 真实轨迹X 与 观测轨迹Z 
    for n=2:N
%         if n == 30
%             X_real(2,n-1) = 1;
%         end
        X_real(:,n) = F * X_real(:,n-1);
    end
    X = X_real + B * noise_u;
    
    for n=1:N
        Z_polar1(1,n) = sqrt((X(1,n)-x1)^2 +(X(3,n)-y1)^2) + noise_R(n);
        Z_polar1(2,n) = atan2((X(3,n)-y1),(X(1,n)-x1)) + noise_beta(n);
        [Z1(1,n),Z1(2,n)]= pol2cart(Z_polar1(2,n), Z_polar1(1,n));
        Z_polar2(1,n) = sqrt((X(1,n)-x2)^2 +(X(3,n)-y2)^2) + noise_R(n);
        Z_polar2(2,n) = atan2((X(3,n)-y2),(X(1,n)-x2)) + noise_beta(n);
        [Z2(1,n),Z2(2,n)]= pol2cart(Z_polar2(2,n), Z_polar2(1,n));
    end
    %这里可以写成function的形式
    P_BC = P1;
    for n=2:N
        x_predict = F * X_EKF1(:,n-1);                       %状态一步预测
        p_predict = F * P1 * F'+ V;                             %协方差一步预测
        r1 = sqrt((x_predict(1)-x1)^2 +(x_predict(3)-y1)^2);
        theta1 = atan2(x_predict(3)-y1,x_predict(1)-x1);
        %求解雅克比矩阵H	
        Hj1 =[(x_predict(1)-x1)/r1,0,(x_predict(3)-y1)/r1,0
             -(x_predict(3)-y1)/r1^2,0,(x_predict(1)-x1)/r1^2,0];
        Fx1 = [r1;theta1];
        S = Hj1 * p_predict * Hj1'+ R1;                       
        K1 = p_predict * Hj1'/ S ;                            
        X_EKF1(:,n) = x_predict + K1 * (Z_polar1(:,n) - Fx1);  
        P1 = (eye(4)-K1*Hj1) * p_predict;  

        x_predict2 = F * X_EKF2(:,n-1);                       
        p_predict2 = F * P2 * F'+ V;                             
        r2 = sqrt((x_predict2(1)-x2)^2 +(x_predict2(3)-y2)^2);
        theta2 = atan2((x_predict2(3)-y2),(x_predict2(1)-x2));
        %求解雅克比矩阵H	
        Hj2 =[(x_predict2(1)-x2)/r2,0,(x_predict2(3)-y2)/r2,0
             -(x_predict2(3)-y2)/r2^2,0,(x_predict2(1)-x2)/r2^2,0];
        Fx2 = [r2;theta2];
        S2 = Hj2 * p_predict2 * Hj2'+ R2;                     
        K2 = p_predict2 * Hj2'/ S2 ;                              
        X_EKF2(:,n) = x_predict2 + K2 * (Z_polar2(:,n) - Fx2); 
        P2 = (eye(4)-K2*Hj2) * p_predict2;  
        
        %Fusion
        P_CC = inv( inv(P1) + inv(P2));
        X_CC(:,n) = P_CC * (P1\X_EKF1(:,n) + P2\X_EKF2(:,n));
        
        P_BC = (eye(4)-K2*Hj2)* F*P_BC*F'*(eye(4)-K1*Hj1)';
        X_BC(:,n) = X_EKF2(:,n)+(P2-P_BC)/(P1+P2-2*P_BC)*(X_EKF1(:,n)-X_EKF2(:,n));
    end
    
    bias(1,:,m) = X(1,:) - X_EKF1(1,:);
    bias(2,:,m) = X(3,:) - X_EKF1(3,:);
    bias(3,:,m) = X(1,:) - X_EKF2(1,:);
    bias(4,:,m) = X(3,:) - X_EKF2(3,:);
    bias(5,:,m) = X(1,:) - X_CC(1,:);
    bias(6,:,m) = X(3,:) - X_CC(3,:);
    bias(7,:,m) = X(1,:) - X_BC(1,:);
    bias(8,:,m) = X(3,:) - X_BC(3,:);
end

Bias = sum(bias, 3) / M;
RMSE = sqrt(sum(bias.^2 , 3)/M);

%% Draw the Result
figure;
hold on;grid on;
plot(X(1,:),X(3,:),'LineWidth',1.5);
plot(X_EKF1(1,:),X_EKF1(3,:),'LineWidth',1.5);
plot(X_EKF2(1,:),X_EKF2(3,:),'LineWidth',1.5);
plot(X(1,1),X(3,1),'pb');
% plot(Z1(1,:),Z1(2,:),'.-');
legend('Real target track','Radar1: Filtered track','Radar2: Filtered track');
xlabel('x ');ylabel('y ');title('\fontsize{12} EKF Track Results of two Radars');

figure;
hold on;grid on;
plot(X(1,:),X(3,:),'LineWidth',1.5);
plot(X_CC(1,:),X_CC(3,:),'LineWidth',1.5);
plot(X_BC(1,:),X_BC(3,:),'LineWidth',1.5);
legend('Real track','CC Fusion track','BC Fusion track');
xlabel('x ');ylabel('y ');title('\fontsize{12} Fusion Tracks of CC & BC');

%----------------bias filtered and fusion----------
% figure;
% grid on;hold on;
% plot(bias(1,:),'LineWidth',1.5);
% plot(bias(2,:),'LineWidth',1.5);
% %ylim([0,0.4]);
% legend('x Bias','y Bias');
% xlabel('t/s');ylabel('Amplitude');title('Bias between Radar1 filtered track and real track');
% 
% figure;
% grid on;hold on;
% plot(bias(3,:),'LineWidth',1.5);
% plot(bias(4,:),'LineWidth',1.5);
% %ylim([0,0.4]);
% legend('x Bias','y Bias');
% xlabel('t/s');ylabel('Amplitude');title('Bias between Radar2 filtered track and real track');
% 
% figure;
% grid on;hold on;
% plot(bias(5,:),'LineWidth',1.5);
% plot(bias(6,:),'LineWidth',1.5);
% %ylim([0,0.4]);
% legend('x Bias','y Bias');
% xlabel('t/s');ylabel('Amplitude');title('Bias between CC Fusion track and real track');
% 
% figure;
% grid on;hold on;
% plot(bias(7,:),'LineWidth',1.5);
% plot(bias(8,:),'LineWidth',1.5);
% %ylim([0,0.4]);
% legend('x Bias','y Bias');
% xlabel('t/s');ylabel('Amplitude');title('Bias between BC Fusion track and real track');

%---------------RMSE fusion-------------
figure;
grid on;hold on;
plot(Bias(5,:),'LineWidth',1.5);
plot(Bias(6,:),'LineWidth',1.5);
legend('x Bias','y Bias');
xlabel('t/s');ylabel('Bias mean value');title('Bias of CC Fusion track and real track');
figure;
grid on;hold on;
plot(Bias(7,:),'LineWidth',1.5);
plot(Bias(8,:),'LineWidth',1.5);
legend('x Bias','y Bias');
xlabel('t/s');ylabel('Bias mean value');title('Bias of BC Fusion track and real track');
figure;
grid on;hold on;
plot(RMSE(5,:),'LineWidth',1.5);
plot(RMSE(6,:),'LineWidth',1.5);
legend('x Bias','y Bias');
xlabel('t/s');ylabel('RMSE');title('RMSE of CC Fusion track and real track');
figure;
grid on;hold on;
plot(RMSE(7,:),'LineWidth',1.5);
plot(RMSE(8,:),'LineWidth',1.5);
legend('x Bias','y Bias');
xlabel('t/s');ylabel('RMSE');title('RMSE of BC Fusion track and real track');

%%
figure;
hold on;grid on;
plot(X_real(2,:),X_real(4,:),'k');
plot(t, X(2,:),'LineWidth',1.5);
plot(t, X_EKF1(2,:),'LineWidth',1.5);
%ylim([0,1.5]);
legend('Real track','Filtered track');
xlabel('t/s');ylabel('v km/s');title('\fontsize{12} Velocity filtering track');



figure;
grid on;hold on;
plot(t(end-50:end),RMSE(1,end-50:end),'LineWidth',1.5);
plot(t(end-50:end),RMSE(2,end-50:end),'LineWidth',1.5);
%ylim([0,1]);legend('Range','Velocity');
xlabel('t/s');ylabel('Amplitude');title('RMSE');


