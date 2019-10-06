
angle = -5:0.001:5;
phi=0.1;%squint angle, o
y1 = sinc(angle - phi);
y2 = sinc(angle + phi);
ysum = y1 + y2;
ydif = y1 - y2;

figure(1);hold on;grid on;
plot(angle,y1);
plot(angle,y2);
xlabel('Angle/^o');ylabel('Amplitude');
title('Squinted patterns');
figure(2);grid on;
subplot(211);plot(angle,ysum);
xlabel('Angle/rad');ylabel('Amplitude');
title('Sum pattern');
subplot(212);plot(angle,ydif);
xlabel('Angle/rad');ylabel('Amplitude');
title('Difference pattern');

dovrs = ydif./ysum;
figure(3);grid on;
plot(angle,dovrs);xlim([-1, 1]);ylim([-3, 3]);
xlabel('Angle/rad');ylabel('Amplitude');
title('difference-to-sum ratio');

%% squint angle change 
angle = -5:0.001:5;
phi=0.1:0.2:0.9;%squint angle, o
for m=1:length(phi)
y1 = sinc(angle - phi(m));
y2 = sinc(angle + phi(m));
ysum = y1 + y2;
ydif = y1 - y2;

dovrs = ydif./ysum;
plot(angle,dovrs);
hold on;
end
xlim([-1, 1]);ylim([-3, 3]);grid on;
xlabel('Angle/rad');ylabel('Amplitude');
title('difference-to-sum ratio');