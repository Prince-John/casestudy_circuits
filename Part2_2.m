
RL = 1e2;% Resistance for L
RC = 1e3;%Resistance C
L = 0.1;%induction
C = 1e-6;% Capacitance
load('container_values.mat');
%where i is C*(dv/dt)
%% GET THE CURRENT FLOWING THRROUGH THE CAPACITOR USING
% use VR = V_in - VC
V_in_C = V(2,:);
VC = V(1,:);
VR = V_in_C - VC;
%find current accross the resistor using
% I = VR/R
IC = VR/RC;
%% CURRNETS AND VOLTAGES OF THROUGH INDUCTOR
IL = container(1,:);
VL = container(2,:);
V_in_L = container(2,:);
%% PLOT
figure;
subplot(4,1,1);
plot(t, VC);
legend ('voltage accross the capacitor');
xlabel('time');
ylabel('voltage');
hold on
subplot(4,1,3);
plot(t, IC);
xlabel('time');
ylabel('current');
legend ('current accross the resistor');
subplot(4,1,2);
plot(t, VL);
xlabel('time');
ylabel('voltage');
legend ('voltage accross the inductor');
subplot(4,1,4);
plot(t, IL);
xlabel('time');
ylabel('current');
legend ('current accross the inductor');
sgtitle('Complementary characteristics between Capacitors and Inductors');