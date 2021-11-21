R = 1e2;% Resistance for L
RC = 1e3;%Resistance C
L = 0.1;%induction
C = 1e-6;% Capacitance
tau = 0.0005;
%where i is C*(dv/dt)
%%
%LOAD THE CONTAINER VALUES.MAT FILE
load('container_values.mat');%loads data for inductor effect on voltage and current
%IL_val = container(1,:);%current from inductor
IC_val = ones(1, length(container));%capacitance current from volt
IC_val(1,:) = V(1,:)/RC;
Il_val = ones(1, length(container));%inductor current from volt
IL_val(1,:) = container(2,:)/R;

%Max_I = max(I_val);
VL_val = container(1,:);%voltage from inductor
VC_val = V(1,:);%voltage from capacitance
Max_Vc = max(VC_val);%max voltage from capacitance
Max_VC = max(VL_val);%max voltage from inductor

IL = ones(1, length(container));
IL(1,:) = (Max_VC/R)+(Max_Vc/RC);%MAX CURRENT
MaxV = ones(1, length(container));
MaxV(1,:) = (Max_Vc/RC)+(Max_VC/R);
%% PLOT
figure;
plot(t, MaxV);
hold on
plot(t,IL_val);
plot(t,IC_val);
legend('max current','current behavior accross inductor', 'current behavior accross capacitor');