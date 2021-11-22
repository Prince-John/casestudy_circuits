%% Part 1 RC Circuit Simulation.
R = 1e3;% Resistance
C = 1e-6;% Capacitance
t0= 0; %Start Time in seconds 
tf= 6*R*C; %End Time in seconds
h = 0.0001;% time step in seconds
V_in_0 = 1; 
V_C_0 = 0;

% Time Series
t = linspace(t0,tf,((tf-t0)/h) +1);
V_in = V_in_0*ones(1,length(t)) ; % input voltage source
V_C = zeros(1, length(t));
V_C(1,1) = V_C_0;



%% Equation to implement: V_c,k+1 = (1 - h/RC)*V_c,k + (h/RC)*V_in,k
%                                 =  a&V_c,k + b*V_in,k
b = h/(R*C);
a = 1-b;

%%
% transformation matrix 
% Keeping all voltages in a Voltage vector V, here we only have 2 voltages,
% V_in and V_C; V = (V_C; V_in)
V = [V_C; V_in];
A = [a b];
for i = 2:length(t)
    V(1,i) = A*V(:,i-1);
end
save ('container_values', 'V','t', '-append');
%% Plotting and Analysis

% *Task 1.2 - 1* 

figure();
plot(t, V);
title('Capacitor Voltage vs Time with constant V_{in} = 1V');
legend('V_C','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');

% *Task 1.2 -2*

% Using capVoltage function to test with a few different h values.
% four test h values, h<<RC, h<RC, h = RC<, h > RC

%%
Vtest = [zeros(1, length(t)); V_in];
Vt1 = capVoltage(0.00001, R, C, Vtest);
Vt2 = capVoltage(0.0001, R, C, Vtest);
Vt3 = capVoltage(0.001, R, C, Vtest);
Vt4 = capVoltage(0.01, R, C, Vtest);

figure();
subplot(2,2,1);
plot(t, Vt1);
title('V_C vs Time, h<<RC');
legend('V_C','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');
subplot(2,2,2);
plot(t, Vt2);
title('V_C vs Time, h<RC');
legend('V_C','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');
subplot(2,2,3);
plot(t, Vt3);
title('V_C vs Time, h=RC');
legend('V_C','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');
subplot(2,2,4);
plot(t, Vt4);
title('V_C vs Time, h>RC');
legend('V_C','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');

% When the value of h is equal to or greater than RC the curve is
% incorrect. The curve is also incorrect when h is very very less than RC.

% Theoretical Vc = 1 - exp(-t/RC)
%%
Vtheoretical = 1 - exp(-t/(R*C));

figure();
hold on;
plot(t, Vt3);
plot(t, Vt2(1,:));
plot(t, Vtheoretical);
hold off;
legend('Incorrect h','V_{in}', 'Correct h', 'Theoretical V_C');
ylabel('Voltage (V)');
xlabel('Time (s)');
title('Various capacitor charging curves vs time');



% *Task 1.2 -3*

% $\tau = RC$ 






