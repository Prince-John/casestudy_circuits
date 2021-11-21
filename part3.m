%% Part 3 RLC Circuit Simulation.
R = 100;% Resistance
C = 1e-6;% Capacitance
L = 10e-3;% Inductance
t0= 0; %Start Time in seconds 
tf= 0.015; %End Time in seconds
h = 1/(192e3);% time step in seconds
V_in_0 = 1; 
V_C_0 = 0;
I_L_0 = 0;


% Time Series for the circuit elements
t = linspace(t0,tf,((tf-t0)/h) +1);
V_in = V_in_0*ones(1,length(t)) ; % input voltage source
V_C = zeros(1, length(t));
I_L = zeros(1, length(t));
V_C(1,1) = V_C_0;
I_L(1,1) = I_L_0;

%State Vector
x = [V_C; I_L];
u = V_in; 

%Transmission matrix and vector

A = [1 (h/C); (-h/L) (1-(h*R/L))];
B = [0; (h/L)];


%% Model loop

for j = 2:length(t)
    x(:,j) = A*x(:,j-1) + B*u(1, j-1);
end

V_R = R*x(2,:);

%% Plots

figure();
plot(t, V_R);
title('Part 3 test');
legend('V_R');
ylabel('Voltage (V)');
xlabel('Time (s)')

%% Exploring VR out with differnt values of RLC,

% using a decade function found online to genrate log scale component test
% values.

R_Values = decade(1,6,3);
L_Values = decade(-4, 1,3);
C_Values = decade(-9, -5,3);
V_R_values = zeros(length(R_Values)*length(L_Values)*length(C_Values), length(t));
V_R_labels = zeros(length(R_Values)*length(L_Values)*length(C_Values), 3);
index = 1;
for i = R_Values
   for j = L_Values
      for k = C_Values
        V_R_values(index,:) = simV_R(i,j,k,V_in,t); 
        V_R_labels(index,:) = [i j k];
        index = index +1;
      end
   end          
end


%%

testV = 500.*t.*exp(-2500.*t);




