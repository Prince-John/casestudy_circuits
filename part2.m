%% Part 2 RL Circuit Simulation.

t0= 0; %Start Time in seconds 
h = 0.0001;% time step in seconds
R = 1e2;% Resistance
L = 0.1;%induction
tf= 0.005; %End Time in seconds
t = linspace(t0,tf,((tf-t0)/h) +1);
V_in = 1;
V_L_0 = 1;
V_in0 = V_in*ones(1,length(t)) ; % input voltage source
V_in0(1,1) = 1; %from the handout; voltage is accross inductor + resistance
V_L = zeros(1, length(t));
V_L(1,1) = V_L_0;

%% equation: (1-hR/L)*i,k + (h/L)V_in,k
B = h/L;
a = R*B;
C = 1-a;
V = [V_L;V_in0];
A = [C B];
for j = 2:length(t)
V(1,j) = A*V(:,j-1);
end
figure();
plot(t, V);
title('Inductance vs Time with constant V_{in} = 1V');
legend('V_L','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');

