%% Part 2 RL Circuit Simulation.

t0= 0; %Start Time in seconds 
h = 0.0001;% time step in seconds
R = 1e2;% Resistance
L = 0.1;%induction
tf= 0.006; %End Time in seconds
tL = linspace(t0,tf,((tf-t0)/h) +1);
V_in_0 = 1;
V_L_0 = V_in_0; % at time zero the current is zero and no voltage is dropped across the resistor.
I_L_0 = 0;


V_in = V_in_0.*ones(1,length(tL)) ; % input voltage source
V_L = zeros(1, length(tL));
V_L(1,1) = V_L_0;
I_L = zeros(1, length(tL));
I_L(1,1) = I_L_0;


%% equation: i(k+1) = (1-hR/L)*i,k + (h/L)V_in,k
% To vectorize the eq, it becomes: i(k+1) = C*i_k + B*V_in_k
b = h/L;
a = R*b;
c = 1-a;

%Container called |container| to house the currents and v_in

container = [I_L; V_L ;V_in];


%%
A = [c b]; 
current_matix = [(-1*R) 1];
for j = 2:length(tL)
container(1,j) =  A*container([1 3], j-1);
container(2,j) = current_matix*container([1 3], j); % Calculates the voltage, V_L = V_I - V_R; V_R = RI
end
%DND..USED FOR PART2_2
save ('container_values', 'container', '-append');
%%
figure();
plot(tL, container([2 3],:));
title('Inductance vs Time with constant V_{in} = 1V');
legend('V_L','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');