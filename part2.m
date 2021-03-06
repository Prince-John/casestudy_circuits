%% Part 2 RL Circuit Simulation.

t0= 0; %Start Time in seconds 
h = 0.0001;% time step in seconds
R = 1e2;% Resistance
L = 0.1;%induction
tf= 0.005; %End Time in seconds
t = linspace(t0,tf,((tf-t0)/h) +1);
V_in_0 = 1;
V_L_0 = V_in_0; % at time zero the current is zero and no voltage is dropped across the resistor.
I_L_0 = 0;
V_in = V_in_0.*ones(1,length(t)) ; % input voltage source
V_L = zeros(1, length(t));
V_L(1,1) = V_L_0;
I_L = zeros(1, length(t));
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
for j = 2:length(t)
container(1,j) =  A*container([1 3], j-1);
container(2,j) = current_matix*container([1 3], j); % Calculates the voltage, V_L = V_I - V_R; V_R = RI
end
%DND..USED FOR PART2_2
save ('container_values', 'container', '-append');
%%

figure();
plot(t, container([2 3],:));
title('Inductance vs Time with constant V_{in} = 1V');
legend('V_L','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');
%% Part 2: RL circuit Simulation

t0= 0; %Start Time in seconds 
R = 1e3;% Resistance
L = 1e-1;% Inductance
tf= 0.005; %End Time in seconds
h = 0.00001;% time step in seconds
V_in_0 = 1; 
V_L_0 = 0;

% Time Series
t = linspace(t0,tf,((tf-t0)/h) +1);
V_L_in = V_in_0*ones(1,length(t)) ; % input voltage source
i_k = zeros(1, length(t));
% V_C(1,1) = V_C_0;

%% Equation to implement: V_c,k+1 = (1 - Rh/L)*V_c,k + (h/L)*V_in,k
%                                 =  a*i_k + b*V_in,k
b = h/L;
a = 1-((R*h)/L);

%% Apply it to each case and graph each case and the result
% make 4 cases
A_container = zeros(1,4);
B_container = zeros(1,4);
h_test_container = [0.00001,0.0001,0.001,0.01];
container = [i_k; V_L_in];
result = [];
for i = 1:4
    A_container(:,i) = 1-(R*h_test_container(i)/L);
    B_container(:,i) = h_test_container(i)/L;
    constant_container = [A_container(:,i) B_container(:,i)];
    Vt = RL_simulation(constant_container, container);
    result = [result; Vt];
end
figure();
plot(result(1,:),result(2,:));

%% Implement the equation
% transformation matrix 
% we only have 1 voltage and 1 current, we can keep them in a vector to do
% the dot product between a vector of coefficients and a vector including
% current and voltage
% V_in; 

function current = RL_simulation(A, container)
    for i = 2:length(container)
        container(1,i) = A*container(:,i-1);
    end
    current = container;
end
