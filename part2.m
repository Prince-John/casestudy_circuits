%% Part 2 RL Circuit Simulation.

t0= 0; %Start Time in seconds 
h = 0.0001;% time step in seconds
R = 1e2;% Resistance
L = 0.1;%induction
tf= 0.006; %End Time in seconds
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


%%
figure();
plot(t, container([2 3],:));
title('Inductance vs Time with constant V_{in} = 1V');
legend('V_L','V_{in}');
ylabel('Voltage (V)');
xlabel('Time (s)');

%% Complementart relationship
% This function return a plot to depict the complementary relationship
% between inductor and capacitor. This function assumes that the initial
% voltage is 1V. It takes in 6 params: capcacitor's resistor, inductor's
% resistor, time_step, final_time, capacitance and inductance
capacitor_inductor_relationship(1e3,1e2,0.0001,0.006,1e-6,0.1);
function capacitor_inductor_relationship(R_c,R_l,time_step,final_time,C,L)
    t0= 0; %Start Time in seconds 
    h = time_step;% time step in seconds
    tf= final_time; %End Time in seconds
    V_in_0 = 1;
    V_L_0 = V_in_0; % at time zero the current is zero and no voltage is dropped across the resistor.
    I_L_0 = 0;
    V_in_0 = 1; 
    % Time Series
    t = linspace(t0,tf,((tf-t0)/h) +1);
    V_in = V_in_0*ones(1,length(t)) ; % input voltage source
    % configure the capacitance circuit
    Vtest = [zeros(1, length(t)); V_in];
    % configure the inductor circuit
    V_L = zeros(1, length(t));
    V_L(1,1) = V_L_0;
    I_L = zeros(1, length(t));
    I_L(1,1) = I_L_0;
    % calculate voltage values of capacitor
    V_C_result = capVoltage(h,R_c,C,Vtest);
    % calculate current values of inductor
    container = [I_L; V_L ;V_in];    
    V_L_result = inductorVoltage(h,R_l,L,container);
    figure();
    hold on;
    plot(t, V_L_result);
    plot(t, V_C_result);
    legend('V_L','V_{in}','V_C');
    ylabel('Voltage (V)');
    xlabel('Time (s)');
    title("Complementary relationship between capacitor and inductor");
end

function voltage_inductor = inductorVoltage(h,R,L,container)
    b = h/L;
    a = R*b;
    c = 1-a;
    A = [c b]; 
    current_matix = [(-1*R) 1];
    for j = 2:length(container(1,:))
        container(1,j) =  A*container([1 3], j-1);
        % Calculates the voltage, V_L = V_I - V_R; V_R = RI
        container(2,j) = current_matix*container([1 3], j);
    end
    voltage_inductor = container([2 3],:);
end
