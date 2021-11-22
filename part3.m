%% Part 3 RLC Circuit Simulation.
R = 100;% Resistance
C = 1e-6;% Capacitance
L = 10e-3;% Inductance
t0= 100*h; %Start Time in seconds 
tf= 0.015; %End Time in seconds
h = 1/(192e3);% time step in seconds
V_in_0 = 1; 
V_C_0 = 0;
I_L_0 = 0;


% Time Series for the circuit elements
t = linspace(0,tf,((tf)/h) +1);
V_in = V_in_0*ones(1,length(t)) ; % input voltage source
V_in(1:t0/h) = 0;
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

R_Values = decade(2,3,0.5);
L_Values = decade(-4, 2,3);
C_Values = decade(-9, -5,3);

% To systematically evaluate the response, we will vary one of the
% parameter values from the decades while keeping the other two constant.

% Varying R with L = 100mH and C = 10nF
V_R_values = zeros(length(R_Values), length(t));
index = 1;
figure();
for i = R_Values
    V_R_values(index,:) = simV_R(i, 100e-3, 10e-9,V_in, t);
    subplot(length(R_Values),1,index);
    plot(t, V_R_values(index,:));
    title("R = "+int2str(i)+" \Omega");
    index = index+1;
end

% From this it can be seen that the flip of oscillation growth to decay starts b/w 500
% \Omega and 550 \Omega. 
%%
% To evaluate that further, evaluating from R = 500 Ohm to 550 Ohm with 10
% Ohm increment.
figure()
R_Values_test = [500 510 520 530 540 550];
index =1;
for i = R_Values_test
    subplot(length(R_Values_test),1,index);
    plot(t, simV_R(i, 100e-3, 10e-9,V_in, t));
    title("R = "+int2str(i)+" \Omega");
    index = index+1;
end

%% Finding relation on L
% Using a similar methodology as above, we will vary L with R = 520 Ohm and
% C = 10 nF


V_L_values = zeros(length(L_Values), length(t));
index = 1;
figure();
L_Values = [0.01 0.05 0.10 0.20 0.50];
for i = L_Values
    V_L_values(index,:) = simV_R(520, i, 10e-9,V_in, t);
    subplot(length(L_Values),1,index);
    plot(t, V_L_values(index,:));
    title("L = "+num2str(i)+" H");
    index = index+1;
end

% Varying L changes the frequency of the oscilation but does not have any
% appretiable affect on the decay or growth of the amplitude. L inversly
% proportional to freq. 


%% Finding relation on C
% Using a similar methodology as above, we will vary C with R = 520 Ohm and
% L = 0.500 H


V_C_values = zeros(length(L_Values), length(t));
index = 1;
figure();
%C_Values = [0.01 0.05 0.10 0.20 0.50];
for i = C_Values
    V_C_values(index,:) = simV_R(520, 0.5, i,V_in, t);
    subplot(length(C_Values),1,index);
    plot(t, V_C_values(index,:));
    title("C = "+num2str(i)+" H");
    index = index+1;
end


% Varying the C value like R changes the amplitude shape, from infinite
% growth to decay, It also changes the frequency higher C corresponding to
% lower freq. Furthemore, a high value of C does not always cause the
% response to decay faster. The decay time gets smaller to a point and then 
% get higher again.

%% Tuning the Circuit Part 3-2.

% Using the relations discovered in the graphs above, we try to get the
% circuit to have a quickly decaying response, sustained response, and
% unstable growing response.

sustained_response = simV_R(521, 100e-2, 10e-9,V_in, t);
decaying_response = simV_R(2000, 100e-2, 10e-9,V_in, t);
growing_response = simV_R(100, 100e-2, 10e-9,V_in, t);

figure();
plot(t,[V_in; sustained_response; decaying_response; growing_response]);
legend('V_in','Sustained R=521\Omega, L = 1 H, C = 10nF','Decaying R=2000\Omega, L = 1 H, C = 10nF','Unstable R=100\Omega, L = 1 H, C = 10nF')
ylabel('Voltage (V)');
xlabel('time (s)');


% Upon playing the signals using soundsc, the higher frequency signals
% sound a note of a higher pitch. The pitch of the sound is inversely
% proportional to the inductor value. 

% The sustained signal had a uniform loudness and played a constant note.
% The growing signal increased in loudness with time.
% The decaying signal decreased in loudness with time. 





