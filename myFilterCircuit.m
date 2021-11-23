%% Case study 3: Circuits as Resonators, Sensors, and Filters
% *ESE 105* 
%
% *Name: Quan Khuc, John Prince, Mordecai Ethapemi
%
% function myFilterCircuit(Vin,h) receives a time-series voltage sequence
% sampled with interval h, and returns the output voltage sequence produced
% by a circuit
%
% inputs:
% Vin - time-series vector representing the voltage input to a circuit
% h - scalar representing the sampling interval of the time series in
% seconds
%
% outputs:
% Vout - time-series vector representing the output voltage of a circuit
function Vout = myFilterCircuit(Vin,h)
vin1 = [Vin(1:length(Vin)/2,:)];
vin2 = [Vin(length(Vin)/2 :end,:)];
Vinl = low(Vin,h) + low2(Vin,h);
Vinh = high(Vin,h) + high2(Vin,h);
Vout = Vin - low(Vin,h)  - high(Vin,h);
end
function V_outl1 = low(Vin1,h)
    bw = 4;
    C = 10e-6;
    L = 1/(C*(2*60*pi)^2); % using the nature response equation
    R = 2*pi*L*bw;
    tf= 5;
    t = linspace(0,tf,((tf)/h) +1);
    V_outl1 = simV_R(R,L,C,Vin1,t,h);
end
function V_outl2 = low2(Vin2,h)
    bw = 4;
    C = 10e-6;
    L = 1/(C*(2*60*pi)^2); % using the nature response equation
    R = 2*pi*L*bw;
    tf= 5;
    t = linspace(0,tf,((tf)/h) +1);
    V_outl2 = simV_R2(R,L,C,Vin2,t,h);
end

function V_outh = high(Vin1,h)
    bw = 5000;
    C = 10e-6;
    L = 1/(C*(10000*pi)^2); % using the nature response equation
    R = 2*pi*L*bw;
    tf= 5;
    t = linspace(0,tf,((tf)/h) +1);
    V_outh = simV_R(R,L,C,Vin1,t,h);
end
function V_outh2 = high2(Vin2,h)
    bw = 5000;
    C = 10e-6;
    L = 1/(C*(10000*pi)^2); % using the nature response equation
    R = 2*pi*L*bw;
    tf= 5;
    t = linspace(0,tf,((tf)/h) +1);
    V_outh2 = simV_R2(R,L,C,Vin2,t,h);
end

function V_R_Out = simV_R(R,L,C,Vin1,t,h)
% h = 1/(192e3);
V_C = zeros(1, length(t));
I_L = zeros(1, length(t));
x = [V_C; I_L];
u = Vin1;
A = [1 (h/C); (-h/L) (1-(h*R/L))];
B = [0; (h/L)];
for j = 2:length(t)
    x(:,j) = A*x(:,j-1) + B*u(1, j-1);
end
V_R_Out = R*x(1,:);
end
function V_R_Out2 = simV_R2(R,L,C,Vin2,t,h)
% h = 1/(192e3);
V_C = zeros(1, length(t));
I_L = zeros(1, length(t));
x = [V_C; I_L];
u = Vin2;
A = [1 (h/C); (-h/L) (1-(h*R/L))];
B = [0; (h/L)];
for j = 2:length(t)
    x(:,j) = A*x(:,j-1) + B*u(1, j-1);
end
V_R_Out2 = R*x(1,:);
end