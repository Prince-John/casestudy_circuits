%% Case study 3: Circuits as Resonators, Sensors, and Filters
% *ESE 105* 
%
% *Name: Quan Khuc, John Prince, Mordecai
%
% function myResonatorCircuit(Vin,h) receives a time-series voltage sequence
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

function Vout = myResonatorCircuit(Vin,h)
% Vout = Vin;
    % for now, I will try random R, L, C values. I guess tuning means
    % finding the good value for each of these values
    R = 500;
    L = 100e-3;
    C = 0.1e-6;
    % generate a time series, assuming the start time is 0
    tf= 0.015; %End Time in seconds
    t = linspace(0,tf,((tf)/h) +1);
    Vout = simV_r(R,L,C,Vin,t, h);
end

function v_r_out = simV_r(R,L,C,V_in,t,h)
V_C = zeros(1, length(t));
I_L = zeros(1, length(t));
x = [V_C; I_L];
u = V_in; 

A = [1 (h/C); (-h/L) (1-(h*R/L))];
B = [0; (h/L)];

for j = 2:length(t)
    % bug at here. How can I multiply a matrix u which is 96000x1 with
    % matrix B which is 1x2881. What would this tell me?
    x(:,j) = A*x(:,j-1) + B*u(1, j-1);
end

v_r_out = R*x(2,:);

end