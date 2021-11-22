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

% if the natural response is given, I can reverse engineer to find the
% right L and C
function Vout = myResonatorCircuit(Vin,h)
    R = 2000;
    L = 5;
    C = 0.1e-6;
    % generate a time series, assuming the start time is 0
    tf= 5; %End Time in seconds
    t = linspace(0,tf,((tf)/h) +1);
    Vout = simV_r(R,L,C,Vin,t, h);
    sound_frequency_response(R,L,C,t,h)
end

function sound_frequency_response(R,L,C,t,h)
    frequencies = decade(1,3,0.5);
    frequency_response = [frequencies' zeros(length(frequencies),1)];
    index = 1;
    for f = frequencies
        V_in = sin(2*pi*f*t);
        V_out = simV_R(R,L,C,V_in,t,h);
        ratio = norm(V_out)/norm(V_in);
        frequency_response(index, 2) = ratio;
        index = index +1;
    end
    figure();
    semilogx(frequency_response(:,1), frequency_response(:,2));
    title('|V_{out}|/|V_{in}| vs frequency of V_{in}');
    ylabel('|V_{out}|/|V_{in}|');
    xlabel('Frequency(Hz)');

end

function v_r_out = simV_r(R,L,C,V_in,t,h)
V_C = zeros(1, length(t));
I_L = zeros(1, length(t));
x = [V_C; I_L];
u = V_in'; 

A = [1 (h/C); (-h/L) (1-(h*R/L))];
B = [0; (h/L)];

for j = 2:length(t)
    x(:,j) = A*x(:,j-1) + B*u(1, j-1);
end

v_r_out = R*x(2,:);

end

function V_R_Out = simV_R(R,L,C,V_in,t,h)
% h = 1/(192e3);
V_C = zeros(1, length(t));
I_L = zeros(1, length(t));
x = [V_C; I_L];
u = V_in; 

A = [1 (h/C); (-h/L) (1-(h*R/L))];
B = [0; (h/L)];

for j = 2:length(t)
    x(:,j) = A*x(:,j-1) + B*u(1, j-1);
end

V_R_Out = R*x(2,:);

end