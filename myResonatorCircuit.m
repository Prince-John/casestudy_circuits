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
    % I pick the frequency to be 440 Hz. However, I have made the function
    % to run any frequencies. Assuming we only run this pitch for 5 secs
    f = 440;
    Vout = resonating_circuit(f,Vin,h);
end

function V_out = resonating_circuit(f,V_in,h)
    bw = 4;
    C = 10e-6;
    L = 1/(C*(2*f*pi)^2); % using the nature response equation
    R = 2*pi*L*bw;
    % generate a time series, assuming the start time is 0
    tf= 5; %End Time in seconds. Simulate in 5 secs
    t = linspace(0,tf,((tf)/h) +1);
    sound_frequency_response(R,L,C,t,h,f)
    V_out = simV_r(R,L,C,V_in,t, h);
end

function sound_frequency_response(R,L,C,t,h,f)
    frequencies = decade(1,ceil(log10(f)),0.5);
    frequency_response = [frequencies' zeros(length(frequencies),1)];
    index = 1;
    for f = frequencies
        V_in = sin(2*pi*f*t);
        V_out = simV_r(R,L,C,V_in,t,h);
        ratio = norm(V_out)/norm(V_in);
        frequency_response(index, 2) = ratio;
        index = index +1;
    end
    figure();
    hold on;
    y = max(frequency_response(:,2))*ones(length(frequency_response(:,2)),1);
    plot(frequency_response(:,1), frequency_response(:,2));
    plot(frequency_response(:,1), y);
    title('|V_{out}|/|V_{in}| vs frequency of V_{in}');
    ylabel('|V_{out}|/|V_{in}|');
    xlabel('Frequency(Hz)');
    legend();

end

function v_r_out = simV_r(R,L,C,V_in,t,h)
    V_C = zeros(1, length(t));
    I_L = zeros(1, length(t));
    x = [V_C; I_L];
    if size(V_in,1) > 1
        u = V_in'; 
    else
        u = V_in;
    end

    A = [1 (h/C); (-h/L) (1-(h*R/L))];
    B = [0; (h/L)];

    for j = 2:length(t)
        x(:,j) = A*x(:,j-1) + B*u(1, j-1);
    end

    v_r_out = R*x(2,:);

end