%% This file is used to tune and validate the mySensor Circuit function.


% Desired freq = 84 Hz

% From previous tests we know that the RLC circuit has a band pass effect
% around its resonent freq. We can use this fact to allow only the rotor
% noise around a narrow band of 84 Hz pass through. 

load('MarsHelicopter_noisy.mat');
Vin = Vsound';
h = 1/Fs;
t = linspace(0,h*length(Vin),length(Vin));

%% RLC Circuit will be simulated using simV function

% resonant freq in an rlc circuit is given by $\frac{1}{2\pi\sqrt{LC}}$
% So solving that for 84 Hz, 84*2*pi = 1/sqrt(LC);  LC = 1/(168pi)^2; 
% Or L = (C(168pi)^2)^-1
% Setting C = 10uF 
% bandwidth of a series RLC circuit is generally given by R/2piL, since we
% want only a narrow band of frequencies to pass through lets choose
% bandwidth starting point as 20 Hz
bw = 4;
C = 10e-6;
L = 1/(C*(168*pi)^2);
R = 2*pi*L*bw;
Vout = simV_R(R,L,C,Vin,t,h);
Vout2 = simV_R(R,L,C,Vout,t,h);
Vout3 = simV_R(2*pi*L*2,L,C,Vout,t,h);

%%
figure();
subplot(2,1,1);
plot(t, Vout);
title('Vout1')
subplot(2,1,2);
plot(t, Vout2);
title('Vout2')
%% FFT analysis
% https://stackoverflow.com/questions/10758315/understanding-matlab-fft-example
NFFT = 2^nextpow2(length(t)); % Next power of 2 from length of y
Y1 = fft(Vin,NFFT)/length(t);
Y2 = fft(Vout,NFFT)/length(t);
Y3 = fft(Vout2,NFFT)/length(t);
Y4 = fft(Vout3,NFFT)/length(t);
f = Fs/2*linspace(0,1,NFFT/2+1);

figure()
% Plot single-sided amplitude spectrum.
subplot(4,1,1)
plot(f,2*abs(Y1(1:NFFT/2+1))) 
xlim([0 200]);
title('FFT of Vin')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
subplot(4,1,2)
plot(f,2*abs(Y2(1:NFFT/2+1))) 
xlim([0 200]);
title('FFT of Vout1')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
subplot(4,1,3)
plot(f,2*abs(Y3(1:NFFT/2+1))) 
xlim([0 200]);
title('FFT of Vout2')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
subplot(4,1,4)
plot(f,2*abs(Y4(1:NFFT/2+1))) 
xlim([0 200]);
title('FFT of Vout3')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
