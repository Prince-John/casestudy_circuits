function V_R_Out = simV_R(R,L,C,V_in,t)
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