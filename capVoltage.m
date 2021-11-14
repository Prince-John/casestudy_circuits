function Vc = capVoltage(h,R,C, V)

b = h/(R*C);
a = 1 - b;
A = [a b];

for i = 2:length(V)
    
    V(1,i) = A*V(:,i-1);
end

Vc = V;

end