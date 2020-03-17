% Zhuofei Xie
% Math 1070
% Professor: Michael Neilan 
% 11/10/2019

function LeastSquareApproximation 
clear;
clc;

a = 1; 
b = -1; 
n = 100; 
h = (b-a)/(n-1); 
xi = linspace(a,b,n); 
fh = @(x)cos(10.*x.^4 - 14.*x).*exp(1)^abs(sin(x)); 
e = ones(n,1); 
A = h/6 * spdiags([e,4*e,e], -1:1, n ,n); 
A(1,1) = h/6 * 2; 
A(n,n) =  h/6 * 2; 

F = ones(n,1); 
for i = 1:n
    F(i) = h * fh(xi(i)); 
end 

F(1) = F(1) * 0.5; 
F(n) = F(n) * 0.5; 

c = A\F;

plot(xi,c)

error = zeros(n,1); 

for i = 1:n
    error(i) = fh(xi(i)) - c(i); 
end 

max(abs(error))




