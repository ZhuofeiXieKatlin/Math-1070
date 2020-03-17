%Zhuofei Xie
%Math 1070
%Professor: Michael Neilan 
%10/14/2019

function calculateSpline 

clear;
clc; 

fh = @(x)1./(1 + 25.*x.^2); 

n = 5; % declear number of nodes 
m = 100; % declear number of grid points for plotting

% declear x nodes with distince entries 
h = 2/(n-1); 
x = linspace(-1,(-1+(n-1)*h),n); % let x be the equal spaced nodes with h = 2/(n-1) 
                                           % Therefore xi+1 - xi = h
% declear y nodes 
y = fh(x); 

% declear a vector X with distinct entries 
t = 2 /(m-1); 
X = linspace(-1, (-1+(m-1)*t), m); % let X also be the equal spaced nodes with t = 2/(m-1) 
                                             % Therefore, Xi+1 - Xi = t

% Call the function mySplineFunction to obtain the value of the spline and
% the points stored in the vector X
Y = mySplineFunction(x,y,X); 

fplot(fh, [-1,1]); % Plot the original function plot with the range [-1,1]
hold on; 
plot(X,Y); % Plot our cublic spline; as realize, if I increase the n from 5 to 10 or 15, 
           %the cublic spline looks much similar to the original function 

end 

function Y = mySplineFunction(x,y,X)

% set up tridiagnal matrix A. A is an (n-2) * (n-2) matrix, 
% where n = length(x) = length(y)

% % % This is for function x 
h = 2/(length(x)-1); 
t = 2/(length(X)-1); 
n = length(x); % n = length(x) = length(y)
m = length(X); 
e = ones(n-2,1); % declear a vector with all index = 1 and length = n-2
A = h/6 * spdiags([e,4*e,e], -1:1,n-2,n-2); % declear a matrix with main diaganol = 4 and the two diaganol around the main diaganol = 1
                                            % the size of the matrix is (n-2) * (n-2)
b = zeros(n-2,1);  % initial a vecor with all index = 1 and the length of the vecot = n-2
for i = 1:n-2
    element = (y(i+2) - 2*y(i+1) + y(i))/h; % bi = (f(xi+2)-2f(xi+1) + f(xi))/h 
    b(i) = element; 
end 


% Set up right-handed vector b (which depends on the values of both x and y)

M = [0;A\b;0]; % AM = b, then M = A\b 

% with the coefficients calculate the values of the cubic spline. Store
% these values in the variable Y  

Y = ones(m,1); % delcear a new vector called Y with 100 elementes all stored 1 at first 

for i = 1:m 
    for j = 2:n 
        if X(i) >= x(j-1) && X(i) <= x(j) % find out whether the the point X(i) inside the range x(j-1) and x(j)
            % if X(i) inside the range, then call the cublic spline funciton 
            Y(i) = (((x(j)-X(i)).^3*M(j-1) + (X(i)-x(j-1)).^3*M(j))/(6*h)) + ((x(j)-X(i))*y(j-1) + (X(i)-x(j-1))*y(j))/h - (1/6)*(h)*((x(j)-X(i))*M(j-1)+(X(i)-x(j-1))*M(j)); 
            % Call the cublic spline function 
        end 
    end 
end
    
end
