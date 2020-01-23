% Skoltech, Fall 2019. Essential Engineering Toolbox
% Matlab tutorial practice excersizes by Prof. Aslan Kasimov.
% 
% Content: 
% 1: Basic elements: vectors, matrices, functions, input/output
% 2: Linear systems: \, pinv, least squares, eigenvalues, SVD, ...
% 3: Graphics
% 4: Nonlinear equations
% 5: Fitting
% 6: ODE: initial/boundary value problems
% 7: Extra, time permitting 

%% 1: Basic elements: vectors, matrices, functions, input/output
clear all   % remove all variables from the workspace
clc         % clear command window

%1. create vectors x = [ -2 -1 0 1 2 3 4] and y = [ -5 -4 -3 -2 -1 0 1]
x = [-2 -1 0 1 2 3 4]'
y = [-5 -4 -3 -2 -1 0 1]'

%2. compute the scalar product of x and y
sp = x'*y

%3. find the outer product of x and y
op = x.*y'

%4. build a teoplitz matrix with column x and row y 
tp = toeplitz(x,y)

%5. make a matrix A with squared x as the first column and 
%   elementwise product of x and y as the second column.
A = [x.^2 x.*y]

%6. make a vector z by concatenating x and y
z = [x y]

%7. make a matrix B by adding to A a column, which is an average of columns of A.   
B = [A, mean(A,2)]

%8. build a block matrix C with zeros in the top left block, ones in the
%   bottom right block, and part of B with the exception of rows 1, 3, and
%   the one before the last in the remaining blocks.
B1 = [B(2,:); B(4,:); B(5,:); B(6,:); B(7,:)]
C = [zeros(5,3), B1; B1, ones(5,3)]

%% 2: Linear systems: \, pinv, least squares, eigenvalues, SVD, ...

%9. Find e-values of A = toeplitz(-5:5) and the e-vector that corresponds
%   to e-values lambda1 = -1 and lambda2 = 0

A = toeplitz(-5:5)
[V,D] = eig(A)
V(:,6)
V(:,end)

%10. Change one number in A to make it full rank and then solve the system 
%   Ax = b with b a uniformly distributed random vector 
%   b = rand(size(A,1),1)

rA = rank(A)
A(1,1) = A(1,1) + 1
rA_ = rank(A)
b = rand(size(A,1),1)
x = linsolve(A, b)

%% 3: Graphics

%% 11. Make a 1D plot of two functions with the same axes
%    add labels on x, y, a title with latex, and a grid
%    add different marker symbols to the graphs
%    add legends at the top left corner
% x = [0 10]
% y1 = 2+tanh(x/10).*sin(x.^2)
% y2 = sin(2*x).*cos(50*x);

x = linspace(0, 10, 1000)
y1 = 2+tanh(x/10).*sin(x.^2)
y2 = sin(2*x).*cos(50*x)

plot(x, y1, '--', x, y2)
legend('$$y1 = 2+tanh(x/10)*sin(x^2)$$', '$$y2 = sin(2*x)*cos(50*x)$$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'NorthWest')
title('2d plots', 'FontSize', 16)
xlabel('x', 'Interpreter','latex')
ylabel('y', 'Interpreter','latex')
grid on

%% 12. Make a 1D plot of two functions with two axes
%    add labels on x, y, a title with latex, and a grid
%    add the functions in the text box at the bottom right
% 
% y = (sin(5x -1))^2*ln(x) 
% z = exp(-1/(2-x)^4)

x = linspace(0, 10, 1000)
y = (sin(5.*x -1)).^2.*log(x)
z = exp(-1./(2-x).^4)
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1])

hold on
subplot(1, 2, 1)
plot(x, y)
title('y(x)', 'FontSize', 16)
xlabel('x', 'Interpreter', 'latex')
ylabel('y', 'Interpreter', 'latex')
legend('$$y = (sin(5x-1))^2*log(x)$$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'SouthEast')
grid on

subplot(1, 2, 2)
plot(x, z)
title('z(x)', 'FontSize', 16)
xlabel('x', 'Interpreter', 'latex')
ylabel('z', 'Interpreter', 'latex')
legend('$$z = exp(-1/(2-x)^4)$$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'SouthEast')
grid on

%% 13. Make two horizontal subplots with different types of figures in each: 
%   1. log-log plot of one line plot of a 1D function, 
%   f(x) = x^(-3)*tanh(x)
%   2. surface plot with a contour underneath of a 2D function
%   g(x,y) = (sin(x+y))^2*(cos(x-y))

x = linspace(0, 10, 1000)
[n, m] = meshgrid(-2:.2:2, -2:.2:2)
f = x.^(-3).*tanh(x)
g = (sin(n+m)).^2.*cos(n-m)
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1])

hold on
subplot(1, 2, 1)
loglog(x, f)
title('f(x)', 'FontSize', 16)
xlabel('x', 'Interpreter', 'latex')
ylabel('f', 'Interpreter', 'latex')
grid on

subplot(1, 2, 2)
surfl(n, m, g)
title('z(x)', 'FontSize', 16)
xlabel('x', 'Interpreter', 'latex')
ylabel('y', 'Interpreter', 'latex')
zlabel('g', 'Interpreter', 'latex')
grid on

%% 4: Nonlinear equations

%14. Find the roots of a polynomial
%skip

%15. Solve an equation by Newton's method
%skip

%% 5: Fitting

%% 16. 
% 1. Take the following data with noise
% N = 30;
% x = linspace(-10,10,N);
% y = tanh(x) + 0.5*rand(N,1)'
% 2. Interpolate these data with a polynomial of degrees 5, 10, and 20 .
% 3. Plot both the original function with noise and the interpolated
% polynomial evaluated on a grid with 10*N points on the same figure. 

N = 30
x = linspace(-10, 10, N)
y = tanh(x) + 0.5*rand(N, 1)'
p1 = polyfit(x, y, 5)
p2 = polyfit(x, y, 10)
p3 = polyfit(x, y, 20)

xx = linspace(-10, 10, N*10)
y1 = polyval(p1, xx)
y2 = polyval(p2, xx)
y3 = polyval(p3, xx)

plot(x, y)
hold on
plot(xx, y1)
plot(xx, y2)
plot(xx, y3)
grid on

xlabel('x', 'Interpreter', 'latex')
ylabel('y', 'Interpreter', 'latex')
legend('$$y$$', '$$y1 (5 degree)$$', '$$y1 (10 degree)$$', '$$y1 (20 degree)$$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'SouthEast')
title('y(x)', 'FontSize', 16)
%% 6: ODE: initial/boundary value problems

%% 17. Solve the IVP for a pendulum with damping and forcing
% Solve the pendulum problem with added nonlinear damping 
% of the type nu*(1-y)*\dot(y)
% plot the phase plane with vector fields as well as the solution
% y as a function of time for three different initial conditions,
% one close to the center, and two somewhat away
% you can try different ODE solvers if ode45 fails

% See the code in the pendulum.m file

