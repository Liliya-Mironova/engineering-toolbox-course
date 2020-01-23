%% 17. Solve the IVP for a pendulum with damping and forcing
% Solve the pendulum problem with added nonlinear damping 
% of the type nu*(1-y)*\dot(y)
% plot the phase plane with vector fields as well as the solution
% y as a function of time for three different initial conditions,
% one close to the center, and two somewhat away
% you can try different ODE solvers if ode45 fails

function pendulum
%Based on Higham's pendulum solver. 
% run this in the command line:
% pendulum, pause(1), end

tspan = [0 30];                       % Solve for 0 <= t <= 10.
yazero = [1; 1];                      % Initial conditions.
ybzero = [-5; 2];                      % Initial conditions.
yczero = [5; -2];                      % Initial conditions.

[ta,ya] = ode45(@pend,tspan,yazero);
[tb,yb] = ode45(@pend,tspan,ybzero);
[tc,yc] = ode45(@pend,tspan,yczero);

[y1, y2] = meshgrid(-5:0.5:5, -3:0.5:3);
Dy1Dt = y2;
nu = 0.05
Dy2Dt = -nu.*(1 - y1).*y2 - sin(y1);
quiver(y1,y2,Dy1Dt,Dy2Dt);
hold on

plot(ya(:,1),ya(:,2),yb(:,1),yb(:,2),yc(:,1),yc(:,2))
xlabel('y_1(t)', 'Interpreter', 'latex')
ylabel('y_2(t)', 'Interpreter', 'latex')
title('Phase plane')
axis equal
axis([-5 5 -3 3])
hold off

figure(2)
subplot(1,3,1)
plot(ta, ya(:,1))
xlabel('xa', 'Interpreter', 'latex')
ylabel('ya', 'Interpreter', 'latex')
grid on

subplot(1,3,2)
plot(tb, yb(:,1))
title('Solutions')
xlabel('xb', 'Interpreter', 'latex')
ylabel('yb', 'Interpreter', 'latex')
grid on

subplot(1,3,3)
plot(tc, yc(:,1))
xlabel('xc', 'Interpreter', 'latex')
ylabel('yc', 'Interpreter', 'latex')
grid on

function yprime = pend(t,y)
%Simple pendulum
    nu = 0.05
    yprime = [y(2); -nu.*(1 - y(1)).*y(2) - sin(y(1))];
    
end

end