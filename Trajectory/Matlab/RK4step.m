function x = RK4step(x,F,t,h)
% Implementation of the classical fourth order Runge-Kutta method for
% an autonomous system of ODEs in vector form

% Input: 
% x - vector of initial conditions
% F - function handle (RHS of ODE to integrate)
% t - time
% h - time step
%
% Output:
% x - solution at time t+h

K1 = F(t,x);
K2 = F(t+0.5*h,x+0.5*h*K1);
K3 = F(t+0.5*h,x+0.5*h*K2);
K4 = F(t,x+h*K3);
x = x + (h/6)*(K1+2*K2+2*K3+K4);

end
