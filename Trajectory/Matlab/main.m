
%% Initialise
close all; clear all
global lam y0 ea0 J thrust exhvel g Cd;
% savepath = sprintf('../../Figures/'); % save figures to this location

%% . . Simulation params:
nt = 1000;                  % Number of time steps
h = 0.01;                   % Time-step
T = nt*h;                   % Length of simulation
TT = linspace(0,T,nt+1);    % Time axis

%% . . Model parameters:

Cd = 1;                   % Quadratic drag coefficient
g = 9.8;                  % Magnitude of gravity
mass = 1000;              % Initial rocket mass

% . . Propulsion terms: thrust = d/dt(fuel mass) * exhaust velocity
thrust = @(t) 10000*t.*exp(-(5*t).^2);
exhvel = @(t) 1;

% . . Rocket moment of inertia tensor (body frame)
J = [1 0 0;
    0 1 0;
    0 0 10];

% . . Particle:
lam = 5; % Aspect ratio of rocket for plot

% . . Initial conditions:
ea0  = [6*pi/180, 0.01, 0]';    % Euler angles
v0   = [0, 0, 0]';              % Linear velocity
om0  = [0, 0, 0]';              % Angular velocity
x0   = [0, 0, 0]';              % Position

% . . Write to output:
q0  = mat2quat(ea2mat(ea0)); Q=quattomat(q0);
y0(1:15,1) = [v0',om0',q0',x0',mass, 0];


%% Integrate with ode45
tic;
fprintf('\nDoing time integration...');
options = odeset('RelTol',1e-13);
[~,yout] = ode45(@VF,TT,y0,options);
fprintf(' time = %5.5f seconds\n',toc);
p   = yout  (:,1:3); om   = yout  (:,4:6); q   = yout  (:,7:10); x   = yout  (:,11:13); mass   = yout  (:,14);



%% Graphs of each solution component of all 3 solutions
fig1 = figure(1);
subplot(3,1,1)
plot(TT  , thrust(TT), 'k'  ,TT  , mass, 'm'  )
grid on
title('Thrust and rocket mass')
xlabel('time')
ylabel('')
legend('Thrust (N)','Rocket mass (kg)')
pbaspect([4 1 1])
subplot(3,1,2)
plot(TT  , p  (:,1), 'r' ,TT  , p  (:,2), 'g' ,TT  , p  (:,3), 'b' )
grid on
title('Linear velocity')
%     axis([0 T min(min(v)) max(max(v))])
xlabel('time')
ylabel('v_i (m/s)')
legend('v_x','v_y','v_z')
pbaspect([4 1 1])
subplot(3,1,3)
plot(TT  , x  (:,1), 'r' ,TT  , x  (:,2), 'g' ,TT  , x  (:,3), 'b' )
grid on
title('Position')
%     axis([0 T min(min(x)) max(max(x))])
xlabel('time')
ylabel('x_i (m)')
legend('x','y','z')
pbaspect([4 1 1])

drawnow

printpdf(fig1,'rocket_results')

fig2 = figure(3);
particle_dynamics(x,q,TT,10)
%%
printpdf(fig2,'rocket_path3d')




