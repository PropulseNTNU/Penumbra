

function [ ydot ] = VF( ~, y )

% . . Get global variables
global J thrust exhvel g Cd

% . . Get variables from input
d = size(y ,1);
v = y(1 :3 ,1);
om = y(4 :6 ,1);
q = y(7 :10,1);
x = y(11:13,1);
mass = y(14,1);
t = y(15   ,1);


% . . Normalize quaternion and convert to rotation matrix
q = real(q)/norm(q);
Q = quattomat(q);

% . . Stokes drag and gravity (stationary frame)
D = -Cd*v*norm(v)/mass ;
G = -[0 0 g]';
T = Q*[0 0 thrust(t)]'; % (stationary frame)

% . . System ODEs for integration
vdot = D + G + T;             % (stationary frame)
omdot = cross(J*om,om);     % (body frame)
qdot = 0.5*R([0;om])*q;     % (body frame)
xdot = v;                   % (stationary frame)
massdot = -thrust(t)/exhvel(t);
tdot = 1;

% . . Write to output
ydot = zeros(d,1);
ydot( 1: 3,1) = vdot;
ydot( 4: 6,1) = omdot;
ydot( 7:10,1) = qdot;
ydot(11:13,1) = xdot;
ydot(14   ,1) = massdot;
ydot(15   ,1) = tdot;

end

