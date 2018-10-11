function Q = ea2mat(theta)
% Converts Euler angles to rotational matrices.

c = cos(theta);
s = sin(theta);
Q1 = [1,0,0;0,c(1),-s(1);0,s(1),c(1)];
Q2 = [c(2),0,s(2);0,1,0;-s(2),0,c(2)];
Q3 = [c(3),-s(3),0;s(3),c(3),0;0,0,1];
Q = Q3*Q2*Q1;
end