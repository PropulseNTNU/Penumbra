function Q=quattomat(y)
% quaternion to rotation matrix. 
S=hat(y(2:4,1));
Q=eye(3)+2*y(1,1)*S+2*S*S;
