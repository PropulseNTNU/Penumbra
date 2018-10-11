function [ ttheta] = quat2ea( q )
%map quaternions to Euler angles
% formulae from Wikipedia http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
q0=q(1,1);
q1=q(2,1);
q2=q(3,1);
q3=q(4,1);
ttheta=[
    atan2(2*(q0*q1+q2*q3),1-2*(q1^2+q2^2))
    asin(2*(q0*q2-q3*q1))
    atan2(2*(q0*q3+q1*q2),1-2*(q2^2+q3^2))
    ];

end

