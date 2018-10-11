function Rq = R( q )
%right multiplication by a quaternion
q0=q(1,1);
qq=q(2:4,1);
Rq=[q0 -qq';
    qq q0*eye(3)-hat(qq)];

end

