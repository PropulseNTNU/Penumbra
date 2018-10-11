function q =  mat2quat(Q)
% Conversion from rotation matrix to quaternions
q = zeros(4,1);
t = trace(Q);
if t > 0
    r = sqrt(1+t);
    s = 0.5/r;
    q(1) = 0.5*r;
    q(2) = (Q(3,2)-Q(2,3))*s;
    q(3) = (Q(1,3)-Q(3,1))*s;
    q(4) = (Q(2,1)-Q(1,2))*s;
else
    [~,i] = max(diag(Q));
    ind = [1:(i-1),(i+1):3];
    if i == 2
        ind = ind([2,1]);
    end
    r = sqrt(1+2*Q(i,i)-t);
    s = 0.5/r;
    q(1) = (Q(ind(2),ind(1))-Q(ind(1),ind(2)))*s;
    q(i+1) = 0.5*r;
    q(ind(1)+1) = (Q(ind(1),i)+Q(i,ind(1)))*s;
    q(ind(2)+1) = (Q(ind(2),i)+Q(i,ind(2)))*s;
end