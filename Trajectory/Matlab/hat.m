function [ Hw ] = hat( w )
% hat map: maps a vector in R^3 into an emelent of so(3) (Lie Algebra)

Hw=[   0    -w(3,1)   w(2,1)
     w(3,1)    0     -w(1,1)
    -w(2,1)  w(1,1)     0   ];

end

