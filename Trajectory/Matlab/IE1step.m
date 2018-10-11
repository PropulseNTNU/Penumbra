function xp = IE1step(x,F,t,h)

% Implicit Euler method

K1 = F(t,x);
K2 = F(t+h,x+h*K1);
x0 = x + h*K2;

objFun = @(xp_) x+h*F(t+h,xp_)-xp_;

options = optimoptions('fsolve','display','none','FunctionTolerance',5e-15,'OptimalityTolerance',5e-15,'StepTolerance',1e-15);

xp = fsolve(objFun,x0,options);

end
