function x = MP2step(x,F,t,h)
% Implicit mid point rule

    K1 = F(t,x);
    K2 = F(t+0.5*h,x+0.5*h*K1);
    x = x + h*K2; 
    
end
