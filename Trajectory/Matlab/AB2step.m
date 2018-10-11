function Y2 = AB2step(Y0,Y1,F,t,h)
    
    % Adam-Bashforth 2 step mehtod
    Y2=Y1+h*(3/2*(F(t+h,Y1))-1/2*(F(t,Y0)));

end


