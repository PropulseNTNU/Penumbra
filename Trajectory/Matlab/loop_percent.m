function loop_percent(it,nt,verbfreq)

% INPUT:
% it        - loop index
% nt        - loop index range
% verbfreq  - how frequently to update

if mod(it,nt/verbfreq)==0
    
    if verbfreq*it/nt>1 % only erase after the first statement
        fprintf(repmat('\b',1,11));
    end
    
    fprintf(': %3.0f%% done',100.*real(it)/real(nt)) % state what percentage of hte loop is complete
    
end

if it==nt % erase the last statement
    fprintf(repmat('\b',1,11));
end

end
