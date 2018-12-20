function dx = dyn(x,u)

    C = 1e-6;
    L = 6.8*1e-6;
    R = 100;
    V_in = 12;
    
    A1 = [0 0; 0 -1/(R*C)];
    A2 = [0 -1/L; 1/C -1/(R*C)];
    K = [V_in/L; 0];
    
    if u == 0 % model 1
        dx = A1*x + K;
    elseif u == 1 % model 2
        dx = A2*x + K;
    else
        error("no such input");
    end
    
end