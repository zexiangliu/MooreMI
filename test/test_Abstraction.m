function tests = test_Abstraction
tests = functiontests(localfunctions);
end

function testAbsContruct(testCase)
    A = [3, -8; -20, 10];
    B = [1;1];
    dyn1 = @(x,u) A*x + B*u;

    K1 = @(x,u) norm(A,2);
    
    % state space [-1,1]x[-1,1], the initial partition has four grids
    P1 = [-1 0;  % [lb, ub] of dim 1
          -1 0]; % [lb, ub] of dim 2
    P2 = [-1 0;
          0 1];
    P3 = [0 1;
          0 1];
    P4 = [0 1;
          -1 0];
    G1 = {P1, P2, P3, P4};
    label1 = ["a","b","c","d"];

    U1 = ["a"];

    % spec []A ^ <>[]B ^ (AND_i []<> (Ri))
    spec1.A = [1,2,3,4];
    spec1.B = [];
    spec1.R = {};
    
    abs = Abstraction(U1, G1, spec1, label1, dyn1, K1);
    testCase.verifyEqual(abs.n, 4);
    testCase.verifyEqual(abs.m, 1);
    testCase.verifyEqual(abs.dim, 2);
    testCase.verifyEqual(abs.A,ones(4));
    testCase.verifyEqual(abs.dyn([1;0],1), A*[1;0]+B*1);
    
    P5 = [-1 -0.2
          -1 1];
    P6 = [-0.2 0.2
           -1 1];
    P7 = [0.2 1
           -1 1];
    G2 = {P5, P6, P7};
    label2 = ["a","b","c"];
    abs = Abstraction(U1, G2, spec1, label2, dyn1, K1);
    testCase.verifyEqual(abs.A,[1 1 0; 1 1 1; 0 1 1]);
end