function tests = test_Abstraction
tests = functiontests(localfunctions);
end

function testAbsContruct(testCase)
    A = [3, -8; -20, 10];
    B = [1;1];
    dyn1 = @(x,u) A*x + B*u;

    K1 = @(x,u) norm(A,'inf');
    
    % state space [-1,1]x[-1,1], the initial partition has four grids
    X = [-1 1; -1 1];
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

    U1 = [1];

    % spec []A ^ <>[]B ^ (AND_i []<> (Ri))
    spec1.A = [1,2,3,4];
    spec1.B = [];
    spec1.R = {};
    
    abs_test = Abstraction(X, U1, G1, spec1, label1, dyn1, K1);
    testCase.verifyEqual(abs_test.n, 4);
    testCase.verifyEqual(abs_test.m, 1);
    testCase.verifyEqual(abs_test.dim, 2);
    testCase.verifyEqual(abs_test.A{1},ones(4));
    testCase.verifyEqual(abs_test.dyn([1;0],1), A*[1;0]+B*1);
    
    P1 = [-1 -0.2
          -1 0];
    P2 = [-1 -0.2
          0 1];
    P3 = [-0.2 0.2
           -1 0];
    P4 = [-0.2 0.2
           0 1];
    P5 = [0.2 1
           -1 0];
    P6 = [0.2 1
           0 1];
    G2 = {P1, P2 ,P3, P4, P5, P6};
    label2 = ["a","b","c","d","e","f"];
    dyn2 = @(x,u) [1;1];
    abs_test = Abstraction(X, U1, G2, spec1, label2, dyn2, @(x,u) 20, 0.1);
    A = ones(6,6);
    A(1:2,5:6) = 0;
    A(5:6,1:2) = 0;
    testCase.verifyEqual(abs_test.A{1},A);
end

function testVerifyTransition(testCase)
    X = [-1 1; -1 1];
    U = [1];
    spec = [];
    P1 = [-1 -0.2
          -1 0];
    P2 = [-1 -0.2
          0 1];
    P3 = [-0.2 0.2
           -1 0];
    P4 = [-0.2 0.2
           0 1];
    P5 = [0.2 1
           -1 0];
    P6 = [0.2 1
           0 1];
    G2 = {P1, P2 ,P3, P4, P5, P6};
    label = ["a","b","c","d","e","f"];
    dyn = @(x,u) [1;1];
    abs_test = Abstraction(X, U, G2, spec, label, dyn, @(x,u) 20, 0.1);
    surface = [-1 1
               -1 1];
           
    % the cover in the center of the grid
    cover = [-0.5, 0.5
             -0.5, 0.5];
    residual = abs_test.surface_diff(surface,cover);
    testCase.verifyEqual(length(residual),8);
    
    cell1 = [-1.0000   -0.5000
            -1.0000   -0.5000];
    cell4 = [-1.0000   -0.5000
             -0.5000    0.5000];
    cell8 = [ 0.5000    1.0000
                0.5000    1.0000];
    testCase.verifyEqual(residual{1},cell1);
    testCase.verifyEqual(residual{4},cell4);
    testCase.verifyEqual(residual{8},cell8);
    
    % the cover completely overlaps the grid
    cover = [-2, 1
             -2, 1];
    residual = abs_test.surface_diff(surface,cover);
    testCase.verifyEqual(length(residual),0);
    
    % the cover overlaps one corner of the grid
    cover = [-1, 0.5
             -1, 0.5];
    residual = abs_test.surface_diff(surface,cover);
    testCase.verifyEqual(length(residual),3);
    testCase.verifyEqual(residual{1},[0.5 1; -1 0.5]);
    testCase.verifyEqual(residual{2},[-1 0.5; 0.5 1]);
    testCase.verifyEqual(residual{3},[0.5 1; 0.5 1]);
    
    % when surface is in lower dimensional
    surface = [0 0
               -1 1];
    cover = [-1, 0.5
             -1, 0.5];
         
    residual = abs_test.surface_diff(surface,cover);
    testCase.verifyEqual(residual{1},[0 0 ;0.5 1]);
    
    x = [0.5;0];
    normal = [1;0];
    u = 1;
    [flow, cover] = abs_test.compute_cover(x,u,normal);
    testCase.verifyTrue(length(flow)==2);
    testCase.verifyTrue(all(cover(:,1)<=cover(:,2)));
    abs_test.verifyTransition(0.001);
    A_ref = [    1     1     1     1     0     0
                 0     1     0     1     0     0
                 0     0     1     1     1     1
                 0     0     0     1     0     1
                 0     0     0     0     1     1
                 0     0     0     0     0     1];
    testCase.verifyEqual(abs_test.A{1},A_ref);
    
    % special case: flow are parallel to some common faces
    dyn = @(x,u) [1;0];
    abs_test = Abstraction(X, U, G2, spec, label, dyn, @(x,u) 20, 0.1);
    abs_test.verifyTransition(0.001);
    A_ref = [1     1     1     0     0     0
             1     1     0     1     0     0
             0     0     1     1     1     0
             0     0     1     1     0     1
             0     0     0     0     1     1
             0     0     0     0     1     1];
    testCase.verifyEqual(abs_test.A{1},A_ref);
end

function testPlot(testCase)
    A = [3, -8; -20, 10];
    B = [1;1];
    dyn1 = @(x,u) A*x + B*u;

    K1 = @(x,u) norm(A,'inf');
    
    % state space [-1,1]x[-1,1], the initial partition has four grids
    X = [-1 1; -1 1];
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

    U1 = [-1, 0 ,1];

    % spec []A ^ <>[]B ^ (AND_i []<> (Ri))
    spec1.A = [1,2,3,4];
    spec1.B = [];
    spec1.R = {};
    
    abs_test = Abstraction(X, U1, G1, spec1, label1, dyn1, K1);
    fig = figure;
    abs_test.phase_portrait(fig);
end