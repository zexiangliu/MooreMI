function tests = test_DFA
tests = functiontests(localfunctions);
end

function testCreateDFA(testCase)
    n = 3;
    m = 2;

    A1 = [0 1 0; 0 0 1;0 0 0];
    A2 = [0 0 0; 1 0 0;0 1 0];

    A = {A1,A2};
    Q_name = ["1","2","3"];
    Q_final = ["1","2"];
    Q_label = ["01","02","03"];
    U_name = ["a","b"];
    % incorrect Q0 assignment
    Q0 = "a";
    
    testCase.verifyError(@()DFA(3,2,A,[],[],[],Q_name,U_name,Q0,...
        Q_final,Q_label),'MATLAB:matrix:singleSubscriptNumelMismatch');
    
    Q0 = "1";
    % incorrect Q_final assignment
    Q_final = ["1","a"];
    testCase.verifyError(@()DFA(3,2,A,[],[],[],Q_name,U_name,Q0,...
        Q_final,Q_label),'MATLAB:matrix:singleSubscriptNumelMismatch');
    
    % test basic functions
    Q_final = ["1","2"];
    G = DFA(3,2,A,[],[],[],Q_name,U_name,Q0,Q_final,Q_label);
    testCase.verifyEqual(G.get_x_name(3),"3");
    testCase.verifyEqual(G.get_x_idx("2"),2);
    testCase.verifyEqual(G.Q0,1);
    testCase.verifyEqual(G.get_u_name(1),"a");
    testCase.verifyEqual(G.get_u_idx("b"),2);
    
    % test pre
    [pre_x,pre_u] = G.pre(2);
    testCase.verifyEqual(pre_x,[1;3]);
    testCase.verifyEqual(pre_u,[1;2]);
    pre_x = G.pre_xu(2,1);
    testCase.verifyEqual(pre_x,1);
    pre_x = G.pre_xu(1,1);
    testCase.verifyTrue(isempty(pre_x));
    
    % test post
    [post_x,post_u] = G.post(2);
    testCase.verifyEqual(post_x,[3;1]);
    testCase.verifyEqual(post_u,[1;2]);
    post_x = G.post_xu(2,2);
    testCase.verifyEqual(post_x,1);
    post_x = G.post_xu(1,2);
    testCase.verifyTrue(isempty(post_x));

    % test run 
    testCase.verifyTrue(G.isaccepting("1"));
    testCase.verifyFalse(G.isaccepting("3"));
    testCase.verifyTrue(G.isaccepting(2));
    % infeasible string
    [status, q] = G.run(["a","a","a"]);
    testCase.verifyEqual(status, -1);
    testCase.verifyEqual(q,3);
    % accepting string
    [status, q] = G.run(["a","a","b","b"]);
    testCase.verifyEqual(status, 1);
    testCase.verifyEqual(q,1);
    % feasible but not accepting string
    [status, q] = G.run(["a","a","b","a"]);
    testCase.verifyEqual(status, 0);
    testCase.verifyEqual(q,3);
end