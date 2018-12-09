function tests = test_MooreMI
tests = functiontests(localfunctions);
end

function testExample1(testCase)
    % ground truth 0 --a--> 1 --c--> 2
    %              |--b--> 3 --d--> 2
    trace_set.x = {["a","c"],["b","d"]};
    trace_set.y = {["0","1","2"],["0","3","2"]};

    DPFA = MooreMI(trace_set,["a";"b";"c";"d"]);
    testCase.verifyEqual(DPFA.n, 4);
    testCase.verifyEqual(DPFA.m, 4);
    testCase.verifyEqual(DPFA.Q_label, ["0","1","2","3"]);
    testCase.verifyEqual(DPFA.Q_final, 2:4);
    A{1} = logical([0 1 0 0
                    0 1 0 0
                    0 0 1 0
                    0 0 0 1]);
    A{2} = logical([0 0 0 1
                    0 1 0 0
                    0 0 1 0
                    0 0 0 1]);
    A{3} = logical([1 0 0 0
                    0 0 1 0
                    0 0 1 0
                    0 0 0 1]);
    A{4} = logical([1 0 0 0
                    0 1 0 0
                    0 0 1 0
                    0 0 1 0]);

    testCase.verifyEqual(DPFA.A,A');
end

function testExample2(testCase)
    % example in Fig 5 in the paper (a)
    trace_set.x = {["b","a"],["b","b"],["a","a"],["a","b"],...
        ["a","b","a"],["a","b","b"]};
    trace_set.y = {["0","1","2"],["0","1","2"],["0","2","0"],["0","2","2"]...
        ["0","2","2","2"],["0","2","2","2"]};

    DPFA = MooreMI(trace_set,["a";"b"]);
    testCase.verifyEqual(DPFA.n, 4);
    testCase.verifyEqual(DPFA.m, 2);
    testCase.verifyEqual(DPFA.Q_label, ["0","1","2","2"]);
    A = {};
    A{1} = logical([  0   0   1   0
                      0   0   1   0
                      1   0   0   0
                      0   0   1   0]);
    A{2} = logical([ 0   1   0   0
                     0   0   1   0
                     0   0   0   1
                     0   0   1   0]);
     testCase.verifyEqual(DPFA.A,A');
end

function testExample3(testCase)
    % example in Fig 5 in the paper (b)

    trace_set.x = {["b","a"],["b","b"],["a","a"],["a","b"],...
        ["a","b","a"],["a","b","b"],["b","a","a"],["b","b","a"]};
    trace_set.y = {["0","1","2"],["0","1","2"],["0","2","0"],["0","2","2"]...
        ["0","2","2","2"],["0","2","2","2"],["0","1","2","2"],["0","1","2","2"]};

    DPFA = MooreMI(trace_set,["a";"b"]);
    testCase.verifyEqual(DPFA.n, 4);
    testCase.verifyEqual(DPFA.m, 2);
    testCase.verifyEqual(DPFA.Q_label, ["0","1","2","2"]);
    A = {};
    A{1} = logical([   0   0   1   0
                       0   0   0   1
                       1   0   0   0
                       0   0   1   0]);
    A{2} = logical([   0   1   0   0
                       0   0   0   1
                       0   0   0   1
                       0   0   1   0]);
     testCase.verifyEqual(DPFA.A,A');
end

