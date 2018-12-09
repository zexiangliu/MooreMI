function tests = test_MooreMI
tests = functiontests(localfunctions);
end

function testMooreMI(testCase)
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