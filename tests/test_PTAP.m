function tests = test_PTAP
tests = functiontests(localfunctions);
end

function testPreprocess(testCase)

    trace.x = {["a","b","c"],["a","b","d"]};
    trace.y = {["0","1","1","1"],["0","1","1","2"]};

    [list_of_pos, list_of_neg, ~] = preprocess_moore_traces(trace);
    
    DFA_list = build_prefix_tree_acceptor_product(list_of_pos,list_of_neg);
    
    % the ground truth is 
    % 1 --a--> 2 --b--> 3 --c--> 5
    %                   | --d--> 4
    DFA1 = DFA_list{1};
    DFA2 = DFA_list{2};
    testCase.verifyEqual(DFA1.Q_final, DFA1.get_x_idx("abd"));
    testCase.verifyEqual(DFA2.Q_final, DFA2.get_x_idx(["a";"ab";"abc"])');
    testCase.verifyEqual(DFA1.A{1}(1,2),true);
    testCase.verifyEqual(DFA1.A{2}(2,3),true);
    testCase.verifyEqual(DFA1.A{3}(3,5),true);
    testCase.verifyEqual(DFA1.A{4}(3,4),true);
end