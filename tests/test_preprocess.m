function tests = test_preprocess
tests = functiontests(localfunctions);
end

function testPreprocess(testCase)

    trace.x = {["a","b","c"],["a","b","d"]};
    trace.y = {["0","1","1","1"],["0","1","1","2"]};

    [list_of_pos, list_of_neg, bits_to_output] = preprocess_moore_traces(trace);
    
    testCase.verifyEqual(size(list_of_pos),[2,1]);
    testCase.verifyEqual(list_of_pos{1}{1},["a","b","d"])
    testCase.verifyEqual(list_of_neg{1}{1},"a");
    testCase.verifyEqual(list_of_neg{1}{2},["a","b"]);
    testCase.verifyEqual(bits_to_output,["0","1","2"]);
    
end
