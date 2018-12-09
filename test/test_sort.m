function tests = test_sort
tests = functiontests(localfunctions);
end

function testSortDFA(testCase)
    str = ["b","a","aa","c","ad"];
    sort_str = sort_DFA(str);
    testCase.verifyEqual(sort_str,["a","b","c","aa","ad"]);
end