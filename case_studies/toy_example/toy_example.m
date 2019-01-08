%% toy example
trace_set.x = {["a","a","a","a","a","a"]};
trace_set.y = {[0,1,2,3,3,3,3]};

% output: 0 --- "0"
%         1 --- "0","1","d"
%         2 --- "0","1","d","3"
%         3 --- "0","1","d","3","2"
%           

DPFA = MooreMI(trace_set,["a"])