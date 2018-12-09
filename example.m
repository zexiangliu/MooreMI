%% example 1
trace_set.x = {["a","c"],["b","d"]};
trace_set.y = {["0","1","2"],["0","3","2"]};

DPFA = MooreMI(trace_set,["a";"b";"c";"d"])

%% example 2 Fig.5 in the paper

trace_set.x = {["b","a"],["b","b"],["a","a"],["a","b"],...
    ["a","b","a"],["a","b","b"]};
trace_set.y = {["0","1","2"],["0","1","2"],["0","2","0"],["0","2","2"]...
    ["0","2","2","2"],["0","2","2","2"]};

DPFA = MooreMI(trace_set,["a";"b";"c";"d"])
