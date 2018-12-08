    trace.x = {["a","b","c"],["a","b","d"]};
    trace.y = {["0","1","1","1"],["0","1","1","2"]};

    [list_of_pos, list_of_neg, bits_to_output] = preprocess_moore_traces(trace);
    
    DFA_list = build_prefix_tree_acceptor_product(list_of_pos);