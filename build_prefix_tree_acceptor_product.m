function DFA_list = build_prefix_tree_acceptor_product(...
    list_of_pos_example_set,list_of_neg_example_set,U)
% Generate PTAP from the positive examples.
% Inputs:     list_of_pos_examples --- returned by preprocess_moore_traces
%             U --- the list of the names of actions

if(nargin < 3)
    U = [];
    for i = 1:length(list_of_pos_example_set)
        example = list_of_pos_example_set{i};
        for j = 1:length(example)
            U = [U, example{j}];
        end
    end
    U = sort(unique(U));
end
    
m = length(U);
A = cell(m,1);
for i = 1:m
    A{i} = 0;
end
Q = "eps0"; % define initial state as eps0
Q0 = 1;
G = DFA(1,length(U),A,[],[],[],Q,U,Q0,[],[]);

num_bits = length(list_of_pos_example_set);
F_list = cell(num_bits,1);

% add positive examples in the tree
for i = 1:num_bits
    % for each bit of output
    example = list_of_pos_example_set{i};
    for j = 1:length(example)
        % for each example in i^th bit
        input_trace = example{j};
        [status,q,idx_u] = G.run(input_trace);
        if status == -1
            if q == 1 % if q is initial state
                x_name = "";
            else
                x_name = G.get_x_name(q);
            end
            
            for k = idx_u:length(input_trace)
                x_name = x_name + input_trace(k);
                G.add_x(x_name,[],false);
                if k == length(input_trace)
                    F_list{i}(end+1) = G.n;
                end
                
                G.add_trans(q,input_trace(k),G.n);
                % update q
                q = G.n;
            end
        elseif status == 0
            F_list{i}(end+1) = q;
        end
    end
end

% make sure all negative examples in the tree
for i = 1:num_bits
    % for each bit of output
    example = list_of_neg_example_set{i};
    for j = 1:length(example)
        % for each example in i^th bit
        input_trace = example{j};
        [status,q,idx_u] = G.run(input_trace);
        if status == -1
            if q == 1 % if q is initial state
                x_name = "";
            else
                x_name = G.get_x_name(q);
            end
            
            for k = idx_u:length(input_trace)
                x_name = x_name + input_trace(k);
                G.add_x(x_name,[],false);
                G.add_trans(q,input_trace(k),G.n);
                % update q
                q = G.n;
            end
        end
    end
end



DFA_list = cell(num_bits,1);
for i = 1:num_bits
    DFA_list{i} = G.copy();
    DFA_list{i}.set_final(unique(F_list{i}));
end
end