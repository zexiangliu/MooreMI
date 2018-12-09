function DFAP = MooreMI(trace_set,U)

[list_of_pos, list_of_neg, bits_to_output] =...
    preprocess_moore_traces(trace_set);

% estimate action alphabet if not exists
if nargin < 2
    U = bits_to_output;
end

DFA_list = build_prefix_tree_acceptor_product(list_of_pos,U); 
N = length(DFA_list);

% construct red and blue cell = {["a1","b1",...], ["a2","b2",...],...}
red = ["eps0"];
[blue,~] = DFA_list{1}.post("eps0");  % get all the successors of the eps
blue = sort_DFA(DFA_list{1}.get_x_name(blue));

% main loop 
while(~isempty(blue))
    q_blue = blue(1);
    blue(1) = [];
 
    merge_accepted = false;  
    
    for i = 1:length(red)
        
        new_DFA_list = copy_DFAs(DFA_list);
        q_red = red(i);
        for j = 1:N
            try
                new_DFA_list{j} = merge(new_DFA_list{j},q_red,q_blue);
            catch 
                keyboard();
            end
        end
        mmflg = true;
        for j = 1:N
            if ~is_consistent(new_DFA_list{j},list_of_neg{j})
                mmflg = false;
                break;
            end
        end 
        if mmflg
            merge_accepted = true;
            break;
        end
    end
    
    
    if merge_accepted
        DFA_list = new_DFA_list; 
        % add more elements into blue
        p_blue = DFA_list{1}.post(red);
        blue(end+1:end+length(p_blue)) = DFA_list{1}.get_x_name(p_blue);
    else 
        red(end+1) = q_blue;
        p_blue = DFA_list{1}.post(q_blue);
        blue(end+1:end+length(p_blue)) = DFA_list{1}.get_x_name(p_blue);
    end
    blue = setdiff(blue,red);
    blue = sort_DFA(unique(blue));
end

DFAP = product(DFA_list,bits_to_output);
DFAP.complete();
end

function new_DFA_list = copy_DFAs(DFA_list)
    N = length(DFA_list);
    new_DFA_list = cell(N,1);
    for j = 1:N
        new_DFA_list{j} = DFA_list{j}.copy(); 
    end
end

function bool = is_consistent(DFA,set_of_neg)
    bool = true;
    % Q0 cannot be final state, because we label Q0 as (0,0,...,0) 
    % by default.
    if ismember(DFA.Q0, DFA.Q_final)
        bool = false;
        return;
    end
    % start from 2 because 1st is eps by default
    for i = 2:length(set_of_neg)
        neg_example = set_of_neg{i};
        status = DFA.run(neg_example);
        % supposed to be feasible and non-accepting
        if status == 1 || status == -1
            bool = false;
            return;
        end
    end
end

function DFAP = product(DFA_list,bits_to_output)
    Q_label = string(DFA_list{1}.Q);
    N = length(DFA_list);
    
    for q = 1:DFA_list{1}.n
        idx = 0;
        for j = 1:N
           if ismember(q,DFA_list{j}.Q_final)
               idx = idx + 2^(N-j);
           end
        end
        Q_label(q) = bits_to_output(idx+1);
    end
    
    DFAP = DFA_list{1}.copy();
    DFAP.set_label(Q_label);
    DFAP.set_final(2:DFAP.n);
end