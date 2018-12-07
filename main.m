[list_of_pos, list_of_neg, bits_to_output,output_total] = preprocess_moore_traces(trace_set);
DFA.list = build_prefix_tree_acceptor_product(); 
N = ceil(log2(length(output_total)));

% construct red and blue cell = {["a1","b1",...], ["a2","b2",...],...}
red = {"eps"};
blue = {};
[p_blue,~] = DFA.list{1}.successor(DFA.list{0},"eps");  % get all the successors of the eps
if(length(p_blue) == length("eps")+1)               % check if they are one-letter successor of eps
    blue = {blue,p_blue};
end

% main loop 
while(length(blue) ~= 0)
    q_blue = blue{1};
    blue(1) = [];
 
    merge_accepted = 0;  
    new_DFA_list = DFA_list; 
    
    for i = 1:length(red)
        for j = 1:N
            q_red = red{j};
            new_DFA_list{j} = merge(DFA_list{j},q_red,q_blue);
        end
          for j = 1:N
            for k = 1:length(list_of_neg)
                if(is_consistent(new_DFA_list{j},list_of_neg{k})==0)
                    mm_flg = 1;
                    break;
                end
            end
            if (mm_flg == 1)
                break;
            end
          end 
          
            if (mm_flg == 0)
                merge_accepted = 1;
            end
            mm_flg = 0;
            
            if (merge_accepted == 1)
                DFA_list = new_DFA_list; 
                for j = 1:N
                    p_blue = DFA_list{1}.successor(DFA_list{1},q_red);
                    if(length(p_blue) == length(q_red) + 1)
                        blue = add_new(blue,p_blue);
                    end
                end
            else 
                red = add_new(red,q_blue);
                for j = 1:N
                    p_blue = DFA_list{1}.successor(DFA_list{1},q_blue);
                    if(length(p_blue) == length(q_blue) + 1)
                        blue = add_new(blue,p_blue);
                    end
                end
            end
    end
end

DFAP = product(DFA_list,bits_to_output);
DFAPF = DFAP.make_complete(DFAP);

    