function [list_of_pos, list_of_neg, bits_to_output] = preprocess_moore_traces(trace_set)
    
    x_list = trace_set.x;
    y_list = trace_set.y;
    
    % compute #output
    output_total = [];
    for i = 1:length(y_list)
        output_total = [output_total,y_list{i}];        
    end
    output_total = unique(output_total);
    % by default "0" is from epsilon string,
    % which is the first element of bits_to_output
    bits_to_output = sort(output_total);
    
    N = ceil(log2(length(output_total)));
    
    list_of_pos = cell(N,1);
    list_of_neg = cell(N,1);
    
    % put epsilon string into negative example
    for i = 1:N
        list_of_pos{i} = {};
        list_of_neg{i} = {};eps
        list_of_neg{i}{1} = "eps";
    end
    
    % pick up a trace
    for i = 1:length(x_list)
        x = x_list{i};
        y = y_list{i};
        
        % pick up a prefix
        for j = 1:length(x)
            y_out = find(bits_to_output == y(j+1));
            
            counter = N;
            while( y_out ~= 0 )
                if(mod(y_out,2))
                    list_of_pos{counter}{end+1} = x(1:j);
                else
                    list_of_neg{counter}{end+1} = x(1:j);
                end
                y_out = floor(y_out/2);
                counter = counter - 1;
            end
        end
    end
    
end