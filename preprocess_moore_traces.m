function [list_of_pos, list_of_neg, bits_to_output] = preprocess_moore_traces(trace_set)
% Input: trace_set --- a struct that contains trace.x and trace.y,
%                      where trace.x = {input 1, ..., input n} and trace.y
%                      = {output 1, ..., output n}. Both input and output
%                      need to be represented by array of strings.
% Output: list_of_pos --- a cell of postive examples. The order is
%                         [1 x x] --- list_of_pos{1}
%                         [x 1 x] --- list_of_pos{2}
%                         [x x 1] --- list_of_pos{3}
%         list_of_neg --- a cell of negative examples. The order is similar
%                         to list_of_pos.
%         bits_to_output --- a mapping between output strings and their
%                            index. e.g. bits_to_output = ["a","b","c"] 
%                            <---> indices [0,1,2].

    x_list = trace_set.x;
    y_list = trace_set.y;
    
    % compute #output
    output_total = [];
    for i = 1:length(y_list)
        output_total = [output_total,y_list{i}];        
    end
    output_total = unique(output_total);
    % make sure that the first element is the output of initial state
    bits_to_output = output_total;
    bits_to_output(2:end) = sort(output_total(2:end));
    
    N = ceil(log2(length(output_total)));
    
    list_of_pos = cell(N,1);
    list_of_neg = cell(N,1);
    
    % put epsilon string into negative example
    for i = 1:N
        list_of_pos{i} = {};
        list_of_neg{i} = {};
    end
    
    % pick up a trace
    for i = 1:length(x_list)
        x = x_list{i};
        y = y_list{i};
        
        % pick up a prefix
        for j = 1:length(x)
            y_out = find(bits_to_output == y(j+1))-1;
            
            counter = N;
            while( counter ~= 0 )
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