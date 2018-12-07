function blue = add_new(blue,q_blue)
% blue is a cell = {["a1","b1",...], ["a2","b2",...],...}
N = length(blue);
flag = "ud";
for i = 1:N
        if (select_small(blue{i},q_blue) == blue{i})
        else
            if ( i == 1)
                blue = [{q_blue}, blue];
            else
                blue = [blue{1:i-1}, {q_blue}, blue{i:end}];
            end
            flag = "d";
            break;
        end
end       
        if (flag == "ud")
            blue{end+1} = q_blue;
        end
end