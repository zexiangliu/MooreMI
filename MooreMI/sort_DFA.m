function sorted_array = sort_DFA(array)
% input a string array, sort the input according to the order
% defined in the paper.

if ~isa(array,"string")
    error("Input array is not a string array.");
end

sorted = sort(array);

str_len = zeros(size(sorted));
for i = 1:length(sorted)
    str_len(i) = strlength(sorted(i));
end

sorted_array = [];
for i = 1:max(str_len)
    idx = str_len == i;
    sorted_array = [sorted_array, sorted(idx)];
end
