function smaller = select_small(a1,a2)
% a1 and a2 are 2 arrays of states
N1 = length(a1);
N2 = length(a2);

if (N1<N2)
    smaller = a1;
elseif (N1>N2)
    smaller = a2;
else
    for (i = 1:1:N1)
        if a1(i)<a2(i)
            smaller = a1;
            break;
        elseif a1(i)>a2(i)
            smaller = a2;
            break;
        else
            smaller = "und";
        end
    end
    if (smaller == "und")
        smaller = a1;
    end
end
    