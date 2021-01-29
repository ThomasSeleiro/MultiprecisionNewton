function Sc = calcSc(S, b)
%calcSc Calculates Sc from Zha and Zhang by removing columns that
%correspond to zero columns of b
    if(size(S,2) ~= size(b,1))
        error("Inputs do not have compatible sizes (S is %dx%d and b"...
            + " is a %d-vector)\n", size(S,1), size(S,2), size(b,1));
    end
    n = size(S,2);
    Sc = S;
    count = 0;
    for i = 1:n*n
        if(b(i) == 0)
            Sc(i-count,:) = [];
            count = count+1;
        end
    end
end


function bc = calcbc(b)
%calcbc Calculates bc from Zha and Zhang by removing zero entries
    n = size(b,1);
    bc = b;
    count = 0;
    for i = 1:n*n
        if(b(i) == 0)
            bc(i-count) = [];
            count = count+1;
        end
    end
end


function Sd = calcSd(Sc)
%calcSd Calculates Sd from Zha and Zhang by removing zero rows of Sc
    n = size(Sc, 1);
    Sd = Sc;
    count = 0;
    for i = 1:n*n
        empty = true;
        for j = 1:n*(n+1)/2
            if(Sc(i,j) ~= 0)
            empty = false;
            break;
            end
        end
        if(empty)
            Sd(i-count, :) = [];
            count = count + 1;
        end
    end
end