function [U, H, sweeps, householder] = maxtracePoldec(A, u, debug)
%maxtracePoldec Computes polar decomp with Givens and Householder matrices
%Computes the polar decomposition of a matrix A by maximising the trace of
%A by applying Givens and Householder transformations.
%The algorithm is outlined in Matthew Smith's PhD thesis (2002 p.84)

    switch nargin
        case 1
            u = float_params("double");
            debug = false;
        case 2
            debug = false;
    end

    n = size(A,1);
    W = eye(n);
    householder = false;

    %Check the input u is not input as a string
    if(isa(u,'string'))
        %Otherwise try converting it using float_params
        u = float_params("double");
    end

    if(debug)
        lastTrace = trace(A);
        fprintf("\nSweep     |A-A'|/|A|     traceDiff \n");
        fprintf("===================================\n");
    end

    %We perform several sweeps to make the A symmetric and maximise its
    %trace
    sweeps = 0;
    symmDist = norm(A - A', inf) / norm(A, inf);
    while(symmDist > u*n)
        %We first make every diagonal element positive.
        for i = 1:n
            if A(i,i) < 0
                A(i, :) = -A(i, :);
                W(i,i) = -W(i,i);
            end
        end

        %We apply Givens rotations/reflexions to A to make the matrix
        %symmetric
        for i = 1:n-1
            for j = i+1:n
                if(det(A([i,j],[i,j])) >= 0)
                    [A, W] = givensMax(A, W, i, j);
                else
                    [A, W] = givensReflMax(A, W, i, j);
                end
            end
        end

        symmDist = norm(A - A', inf) / norm(A, inf);
        sweeps = sweeps + 1;
        if(debug)
           traceDiff = abs((trace(A) - lastTrace) / trace(A));
           %fprintf("Sweep     |A-A'|/|A|     traceDiff \n");
            fprintf(" %3d      %.4e"   +"     %.4e\n", sweeps, ...
                symmDist, traceDiff);
            lastTrace = trace(A);
        end
    end

    %If the resulting matrix is not positive definite apply a householder
    %transformation
    [V, D] = eig(A);
    [lambdaMin, indexMin] = min(diag(D));
    if(lambdaMin < 0)
        x = V(:, indexMin);
        G = eye(n) - 2*(x*x');
        A = G' * A;
        W = W  * G;
        householder = true;
    end


    %Finally from U and H
    U = W;
    H = A;
end


function [Anew, Wnew] = givensMax(A, W, i, j)
%givensMax Returns the max trace symm submatrix Aij using a Givens

    %Calculate the optimal angle of rotation
    Aij = A([i j], [i j]);
    x = [Aij(2,1) - Aij(1,2)   Aij(1,1) + Aij(2,2)];
    c =  x(2) / norm(x,2);
    s = -x(1) / norm(x,2);

    %Form the new matrices
    Anew = A;
    Wnew = W;

    %Update the matrices by performing the matrix product
    Anew(i, :) = c*A(i, :) - s*A(j, :);
    Anew(j, :) = s*A(i, :) + c*A(j, :);
    Wnew(:, i) = c*W(:, i) - s*W(:, j);
    Wnew(:, j) = s*W(:, i) + c*W(:, j);
end



function [Anew, Wnew] = givensReflMax(A, W, i, j)
%givensReflMax Returns the max trace symm submatrix Aij using a Givens refl

    %Calculate the optimal angle of rotation
    Aij = A([i j], [i j]);
    x = [Aij(2,2) - Aij(1,1)   Aij(1,2) + Aij(2,1)];
    c = -x(1) / norm(x,2);
    s =  x(2) / norm(x,2);

    %Form the new matrices
    Anew = A;
    Wnew = W;

    %Update the matrices by performing the matrix product
    Anew(i, :) = c*A(i, :) + s*A(j, :);
    Anew(j, :) = s*A(i, :) - c*A(j, :);
    Wnew(:, i) = c*W(:, i) + s*W(:, j);
    Wnew(:, j) = s*W(:, i) - c*W(:, j);
end
