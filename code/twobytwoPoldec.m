function [U, H, sweeps, householder] = twobytwoPoldec(A, u, debug)
%twobytwoPoldec Computes the polar decomposition by using 2x2 poldec

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
        %Loop through the principal submatrices and calculate their polar
        %decompositions
        for i = 1:n-1
            for j = i+1:n
                [A, W] = submatrixPoldec(A, W, i, j);
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
        G = eye(n) - 2*x*x';
        A = G' * A;
        W = W  * G;
        householder = true;
    end
    
    %Finally from U and H
    U = W;
    H = A;
end


function [Anew, Wnew] = submatrixPoldec(A, W, i, j)
%submatrixPoldec Computes the 2x2 poldec of the i,jth submatrix

    %Calculate the polar decomposition of the 2x2 principal submatrix
    Aij = A([i j], [i j]);
    detAij = det(Aij);
    detXinvAijT = sign(detAij) * [A(j,j) -A(j,i); -A(i,j) A(i,i)];
    lambda = abs(det(Aij + detXinvAijT))^(-0.5);
    Uij = lambda * (Aij + detXinvAijT);
    % Hij = lambda * (Aij'*Aij + abs(detAij)* eye(2));

    %Form the new matrices
    Anew = A;
    Wnew = W;

    %Update the matrices by performing the matrix product
    Anew(i, :) = Uij(1,1) * A(i, :) + Uij(2,1) * A(j, :);
    Anew(j, :) = Uij(1,2) * A(i, :) + Uij(2,2) * A(j, :);
    Wnew(:, i) = Uij(1,1) * W(:, i) + Uij(2,1) * W(:, j);
    Wnew(:, j) = Uij(1,2) * W(:, i) + Uij(2,2) * W(:, j);

end
