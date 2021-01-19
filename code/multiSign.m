function [singleS, S, its, devArray, commutArray, E] = multiSign(A, debug)
%multiSign Computes matrix sign function in mixed single-double precision.
%   Computes the matrix sign function of input matrix A using a norm scaled
%   Newton iteration in single precision. We then add a correction to the
%   result to increase its commutativity with A in double precision, and
%   perform a final double precision iteration.

    %Check inputs
    if(nargin == 1) debug = false; end
    
    n = size(A, 2);
    k = 0;
    muFunc = @(B, invB, p) sqrt(norm(invB, p) / norm(B, p));
    
    singleA = cast(A, "single");
    singleX = cast(A, "single");
    X = cast(A, "double");
    singleXnew = zeros(n);
    Xnew = zeros(n);
    
    %Array to store the single and double values for |X_{k+1} - X_k|/|X_k|
    iterDist    = ones(1, 2);
    %Array to store the single and double values for |X_k^2 - I| (ie the
    %distance to an involutory matrix)
    involDist   = ones(1, 2);
    %List to store the deviation between single and double matrices
    devArray    = [];
    %List to store the commutativity of the iterates with A
    commutArray = [];
    
    %Print out debug headers if debugging
    if(debug)
        fprintf("\nk  \t |X_k-X_{k-1}|/|X_k| \t     |I - X_k^2|     \n");
        fprintf("---\t---------------------\t---------------------\n");
        fprintf("   \t  Single     Double  \t  Single     Double  \n");
        fprintf("===\t=====================\t=====================\n");
    end
    
    %While loop that calculates the Newton iterates for the sign function
    %in single and double precision, until the single precision iterate
    %converges.
    while(iterDist(1) >= 1e-8 * n && involDist(1) >= 1e-8 * n && k <= 100)
        %Store the inverses X_k to only have to calculate them once.
        singleInvX = inv(singleX);
        invX = inv(X);
        
        %Calculate the scaling factors mu
        singleMu = muFunc(singleX, singleInvX, inf);
        mu = muFunc(X, invX, inf);
        
        %Calculate the new iterates
        singleXnew = (singleMu * singleX + singleInvX / singleMu) / 2;
        Xnew = (mu * X + invX / mu) / 2;
        
        %Calculate various metrics from the iterates
        iterDist(1)  = norm(singleXnew - singleX, inf) / norm(singleXnew, inf);
        iterDist(2)  = norm(Xnew - X, inf) / norm(Xnew, inf);
        involDist(1) = norm(singleXnew*singleXnew - eye(n), inf);
        involDist(2) = norm(Xnew*Xnew - eye(n), inf);
        devArray     = [devArray norm(singleXnew - cast(Xnew, "single"), inf)];
        commutArray  = [commutArray; norm(singleA*singleXnew - singleXnew * singleA, inf) norm(A*Xnew - Xnew*A, inf)];
        
        %Move Xnew into X and increment k
        singleX = singleXnew;
        X = Xnew;
        k = k+1;

        %Print the metrics if debugging
        if(debug)
            fprintf("%3d\t%.4e %.4e\t%.4e %.4e\n", k, iterDist(1), iterDist(2), involDist(1), involDist(2));
        end
    end
    
    %We calculate the final iteration of the algorithm by converting the
    %single precision iterate to double precision and running the algorithm
    %one more time
    castX = cast(singleX, "double");
    
    %Print the distance from commutativity of the converted iterate
    if (debug)
        fprintf("\n----------CONVERTING TO DOUBLE PRECISION---------\n");
        fprintf("Norm of the commutator:\t%e\n", norm(castX*A - A*castX, inf));
    end
    
    %We add a correction to the iterate to improve its commutativity with A
    [castX, E] = commutLSQ(castX, A, debug);
    %[castX, E] = randomCommut(castX, A, 10, debug);
    
    %We repeat one final iteration of the algorithm in double precision
    castIts = 0;
    castIterDist = 1;
    castInvolDist = 1;
    castXnew = zeros(n);
    while(castIterDist >= 1e-16 * n && castInvolDist >= 1e-16 *n && k <= 100)
        %Store the inverses X_k to only have to calculate them once.
        castInvX = inv(castX);

        %Calculate the scaling factors mu
        castMu = muFunc(castX, castInvX, inf);

        %Calculate the new iterates
        castXnew = (castMu * castX + castInvX / castMu) / 2;

        %Calculate various metrics from the iterates
        castIterDist  = norm(castXnew - castX, inf) / norm(castXnew, inf);
        castInvolDist = norm(castXnew*castXnew - eye(n), inf);

        %Move castXnew into castX and increment castIts
        castX = castXnew;
        castIts = castIts+1;
        
        %Print the metrics if debugging
        if(debug)
            fprintf("%3d\t%.4e           \t%.4e           \n", k+castIts, castIterDist, castInvolDist);
            fprintf("Norm of the commutator:\t%e\n", norm(castX*A - A*castX, inf));
        end
    end
    
    singleS = castXnew;
    S = X;
    its = k+castIts;
end


function [Xnew, E] = commutLSQ(X, A, debug)
%commutLSQ Computes the nearest matrix to X that commutes with A by LSQ
%   Computes the lowest norm matrix E such that X + E commutes with A by
%   solving a least squares problem in a system of n^2 equations (for X and
%   A nxn matrices). The function returns the matrix X + E.
    
    %We check the input matrices have the same size
    n = size(X,1);
    if (size(A,1) ~= n)
        error("Matrices X (%dx%d) and A (%dx%d) have incompatible sizes\n", size(X,1), size(X,2), size(A,1), size(A,2));
    end
    
    %We form the matrix of the n^2 system and the target vector b
    S = kron(eye(n), A) - kron(A', eye(n));
    b = -S * X(:);
    
    %We solve the least squares problem using the SVD to find the minimal
    %norm solution
    [U, S, V] = svd(S);
    vE = 0;
    for i = 1:rank(S)
        vE = vE + ((U(:,i)' * b) * V(:,i)) / S(i,i);
    end
    
    %Form the corresponding matrix E and add it to X
    E = reshape(vE, [n,n]);
    Xnew = X + E;
    
    if(debug)
        fprintf("         --- Commutativity correction ---        \n");
        fprintf("              Norm(E) = %.4e                     \n", norm(E, inf));
        fprintf("       Norm of the commutator:\t%e\n", norm(Xnew*A - A*Xnew, inf));
    end
end

function [Xnew, E] = randomCommut(X, A, numTrials, debug)

    %We check the input matrices have the same size
    n = size(X,1);
    if (size(A,1) ~= n)
        error("Matrices X (%dx%d) and A (%dx%d) have incompatible sizes\n", size(X,1), size(X,2), size(A,1), size(A,2));
    end
    
    E = zeros(n);
    commut = norm((X+E)*A - A*(X+E), inf);
    for i=1:numTrials
        temp = rand(n);
        temp = temp * (5e-7 / norm(temp, inf));
        tempX = X+temp;
        tempNorm = norm(tempX*A - A*tempX, inf);
        if tempNorm < commut
            E = temp;
            commut = tempNorm;
        end
    end
    Xnew = X + E;
    
    if(debug)
        fprintf("         --- Commutativity correction ---        \n");
        fprintf("              Norm(E) = %.4e                     \n", norm(E, inf));
        fprintf("       Norm of the commutator:\t%e\n", norm(Xnew*A - A*Xnew, inf));
    end
end
