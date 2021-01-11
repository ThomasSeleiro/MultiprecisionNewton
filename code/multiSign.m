function [singleS, S, its, devArray, commutArray] = multiSign(A, debug)
%MULTISIGN Computes matrix sign function in mixed single-double precision.
%   Computes the matrix sign function of input matrix A using a norm scaled
%   Newton iteration in single precision. We then add a correction to the
%   result to increase its commutativity with A in double precision, and
%   perform a final double precision iteration.

    %Check inputs
    if(nargin == 1) debug = false; end
    
    n = size(A, 2);
    k = 1;
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
        fprintf("k  \t |X_k-X_{k-1}|/|X_k| \t     |I - X_k^2|     \n");
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
        iterDist(1)  = norm(singleXnew - singleX, inf) / norm(singleX, inf);
        iterDist(2)  = norm(Xnew - X, inf) / norm(X, inf);
        involDist(1) = norm(singleX*singleX - eye(n), inf);
        involDist(2) = norm(X*X - eye(n), inf);
        devArray     = [devArray norm(singleXnew - cast(Xnew, "single"), inf)];
        commutArray  = [commutArray; norm(singleA*singleXnew - singleXnew * singleA, inf) norm(A*Xnew - Xnew*A, inf)];

        %Print the metrics if debugging
        if(debug)
            fprintf("%3d\t%.4e %.4e\t%.4e %.4e\n", k, iterDist(1), iterDist(2), involDist(1), involDist(2));
        end
        
        %Move newX into X and increment k
        singleX = singleXnew;
        X = Xnew;
        k = k+1;
    end
    
    %We calculate the final iteration of the algorithm by converting the
    %single precision iterate to double precision and running the algorithm
    %one more time
    singleX = cast(singleX, "double");
    
    %Print the distance from commutativity of the converted iterate
    if (debug)
        fprintf("\n----------CONVERTING TO DOUBLE PRECISION---------\n");
        fprintf("Norm of the commutator:\t%e\n", norm(singleX*A - A*singleX, inf));
    end
    
    %We add a correction to the iterate to improve its commutativity with A
    singleX = singleX;
    
    %We repeat one final iteration of the algorithm in double precision
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
    iterDist(1)  = norm(singleXnew - singleX, inf) / norm(singleX, inf);
    iterDist(2)  = norm(Xnew - X, inf) / norm(X, inf);
    involDist(1) = norm(singleX*singleX - eye(n), inf);
    involDist(2) = norm(X*X - eye(n), inf);
    devArray     = [devArray norm(singleXnew - cast(Xnew, "single"), inf)];
    commutArray  = [commutArray; norm(singleA*singleXnew - singleXnew * singleA, inf) norm(A*Xnew - Xnew*A, inf)];

    %Print the metrics if debugging
    if(debug)
        fprintf("%3d\t%.4e %.4e\t%.4e %.4e\n", k, iterDist(1), iterDist(2), involDist(1), involDist(2));
    end

    %Move newX into X and increment k
    singleX = singleXnew;
    X = Xnew;
    k = k+1;
    
    singleS = singleX;
    S = X;
    its = k-1;
end

