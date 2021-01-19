function [S, its, E] ... 
    = multiSign(A, type, debug)
%multiSign Computes matrix sign function in mixed single-double precision.
%   Computes the matrix sign function of input matrix A using a norm scaled
%   Newton iteration in single precision. We then add a correction to the
%   result to increase its commutativity with A in double precision, and
%   perform a final double precision iteration.

    %%Process the input argument
    switch nargin
        case 1
            debug = false;
            type = "single";
            switchPrecision = true;
            typeRoundoff = 1e-8;
        case 2
            debug = false;
            [type, switchPrecision, typeRoundoff] = processTypeInput(type);
        case 3
            debug = debug;
            [type, switchPrecision, typeRoundoff] = processTypeInput(type);
    end
    
    %We declare variables to be used in the function
    n = size(A, 1);
    X = cast(A, type); Xnew = zeros(n);
    iterDist = 1; involDist = norm(A*A - eye(n));
    
    %Print out debug headers if debugging
    if(debug)
        fprintf("\nk  \t |X_k-X_{k-1}|/|X_k| \t     |I - X_k^2|     \n");
        fprintf("---\t---------------------\t---------------------\n");
    end
    
    %Begin the loop to calculate the iterates
    k = 0;
    while(iterDist >= n*typeRoundoff && involDist >= n*typeRoundoff ...
            && k <= 100)
        %Calculate the next Newton iterate
        [Xnew, iterDist, involDist] = signNewtonStep(X, k, debug);
        
        %MORE METRICS FOR ANALYSIS CAN BE PLACED HERE
        
        %Move Xnew into X and increment k
        X = Xnew;
        k = k+1;
    end
    its = k;
    
    %If specified in the function inputs, we switch precision and iterate
    %again
    if(switchPrecision && type == "single")
        X = cast(X, "double");
        typeRoundoff = 1e-16;
    
        %Print the distance from commutativity of the converted iterate
        if (debug)
            fprintf("\n----------CONVERTING TO DOUBLE PRECISION---------\n");
            fprintf("Norm of the commutator:\t%e\n", norm(X*A - A*X, inf));
        end

        %We add a correction to the iterate to improve its commutativity with A
        [X, E] = commutLSQ(X, A, debug);

        %We calculate at most 3 iterations, since quadratic convergence
        %should ensure that the convergence is unitary after only one
        %iteration (in practice this might not be guaranteed immediately
        %and the accuracy can be improved by one or two more iterations)
        newK = 0;
        while(iterDist >= n*typeRoundoff && involDist >= n*typeRoundoff ...
            && newK <= 3)
            %Calculate the next Newton iterate
            [Xnew, iterDist, involDist] = signNewtonStep(X,k+newK,debug);
            
            %MORE METRICS FOR ANALYSIS CAN BE PLACED HERE
            if(debug)
                fprintf("Norm of the commutator:\t%e\n", norm(X*A - A*X, inf));
            end
        
            %Move Xnew into X and increment newK
            X = Xnew;
            newK = newK+1;
        end
        its = its + newK;
    end

    %Form S to return it
    S = X;
end



function [type, switchPrecision, typeRoundoff] = processTypeInput(inputType)
%processTypeInput Processes the type argument for multiPoldec
    if ismember(inputType, ["single", "s"])
        type = "single";
        switchPrecision = true;
        typeRoundoff = 1e-8;
    elseif ismember(inputType, ["double", "d"])
        type = "double";
        switchPrecision = false;
        typeRoundoff = 1e-16;
    elseif ismember(inputType, ["singleOnly", "so"])
        type = "single";
        switchPrecision = false;
        typeRoundoff = 1e-8;
    else
        error('Input argument %s not recognised. Use "single", '...
            + '"double" or "singleOnly"', inputType);
    end
end



function [Xnew, iterDist, involDist] = signNewtonStep(X, k, debug)
%poldecNewtonStep Computes one step of the polar decomp Newton iteration
%   Function that computes the next iterate in the scaled Newton iteration
%   using the inf norm scaling (to reduce the computational overhead from
%   using the determinant and spectral scaling factor).
%   Note that k corresponds to the index of the iterate X NOT Xnew

    %Store the inverse to only calculate it once
    invX = inv(X);
    %Calculate the 1,inf norm scaling factor
    mu = sqrt(norm(invX, inf) / norm(X, inf));
    
    %Calculate the new iterate
    Xnew = (mu * X + invX / mu) / 2;
    
    %Calculate various metrics from the iterates
    iterDist = norm(Xnew - X, inf) / norm(Xnew, inf);
    involDist = norm(eye(n) - Xnew * Xnew, inf);
    
    %Print the metrics if debugging
    if(debug)
        fprintf("%3d\t  %.11e  \t  %.11e  \n", k+1, iterDist, involDist);
    end
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
