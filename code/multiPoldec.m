function [U, H, its] = multiPoldec(A, type, debug)
%multiPoldec Computes the polar decomposition in multiprecision
    
    %Process the input argument
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
    iterDist = 1; unitDist = norm(A'*A - eye(n));
    
    %Print out debug headers if debugging
    if(debug)
        fprintf("\nk  \t |X_k-X_{k-1}|/|X_k| \t    |I - X_k*X_k|    \n");
        fprintf("---\t---------------------\t---------------------\n");
    end
    
    %Begin the loop to calculate the iterates
    k = 0;
    while(iterDist >= n*typeRoundoff && unitDist >= n*typeRoundoff ...
            && k <= 100)
        %Calculate the next Newton iterate
        [Xnew, iterDist, unitDist] = poldecNewtonStep(X, k, debug);
        
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
        
        %We calculate at most 3 iterations, since quadratic convergence
        %should ensure that the convergence is unitary after only one
        %iteration (in practice this might not be guaranteed immediately
        %and the accuracy can be improved by one or two more iterations)
        newK = 0;
        while(iterDist >= n*typeRoundoff && unitDist >= n*typeRoundoff ...
            && newK < 3)
            %Calculate the next Newton iterate
            [Xnew, iterDist, unitDist] = poldecNewtonStep(X,k+newK,debug);
            
            %MORE METRICS FOR ANALYSIS CAN BE PLACED HERE
        
            %Move Xnew into X and increment newK
            X = Xnew;
            newK = newK+1;
        end
        its = its + newK;
        
        %Set U and calculate H
        U = X;
        H = U' * cast(A, "double");
    else
        %Set U and calculate H
        U = X;
        H = U' * cast(A, type);
    end
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



function [Xnew, iterDist, unitDist] = poldecNewtonStep(X, k, debug)
%poldecNewtonStep Computes one step of the polar decomp Newton iteration
%   Function that computes the next iterate in the scaled Newton iteration
%   using the 1,inf norm scaling (to reduce the computational overhead from
%   using the optimal scaling factor using the 2 norm).
%   Note that k corresponds to the index of the iterate X NOT Xnew

    %Store the inverse to only calculate it once
    invX = inv(X);
    %Calculate the 1,inf norm scaling factor
    mu = (norm(invX,1) * norm(invX,inf)/(norm(X,1) * norm(X,inf)))^(0.25);
    
    %Calculate the new iterate
    Xnew = (mu * X + invX' / mu) / 2;
    
    %Calculate various metrics from the iterates
    iterDist = norm(Xnew - X, inf) / norm(Xnew, inf);
    unitDist = norm(eye(n) - Xnew' * Xnew, inf);
    
    %Print the metrics if debugging
    if(debug)
        fprintf("%3d\t  %.11e  \t  %.11e  \n", k+1, iterDist, unitDist);
    end
end