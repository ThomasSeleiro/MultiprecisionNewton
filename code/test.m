NUM_TRIALS = 10; N = 100:100:1000; k = size(N,2);

%Time to compute
t_givens = zeros(1, k);
t_newton = zeros(1, k);

%Number of sweeps (givens only)
s_givens = zeros(1, k);

%Unitary polar factor distance
unitDist_givens = zeros(1, k);
unitDist_newton = zeros(1, k);

%Symmetry of the Hermitian polar factor
hermDist_givens = zeros(1, k);
hermDist_newton = zeros(1, k);

%Accuracy of the polar decomposition
accuracy_givens = zeros(1, k);
accuracy_newton = zeros(1, k);

%Difference between the polar factors
unitFactDiff = zeros(1, k);
hermFactDiff = zeros(1, k);


for i = 1:k
    n = N(i);
    disp(n);

    for j = 1:NUM_TRIALS
        A = rand(n);

        t = tic();
        [Ug, Hg, sweep] = maxtracePoldec(A);
        t_givens(i) = t_givens(i) + toc(t)/NUM_TRIALS;

        t = tic();
        [Un, Hn, ~] = multiPoldec(A, "double", false);
        t_newton(i) = t_newton(i) + toc(t)/NUM_TRIALS;

        s_givens(i) = s_givens(i) + sweep/NUM_TRIALS;

        unitDist_givens(i) = unitDist_givens(i) + norm(eye(n) - Ug'*Ug, inf) / NUM_TRIALS;
        unitDist_newton(i) = unitDist_newton(i) + norm(eye(n) - Un'*Un, inf) / NUM_TRIALS;

        hermDist_givens(i) = hermDist_givens(i) + norm(Hg - Hg', inf) / (2 * NUM_TRIALS);
        hermDist_newton(i) = hermDist_newton(i) + norm(Hn - Hn', inf) / (2 * NUM_TRIALS);

        accuracy_givens(i) = accuracy_givens(i) + norm(A - Ug*Hg, inf) / NUM_TRIALS;
        accuracy_newton(i) = accuracy_newton(i) + norm(A - Un*Hn, inf) / NUM_TRIALS;

        unitFactDiff(i) = unitFactDiff(i) + norm(Un - Ug, inf) / NUM_TRIALS;
        hermFactDiff(i) = hermFactDiff(i) + norm(Hn - Hg, inf) / NUM_TRIALS;
    end
end
