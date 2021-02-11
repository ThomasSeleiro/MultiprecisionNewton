NUM_TRIALS = 4; N = 10:5:100; k = size(N,2);

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
    fprintf("n = %4d\n", n);

    temp_t_givens = 0;
    temp_t_newton = 0;
    temp_s_givens = 0;
    temp_unitDist_givens = 0;
    temp_unitDist_newton = 0;
    temp_hermDist_givens = 0;
    temp_hermDist_newton = 0;
    temp_accuracy_givens = 0;
    temp_accuracy_newton = 0;
    temp_unitFactDiff = 0;
    temp_hermFactDiff = 0;

    parfor j = 1:NUM_TRIALS
        A = rand(n);

        t = tic();
        [Ug, Hg, sweep] = maxtracePoldec(A);
        temp_t_givens = temp_t_givens + toc(t)/NUM_TRIALS;

        t = tic();
        [Un, Hn, ~] = multiPoldec(A, "double", false, false);
        temp_t_newton = temp_t_newton + toc(t)/NUM_TRIALS;

        temp_s_givens = temp_s_givens + sweep/NUM_TRIALS;

        temp_unitDist_givens = temp_unitDist_givens + norm(eye(n) - Ug'*Ug, inf) / NUM_TRIALS;
        temp_unitDist_newton = temp_unitDist_newton + norm(eye(n) - Un'*Un, inf) / NUM_TRIALS;

        temp_hermDist_givens = temp_hermDist_givens + norm(Hg - Hg', inf) / (2 * NUM_TRIALS);
        temp_hermDist_newton = temp_hermDist_newton + norm(Hn - Hn', inf) / (2 * NUM_TRIALS);

        temp_accuracy_givens = temp_accuracy_givens + norm(A - Ug*Hg, inf) / NUM_TRIALS;
        temp_accuracy_newton = temp_accuracy_newton + norm(A - Un*Hn, inf) / NUM_TRIALS;

        temp_unitFactDiff = temp_unitFactDiff + norm(Un - Ug, inf) / NUM_TRIALS;
        temp_hermFactDiff = temp_hermFactDiff + norm(Hn - Hg, inf) / NUM_TRIALS;
    end

    t_givens(i) = temp_t_givens;
    t_newton(i) = temp_t_newton;
    s_givens(i) = temp_s_givens;
    unitDist_givens(i) = temp_unitDist_givens;
    unitDist_newton(i) = temp_unitDist_newton;
    hermDist_givens(i) = temp_hermDist_givens;
    hermDist_newton(i) = temp_hermDist_newton;
    accuracy_givens(i) = temp_accuracy_givens;
    accuracy_newton(i) = temp_accuracy_newton;
    unitFactDiff(i) = temp_unitFactDiff;
    hermFactDiff(i) = temp_hermFactDiff;

end
