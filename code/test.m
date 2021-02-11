NUM_TRIALS = 8; N = 20:20:100; k = size(N,2);

%Time to compute
t_givens = zeros(1, k);
t_poldec = zeros(1, k);
t_newton = zeros(1, k);

%Number of sweeps (givens only)
s_givens = zeros(1, k);
s_poldec = zeros(1, k);

%Symmetry of the Hermitian polar factor
hermDist_givens = zeros(1, k);
hermDist_poldec = zeros(1, k);
hermDist_newton = zeros(1, k);



for i = 1:k
    n = N(i);
    fprintf("n = %4d\n", n);

    temp_t_givens = 0;
    temp_t_poldec = 0;
    temp_t_newton = 0;
    temp_s_givens = 0;
    temp_s_poldec = 0;
    temp_hermDist_givens = 0;
    temp_hermDist_poldec = 0;
    temp_hermDist_newton = 0;

    parfor j = 1:NUM_TRIALS
        A = rand(n);

        [U1, H1, ~] = multiPoldec(A, "single", false, false);

        t = tic();
        [~, Hg, sg, gFlag] = maxtracePoldec(H1, "double", false);
        temp_t_givens = temp_t_givens + toc(t)/NUM_TRIALS;

        t = tic();
        [~, Hp, sp, pFlag] = twobytwoPoldec(H1, "double", false);
        temp_t_poldec = temp_t_poldec + toc(t)/NUM_TRIALS;

        t = tic();
        [~, Hn, ~] = multiPoldec(H1, "double", false, false);
        temp_t_newton = temp_t_newton + toc(t)/NUM_TRIALS;

        temp_s_givens = temp_s_givens + sg/NUM_TRIALS;
        temp_s_poldec = temp_s_poldec + sp/NUM_TRIALS;

        temp_hermDist_givens = temp_hermDist_givens + norm(Hg - Hg', inf) / (2 * NUM_TRIALS);
        temp_hermDist_poldec = temp_hermDist_poldec + norm(Hp - Hp', inf) / (2 * NUM_TRIALS);
        temp_hermDist_newton = temp_hermDist_newton + norm(Hn - Hn', inf) / (2 * NUM_TRIALS);

        if(gFlag || pFlag) fprintf("Householder was used\n"); end
    end

    t_givens(i) = temp_t_givens;
    t_poldec(i) = temp_t_poldec;
    t_newton(i) = temp_t_newton;
    s_givens(i) = temp_s_givens;
    s_poldec(i) = temp_s_poldec;
    hermDist_givens(i) = temp_hermDist_givens;
    hermDist_poldec(i) = temp_hermDist_poldec;
    hermDist_newton(i) = temp_hermDist_newton;

end
