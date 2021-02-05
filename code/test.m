NUM_TRIALS = 1000; n = 16;
inputMatrices = {eye(8), hilb(6), magic(6), hadamard(8)}

testPoldec(n, NUM_TRIALS);

fprintf("\n\n ========     POLDEC    ======== \n\n");
[singleTime, doubleTime, singleDoubleDiff, singleSkew, singleAcc, ...
    doubleSkew, doubleAcc, singleAvgIts, doubleAvgIts] ...
    = testInputPoldec(inputMatrices);
singleTime
doubleTime
singleDoubleDiff
singleSkew
singleAcc
doubleSkew
doubleAcc
singleAvgIts
doubleAvgIts


fprintf("\n\n ======== SIGN FUNCTION ======== \n\n");
[singleDist, doubleDist, singleDoubleDist, ...
    correctedMultiprecisionTime, noCorrectionTime, doubleTime, ...
    correctedIts, noCorrectionIts, doubleIts] ...
    = testInputSign(inputMatrices);
singleDist
doubleDist
singleDoubleDist
correctedMultiprecisionTime
noCorrectionTime
doubleTime
correctedIts
noCorrectionIts
doubleIts



function [singleTime, doubleTime, singleDoubleDiff, singleSkew, ...
    singleAcc, doubleSkew, doubleAcc, singleAvgIts, doubleAvgIts] ...
    = testInputPoldec(matrices)
    
    k = size(matrices,2);
   singleTime = zeros(1,k); doubleTime = zeros(1,k);
   singleDoubleDiff = zeros(1,k);
   singleSkew = zeros(1,k); singleAcc = zeros(1,k);
   doubleSkew = zeros(1,k); doubleAcc = zeros(1,k);
   singleAvgIts = zeros(1,k); doubleAvgIts = zeros(1,k);
    
   for i = 1:k
        A = matrices{i};
        time = tic();
        [U_s, H_s, its_s] = multiPoldec(A, "single", false);
        singleTime(i) = toc(time);
        time = tic();
        [U_d, H_d, its_d] = multiPoldec(A, "double", false);
        doubleTime(i) = toc(time);

        singleDoubleDiff(i) = norm(U_s - U_d, inf);

        singleSkew(i) = norm(H_s - H_s', inf)/2;
        singleAcc(i) = norm(A - U_s*H_s);
        doubleSkew(i) = norm(H_d - H_d', inf)/2;
        doubleAcc(i) = norm(A - U_d*H_d);

        singleAvgIts(i) = its_s;
        doubleAvgIts(i) = its_d;
   end
end


function [singleTime, doubleTime, singleDoubleDiff, singleSkew, ...
    singleAcc, doubleSkew, doubleAcc, singleAvgIts, doubleAvgIts] ...
    = testPoldec(n, numTrials)
    %Initialise Variables
    singleTime = 0; doubleTime = 0;
    singleDoubleDiff = 0;
    singleSkew = 0; singleAcc = 0;
    doubleSkew = 0; doubleAcc = 0;
    singleAvgIts = 0; doubleAvgIts = 0;
    
    %Compute polar decomposition numTrials number of times
    for i = 1:numTrials
        A = rand(n);
        time = tic();
        [U_s, H_s, its_s] = multiPoldec(A, "single", false);
        singleTime = singleTime + toc(time);
        time = tic();
        [U_d, H_d, its_d] = multiPoldec(A, "double", false);
        doubleTime = doubleTime + toc(time);
        
        singleDoubleDiff = singleDoubleDiff + norm(U_s - U_d, inf);
        
        singleSkew = singleSkew + norm(H_s - H_s', inf)/2;
        singleAcc = singleAcc + norm(A - U_s*H_s);
        doubleSkew = doubleSkew + norm(H_d - H_d', inf)/2;
        doubleAcc = doubleAcc + norm(A - U_d*H_d);
        
        singleAvgIts = singleAvgIts + its_s;
        doubleAvgIts = doubleAvgIts + its_d;
    end

    %Form averages of the quantities
    singleTime = singleTime / numTrials;
    doubleTime = doubleTime / numTrials;
    singleDoubleDiff = singleDoubleDiff / numTrials;
    singleSkew = singleSkew / numTrials;
    singleAcc = singleAcc / numTrials;
    doubleSkew = doubleSkew / numTrials;
    doubleAcc = doubleAcc / numTrials;
    singleAvgIts = singleAvgIts / numTrials;
    doubleAvgIts = doubleAvgIts / numTrials;
end


function [singleDist, doubleDist, singleDoubleDist, ...
    correctedMultiprecisionTime, noCorrectionTime, doubleTime, ...
    correctedIts, noCorrectionIts, doubleIts] = testInputSign(matrices)

    k = size(matrices, 2);
    singleDist = zeros(1,k); doubleDist = zeros(1,k);
    singleDoubleDist = zeros(1,k);

    correctedMultiprecisionTime = zeros(1,k);
    noCorrectionTime = zeros(1,k); doubleTime = zeros(1,k);

    correctedIts = zeros(1,k); noCorrectionIts = zeros(1,k);
    doubleIts = zeros(1,k);
    
    for i = 1:k
        A = matrices{i};
        time = tic();
        [S_s, its_s, E_s] = multiSign(A, "single", true, false);
        correctedMultiprecisionTime(i) = toc(time);
        time = tic();
        [S_noE, its_noE, E_noE] = multiSign(A, "single", false, false);
        noCorrectionTime(i) = toc(time);
        time = tic();
        [S_d, its_d, E_d] = multiSign(A, "double", false, false);
        doubleTime(i) = toc(time);

        singleDist(i) = norm(S_s - S_noE);
        doubleDist(i) = norm(S_d - S_noE);
        singleDoubleDist(i) = norm(S_d - S_s);

        correctedIts(i) = its_s;
        noCorrectionIts(i) = its_noE;
        doubleIts(i) = its_d;
    end
end



function [singleDist, doubleDist, singleDoubleDist, ...
    correctedMultiprecisionTime, noCorrectionTime, doubleTime, ...
    correctedIts, noCorrectionIts, doubleIts] = testSign(n, numTrials)
    %initialise variables
    singleDist = 0; doubleDist = 0; singleDoubleDist = 0;
    correctedMultiprecisionTime = 0; noCorrectionTime = 0; doubleTime = 0;
    correctedIts = 0; noCorrectionIts = 0; doubleIts = 0;
    
    %Compute sign function numTrials number of times
    for i = 1:numTrials
        A = rand(n);
        time = tic();
        [S_s, its_s, E_s] = multiSign(A, "single", true, false);
        correctedMultiprecisionTime = correctedMultiprecisionTime ...
            + toc(time);
        time = tic();
        [S_noE, its_noE, E_noE] = multiSign(A, "single", false, false);
        noCorrectionTime = noCorrectionTime + toc(time);
        time = tic();
        [S_d, its_d, E_d] = multiSign(A, "double", false, false);
        doubleTime = doubleTime + toc(time);

        singleDist = singleDist + norm(S_s - S_noE);
        doubleDist = doubleDist + norm(S_d - S_noE);
        singleDoubleDist = singleDoubleDist + norm(S_d - S_s);

        correctedIts = correctedIts + its_s;
        noCorrectionIts = noCorrectionIts + its_noE;
        doubleIts = doubleIts + its_d;
    end

    %Form averages of the quantities
    singleDist = singleDist / numTrials;
    correctedMultiprecisionTime = correctedMultiprecisionTime / numTrials;
    doubleDist = doubleDist / numTrials;
    noCorrectionTime = noCorrectionTime / numTrials;
    singleDoubleDist = singleDoubleDist / numTrials;
    doubleTime = doubleTime / numTrials;

    correctedIts = correctedIts / numTrials;
    noCorrectionIts = noCorrectionIts / numTrials;
    doubleIts = doubleIts / numTrials;

end