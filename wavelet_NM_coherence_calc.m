function [NMC, NMP] = wavelet_NM_coherence_calc(X, Y, in, out)
    %WAVELET_NM_COHERENCE_CALC Compute trialwise n:m coherence and phase
    %   [NMC, NMP] = WAVELET_NM_COHERENCE_CALC(X, Y, in, out) computes the
    %   generalized n:m coherence magnitude-squared (NMC) and phase in degrees (NMP)
    %   for each trial in the 3D time-frequency arrays X and Y.
    %
    %   Inputs:
    %       X   - [M_X × T × N] complex, TF coefficients of signal X
    %       Y   - [M_Y × T × N] complex, TF coefficients of signal Y
    %       in  - [1 × M_X] integer, input frequencies
    %       out - [1 × M_Y] integer, output frequencies
    %
    %   Outputs:
    %       NMC - [M_X × M_Y × N] double, coherence magnitude-squared per trial
    %       NMP - [M_X × M_Y × N] double, phase angles (degrees) per trial
    %
    %   Example:
    %       % Suppose X and Y are 3D wavelet outputs (freq × time × trials)
    %       m = 1:size(X,1);
    %       n = 1:size(Y,1);
    %       [C_all, P_all] = wavelet_NM_coherence_calc(X, Y, m, n);
    
    % Preallocate outputs: freqs_X × freqs_Y × trials
    [Mx, ~, nTrials] = size(X);
    My = size(Y,1);
    NMC = zeros(Mx, My, nTrials);
    NMP = zeros(Mx, My, nTrials);
    
    % Loop (parallel) over trials
    % Requires Parallel Computing Toolbox for parfor speed-up
    parfor t = 1:nTrials
        % Extract single-trial TF matrices
        Xt = X(:,:,t);  % [Mx × T]
        Yt = Y(:,:,t);  % [My × T]
    
        % Compute coherence and phase for this trial
        [Ct, Pt] = computeNMCoherence(Xt, Yt, in, out);
    
        % Store results
        NMC(:,:,t) = Ct;
        NMP(:,:,t) = Pt;
        % Optional debug: display progress (uncomment if needed)
        % fprintf('Processed trial %d of %d\n', t, nTrials);
    end
    
    % Convert phase from radians to degrees for interpretability
    NMP = rad2deg(NMP);
end
