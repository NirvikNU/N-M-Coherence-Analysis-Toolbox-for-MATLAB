function [NMcoh, NMp] = computeNMCoherence(X,Y,in,out)
    %COMPUTENMCOHERENCE Compute generalized n:m coherence and phase between two signals
    %   [NMCOH, NMP] = COMPUTENMCOHERENCE(X, Y, in, out) returns the magnitude-squared
    %   coherence (NMCOH) and phase angle (NMP) for all n:m frequency couplings between
    %   wavelet tranforms of X and Y.
    %
    %   Inputs:
    %       X   - [frequency (M_X) × trials] complex, wavelet tranform of signals at a particular time  
    %       Y   - [frequency (M_Y) × trials] complex, wavelet tranform of output signals at the same time
    %       in  - [1 × M_X] integer, input frequencies
    %       out - [1 × M_Y] integer, output frequencies
    %
    %   Outputs:
    %       NMcoh - [M_X × M_Y] double, coherence magnitude squared for each n:m pair
    %       NMp   - [M_X × M_Y] double, phase angles (radians) of the cross-spectrum
    
    % Ensure vectors 'in' and 'out' are row and column
    in  = in(:)';    % row vector [1 × M_X]
    out = out(:);    % column [M_Y × 1]
    
    % Compute least common multiple for each m:n pair
    In  = repmat(in,  length(out), 1);  % [M_Y × M_X]
    Out = repmat(out, 1, length(in));   % [M_Y × M_X]
    LCM = lcm(In, Out);
    
    % Mask for non-integer and non-isofrequency couplings
    mask = ~((In == Out) | mod(In,Out) == 0 | mod(Out,In) == 0);
    
    % Compute exponent N and M matrices to raise X and Y
    N = LCM ./ In;      % [M_Y × M_X]
    M = LCM ./ Out;     % [M_Y × M_X]

    % Replicate across trial dimension
    N = repmat(N(:), 1, size(X,2));    % [(M_X·M_Y) x trials]
    M = repmat(M(:), 1, size(Y,2));    % same size  

    % % Normalize by magnitude (comment-in to ignore signal amplitude)
    % X = X ./ abs(X); Y = Y ./ abs(Y);
    
    % Reshape wavelet tranform matrices into long form for exponentiation
    Xr = repelem(X, size(Y,1), 1);     % [(M_X·M_Y) × trials]
    Yr = repmat(Y, size(X,1), 1);      % [(M_X·M_Y) × trials]
    
    % Compute autospectra: E[|X|^(2N)] and E[|Y|^(2M)]
    ASD_X  = mean(abs(Xr).^(2 * N), 2);  % [(M_X·M_Y)×1]
    ASD_Y  = mean(abs(Yr).^(2 * M), 2);  % [(M_X·M_Y)×1]
    
    % Compute cross-spectrum: E[X^N · conj(Y^M)]
    CSD_XY = mean((Xr.^N) .* conj(Yr.^M), 2);  % [(M_X·M_Y)×1]
    
    % Coherence magnitude squared
    NMcoh = abs(CSD_XY).^2 ./ (ASD_X .* ASD_Y);
    % Reshape back to [M_X×M_Y]
    NMcoh = reshape(NMcoh,[length(out),length(in)])';
    
    % Phase angles
    NMp = angle(CSD_XY);
    NMp = reshape(NMp,[length(out),length(in)])';    
        
    % Optional: Bonferroni correction for multiple comparisons
    % Uncomment the following lines to apply (remove if not required)
    % threshold = 1.0 - (0.05/length(N(:)))^(1/(size(X,2) - 1));
    % NMcoh(NMcoh < threshold) = 0;
    
    % Replace NaNs with zeros for coherence
    NMcoh(isnan(NMcoh)) = 0;
    
    % Mask out non-integer couplings
    NMcoh(mask') = NaN;
    NMp(mask')   = NaN;    
    
end



