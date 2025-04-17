function demo()
    %DEMO Demonstration of wavelet-based n:m coherence analysis on sample_data.mat
    %   DEMO loads sample_data.mat (el1, el2, movement_onset, sampling_rate), computes
    %   Morlet wavelet transforms over 3–100 Hz, epochs signals ±500 ms around movement
    %   onsets, plots trial-averaged LFPs and spectrograms (Figure 1), and animates
    %   the N:M coherence and phase maps across time, ultimately showing their averages
    %   (Figure 2).
    %
    %   Requirements:
    %     • wavelet_transform.m
    %     • computeNMCoherence.m
    %     • sample_data.mat in current folder
    %
    %   Usage:
    %       >> demo
    
    %% Load data
    D = load('sample_data.mat');
    ColorMaps = load('colormaps.mat');
    el1 = D.el1;      % [1×T] LFP from electrode 1
    el2 = D.el2;      % [1×T] LFP from electrode 2
    movement_onset = D.movement_onset;  % [1×N] onset times (sec)
    Fs = D.sampling_rate;               % sampling frequency (Hz)
    
    %% Define frequencies and wavelet cycles
    F = 3:1:100;   % frequencies of interest (Hz)
    n_cyc = logspace(log10(3), log10(12), length(F));
    
    %% Compute full-signal wavelet transforms
    wt1 = wavelet_transform(el1(:), Fs, F(:), n_cyc(:));  % [freq×time]
    wt2 = wavelet_transform(el2(:), Fs, F(:), n_cyc(:));
    
    %% Epoch parameters (±500 ms)
    win_ms = -500:500;                           % window in milliseconds
    win_samples = round(win_ms * Fs / 1000);     % convert to samples
    nTime = numel(win_ms);
    onset_idx = round(movement_onset * Fs);      % convert onsets to sample indices
    nTrials = length(onset_idx);
    
    % Preallocate epoch arrays
    lfp1_ep = zeros(nTime, nTrials);
    lfp2_ep = zeros(nTime, nTrials);
    spec1_ep = zeros(length(F), nTime, nTrials);
    spec2_ep = zeros(length(F), nTime, nTrials);
    
    %% Extract epochs for each trial
    for t = 1:nTrials
        idx = onset_idx(t) + win_samples;
        % Time-domain LFP epochs
        lfp1_ep(:,t) = el1(idx);
        lfp2_ep(:,t) = el2(idx);
        % Time-frequency epochs
        spec1_ep(:,:,t) = wt1(:, idx);
        spec2_ep(:,:,t) = wt2(:, idx);
    end
    
    %% Figure 1: Trial-averaged LFPs and Spectrograms
    figure(1); clf;
    % Time vector in ms
    
    subplot(2,2,1);
    plot(win_ms, mean(lfp1_ep,2),'LineWidth',3);
    ylabel('LFP (\muV)');
    title('Avg LFP - Electrode 1');
    set(gca,'FontSize',16,'TickDir','out');
    
    subplot(2,2,3);
    plot(win_ms, mean(lfp2_ep,2),'LineWidth',3);
    xlabel('Time (ms)'); ylabel('LFP (\muV)');
    title('Avg LFP - Electrode 2');
    set(gca,'FontSize',16,'TickDir','out');
    
    % Spectrograms: mean power across trials
    spec1_avg = mean(10*log10(abs(spec1_ep).^2), 3);
    subplot(2,2,2);
    imagesc(win_ms, F, spec1_avg);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)');
    title('Spectrogram - Electrode 1'); 
    set(gca,'FontSize',16,'TickDir','out');
    set(gca,'YDir','normal');
    cb = colorbar;             % get colorbar handle
    cb.Label.String = 'Power (dB)';  % set its label
    
    spec2_avg = mean(10*log10(abs(spec1_ep).^2), 3);
    subplot(2,2,4);
    imagesc(win_ms, F, spec2_avg);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)');
    title('Spectrogram - Electrode 2'); 
    set(gca,'FontSize',16,'TickDir','out');
    set(gca,'YDir','normal');
    cb = colorbar;             % get colorbar handle
    cb.Label.String = 'Power (dB)';  % set its label
    colormap(ColorMaps.colormaps(1).Colors{3});
    
    %% Figure 2: Animated N:M coherence & phase maps
    figure(2); clf;
    in = 1:length(F);   % m-values (freq index)
    out = 1:length(F); % n-values
    
    % Preallocate storage for maps across timepoints
    Cmap = zeros(length(in), length(out), nTime);
    Pmap = zeros(length(in), length(out), nTime);
    
    for ti = 1:nTime
        % Extract TF coefficients across trials at this latency
        X_tp = squeeze(spec1_ep(:, ti, :));  % [freq×trials]
        Y_tp = squeeze(spec2_ep(:, ti, :));
    
        % Compute NM coherence & phase (radians)
        [C, P_rad] = computeNMCoherence(X_tp, Y_tp, in, out);
        Cmap(:,:,ti) = C;
        Pmap(:,:,ti) = P_rad;
    
        % Plot instantaneous map
        ax(1) = subplot(1,2,1);
        imagesc(in, out, C);
        axis xy; clim([0 1]); 
        title(sprintf('Coherence at %d ms', win_ms(ti)));
        xlabel('Electrode-1 frequencies'); ylabel('Electrode-2 frequencies');
        set(gca,'FontSize',16,'TickDir','out');
        cb = colorbar;             % get colorbar handle
        cb.Label.String = 'Coherence';  % set its label  
        set(gca,'FontSize',16,'TickDir','out');
        axis square;
    
        ax(2) = subplot(1,2,2);
        imagesc(in, out, rad2deg(P_rad), 'AlphaData', ~isnan(P_rad));
        axis xy; clim([-180 180]); colorbar;
        title(sprintf('Phase at %d ms', win_ms(ti)));
        xlabel('Electrode-1 frequencies'); ylabel('Electrode-2 frequencies');
        cb = colorbar;             % get colorbar handle
        cb.Label.String = 'Phase';  % set its label    
        set(gca,'FontSize',16,'TickDir','out');
        axis square;
    
        colormap(ax(1),'hot');
        colormap(ax(2),ColorMaps.colormaps(6).Colors{6});
        drawnow;
        pause(0.001);
    end
    
    % After animation, show average maps
    C_avg = mean(Cmap, 3);
    P_avg = mean(Pmap, 3);
    ax(1) = subplot(1,2,1);
    imagesc(in, out, C_avg);
    axis xy; clim([0 1]); 
    title('Avg N:M Coherence'); 
    xlabel('Electrode-1 frequencies'); ylabel('Electrode-2 frequencies');
    cb = colorbar;             % get colorbar handle
    cb.Label.String = 'Phase';  % set its label    
    set(gca,'FontSize',16,'TickDir','out');
    axis square;
    
    ax(2) = subplot(1,2,2);
    imagesc(in, out, rad2deg(P_avg));
    axis xy; clim([-180 180]); 
    title('Avg N:M Phase'); 
    xlabel('Electrode-1 frequencies'); ylabel('Electrode-2 frequencies');
    cb = colorbar;             % get colorbar handle
    cb.Label.String = 'Phase';  % set its label    
    set(gca,'FontSize',16,'TickDir','out');
    axis square;

    colormap(ax(1),'hot');
    colormap(ax(2),ColorMaps.colormaps(6).Colors{6});    
end
