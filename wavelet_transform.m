function wtresult = wavelet_transform(input, Fs, F, n_cyc)
%WAVELET_TRANSFORM Compute complex Morlet wavelet transform (scalogram).
%   WTRESULT = WAVELET_TRANSFORM(INPUT, Fs, F, n_cyc) returns the complex-valued
%   wavelet coefficients for the 1-D signal INPUT sampled at Fs (Hz), evaluated
%   at each frequency in vector F (Hz) using the corresponding cycle counts in
%   n_cyc.
%
%   Inputs:
%       INPUT   - [N×1] double, column vector of time-domain samples
%       Fs      - double, sampling frequency (Hz)
%       F       - [M×1] vector of center (carrier) frequencies (Hz)
%       n_cyc   - [M×1] vector specifying the number of cycles per wavelet
%
%   Output:
%       WTRESULT - [M×N] complex, wavelet coefficients (frequency × time)
%
%   Example:
%       t = (0:1/1000:1-1/1000)';           % 1 second at 1 kHz
%       x = sin(2*pi*10*t);                % 10 Hz sine wave
%       F = [5; 10; 20];                    % frequencies of interest
%       n_cyc = [4; 7; 10];                 % cycles for each frequency
%       WT = wavelet_transform(x, 1000, F, n_cyc);
%       
%       Notes on defining n_cyc: 
%       The vector n_cyc specifies, for each center frequency in F, 
%       how many oscillation cycles your Morlet wavelet should contain. In practice:
%       A wavelet with n_cyc cycles at frequency f lasts about n_cyc / f seconds, 
%       so more cycles ⇒ a longer wavelet.
%       Longer wavelets (higher n_cyc) give you better frequency resolution
%       (narrower bandwidth) but poorer temporal precision.
%       Conversely, fewer cycles yield sharper time localization 
%       at the expense of frequency specificity.
%       You supply one n_cyc value per target frequency so you can 
%       balance the time–frequency trade‑off differently across your band of interest.
%       For example, you can define n_cyc =
%       logspace(log10(3),log10(12),length(F)) so
%       you end up with ~3 cycles at your lowest F smoothly ramping up to
%       ~12 cycles at your highest F (very sharp frequency precision).
%       Low freqs (3 cycles): if your lowest F is, say, 5 Hz, 3 cycles
%       span 0.6 s—enough to localize a slow rhythm in time.
%       High freqs (12 cycles): if your highest F is 100 Hz, 12 cycles span
%       0.12 s—short enough to track fast events but narrow enough in frequency.
%       You can tweak the 3→12 endpoints to shift the time–frequency trade‑off:
%       lower minimum cycles for crisper timing, higher maximum for crisper frequency. 
%       Adjust based on your signal's bandwidth and the durations you need to resolve.

% Validate inputs
if size(input,2) > 1
    error('wavelet_transform:InvalidInput', 'INPUT must be a 1D column vector');
end
if isempty(F)
    error('wavelet_transform:NoFrequencies', 'Frequency vector F is empty');
end
if isempty(input)
    error('wavelet_transform:EmptySignal', 'Input signal is empty');
end
if numel(F) ~= numel(n_cyc)
    error('wavelet_transform:LengthMismatch', ...
        'Length of F and n_cyc must match');
end

% Number of time points in input
Npoints = length(input);

% Precompute Morlet wavelets for each frequency
padding  = 0;                  % half-length of largest wavelet for zero-padding
wavelets = cell(numel(F),1);
for i = 1:numel(F)
    % Generate complex Morlet wavelet for center freq F(i) with n_cyc(i)
    wavelets{i} = cxmorlet(F(i), n_cyc(i), Fs);
    % Track maximum half-length for buffer padding
    padding = max(padding, floor(numel(wavelets{i})/2));
end

% Prepare zero-padded input buffer to avoid edge artifacts
buffer = zeros(Npoints + 2*padding, 1);
% Allocate output matrix: rows=frequencies, cols=time points
wtresult = zeros(numel(F), Npoints);

% Indices delimiting the original signal within the buffer
bufStart = padding + 1;
bufEnd   = padding + Npoints;

% Convolve each wavelet with the padded signal
for i = 1:numel(F)
    % Place the raw input into central buffer region
    buffer(bufStart:bufEnd) = input;
    % Convolution with 'same' returns central part of full convolution
    convResult = conv(buffer, wavelets{i}, 'same');
    % Extract region corresponding to original signal length
    wtresult(i,:) = convResult(bufStart:bufEnd)';
end
end


%% Subfunctions
function w = cxmorlet(Fc, Nc, Fs)
% CXMORLET Generate complex Morlet wavelet
%   W = CXMORLET(Fc, Nc, Fs) returns a column vector W containing a
%   zero-mean, unit-norm complex Morlet wavelet centered at frequency Fc
%   (Hz) with Nc cycles, sampled at Fs (Hz).

% Temporal standard deviation: ensure wavelet spans ±2.5 SD
sd = (Nc/(2*Fc)) / 2.5;
% Total wavelet length: cover ±3 SD (~6 SD total), ensure odd length
wl = 2*floor((6*sd*Fs)/2) + 1;
w  = zeros(wl,1);
half = floor(wl/2);

% Normalize factor accumulator
normFactor = 0;
for k = 1:wl
    t = (k-1-half)/Fs;
    % Gaussian-windowed complex sinusoid
    w(k) = bw_cf(t, sd, Fc);
    % For normalization, accumulate real Gaussian envelope
    normFactor = normFactor + gauss(t, sd);
end
% Normalize wavelet so area under envelope = 1
w = w / normFactor;
end

function res = bw_cf(t, bw, cf)
% BW_CF Compute bandwidth-limited complex sinusoid
%   RES = BW_CF(T, BW, CF) = exp(2i*pi*CF*T) .* gaussian(T, BW)
%   scaled by 1/(BW*sqrt(2*pi)).
cnorm = 1/(bw*sqrt(2*pi));
res   = cnorm * exp(-(t.^2)/(2*bw^2)) .* exp(2i*pi*cf*t);
end

function res = gauss(t, sd)
% GAUSS Compute real Gaussian envelope
%   RES = GAUSS(T, SD) = exp(-(T^2)/(2*SD^2)) scaled by 1/(SD*sqrt(2*pi)).
res = (1/(sd*sqrt(2*pi))) * exp(-(t.^2)/(2*sd^2));
end
