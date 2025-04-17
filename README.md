# N:M Coherence Analysis Toolbox for MATLAB

This repository provides MATLAB functions to compute complex Morlet wavelet transforms and generalized n:m coherence between two signals (e.g., simultaneous LFP recordings). A demonstration script shows how to load sample data, perform time–frequency analyses, and visualize trial-averaged and time-resolved coherence/phase maps.

---
## Repository Contents

- **wavelet_transform.m**  
  Compute complex Morlet wavelet transform (scalogram) for a 1-D signal.

- **computeNMCoherence.m**  
  Calculate generalized n:m coherence magnitude-squared and phase between two time–frequency representations.

- **wavelet_NM_coherence_calc.m**  
  Compute trial-by-trial n:m coherence and phase (in degrees) for 3-D wavelet outputs using `parfor`.

- **demo.m**  
  Load `sample_data.mat` and `colormaps.mat`, perform full-signal wavelet transforms, epoch signals ±500 ms around movement onsets, and plot:
  1. **Figure 1**: Trial-averaged LFPs and spectrograms (dB scale).
  2. **Figure 2**: Animated N:M coherence & phase maps over time and their averages.

- **sample_data.mat**  
  Contains example variables: `el1`, `el2` (LFPs), `movement_onset`, `sampling_rate`.

- **colormaps.mat**  
  Custom colormap definitions used by `demo.m`.

---
## Requirements

- MATLAB R2025a or later
- Signal Processing Toolbox (for `conv`)  
- Parallel Computing Toolbox (for `parfor` in `wavelet_NM_coherence_calc`)  

---
## Usage

1. **Add this folder to your MATLAB path**:
   ```matlab
   addpath(genpath('path/to/this/repo'));
   ```

2. **Run the demo**:
   ```matlab
   demo;
   ```
   This will open two figures illustrating LFP averages, spectrograms, and animated N:M coherence/phase maps.

3. **Call individual functions** for custom analyses:
   ```matlab
   % Compute wavelet transform
   WT = wavelet_transform(signal, Fs, freqs, nCycles);

   % Compute coherence between two TF matrices
   [Coh, Phase] = computeNMCoherence(Xmat, Ymat, mVals, nVals);

   % Compute trialwise coherence
   [Ctrial, Ptrial] = wavelet_NM_coherence_calc(X3D, Y3D, mVals, nVals);
   ```

---
## Function Summaries

### wavelet_transform
```matlab
WT = wavelet_transform(x, Fs, F, n_cyc);
```
- **Inputs**:  
  - `x`: [N×1] signal  
  - `Fs`: sampling rate (Hz)  
  - `F`: [M×1] center frequencies (Hz)  
  - `n_cyc`: [M×1] cycles per wavelet  
- **Output**:  
  - `WT`: [M×N] complex wavelet coefficients.

### computeNMCoherence
```matlab
[NMcoh, NMp] = computeNMCoherence(X, Y, in, out);
```
- **Inputs**:  
  - `X`: [M_X×T] TF matrix for signal X  
  - `Y`: [M_Y×T] TF matrix for signal Y  
  - `in`: 1×M_X integer vector  
  - `out`: 1×M_Y integer vector  
- **Outputs**:  
  - `NMcoh`: [M_X×M_Y] coherence magnitude²  
  - `NMp`: [M_X×M_Y] phase angles (rad).

### wavelet_NM_coherence_calc
```matlab
[NMC, NMP] = wavelet_NM_coherence_calc(X, Y, in, out);
```
- **Inputs**:  
  - `X`: [M_X×T×N] TF array  
  - `Y`: [M_Y×T×N] TF array  
  - `in`, `out`: integer vectors  
- **Outputs**:  
  - `NMC`: Coherence [M_X×M_Y×N]  
  - `NMP`: Phase [M_X×M_Y×N] (deg).

---
## Customization
- **Frequencies & cycles**: Adjust `F` and `n_cyc` in `demo.m` to target different bands.  
- **Epoch window**: Change `win_ms` for longer/shorter peri-event epochs.  
- **Bonferroni correction**: Uncomment in `computeNMCoherence.m` to enforce significance thresholds.

---
## License
This toolbox is released under the MIT License. Feel free to adapt and redistribute.
