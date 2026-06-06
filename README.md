# Iterative Dispersive Vold-Kalman Filter (IDVKF)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-blue.svg)](https://www.mathworks.com/)

Official MATLAB implementation of the paper:

> **Iterative Dispersive Vold-Kalman Filter with Local Adaptive Bandwidth for Dispersive Signal Decomposition in Structural Health Monitoring**
>
> *IEEE Transactions on Industrial Informatics*, Accepted June 2026.

## Abstract

Dispersive signals, with frequency-dependent propagation velocities, are omnipresent in structural health monitoring (SHM) and non-destructive testing. This paper proposes a novel decomposition algorithm, **Iterative Dispersive Vold-Kalman Filter (IDVKF)**, for accurate group delay (GD) estimation and robust damage feature detection from multi-component broadband dispersive signals. Specifically:

- A **dispersive VKF (DVKF)** framework is developed by formulating the decomposition problem in the **frequency domain**, enabling effective demodulation of highly dispersive signals.
- **Variable filtering bandwidths** and a **joint decomposition scheme** are integrated to handle overlapping components and non-uniform bandwidths.
- A **GD refinement strategy** iteratively updates GD estimates via reconstructed envelopes, improving time–frequency resolution.
- An **iterative local bandwidth adaptation** mechanism, inspired by windowed signal orthogonality, suppresses noise and mitigates feature leakage under frequency-dependent bandwidth variations.



## File Structure

```
IDVKF_release/
├── main_IDVKF.m            # Main demo script (runs the full pipeline)
├── IDVKF.m                 # Core IDVKF algorithm
├── ridgeDetectMult_F.m     # Adaptive TF ridge detection (frequency domain)
├── RPRG.m                  # Ridge Path ReGrouping
├── IFsmooth.m              # Smoothed group delay estimation from ridge indices
├── DFRE_F.m                # Frequency-domain ridge extraction
├── mySFFT.m                # Short-frequency Fourier transform
├── Differ.m                # Numerical differentiation
├── curvesmooth.m           # Curve smoothing utility
├── low_filter.m            # FIR low-pass filter (frequency domain)
├── sigshow.m               # Time-domain signal visualization
├── stftshow.m              # TF spectrum visualization
├── data_generation_origin.m  # Generate clean 3-component synthetic signal
├── data_generation_noise.m   # Generate noisy synthetic signal (5 dB SNR)
├── data_origin.mat           # Pre-generated clean signal (by data_generation_origin.m)
├── data_5dB.mat              # Pre-generated noisy signal (by data_generation_noise.m)
└── LICENSE
```

## Quick Start

### Run the numerical example

The repository includes a pre-generated 3-component synthetic dispersive signal at 5 dB SNR (sampling frequency 2000 Hz, duration 10 s). To run the full decomposition demo:

```matlab
main_IDVKF
```

This script will:
1. Extract TF ridge curves from the frequency-domain signal
2. Apply RPRG to resolve ridge crossings
3. Run IDVKF to decompose the three dispersive components
4. Plot the estimated GDs, frequency-domain components, and time-domain waveforms

### Generate synthetic data

To regenerate the synthetic signals from scratch:

```matlab
data_generation_origin   % clean signal -> data_origin.mat
data_generation_noise    % noisy signal (5 dB) -> data_5dB.mat
```

### Use IDVKF on your own signal

```matlab
% Inputs
%   Sig    : frequency-domain signal (length-N vector)
%   T      : signal duration (s)
%   iniGD  : initial group delays (M x N matrix, one row per component)
%   r0     : initial bandwidth controller (scalar or M x 1 vector)
%   betaC  : GD smoothness controller
%   winLen : local bandwidth adaptation window length
%   tol    : convergence threshold
%   maxit  : maximum iterations

[SigFest, GDest, IAest] = IDVKF(Sig, T, iniGD, r0, betaC, winLen, tol, maxit);
```

### Hyperparameter guidance

| Parameter | Recommended range | Effect |
|-----------|-------------------|--------|
| `r0` | `5e2`–`1e4` | Initial bandwidth controller; larger → narrower bandwidth |
| `betaC` (µ) | `1e-5`–`1e-11` | GD smoothness; smaller → smoother GD estimates |
| `winLen` (W) | ~10% of signal length | Bandwidth adaptation window; smaller → more local but noisier |
| `tol` | `1e-5` | Convergence threshold |

> For large initial GD errors (>10%), use smaller `r0` (e.g., `1e2`) and smaller `betaC` (e.g., `1e-9`) to allow wider initial bandwidth and stronger smoothing.

## Results

### Numerical Simulation (5 dB SNR)

Quantitative comparison over 50 independent runs on a 3-component synthetic dispersive signal:

| Method | Output SNR C1 (dB) | Output SNR C2 (dB) | Output SNR C3 (dB) | Runtime (s) |
|--------|--------------------|--------------------|--------------------|----|
| **IDVKF** | **21.85±0.23** | **16.66±0.11** | **22.56±0.22** | 1.34 |
| GHST   | 12.04±0.19 | 13.13±0.18 | 13.81±0.26 | 27.01 |
| IFETF  | 20.76±0.35 | 13.86±0.08 | 21.48±0.41 | 0.56 |
| GDMD   | 19.88±0.34 | 15.72±0.19 | 20.31±0.36 | 1.09 |

### Experimental Validation

Experimental data are not released.

- **Lamb wave decomposition** (aluminum alloy plate, 500×500×3.7 mm): IDVKF achieves at least 2 dB improvement in GD estimation accuracy over all baselines while maintaining comparable runtime (~0.16 s).
- **Rail damage detection** (100% low-floor tram, 60-s axle box vibration): IDVKF successfully extracts all impulsive damage features, completing detection in 2.1 s — 10× faster than TMSST, and outperforming IFETF and GDMD in low-amplitude impulse recovery.

## Citation

The paper is accepted and currently in early access. BibTeX will be updated upon publication. For now, please cite as:

```bibtex
@article{jiang2026idvkf,
  author  = {Jiang, Yuan and Jiang, Yiyue and Liu, Zheng and Chen, Yuejian and Wang, Pingfeng},
  title   = {Iterative Dispersive Vold-Kalman Filter with Local Adaptive Bandwidth
             for Dispersive Signal Decomposition in Structural Health Monitoring},
  journal = {IEEE Transactions on Industrial Informatics},
  year    = {2026},
  note    = {Early Access, Paper No. TII-26-2374},
  doi     = {to be updated upon publication}
}

@article{jiang2024iterative,
  title={An iterative adaptive Vold--Kalman filter for nonstationary signal decomposition in mechatronic transmission fault diagnosis under variable speed conditions},
  author={Jiang, Yuan and Chen, Yuejian and Wang, Pingfeng},
  journal={IEEE Transactions on Industrial Informatics},
  volume={20},
  number={8},
  pages={10510--10519},
  year={2024},
  publisher={IEEE}
}

@article{jiang2022iterative,
  title={An iterative frequency-domain envelope-tracking filter for dispersive signal decomposition in structural health monitoring},
  author={Jiang, Yuan and Niu, Gang},
  journal={Mechanical Systems and Signal Processing},
  volume={179},
  pages={109329},
  year={2022},
  publisher={Elsevier}
}
```

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

Yuan Jiang — yuanjiang.phm@gmail.com

Department of Industrial and Enterprise Systems Engineering

University of Illinois Urbana-Champaign
