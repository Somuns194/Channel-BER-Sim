# Channel-BER-Sim
Lab 1 MATLAB Simulation: BER vs SNR for Multiple Modulation under AWGN & Rayleigh Channel

## Project Introduction
This lab implements Monte Carlo simulation to calculate Bit Error Rate (BER) of 5 digital modulation schemes:
QPSK, 8PSK, 16PSK, 16QAM, 32QAM.
Two classic communication channel models are built:
1. AWGN additive white Gaussian noise channel
2. Rayleigh flat multipath fading channel
Simulated BER curves are compared with theoretical analytical formula values to verify communication theory.

## Tech Stack
MATLAB R2023b
- AWGN noise generation
- Rayleigh fading channel modeling
- Multi-modulation signal mapping & demodulation
- Monte Carlo BER statistics
- Theoretical BER formula calculation & curve drawing

## File Directory
- AWGN.m: AWGN channel noise generation function
- Rayleigh.m: Rayleigh flat fading channel generation function
- *.docx / *.pdf: Complete experiment report, homework and theoretical formula comparison analysis document
- Reference papers: Digital modulation & fading channel communication theory

## Experiment Content
1. Build independent AWGN and Rayleigh channel simulation modules.
2. Implement signal modulation & demodulation for QPSK/8PSK/16PSK/16QAM/32QAM.
3. Traverse SNR range, count error bits to draw BER-SNR simulation curve.
4. Calculate theoretical BER formula values, overlay two sets of curves for contrast analysis.
