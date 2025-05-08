# Simulation Code – TUMA Achievability Bound for AWGN Channel

## Overview
This repository contains  MATLAB codes for computing the TUMA achievability bound for an AWGN MAC with a Zipf-like type profile. The code corresponds to the bound presented in [arXiv:2504.19916](https://arxiv.org/abs/2504.19916).

## How to Use
1. Run `generateSubsetData.m` to generate subset structures (`.mat` file).
2. Run `main_compute_tv_bound_zipf_type.m` to evaluate the bound.

## Dependencies
- MATLAB R2022b or later
- No explicit toolboxes required

## Notes
- Modify parameters like `Ma`, `Ka`, `n`, `k`, `eps_target` in the main script.
- Ensure the `.mat` file generated matches the filename expected in the main script.

## Contact
Deekshith Pathayappilly Krishnan  
Chalmers University of Technology, Gothenburg, Sweden  
deepat@chalmers.se
