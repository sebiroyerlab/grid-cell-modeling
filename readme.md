# grid_interference_model

Minimal MATLAB function to simulate a grid-cell model and extract per-cell **theta frequency**, **grid score**, **grid size**.

## Files

* `grid_interference_model.m` — main function

## Requirements

* MATLAB (R2014b+ recommended)
* Signal Processing Toolbox (`findpeaks`)

## Function

```matlab
[thetaF, grid_score, grid_size] = grid_interference_model()
```

### Output

* thetaF: Cells' theta frequencies based on inter-peak-intervals within 50-300ms time range.
* grid_score: Cells' grid scores based on peak/through-ratio of ACGs (using ACGs of firing rate envelops)
* grid_size: Cells' grid scores based on inter-peak-intervals of firing rate envelops (beat size)

* Figure: 
* For 3 cell examples, firing rate (blue) and firing rate envelop (red)
* Firing rate envelops for grid cells sorted by peak position
* Autocorrelograms of envelops for grid cells sorted by grid size
* Cells' theta frequency (mean±s.e.m and std) 
* Cells' grid size (mean±s.e.m)


## Quickstart

```matlab

% To select frequency settings for PI, Cue or OF, uncomment the appropriate line in the 'parameters' section of the script
% then simply type in the command:

[thetaF, grid_score, grid_size] = grid_interference_model();


% Notice: each model run will output a slightly different figure because of the varying random seeds 
