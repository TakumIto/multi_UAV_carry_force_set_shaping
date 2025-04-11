# multi_UAV_carry_force_set_shaping
Basic code associated with the paper **Feasible Force Set Shaping for a Payload-Carrying Platform Consisting of Tiltable Multiple UAVs Connected Via Passive Hinge Joints**.
[arxiv version](https://arxiv.org/abs/2503.00341) 

---
## Requirement
- MATLAB R2023b or more

## Initialize
MATLAB Command window
```matlab
>> run init.m
```

## Physical parameter settings
Edit `settings/physical_parameters.m`


## Programs
### Force set shaping
Set required force set (RFS) on `tilt_angle_optimization.m > %% dedine the parameters` and run
```matlab
>> run tilt_angle_optimization.m
```

### Simulation with pre-optimized tilt angles
Run
```matlab
>> run multi__UAV_simulation_main.m
```

Notes
- Pre-optimized tilt angles are at `data/optimal_tilt_angles.m`.
- Trajectory states and inputs are at `data/Fp_nom.mat` and `data/x_nom.mat`, respectively.
- You can set log fine name on `multi__UAV_simulation_main.m > %% Simulation parameters`.

