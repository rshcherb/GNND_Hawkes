## Generalized Nearest Neighbour Distance (GNND) method and Hawkes point process modelling

The scripts were written for the following publication: M. Sedghizadeh, R. Shcherbakov, M. van den Berghe (2025) Generalized Nearest-Neighbor Distance and Hawkes point process modeling applied to mining induced seismicity, submitted to The Seismic Record.

### Directory structure
**MATLAB/** folder contains the Matlab scripts that perform specific tasks. The following subfolderes contain functions:
- **etas/** to fit the ETAS model.
- **hawkes/** to fit the Hawkes point process.
- **fm/** to analyze the frequency-magnitude statistics.
- **nnd/** to perfomr GNND analysis.

### The output of the fits are written into the following folders
- **eq_rate_etas/**
- **eq_rate_hawkes/** 
- **fm/**
- **nnd/**

**common_parameters.m** script contains information about the geographical boundaries of the full and target regions, time intervals, initial conditions, and other parameters.

### To perform analsyis:
1. run **fm_analysis.m** to model the ferquency-magnitude distribution of the seismic catalog.
2. run **nnd_analysis.m** to perform the GNND analysis.
3. run **rate_analysis.m** to fit the Hawkes point process or ETAS model.

These entry scripts specify all the initial parameters needed to perform the corresponding tasks. **Model** structure is used to pass the information into most functions. 
