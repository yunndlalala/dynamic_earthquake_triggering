# DynTriPy 
## A package for detecting dynamic triggering.
-------
### Install
Run the following code in the terminal:

```
pip install dyntripy
```

### Usage 
#### *Triggering* class
*Triggering* is the implement of HiFi [Yun et al., 2019](https://doi.org/10.1029/2019GL083913) to detect dynamic triggering. An instance should be defined firstly by passing an input file containing all parameters used in the detection process.    
We provide some test data, examples of input file and performing script of *dyntripy* on [our github](https://github.com/yunndlalala/dynamic_earthquake_triggering/blob/master/test.zip). Please download and run them to test whether the *dyntripy* package can perform normally.   
Detailed descriptions of the key words in the input file are as follows:
- "data_source"
	+ "station_file": file containing station names
	+ "waveform_path": folder of raw sac/mseed data
	+ "response_file": file of PZ files summary 
	+ "remote_earthquake_catalog": catalog file of remote earthquakes
- "net_database"
    + "days": [M_1, M_2]; power integrals of M1 days before and M2 days after the dates of distant earthquakes are computed
	+ "time_segment": length of time segment
	+ "frequency_segment": [min, step, max]; minimum (min), segment length (step), and maximum (max) of the frequency segments 
    + "database_path": path to store the power integral database
- "net_ratio"
	+ "RE_path": folder to store the logarithmic ratio values of the remote earthquakes
	+ "background_days": [N_1, N_2]; select N_1 days before and N_2 days after the dates of distant earthquakes as background days
    + "background_catalog": catalog file of the virtual events in background days
    + "RB_path": folder to store the logarithmic ratio values of the background days
- "net_cl"
	+ "matched_ratio_path": folder to store the matched ratio of remote earthquakes and background days
	+ "cl_path": folder to store the confidence levels
	+ "threshold": threshold of confidence level to identify triggering
	+ "figure_path": folder to store the figures about the distribution of the logarithmic ratio values for each remote earthquake at each station. If don't want to plot these figures, just set this parameter as "None".
	
#### Utility
- dyntripy.utils.**gen_time_windows**(*catalog_file=None, reference_lat=None, reference_lon=None, tb=18000, te_b_vel=5.0, te_e_vel=2.0, out__file=None*)  
Generate a class of $T_b$ and $T_e$, and output a catalog of remote earthquakes with origin time, start and end of $T_b$ and $T_e$.   
    **Parameters:**   
    + **catalog_file**: str  
    The raw catalog of remote earthquakes with source locations.
    + **reference_lat**: float  
    The latitude of the reference point in the study area.
    + **reference_lon**: float  
    The longitude of the reference point in the study area.
    + **tb**: float  
    The length of $T_b$ in seconds.
    + **te_b_vel**: float  
    The velocity of a surface wave whose arrival is the start of $T_e$.
    + **te_e_vel**: float  
    The velocity of a surface wave whose arrival is the end of $T_e$.
    + **out_file**: str
    The catalog of remote earthquakes with the origin time, start and end time of $T_b$ and $T_e$.  

    **Returns:**
    + None