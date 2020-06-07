# DynTriPy 

A package for detecting dynamic triggering.
=======

Download all files and run setup.py to install this package.

We provide some test data and an example script in the folder of *test*. 

### Input file 
Before calculating, the parameters and data sources are needed to be configured in the file *input.json*. 
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
