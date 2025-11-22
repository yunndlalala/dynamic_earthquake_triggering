# DynTriPy 
## A package for detecting dynamic triggering.
-------
### Install
Run the following code in the terminal in the source code package:

```
python setup.py install
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
	+ "nperseg": number of points for each time segment
	+ "noverlap": overlap of points for two time segments
	+ "nfft": point number to calculate fft
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
	
#### Notes
example_data_prepare.py gives some tools to preprocess the data. Only function of catalog_processing() is suitable for most of tasks and data, so we give its using description. For other functions, use according to yourself's task and data.