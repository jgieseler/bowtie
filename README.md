# The bowtie analysis tool
---

This analysis tool runs bowtie analysis for particle instruments' energy channels.
The input is a csv table of channel responses indexed by the incident energy.
The results of analysis are the geometric factor (with errors) and the
effective energy of the channel.

The package operates with two main classes, which are called `Bowtie` and `Spectra`. `Bowtie` stores 
response functions and contains the methods to run bowtie analysis, while `Spectra` contains the 
information of the spectral indices and the amount of different spectra that are used in the bowtie calculation.

Bowtie
---
The `Bowtie` class contains the data that the bowtie analysis is applied on, and the energy range to be
considered in the calculations. Its methods make running analysis easy and straightforward.

Methods:
#
	set_energy_range(energy_min, energy_max):
 		energy_min : {float} The minimum energy in MeV to consider
   		energy_max : {float} See energy_min.
#
 	bowtie_analysis(channel, spectra, plot):
  		channel : {str} The channel name as it appears in the csv table.
		spectra : {Spectra} The Spectra class object, introduced in this package. Contains the 
  							spectral indices and the power law spectra used in the bowtie analysis.
  		plot : {bool} A boolean switch to produce a plot visualizing the analysis.
#
	bowtie_analysis_full_stack(spectra, plot):
 		A wrapper for bowtie_analysis(). Runs the analysis on all channels that appear in the input file.
 		spectra : {Spectra} See bowtie_analysis().
   		plot : {bool} See bowtie_analysis()

Spectra
---
The `Spectra` class contains the range of spectra that are applied on the response function to run 
bowtie analysis.

Methods:
#
	set_spectral_indices(gamma_min, gamma_max):
 		gamma_min : {float} The minimum spectral index to consider in the calculation.
   		gamma_max : {float} See gamma_min.
#
 	produce_power_law_spectra(response_df):
  		response_df : {pandas.DataFrame} The input csv table read in to a pandas DataFrame. Contains
										 the channel responses as a function of incident energy.
