# GeCalibration

Lists information about about the High Voltage and Amplifier and other settings

Reads in a text file of the ASCII data from Maestro and bins the data to a TH1D histogram

Searches the histogram for gamma peaks.

Peak locations are modified by fitting a gaussian + polynomial fit, which gives a new mean and standard deviation.

Linear fit applied to peak locations, which converts ADC channels to energies.

Use ConvertMaestroToHistogram to convert from Maestro txt file to a Histogram. Inputs are file location on machine and name of the Histogram.

Use CalibrationFunction to find the linear fit of the ADC channels to energies. Takes an input of the file location and the name of the function.

Use the function that was found with CalibrationFunction as an input for energyCalibratedSpectrum as well as the ADC histogram. Will output a calibrated spectrum as a histogram.
