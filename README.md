Splinectomy: A small suite of statistical analysis tools
=======
Tools for messy longitudinal data analysis in R, using splines and non-parametric comparisons

## Installation
Clone the git repo, then add the bin directory to your path for easiest execution. Or, run the scripts from the bin folder by typing
```
[path_to_repo]/bin/script.R [arguments]
```
If the message `Permission denied` is returned in error, ensure that the scripts are executable as follows (from the `bin` directory):
```
chmod +x permusplinectomy.R
chmod +x sliding_spline_test.R
```
### Dependencies
These scripts have been developed and tested in R version 3.3.1. Currently, the following R packages are required (future versions will reduce the number of large dependencies):
```
# permusplinectomy.R
dplyr
optparse
# sliding_spline_test.R
dplyr
ggplot2
reshape2
cowplot
optparse
```
## How to use
### permusplinectomy.R
Given a longitudinal dataset, permusplinectomy reports a p-value for a binary categorical variable (e.g. Healthy vs Disease). It is robust to differing number and frequency of sampling between individuals ("units" in the code), and allows the user to alter key parameters controlling the stringency of the pre-test filtering. The input dataset should be a tab-delimited file in long format and have a column for the patient/individual ID, the independent and response variables (e.g. time and continuous measured variable), and the categorical variable. Before running, ensure there are no missing values in the measured variable or the categorical variable columns.

The p-value is calculated by permutation of the categorical label across the individuals, and non-parametrically measuring the likelihood of the true area between the categories being a result of random chance.
```
# Sample run command
permusplinectomy.R -i data_table.txt -x Years -y Blood_glucose -c Disease_status -p PATIENT_ID --perms 999
# To see usage and all commandline options
permusplinectomy.R --help
```
### sliding_spline_test.R
Given the same dataset, the sliding_spline_test treats each individual separately, effectively converting their time series into a dense time series through extrapolation of a spline. Each interval (default = 100) is then tested for non-parametric significance between the two groups, provided there are enough data points (default = 3+ per group). The script produces three output files: a plot (png) of the splines for each individual, a plot of the negative log of p-values over the independent (x) variable where the size of the points is scaled by the number of data points contributing to that test (thicker line = greater n at that x), and a table of p-value at each interval. The user may provide a file prefix for each of these files, which are saved in the current working directory.
Note: It is worth highlighting that in the p-value plot, the *negative log* is reported; thus, 0.05 = ~1.301, and lower p-values appear as larger positive values (greater than 1.3).
```
# Sample run command
sliding_spline_test.R -i data_table.txt -x Years -y Blood_glucose -c Disease_status -p PATIENT_ID --spline_intervals=200 --prefix=Foo_
# To see usage and all command line options
sliding_spline_test.R --help
```
A good practice set is the `ChickWeight` dataset in R's datasets package. To export this data for use in the `splinectomy` scripts run the following from the R console:
```
> write.table(ChickWeight, file = 'ChickWeight.txt', sep='\t', quote =F, row.names = F)
```
Try it out! A few examples to get you started:
```
# Compare diets 1 and 2 overall:

# Should return a non-significant pvalue. But the separation isn't consistent across the time series... Try a sliding spline:

# Look at the pvalues plot...
# Cool! So, this suggests that the diets may lead to significantly different 
# chick weights at early timepoints, but this difference is not maintained 
# through the end of the experiment.
```
