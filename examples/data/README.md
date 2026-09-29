# Example data

`climate.csv` is daily rainfall and potential evaporation, in metres per day,
from 31 January 1990 to 5 January 2000: 3,627 days.

It is the forcing of the LUMPREP example shipped with LUMPREM
(`lumprem/lumprep_example/rain.dat` and `epot.dat`), relabelled from day
numbers to dates using that example's `START_DATE` and written as CSV. The
series come from SILO patched-point data for station 39104, Monto Township,
Queensland (`lumprem/lumprep_example/39104pp.txt`), extracted in November
2017. That file states the data are copyright to the Queensland Government
(then DSITI) and the Bureau of Meteorology, and that the supply-to-licensee
restriction on SILO data does not apply to Queensland patched-point data.
Evaporation is read at 9 am and has been shifted to the day before.
