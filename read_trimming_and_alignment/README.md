# Bicyclus anynana SRR read cleaning and alignment#

updated 6 Dec 2021

Steward, R.A., de Jong, M., Oostra, V., Wheat, C.W. 

***************************************

### Contents ###

+ A. Cleaning
+ B. STAR 2-pass alignment

***************************************

### A. Cleaning ###
We used an in house script `getlog_noclonefilter_Q20.py` to filter adapters and trim reads by quality (rl, PHRED>=20)

***************************************

### B. Alignment ###
We used the STAR aligner to generate a genome directory and map reads using default settings and a two-pass approach. 
