## Overview
Creates halfsat curves as reference for additional sequencing. Scrapes web_summary.html
to obtain reference information for different tables. The package contains two main files,
'webSum_10x_halfSat' and 'metrics'.

**webSum_10x_halfSat** will take in a web_summary.html file and output half-saturation plots
and relevant information such as sequencing saturation half-sat point, current sequencing
saturation level, current reads per cell, and current genes per cell. Additionally, the
user can input the desired reads per cell for the samples and see the corresponding
sequencing saturation and unique genes per cell based on the curves generated.

**metrics** will scrape a list of web_summary files for all relevant information depending on
the table-type provided: full, delivery doc, or repooling. Repooling allows for the
additional parameter of 'readsDesired', which is used in the repooling calculation.


## Installation
### 1. Create a new directory and clone repository
```
(your_env) $ git clone https://github.com/jpreall/halfsat.git
(your_env) $ pip install -e halfsat
```

## Example
### Import and load file path
```
from halfsat import webSum_10x_halfSat, metrics
web_summary_list = ['Brain_3p_web_summary.html',
                    'Brain_3p_LT_web_summary.html',
                    'Breast_Cancer_3p_LT_web_summary.html']
```

### 1. webSum_10x_halfSat
```
webSum_10x_halfSat.satcurves(web_summary_list[0], readsDesired=100000)
webSum_10x_halfSat.satcurves(web_summary_list[1], readsDesired=100000)
```

### 2. metrics
```
metrics.tableGenerator(web_summary_list, 'full')
```
