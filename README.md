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
web_summary_list = ['500_PBMC_3p_LT_Chromium_X_web_summary.html']
```

### 1. metrics
#### The 'delivery doc' setting will provide an abbreviated metrics table from your sample(s) web_summary.html files
```
metrics.tableGenerator(web_summary_list, 'delivery doc')
```
![delivery_doc](https://user-images.githubusercontent.com/70353129/137335098-984b5f96-07e3-4bc9-8dca-c5bbccaed7c6.JPG)

#### The 'full' setting will provide a metrics table with many different attributes scrapped from your sample(s) web_summary.html files
```
metrics.tableGenerator(web_summary_list, 'full')
```
![full](https://user-images.githubusercontent.com/70353129/137335398-609ff8b8-84b0-48b8-ad46-adf8c09853f1.JPG)

#### The 'repooling' setting will allow the user to see how many reads are necessary to reach a desired read depth
```
metrics.tableGenerator(web_summary_list, 'repooling', readsDesired=120000)
```
![repooling](https://user-images.githubusercontent.com/70353129/137335634-f031c261-7b8c-4d84-847a-cb9df86a6ef5.JPG)

#### Here, we can see that approximately a billion reads are necessary for both samples to reach 120,000 reads per cell


### 2. webSum_10x_halfSat
#### We can look at the saturation curves for our samples and extrapolate sequencing saturation and median reads per cell for a given number of reads per cell
```
webSum_10x_halfSat.satcurves(web_summary_list[0], readsDesired=120000)
webSum_10x_halfSat.satcurves(web_summary_list[1], readsDesired=120000)
```
![10k_PBMC_3p_nextgem_Chromium_X](https://user-images.githubusercontent.com/70353129/137337258-425dab32-d4e9-47e2-af20-555758ed2663.png)\
Sequencing saturation half-saturation point: 27,304 reads\
Current sequencing saturation level: 61.1%\
Current reads per cell: 41,379\
Current genes per cell: 2,049

Desired reads per cell: 120,000\
Sequencing saturation for desired reads per cell: 81.5%\
Uniques genes per cell for desired reads per cell: 2,395

![500_PBMC_3p_LT_Chromium_Controller](https://user-images.githubusercontent.com/70353129/137337289-e5442a21-4552-4a89-acd4-b269dedefd95.png)\

Sequencing saturation half-saturation point: 35,879 reads\
Current sequencing saturation level: 73.8%\
Current reads per cell: 99,250\
Current genes per cell: 2,313

Desired reads per cell: 120,000\
Sequencing saturation for desired reads per cell: 77.0%\
Uniques genes per cell for desired reads per cell: 2,358
