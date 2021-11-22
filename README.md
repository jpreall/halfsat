## Overview
Creates halfsat curves as reference for additional sequencing. Scrapes web_summary.html
to obtain reference information for different tables. The package contains two main files,
'webSum_10x_halfSat' and 'metrics'

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
web_summary_arc = 'human_brain_3k_web_summary.html'
```

### 1. metrics
#### The 'delivery doc' setting will provide an abbreviated metrics table from your sample(s) web_summary.html files
```
metrics.tableGenerator(web_summary_list, 'GEX', 'delivery doc')
```
![image](https://user-images.githubusercontent.com/70353129/142926083-09bea16e-5a9c-4aa7-b12f-5d769c1df19c.png)

#### The 'full' setting will provide a metrics table with many different attributes scrapped from your sample(s) web_summary.html files
```
metrics.tableGenerator(web_summary_list, 'GEX', 'full')
```
![image](https://user-images.githubusercontent.com/70353129/142926261-ad602faa-4f35-4da6-805e-40083d4cbd1c.png)

#### The 'full' setting can also be obtained for ARC web_summary.html files
```
atacTable, gexTable = metrics.tableGenerator([web_summary_arc], webSummaryType='ARC', tableType='full')
```
![image](https://user-images.githubusercontent.com/70353129/142925906-f411d27c-2314-48d5-957d-ded14bfe17e1.png)
![image](https://user-images.githubusercontent.com/70353129/142926309-0358c90d-b130-4e2f-8651-3805a9ce50f2.png)

#### The 'repooling' setting will allow the user to see how many reads are necessary to reach a desired read depth
```
metrics.tableGenerator(web_summary_list, 'repooling', readsDesired=120000)
```
![image](https://user-images.githubusercontent.com/70353129/142926398-7af4d07b-7a12-43d0-aad5-0fce4af13881.png)

#### Here, we can see that approximately a billion reads are necessary for both samples to reach 120,000 reads per cell


### 2. webSum_10x_halfSat
#### We can look at the saturation curves for our samples and extrapolate sequencing saturation and median reads per cell for a given number of reads per cell
```
webSum_10x_halfSat.satcurves(web_summary_list[0], readmax=100000, readsDesired=80000)
webSum_10x_halfSat.satcurves(web_summary_list[1], readmax=100000, readsDesired=80000)
```
![image](https://user-images.githubusercontent.com/70353129/142926473-1d2f0392-8b68-4aed-b36e-889fb66fe9dd.png)

Sequencing saturation half-saturation point: 40,018 reads per cell\
Current sequencing saturation level: 50.3%\
Current reads per cell: 39,865\
Current genes per cell: 3,358

Desired reads per cell: 80,000\
Sequencing saturation for desired reads per cell: 66.7%\
Uniques genes per cell for desired reads per cell: 3,842

![image](https://user-images.githubusercontent.com/70353129/142926774-d574de19-2439-4e5b-8ccd-725b7f24f9c0.png)

Sequencing saturation half-saturation point: 43,152 reads per cell
Current sequencing saturation level: 58.5%\
Current reads per cell: 59,264\
Current genes per cell: 3,764

Desired reads per cell: 80,000\
Sequencing saturation for desired reads per cell: 65.0%\
Uniques genes per cell for desired reads per cell: 3,965

```
webSum_10x_halfSat.satcurves(web_summary_arc, webSummaryType='ARC', readmax=180000, readsDesired=100000)
```
![image](https://user-images.githubusercontent.com/70353129/142926999-be740fe8-52b5-4261-af06-10d5a2c425a2.png)

Sequencing saturation half-saturation point: 18,409 reads per cell\
Current sequencing saturation level: 86.1%\
Current reads per cell: 122,334.93\
Current genes per cell: 2,600

Desired reads per cell: 100,000\
Sequencing saturation for desired reads per cell: 84.5%\
Uniques genes per cell for desired reads per cell: 2,526
Fragments per cell for desired reads per cell: 23,489
