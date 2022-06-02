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

**predictUMIs** will take in a metrics_summary_json.json file and output UMIs/cell predictions
and half-saturation based a michaelis menten curve fit. 


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
webSum_10x_halfSat.satcurves(web_summary_list[0], readMax=100000, readsDesired=80000)
```
![image](https://user-images.githubusercontent.com/70353129/142926473-1d2f0392-8b68-4aed-b36e-889fb66fe9dd.png)

Sequencing saturation half-saturation point: 40,018 reads per cell
Median genes per cell half-saturation point: 2,260 genes per cell ; 14,091 reads per cell
Current sequencing saturation level: 50.3%
Current reads per cell: 39,865
Current genes per cell: 3,358

Desired reads per cell: 80,000
Sequencing saturation for desired reads per cell: 66.7%
Median genes per cell for desired reads per cell: 3,842

<img width="903" alt="Screen Shot 2021-12-30 at 9 48 32 PM" src="https://user-images.githubusercontent.com/70353129/147800243-53fc98f1-1c68-4f96-baac-a68438df2845.png">

```
webSum_10x_halfSat.satcurves(web_summary_arc, webSummaryType='ARC', readMax=180000, readsDesired=100000, readPairsDesired=150000)
```
![image](https://user-images.githubusercontent.com/70353129/142926999-be740fe8-52b5-4261-af06-10d5a2c425a2.png)

Sequencing saturation half-saturation point: 18,409 reads per cell
Median genes per cell half-saturation point: 1,378 genes per cell ; 9,163 reads per cell
Median fragments per cell half-saturation point: 26,618 fragments per cell ; 126,635 read pairs per cell
Current sequencing saturation level: 86.1%
Current reads per cell: 122,334.93
Current genes per cell: 2,600
Current median fragments per cell: 22,881
Current mean raw read pairs per cell 95,897.34

Desired reads per cell: 100,000
Sequencing saturation for desired reads per cell: 84.5%
Median genes per cell for desired reads per cell: 2,526
Desired read pairs per cell: 150,000
Fragments per cell for desired reads per cell: 28,866

<img width="900" alt="Screen Shot 2021-12-30 at 9 56 05 PM" src="https://user-images.githubusercontent.com/70353129/147800450-4cf20f95-cc17-4738-bb6e-0a26a802447f.png">

### 3. predictUMIs
#### Plot UMIs versus reads graph and provide dataframe with prediction of UMIs per cell given a specified number of reads.
```
# metrics_summary_json was generated from 10x's public 500_PBMC_3p_LT_Chromium_X dataset
predictUMIs.plotUMIcurve(metrics_summary_json, readMax=175000, readsDesired=150000)
```
![image](https://user-images.githubusercontent.com/70353129/147623094-7dd5396b-8ac5-4257-952d-b837cd28b7ce.png)
![image](https://user-images.githubusercontent.com/70353129/147623125-53742b65-c75b-401a-a3e4-88ba1f608e84.png)

## Mathematical Background
### Halfsat using lander-waterman modeling

From Picard tools, we have the Lander-Waterman equation as follows:

$$
\begin{equation} \frac{C}{X} = 1 - e^\frac{-N}{X} \end{equation} \qquad (1)
$$

where

$X = number \enspace of \enspace distinct \enspace molecules \enspace in \enspace library$

$N = number \enspace  of \enspace  read \enspace  pairs$

$C = number\enspace of\enspace distinct\enspace fragments\enspace observed\enspace in\enspace read\enspace pairs$

From 10x genomics’ webpage, sequencing saturation is calculated as 

$s = 1 - \frac{n\textunderscore deduped\textunderscore reads}{n\textunderscore reads}$ where $s$ is the sequencing saturation, $n\textunderscore deduped\textunderscore reads$ is the number of unique (valid cell-barcode, valid-UMI, gene) combinations among confidently mapped reads, and $n\textunderscore reads$ is the total number of confidently mapped, valid cell-barcode, valid-UMI reads.

Using notation from the Lander-Waterman equation, we can rewrite the sequencing saturation equation as follows:

$$
\begin{equation}
s = 1 - \frac{C}{N} \qquad (2)
\end{equation}
$$

Our goal is to fit a model to predict sequencing saturation $s$ as a function of read pairs $N$ and number of distinct molecules in library $X$, optimizing for $X$.

In other words, $s = f(N)$. 

We can rewrite equation (2) as $C = (1-s)*N$ and use systems of equations to solve for $s$.

$$
\frac{(1-s)*N}{X} = 1-e^\frac{-N}{X}
$$

$$
\begin{equation}
s = 1-\frac{(1-e^\frac{-N}{X})*X}{N} \qquad (3)
\end{equation}
$$

where $s$ is a function of read pairs $N$.

From here, we can create a plot of sequencing saturation versus mean reads per cell (total reads divided by the number of cells). In order to find the half saturation point, we must use a root-finding algorithm such as Brent’s method to obtain a number $N$ such that $f(N) = 0.5$

$$
\begin{equation}
s = 1-\frac{(1-e^\frac{-N}{X})*X}{N}-0.5 \qquad (4)
\end{equation}
$$

 <br/><br/> 
### Halfsat using michaelis-menten modeling
The Michaelis-Menten equation is as follows:

$$
\begin{equation}
v = \frac{\mathrm{d}[P]}{\mathrm{d}t} = V_{max}\frac{[S]}{K_M+[S]} \qquad (5)
\end{equation}
$$

where

$v = rate\enspace of\enspace product\enspace formation$

$P = concentration\enspace of\enspace product$

$S = concentration\enspace of\enspace substrate$

$V_{max} = max\enspace rate\enspace achieved\enspace by\enspace system$

$K_M = michaelis\enspace constant$

When $K_M$ is  equivalent to the substrate concentration the reaction rate $v$ is half of $V_{max}$

Analogously, we can model the following:

$$
\begin{equation}
s = S_{max}\frac{R}{K_M+R} \qquad (6)
\end{equation}
$$

where

$s = sequencing\enspace saturation\enspace [0,1]$

$S_{max} = maximum\:sequencing\:saturation$

$R = reads\enspace per\enspace cell$

$K_M = michaelis\enspace constant$

Since $S_{max}=1$, we can simplify equation (5)

$$
\begin{equation}
s = \frac{R}{K_M+R} \qquad (7)
\end{equation}
$$

We can also use the Michaelis-Menten equation (5) to model genes per cell and UMIs per cell.

$s = median\enspace genes\enspace per\enspace cell$

$s = median\enspace UMIs\enspace per\enspace cell$
