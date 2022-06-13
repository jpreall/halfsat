## Overview
Creates halfsat curves as reference for additional sequencing. Scrapes web_summary.html
to obtain reference information for different tables. The package contains two main files,
'metrics', 'predictUMIs', and 'webSum_10x_halfsat'

**metrics** will scrape a list of web_summary files for all relevant information depending on
the table-type provided: full, delivery doc, or repooling. Repooling allows for the
additional parameter of 'readsDesired', which is used in the repooling calculation.

**predictUMIs** will take in a metrics_summary_json.json file and output UMIs/cell predictions
and half-saturation based a michaelis menten curve fit. 

**webSum_10x_halfSat** will take in a web_summary.html file and output half-saturation plots
and relevant information such as sequencing saturation half-sat point, current sequencing
saturation level, current reads per cell, and current genes per cell. Additionally, the
user can input the desired reads per cell for the samples and see the corresponding
sequencing saturation and unique genes per cell based on the curves generated.

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


### 2. predictUMIs
#### Build a UMI_model object with various attributes scraped from the 10x json file
```
# metrics_summary_json was generated from 10x's public 500_PBMC_3p_LT_Chromium_X dataset
# build UMI_model first
my_UMI_model = predictUMIs.UMI_model(metrics_json)
# check out the model attributes
my_UMI_model.current_reads_per_cell, my_UMI_model.current_UMIs
```
Out: (128910.867120954, 8935.0)

#### Scrape reads and UMI information from the json file exclusively
```
my_UMI_model_reads, my_UMI_model_UMIs = predictUMIs.get_reads_and_UMIs_from_json(metrics_json)
print('reads: ', my_UMI_model_reads) 
print('UMIs: ', my_UMI_model_UMIs)
```
reads:  [5000.  10000.  12891.  20000.  25782.  30000.  38673.  50000.  51564. 64455.  77346.  90237. 103128. 116019. 128910.]

UMIs:  [1596.0, 2827.0, 3380.0, 4551.0, 5253.0, 5705.0, 6388.0, 7049.0, 7128.0, 7638.0, 8021.0, 8287.0, 8552.0, 8731.0, 8935.0]


#### fit UMI model with read and UMI data from json
```
my_UMI_model.fit_UMIs()
my_UMI_model.make_UMI_table()
```
<img width="881" alt="image" src="https://user-images.githubusercontent.com/70353129/173153046-84e8bcb9-5dd4-4844-a9e1-36eb965694eb.png">

```
my_UMI_model.plot_UMIs(readMax=200000)
```
![image](https://user-images.githubusercontent.com/70353129/173153262-8af7dc37-da84-4a5e-909d-b4e5d41a37ee.png)

#### Plot new data on a model's plot
```
my_UMI_model2.plot_UMIs(reads_test=test_array_reads, UMIs_test=test_array_saturations)
```
<img width="401" alt="image" src="https://user-images.githubusercontent.com/70353129/173390201-deb4f134-f19d-466f-b922-52851a40cead.png">

#### Predict UMI values based on model and determine the goodness of fit
```
my_UMI_model.predict(test_array_reads)
my_UMI_model.score(test_array_reads, test_array_saturations)
```

### 3. webSum_10x_halfSat
#### Build web summary model
```
Brain_3p = webSum_10x_halfSat.webSum_model(web_summary_list[0])
```

#### Scrape reads, saturations, and genes information from web summary test file exclusively
```
reads_test, sat_test, genes_test = webSum_10x_halfSat.get_reads_sats_genes_from_web(web_summary[1])
```

#### Fit model and plot train and test data on it
```
Brain_3p.fit_model(seqSatModel='mm')
Brain_3p.plot(readMax=250000, reads_test=reads_test, saturations_test=sat_test, genes_test=genes_test)
```
<img width="725" alt="image" src="https://user-images.githubusercontent.com/70353129/173396430-443b5631-2ab3-4ce9-9a21-b78d5c0d378a.png">

#### Use models for prediction and test the goodness of fit
```
Brain_3p.predict_genes(genes_test)
Brain_3p.score_genes(reads_test, genes_test)
```
Out: [350, 389, 524, 570, 682, 712, 756, 812, 816, 856, 888, 913, 923, 934, 953]

Out: RSS: 50334.5  ymean: 2794.3333333333335  TSS: 10176173.833333332  Rsquared: 0.9950536910213618


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
