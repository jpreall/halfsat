from bs4 import BeautifulSoup
import re
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
import os

def scrape_saturation_stats(web_summary_html_file):
    """
    Scrape saturation stats from web_summary.html produced by 10x Cellranger.
    Only works with Cellranger version 3.x and up.
    Args:
        web_summary_html_file (string): web_summary.html file to scrape.

    Returns:
        list of int: The number of reads for each point of seq_saturation_plot.
        list of int: The number of genes for each point of median_gene_plot.
        list of int: The sequencing saturation for each point of seq_saturation_plot.
        string: The sequencing saturation percent for that sample.
        string: The sample name.

    Raises:
        None.

    """

    f = open(web_summary_html_file, encoding="utf8")
    soup = BeautifulSoup(f, "html.parser")
    for line in soup.find('script'):
        if 'const data' in line:
            const_data = line
    f.close()

    constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n')+1])

    seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
    current_sat = [name for name in seq_summary_table if 'Sequencing Saturation' in name][0][1]
    sampname = constant_data['summary']['sample']['id']

    reads = np.array(constant_data['summary']['analysis_tab']['seq_saturation_plot']['plot']['data'][0]['x'])
    genes = np.array(constant_data['summary']['analysis_tab']['median_gene_plot']['plot']['data'][0]['y'])
    saturations = np.array(constant_data['summary']['analysis_tab']['seq_saturation_plot']['plot']['data'][0]['y'])

    cells_table = constant_data['summary']['summary_tab']['cells']['table']['rows']
    for entry in cells_table:
        if entry[0] == 'Median Genes per Cell':
            current_median_genes = entry[1]
        if entry[0] == 'Mean Reads per Cell':
            current_mean_reads = entry[1]


    return reads, genes, saturations, current_sat, sampname, current_median_genes, current_mean_reads


def satcurves(web_summary_html_file, readmax=250000, title=None, readsDesired=40000):
    """
    Plot saturation curves from read/gene/sat data scraped from web_summary.html file.

    Args:
        web_summary_html_file (string): web_summary.html file to scrape.
        readmax (int): The end of the plot.
        title (string): The title  of the plot.
        readsDesired (int): Mean reads/cell desired.

    Returns:
        None.

    Raises:
        None.

    """

    pl.rcParams['figure.dpi'] = 120

    reads, genes, saturations, current_sat, sampname, current_median_genes, current_mean_reads = scrape_saturation_stats(web_summary_html_file)
    #a is the "kd"
    #b is the max genes per cell, or fraction saturation
    # here b=1

    ## initiate the figure
    fig, ax = pl.subplots(1,2, figsize=[6.75,2.75])

    ## Fit the curves
    ## Step 1: Saturation curve.  b = 1
    def f(x, a):
        return(x /(x+a))

    popt, pcov = curve_fit(f, reads, saturations)
    halfsat_sat = np.round(popt[0],0).astype('int')
    #ymax = np.round(popt[1],0).astype('int')
    ymax_sat = 1

    assert readsDesired > 0, "Cannot have negative mean reads/cell"
    desiredSeqSat = np.round(f(readsDesired, halfsat_sat),3)

    ## Plot
    ax[0].scatter(reads, saturations)
    ax[0].plot(reads, f(reads, *popt), 'r-')

    # extrapolate with the curve fit
    x_more = np.linspace(0,readmax, 100)
    xmax = np.max(x_more)
    ax[0].plot(x_more, f(x_more, *popt), 'k-')

    #pl.text(0.1*xmax,ymax*0.9,str(ymax) + ' genes', size=8)

    # set plot boundaries and add asymptotes and lines
    color_sat='purple'
    ax[0].set_ylim(0,1.1)
    ax[0].axhline(ymax_sat, linestyle='--', color=color_sat)
    ax[0].vlines(x= halfsat_sat, ymin=0, ymax=f(halfsat_sat, halfsat_sat), linestyle=':', color=color_sat)
    ax[0].hlines(y=f(halfsat_sat, halfsat_sat), xmin=0, xmax = halfsat_sat, linestyle=':',color=color_sat)
    ax[0].text(halfsat_sat + xmax/20,0.4*f(halfsat_sat, halfsat_sat),'half-saturation point: \n' + format(halfsat_sat,',') + ' reads/cell', size=8)
    #ax[0].text(halfsat_sat, f(halfsat_sat, halfsat_sat),'half-saturation point:' + str(halfsat_sat) + ' reads', size=6)
    #ax[0].text(xmax*0.065,ymax_sat*0.05,'current saturation:' + str(current_sat) + ', ' + str(current_mean_reads) +' reads/cell', size=7)

    #pl.show()

    #label the axes
    ax[0].set_xlabel('Reads per cell')
    ax[0].set_ylabel('Sequencing Saturations')

    #print('popt:',popt)
    #print('pcov:',pcov)

    print()
    #print('max genes',ymax_sat,'genes')
    print('half-saturation point:',halfsat_sat,'reads')
    print('current saturation level:', current_sat)

    ## Gene saturation curve.
    def f(x, a, b):
        return(b * x /(x+a))

    popt, pcov = curve_fit(f, reads, genes, p0=[22000,1000])
    halfsat_genes = np.round(popt[0],0).astype('int')
    ymax_genes = np.round(popt[1],0).astype('int')
    halfsat_genes_counts = np.round(f(halfsat_genes, halfsat_genes, ymax_genes)).astype(int)

    desiredUniqueGenes = np.round(f(readsDesired, halfsat_genes, ymax_genes)).astype(int)
    #ymax_genes = 1

    ## Plot
    ax[1].scatter(reads, genes, color='red')
    ax[1].plot(reads, f(reads, *popt), 'r-')

    # extrapolate with the curve fit
    x_more = np.linspace(0,readmax, 100)
    xmax = np.max(x_more)
    ax[1].plot(x_more, f(x_more, *popt), 'k-')

    # set plot boundaries and add asymptotes and lines
    color_genes='r'
    ax[1].set_ylim(0,ymax_genes*1.1)
    ax[1].axhline(ymax_genes, linestyle='--', color=color_genes)
    ax[1].vlines(x= halfsat_genes, ymin=0, ymax=f(halfsat_genes, halfsat_genes, ymax_genes), linestyle=':', color=color_genes)
    ax[1].hlines(y=f(halfsat_genes, halfsat_genes, ymax_genes), xmin=0, xmax = halfsat_genes, linestyle=':', color=color_genes)
    ax[1].text(halfsat_genes + xmax/20,0.65*f(halfsat_genes, halfsat_genes, ymax_genes),'half-saturation point: \n' +
               format(halfsat_genes,',') + ' reads/cell, ' + format(halfsat_genes_counts, ',') + ' genes/cell', size=8)
    ax[1].text(0.1*xmax,ymax_genes*1.01,format(ymax_genes,',') + ' genes max', size=8)
    #ax[1].text(xmax*0.08,ymax_genes*0.05,'current saturation: \n' + str(current_mean_reads) + ' reads/cell, ' +
    #           str(current_median_genes) + ' genes/cell', size=7)

    #label the axes
    ax[1].set_xlabel('Reads per cell')
    ax[1].set_ylabel('Unique Genes Detected')

    ## figure title
    if not title:
        if isinstance(sampname, str):
            title = sampname
    fig.suptitle(title, y=1, size=10)

    pl.subplots_adjust(wspace=0.5)
    pl.tight_layout()
    pl.show()
    #print('popt:',popt)
    #print('pcov:',pcov)

    print()
    #print('max genes',ymax_sat,'genes')
    print('Sequencing saturation half-saturation point:',format(halfsat_sat, ','),'reads')
    print('Current sequencing saturation level:', current_sat)
    print('Current reads per cell:', current_mean_reads)
    print('Current genes per cell:', current_median_genes)
    print()
    print('Desired reads per cell:', readsDesired)
    print('Sequencing saturation for desired reads per cell:', desiredSeqSat)
    print('Uniques genes per cell for desired reads per cell:', desiredUniqueGenes)

def find_satcurves(folder):
    """
    Walk through a folder of Cellranger outputs to find a bunch of web_summary files and run them all.

    Args:
        folder (string): Path to folder containing web_summary files.

    Returns:
        None.

    Raises:
        None.

    """
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(".html"):
                file_path = os.path.join(root, file)
                path_term_list = file_path.split('/')
                sample_name_index_in_path_term_list = path_term_list.index('outs') - 1
                title = path_term_list[sample_name_index_in_path_term_list]
                print(file_path)
                #reads, genes, saturations, current_sat, sampname = scrape_saturation_stats(file_path)
                satcurves(file_path)
                #satcurves(reads, genes, saturations)
