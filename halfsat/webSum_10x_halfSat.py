from bs4 import BeautifulSoup
import re
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import pandas as pd
import os

__all__ = ['satcurves', 'find_satcurves']


def __scrape_saturation_stats(web_summary_html_file, webSummaryType='GEX'):
    """
    Scrape saturation stats from web_summary.html produced by 10x Cellranger.
    Only works with Cellranger version 3.x and up.
    Args:
        web_summary_html_file (string): web_summary.html file to scrape.
        webSummaryType (string): either 'GEX' or 'ARC' based on Cellranger pipeline used

    Returns:
        list of int: The number of reads for each point of seq_saturation_plot.
        list of int: The number of genes for each point of median_gene_plot.
        list of int: The sequencing saturation for each point of seq_saturation_plot.
        string: The sequencing saturation percent for that sample.
        string: The sample name.

    Raises:
        KeyError: when attributes in specified web sumamry html file does not match
                  the supplied format.
        ValueError: when neither 'GEX' nor 'ARC' is supplied.

    """

    f = open(web_summary_html_file, encoding="utf8")
    soup = BeautifulSoup(f, "html.parser")
    for line in soup.find('script'):
        if 'const data' in line:
            const_data = line
    f.close()

    constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n')+1])

    if webSummaryType == 'GEX':
        try:
            seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
            current_sat = [
                name for name in seq_summary_table if 'Sequencing Saturation' in name][0][1]
            sampname = constant_data['summary']['sample']['id']
        except KeyError:
            print('check if GEX web summary html file because attributes not found')
            raise

        reads = np.array(constant_data['summary']['analysis_tab']
                         ['seq_saturation_plot']['plot']['data'][0]['x'])
        genes = np.array(constant_data['summary']['analysis_tab']
                         ['median_gene_plot']['plot']['data'][0]['y'])
        saturations = np.array(constant_data['summary']['analysis_tab']
                               ['seq_saturation_plot']['plot']['data'][0]['y'])

        cells_table = constant_data['summary']['summary_tab']['cells']['table']['rows']
        for entry in cells_table:
            if entry[0] == 'Median Genes per Cell':
                current_median_genes = entry[1]
            if entry[0] == 'Mean Reads per Cell':
                current_mean_reads = entry[1]
        return reads, genes, saturations, current_sat, sampname, current_median_genes, current_mean_reads

    elif webSummaryType == 'ARC':
        try:
            sampleName = constant_data['sample']['id']
        except KeyError:
            print('check if ARC web summary html file because attributes not found')
            raise

        meanReadPairsArray = constant_data['atac_seq_saturation_plot']['data'][0]['x']
        medianUniqFragArray = constant_data['atac_seq_saturation_plot']['data'][0]['y']
        gex_percentDupArray = constant_data['gex_seq_saturation_plot']['data'][0]['y']
        gex_meanReadsArray = constant_data['gex_seq_saturation_plot']['data'][0]['x']
        gex_medianGenesArray = constant_data['gex_genes_per_cell_plot']['data'][0]['y']

        gex_cells_table = constant_data['gex_cells_table']['rows']
        for entry in gex_cells_table:
            if entry[0] == 'Mean raw reads per cell':
                current_mean_raw_reads_gex = entry[1]
            if entry[0] == 'Median genes per cell':
                current_median_genes = entry[1]

        gex_sequencing_table = constant_data['gex_sequencing_table']['rows']
        for entry in gex_sequencing_table:
            if entry[0] == 'Percent duplicates':
                gex_percentDup = entry[1]

        atac_cells_table = constant_data['atac_cells_table']['rows']
        for entry in atac_cells_table:
            if entry[0] == 'Mean raw read pairs per cell':
                current_mean_raw_reads_atac = entry[1]
            if entry[0] == 'Median high-quality fragments per cell':
                current_median_frags = entry[1]

        atac_sequencing_table = constant_data['atac_sequencing_table']['rows']
        for entry in atac_sequencing_table:
            if entry[0] == 'Percent duplicates':
                atac_percentDup = entry[1]

        arc_dictionary = {'atac_dictionary': {'sample name': sampleName,
                          'percent duplicates': atac_percentDup,
                                              'mean read pairs array': meanReadPairsArray,
                                              'median unique fragments array': medianUniqFragArray,
                                              'current mean raw reads': current_mean_raw_reads_atac,
                                              'current median fragments': current_median_frags},
                          'gex_dictionary': {'sample name': sampleName,
                                             'percent duplicates': gex_percentDup,
                                             'percent duplicates array': gex_percentDupArray,
                                             'mean reads array': gex_meanReadsArray,
                                             'median genes array': gex_medianGenesArray,
                                             'current mean raw reads': current_mean_raw_reads_gex,
                                             'current median genes': current_median_genes}
                          }

        return arc_dictionary

    else:
        raise ValueError('webSummaryType neither GEX nor ARC')


def satcurves(web_summary_html_file, webSummaryType='GEX', readMax=150000, title=None, readsDesired=40000, readPairsDesired=40000, verbose=True):
    """
    Plot saturation curves from data scraped from web_summary.html file.
    Only works with Cellranger version 3.x and up.

    Args:
        web_summary_html_file (string): web_summary.html file to scrape.
        webSummaryType (string): either 'GEX' or 'ARC' based on Cellranger pipeline used
        readMax (int): The limit of the plot.
        title (string): The title  of the plot.
        readsDesired (int): Mean reads/cell desired.
        readPairsDesired (int): Mean read pairs/cell desired (for atac portion).
        verbose (boolean): True to print information pertaining to plots.

    Returns:
        df (pandas DataFrame): DataFrame containing sample information and michaelis
        menten fit predictions.

    Raises:
        ValueError: when neither 'GEX' nor 'ARC' is supplied.

    """

    pl.rcParams['figure.dpi'] = 120

    if webSummaryType == 'GEX':
        reads, genes, saturations, current_sat, sampname, current_median_genes, current_mean_reads = __scrape_saturation_stats(
            web_summary_html_file, webSummaryType='GEX')
        fig, ax = pl.subplots(1, 2, figsize=[6.75, 2.75])

    elif webSummaryType == 'ARC':
        arc_dictionary = __scrape_saturation_stats(
            web_summary_html_file, webSummaryType='ARC')
        reads = np.array(arc_dictionary['gex_dictionary']['mean reads array'])
        genes = np.array(arc_dictionary['gex_dictionary']['median genes array'])
        saturations = np.array(arc_dictionary['gex_dictionary']['percent duplicates array'])
        current_sat = arc_dictionary['gex_dictionary']['percent duplicates']
        sampname = arc_dictionary['gex_dictionary']['sample name']
        current_median_genes = arc_dictionary['gex_dictionary']['current median genes']
        current_mean_reads = arc_dictionary['gex_dictionary']['current mean raw reads']

        atac_reads = np.array(arc_dictionary['atac_dictionary']['mean read pairs array'])
        fragments = np.array(arc_dictionary['atac_dictionary']['median unique fragments array'])
        current_mean_atac_reads = arc_dictionary['atac_dictionary']['current mean raw reads']
        current_median_frag = arc_dictionary['atac_dictionary']['current median fragments']
        fig, ax = pl.subplots(1, 3, figsize=[15, 4])

    else:
        raise ValueError('webSummaryType neither GEX nor ARC')

    # constraining to only 3 bins in x-axis for visiblity purposes
    ax[0].locator_params(axis='x', nbins=4)
    ax[1].locator_params(axis='x', nbins=4)
    if webSummaryType == 'ARC':
        ax[2].locator_params(axis='x', nbins=4)

    # Fit the curves
    # Step 1: Saturation curve.  b = 1

    def f(x, a):
        return(x / (x+a))

    popt, pcov = curve_fit(f, reads, saturations)
    halfsat_sat = np.round(popt[0], 0).astype('int')
    # ymax = np.round(popt[1],0).astype('int')
    ymax_sat = 1

    assert readsDesired > 0, "Cannot have negative mean reads/cell"
    desiredSeqSat = np.round(f(readsDesired, halfsat_sat), 3)

    # Plot
    ax[0].scatter(reads, saturations)
    ax[0].plot(reads, f(reads, *popt), 'r-')

    # extrapolate with the curve fit
    x_more = np.linspace(0, readMax, 100)
    xmax = np.max(x_more)
    ax[0].plot(x_more, f(x_more, *popt), 'k-')

    # pl.text(0.1*xmax,ymax*0.9,str(ymax) + ' genes', size=8)

    # set plot boundaries and add asymptotes and lines
    color_sat = 'purple'
    ax[0].set_ylim(0, 1.1)
    ax[0].axhline(ymax_sat, linestyle='--', color=color_sat)
    ax[0].vlines(x=halfsat_sat, ymin=0, ymax=f(
        halfsat_sat, halfsat_sat), linestyle=':', color=color_sat)
    ax[0].hlines(y=f(halfsat_sat, halfsat_sat), xmin=0,
                 xmax=halfsat_sat, linestyle=':', color=color_sat)
    ax[0].text(halfsat_sat + xmax/20, 0.4*f(halfsat_sat, halfsat_sat),
               'half-saturation point: \n' + format(halfsat_sat, ',') + ' reads/cell', size=8)

    # pl.show()

    # label the axes
    ax[0].set_xlabel('Reads per cell')
    ax[0].set_ylabel('Sequencing Saturations')

    # print('popt:',popt)
    # print('pcov:',pcov)

    # Gene saturation curve.
    # a is the "kd"
    # b is the max genes per cell, or fraction saturation
    # here b=1

    def f(x, a, b):
        return(b * x / (x+a))

    popt, pcov = curve_fit(f, reads, genes, p0=[22000, 1000])
    halfsat_genes = np.round(popt[0], 0).astype('int')
    ymax_genes = np.round(popt[1], 0).astype('int')
    halfsat_genes_counts = np.round(f(halfsat_genes, halfsat_genes, ymax_genes)).astype(int)

    desiredUniqueGenes = np.round(f(readsDesired, halfsat_genes, ymax_genes)).astype(int)
    # ymax_genes = 1

    # Plot
    ax[1].scatter(reads, genes, color='red')
    ax[1].plot(reads, f(reads, *popt), 'r-')

    # extrapolate with the curve fit
    x_more = np.linspace(0, readMax, 100)
    xmax = np.max(x_more)
    ax[1].plot(x_more, f(x_more, *popt), 'k-')

    # set plot boundaries and add asymptotes and lines
    color_genes = 'r'
    ax[1].set_ylim(0, ymax_genes*1.1)
    ax[1].axhline(ymax_genes, linestyle='--', color=color_genes)
    ax[1].vlines(x=halfsat_genes, ymin=0, ymax=f(
        halfsat_genes, halfsat_genes, ymax_genes), linestyle=':', color=color_genes)
    ax[1].hlines(y=f(halfsat_genes, halfsat_genes, ymax_genes), xmin=0,
                 xmax=halfsat_genes, linestyle=':', color=color_genes)
    ax[1].text(halfsat_genes + xmax/20, 0.65*f(halfsat_genes, halfsat_genes, ymax_genes), 'half-saturation point: \n' +
               format(halfsat_genes, ',') + ' reads/cell, ' + '\n' + format(halfsat_genes_counts, ',') + ' genes/cell', size=8)
    ax[1].text(0.1*xmax, ymax_genes*1.01, format(ymax_genes, ',') + ' genes max', size=8)

    # label the axes
    ax[1].set_xlabel('Reads per cell')
    ax[1].set_ylabel('Median genes per cell')

    # ATAC saturation curve.
    if webSummaryType == 'ARC':
        def f(x, a, b):
            return(b * x / (x+a))

        popt, pcov = curve_fit(f, atac_reads, fragments, p0=[22000, 1000])
        halfsat_frags = np.round(popt[0], 0).astype('int')
        ymax_frags = np.round(popt[1], 0).astype('int')
        halfsat_frags_counts = np.round(f(halfsat_frags, halfsat_frags, ymax_frags)).astype(int)

        desiredUniqueFrags = np.round(f(readPairsDesired, halfsat_frags, ymax_frags)).astype(int)
        # ymax_genes = 1

        # Plot
        ax[2].scatter(atac_reads, fragments, color='green')
        ax[2].plot(atac_reads, f(atac_reads, *popt), 'r-')

        # extrapolate with the curve fit
        x_more = np.linspace(0, readMax, 100)
        xmax = np.max(x_more)
        ax[2].plot(x_more, f(x_more, *popt), 'k-')

        # set plot boundaries and add asymptotes and lines
        color_frags = 'g'
        ax[2].set_ylim(0, ymax_frags*1.1)
        ax[2].axhline(ymax_frags, linestyle='--', color=color_frags)
        ax[2].vlines(x=halfsat_frags, ymin=0, ymax=f(
            halfsat_frags, halfsat_frags, ymax_frags), linestyle=':', color=color_frags)
        ax[2].hlines(y=f(halfsat_frags, halfsat_frags, ymax_frags), xmin=0,
                     xmax=halfsat_frags, linestyle=':', color=color_frags)
        ax[2].text(halfsat_frags + xmax/20, 0.65*f(halfsat_frags, halfsat_frags, ymax_frags), 'half-saturation point: \n' +
                   format(halfsat_frags, ',') + ' read pairs/cell, ' + '\n' + format(halfsat_frags_counts, ',') + ' fragments/cell', size=8)
        ax[2].text(0.1*xmax, ymax_frags*1.01, format(ymax_frags, ',') + ' fragments max', size=8)
        # ax[1].text(xmax*0.08,ymax_genes*0.05,'current saturation: \n' + str(current_mean_reads) + ' reads/cell, ' +
        #           str(current_median_genes) + ' genes/cell', size=7)

        # label the axes
        ax[2].set_xlabel('Mean read pairs per cell')
        ax[2].set_ylabel('Median fragments per cell')

    # figure title
    if not title:
        if isinstance(sampname, str):
            title = sampname
    fig.suptitle(title, y=1, size=10)

    pl.subplots_adjust(wspace=0.5)
    pl.tight_layout()
    pl.show()
    # print('popt:',popt)
    # print('pcov:',pcov)
    if verbose is True:
        print()
        print('Sequencing saturation half-saturation point:',
              format(halfsat_sat, ','), 'reads per cell')
        print('Median genes per cell half-saturation point:',
              format(halfsat_genes_counts, ','), 'genes per cell', ';',
              format(halfsat_genes, ','), 'reads per cell')
        if webSummaryType == 'ARC':
            print('Median fragments per cell half-saturation point:',
                  format(halfsat_frags_counts, ','), 'fragments per cell', ';',
                  format(halfsat_frags, ','), 'read pairs per cell')
        print('Current sequencing saturation level:', current_sat)
        print('Current reads per cell:', current_mean_reads)
        print('Current genes per cell:', current_median_genes)
        if webSummaryType == 'ARC':
            print('Current median fragments per cell:', current_median_frag)
            print('Current mean raw read pairs per cell', current_mean_atac_reads)
        print()
        print('Desired reads per cell:', format(readsDesired, ','))
        print('Sequencing saturation for desired reads per cell:', "{:.1%}".format(desiredSeqSat))
        print('Median genes per cell for desired reads per cell:', format(desiredUniqueGenes, ','))
        if webSummaryType == 'ARC':
            print('Desired read pairs per cell:', format(readPairsDesired, ','))
            print('Fragments per cell for desired reads per cell:', format(desiredUniqueFrags, ','))

    if webSummaryType == 'ARC':
        ARC_table_data = {'sample_name': [sampname],
                          'sequencing saturation halfsat point (reads per cell)': [halfsat_sat],
                          'median genes halfsat (genes per cell)': [halfsat_genes_counts],
                          'median genes halfsat (reads per cell)': [halfsat_genes],
                          'median frags halfsat (frags per cell)': [halfsat_frags_counts],
                          'median frags halfsat (reads per cell)': [halfsat_frags],
                          'current sequencing saturation': [current_sat],
                          'current reads per cell': [current_mean_reads],
                          'current genes per cell': [current_median_genes],
                          'current median frags per cell': [current_median_frag],
                          'current mean raw read pairs per cell': [current_mean_atac_reads],
                          'desired reads per cell': [readsDesired],
                          'predicted sequencing sat for ' + str(readsDesired) +
                          ' reads per cell': ["{:.1%}".format(desiredSeqSat)],
                          'predicted median genes for ' + str(readsDesired) +
                          ' reads per cell': [desiredUniqueGenes],
                          'desired read pairs per cell': [readPairsDesired],
                          'predicted frags for ' + str(readPairsDesired) +
                          ' read pairs per cell': [desiredUniqueFrags]
                          }
        df = pd.DataFrame.from_dict(ARC_table_data)
    else:
        GEX_table_data = {'sample_name': [sampname],
                          'sequencing saturation halfsat point (reads per cell)': [halfsat_sat],
                          'median genes halfsat (genes per cell)': [halfsat_genes_counts],
                          'median genes halfsat (reads per cell)': [halfsat_genes],
                          'current sequencing saturation': [current_sat],
                          'current reads per cell': [current_mean_reads],
                          'current genes per cell': [current_median_genes],
                          'desired reads per cell': [readsDesired],
                          'predicted sequencing sat for ' + str(readsDesired) +
                          ' reads per cell': ["{:.1%}".format(desiredSeqSat)],
                          'predicted median genes for ' + str(readsDesired) +
                          ' reads per cell': [desiredUniqueGenes],
                          }
        df = pd.DataFrame.from_dict(GEX_table_data)
    return df


def find_satcurves(folder):
    """
    ***Needs updating*** currently deprecated
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
                # reads, genes, saturations, current_sat, sampname = scrape_saturation_stats(file_path)
                satcurves(file_path)
                # satcurves(reads, genes, saturations)
