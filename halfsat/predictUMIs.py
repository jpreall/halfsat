import re
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os

__all__ = ['plot_UMI_curve']


def __scrapeUMIinfo(jsonPath):
    """
    Scrape information relating to UMI of a sample's metrics summary json file.
    Currently only tested on Cellranger 6.0.0+ metrics_summary_json.json files.

    Args:
        jsonPath (string): path to metrics_summary_json.json file

    Returns:
        dict: dictionary of reads_per_cell key and UMI values
        string: sample_id of sample
        float: current UMI count for sample
        float: current read per cell for sample
    Raises:
        ValueError: Raise when JSON file not provided.
        IndexError: Raise when unable to index JSON file properly.

    """
    with open(jsonPath) as jsonFile:
        try:
            metrics_summary = json.load(jsonFile)
            # print(metrics_summary)
        except ValueError:
            print('Decoding JSON has failed. JSON file likely not provided')
            raise

    sample_id = metrics_summary['sample_id']
    try:
        current_UMIs = [metrics_summary[key]
                        for key in metrics_summary if 'filtered_bcs_median_counts' in key and 'subsampled' not in key][0]
    except IndexError:
        print('Unable to index properly. Try again with different json.')
        raise
    current_reads_per_cell = metrics_summary['reads_per_cell']

    UMI_reads_per_cell_list = ['raw_rpc', 'subsampled_filtered_bcs_median_counts']
    subsampled_reads_UMIs = [
        # https://stackoverflow.com/questions/58045147/simultaneously-checking-for-multiple-substrings-using-list-comprehension
        x for x in metrics_summary if all(entry in x for entry in UMI_reads_per_cell_list)]

    UMI_dict = {}
    for i in subsampled_reads_UMIs:
        if i not in UMI_dict:
            UMI_dict[i] = metrics_summary[i]
    # print('UMI_dict: ', UMI_dict)

    UMI_dict_abbr = {}
    for key in UMI_dict:
        # slice out the genome identifier before raw_rpc
        sliced_key = __slicer(key, 'raw_rpc')
        abbr_key = re.sub("raw_rpc_", "", re.sub(
            "_subsampled_filtered_bcs_median_counts", "", sliced_key))
        UMI_dict_abbr[abbr_key] = UMI_dict[key]

    # remove the impossible values for 0 UMIs
    for k, v in list(UMI_dict_abbr.items()):
        if v == 0:
            del UMI_dict_abbr[k]
    # print('UMI_dict_abbr: ', UMI_dict_abbr)
    return UMI_dict_abbr, sample_id, current_UMIs, current_reads_per_cell


def __slicer(my_str, sub):
    """
    Remove everything in a string before a specified substring is found.
    Throw exception if substring is not found in string

    https://stackoverflow.com/questions/33141595/how-can-i-remove-everything-in-a-string-until-a-characters-are-seen-in-python

    Args:
        my_str (string): the string to slice.
        sub (string): the substring to stop slicing at.

    Returns:
        str: substring of my_str, without everything before sub.

    Raises:
        Exception: Sub string specified is not found in my_str.

    """

    index = my_str.find(sub)
    if index != -1:
        return my_str[index:]
    else:
        # raise Exception('Sub string not found!')
        return my_str


def plot_UMI_curve(jsonPath, readMax=80000, readsDesired=40000, showPlot=True):
    """
    Plot Unique UMIs detected per cell versus reads per cell.

    Plot UMIs versus reads graph and provide dataframe with prediction of UMIs
    per cell given a specified number of reads.

    Args:
        jsonPath (string): path to metrics_summary_json.json file.
        readMax (int): length of the plot.
        readsDesired: number of reads to predict UMIs per cell.
        showPlot: whether or not to plot the curve or only return dataframe

    Returns:
        Pandas dataframe: Table with sample id, current metrics, and predicted UMIs.

    Raises:
        Exception: description

    """

    UMI_dict, my_sample, my_UMIs, my_reads_per_cell = __scrapeUMIinfo(jsonPath)
    # print('UMI_dict_abbr: ', UMI_dict)

    # sort reads (and associated UMIs) so that ax.plot(reads, f(reads, *popt), 'r-')
    # connects the points with a line in order
    reads = np.sort(np.array(list(UMI_dict.keys())).astype(float))
    # print('reads: ', reads)
    UMIs = [UMI_dict[str(int(key))] for key in reads]

    # print('reads: ', reads)
    # print('UMIs: ', UMIs)

    # UMI saturation curve. b is vmax (max umis)
    # a is km the michaelis menten constant, reads per cell at half vmax
    # [reads, UMIs]
    # better way to optimize it. Might break with low seq depth (2 mill reads)
    # KC's myseq, or JG, DTPark
    def f(x, a, b):
        return(b * x / (x+a))

    popt, pcov = curve_fit(f, reads, UMIs, p0=[22000, 1000])
    halfsat_UMIs = np.round(popt[0], 0).astype('int')
    ymax_UMIs = np.round(popt[1], 0).astype('int')
    halfsat_UMI_counts = np.round(f(halfsat_UMIs, halfsat_UMIs, ymax_UMIs)).astype(int)

    desiredUniqueUMIs = np.round(f(readsDesired, halfsat_UMIs, ymax_UMIs)).astype(int)
    # ymax_genes = 1

    if showPlot is True:
        # Plot
        plt.rcParams['figure.dpi'] = 120
        fig, ax = plt.subplots(1, 1, figsize=[4, 2.75])
        ax.scatter(reads, UMIs, color='red')
        ax.plot(reads, f(reads, *popt), 'r-')

        # extrapolate with the curve fit
        x_more = np.linspace(0, readMax, 100)
        xmax = np.max(x_more)
        ax.plot(x_more, f(x_more, *popt), 'k-')

        # set plot boundaries and add asymptotes and lines
        color_UMIs = 'r'
        ax.set_ylim(0, ymax_UMIs*1.1)
        ax.axhline(ymax_UMIs, linestyle='--', color=color_UMIs)
        ax.vlines(x=halfsat_UMIs, ymin=0, ymax=f(
            halfsat_UMIs, halfsat_UMIs, ymax_UMIs), linestyle=':', color=color_UMIs)
        ax.hlines(y=f(halfsat_UMIs, halfsat_UMIs, ymax_UMIs), xmin=0,
                  xmax=halfsat_UMIs, linestyle=':', color=color_UMIs)
        ax.text(halfsat_UMIs + xmax/20, 0.65*f(halfsat_UMIs, halfsat_UMIs, ymax_UMIs), 'half-saturation point: \n' +
                format(halfsat_UMIs, ',') + ' reads/cell, ' + '\n' + format(halfsat_UMI_counts, ',') + ' UMIs/cell', size=8)
        ax.text(0.1*xmax, ymax_UMIs*1.01, format(ymax_UMIs, ',') + ' UMIs max', size=8)

        # label the axes
        ax.set_xlabel('Reads per cell')
        ax.set_ylabel('Unique UMIs detected per cell')

        ax.set_title(my_sample)

        plt.subplots_adjust(wspace=0.5)
        plt.tight_layout()
        plt.show()

    UMI_table_data = {'sample_id': [my_sample],
                      'reads_per_cell': [int(my_reads_per_cell)],
                      'current_UMIs': [int(my_UMIs)],
                      'halfsat_UMIs_reads_per_cell': [halfsat_UMIs],
                      'halfsat_UMIs': [halfsat_UMI_counts],
                      'max_UMIs': [ymax_UMIs],
                      'predicted UMIs for ' + str(readsDesired) +
                      ' reads per cell': [desiredUniqueUMIs]
                      }

    # pass column names in the columns parameter
    df = pd.DataFrame.from_dict(UMI_table_data)
    return df


def aggr_UMI_table(jsonPathList, readsDesired=40000, showPlot=False):
    """
    Create table of UMI halfsat values and UMI predictions for multiple samples.

    Args:
        jsonPathList (list): a list of json pathes to samples.
        readsDesired (int): number of reads to estimate the UMIs/cell for.
        showPlot (boolean): whether to print halfsat UMI plots for each sample

    Returns:
        sample_df (pandas DataFrame): DataFrame containing UMI numbers for each sample.

    Raises:
        Exception: description

    """

    initial = True
    for i in jsonPathList:
        if initial is True:
            sample_df = plot_UMI_curve(i, readsDesired=readsDesired, showPlot=showPlot)
            initial = False
        else:
            tmp = plot_UMI_curve(i, readsDesired=readsDesired, showPlot=showPlot)
            sample_df = sample_df.append(tmp)
    return sample_df
