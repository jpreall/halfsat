import re
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os

#__all__ = ['plot_UMI_curve']


class UMI_model:
    """Class that contains UMIs per cell versus reads per cell model. Based on michaelis-menten equation.

    Attributes:
        jsonPath (str): path to 10x CellRanger count metrics summary json file
        UMI_dictionary (dict): dictionary of reads and corresponding UMI amounts
        sample_id (str): array of genes corresponding to 10x CellRanger count genes plot
        current_UMIs (str): current UMIs of sample
        current_reads_per_cell (str): current reads per cell of sample
        reads (array_like): array of reads for sample
        UMIs (array_like): array of UMIs for given reads of sample.
        halfsat_reads_per_cell (int): estimated reads per cell for half saturation of UMIs
        halfsat_UMI_counts (int): estimated UMIs per cell for half saturation of UMIs
        ymax_UMIs (int): max number of UMIs per cell based on model
        popt (array): optimal values for the parameters so that the ssr of UMIs function is minimized
        pcov (2d array): the estimated covariance of popt

    """
    def __init__(self, jsonPath):
        self.jsonPath = jsonPath
        self.UMI_dictionary, self.sample_id, \
            self.current_UMIs, self.current_reads_per_cell = UMI_model.scrapeUMIinfo(jsonPath)
        # sort reads (and associated UMIs) so that ax.plot(reads, f(reads, *popt), 'r-')
        # connects the points with a line in order
        self.reads = np.sort(np.array(list(self.UMI_dictionary.keys())).astype(float))
        self.UMIs = [self.UMI_dictionary[str(int(key))] for key in self.reads]
        self.halfsat_reads_per_cell = None  # default
        self.halfsat_UMI_counts = None  # default
        self.ymax_UMIs = None  # default
        self.popt = None  # default
        self.pcov = None  # default

    @staticmethod
    def scrapeUMIinfo(jsonPath):
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
                            for key in metrics_summary if
                            'filtered_bcs_median_counts' in key and 'subsampled' not in key][0]
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
            sliced_key = UMI_model.__slicer(key, 'raw_rpc')
            abbr_key = re.sub("raw_rpc_", "", re.sub(
                "_subsampled_filtered_bcs_median_counts", "", sliced_key))
            UMI_dict_abbr[abbr_key] = UMI_dict[key]

        # remove the impossible values for 0 UMIs
        for k, v in list(UMI_dict_abbr.items()):
            if v == 0:
                del UMI_dict_abbr[k]
        # print('UMI_dict_abbr: ', UMI_dict_abbr)
        return UMI_dict_abbr, sample_id, current_UMIs, current_reads_per_cell

    @staticmethod
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


    @staticmethod
    def f(x, a, b):
        """
        Michaelis menten model.
        Args:
            x: reads per cell
            a: a is the "kd"
            b: b is the max UMIs per cell

        Returns:
            Michaelis-menten model.

        """
        return (b * x / (x + a))


    def fit_UMIs(self, verbose=True):
        """
        Fit michaelis-menten model with given reads and UMIs
        Args:
            verbose: bool, whether to print information for user.
        """

        self.popt, self.pcov = curve_fit(UMI_model.f, self.reads, self.UMIs, p0=[22000, 1000])
        self.halfsat_reads_per_cell = np.round(self.popt[0], 0).astype('int')
        self.ymax_UMIs = np.round(self.popt[1], 0).astype('int')
        self.halfsat_UMI_counts = np.round(UMI_model.f(self.halfsat_reads_per_cell,
                                             self.halfsat_reads_per_cell, self.ymax_UMIs)).astype(int)

        if verbose is True:
            print('Here is the information pertaining to the model:')
            print('reads from JSON path: ', self.reads)
            print('UMIs from JSON path: ', self.UMIs)
            print('popt: {} and pcov: {}'.format(self.popt, self.pcov))
            print('halfsat_reads_per_cell: {}  ymax_UMIs: {}  halfsat_UMI_counts: {}'.format(
                self.halfsat_reads_per_cell, self.ymax_UMIs, self.halfsat_UMI_counts))

    def plot_UMIs(self, readMax=80000, reads_test=None, UMIs_test=None):
        """
        Plot Unique UMIs detected per cell versus reads per cell.

        Plot UMIs versus reads graph and provide dataframe with prediction of UMIs
        per cell given a specified number of reads.

        Args:
            readMax (int): length of the plot.
            readsDesired: number of reads to predict UMIs per cell.
            reads_test (array_like): array for read values to test model against.
            UMIs_test (array_like): array of UMI values to test model against.

        Returns:
            None

        Raises:
            AttributeError: run fit_UMIs() first to set necessary attributes

        """
        fit_attributes = [self.popt, self.pcov, self.halfsat_reads_per_cell, self.ymax_UMIs,
                          self.halfsat_UMI_counts]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_UMIs() first to set necessary attributes.')
        else:
            # Plot
            plt.rcParams['figure.dpi'] = 120
            fig, ax = plt.subplots(1, 1, figsize=[4, 2.75])
            ax.scatter(self.reads, self.UMIs, color='red', label='train')
            ax.plot(self.reads, UMI_model.f(self.reads, *self.popt), 'r-')

            # extrapolate with the curve fit
            x_more = np.linspace(0, readMax, 100)
            xmax = np.max(x_more)
            ax.plot(x_more, UMI_model.f(x_more, *self.popt), 'k-')

            # set plot boundaries and add asymptotes and lines
            color_UMIs = 'r'
            ax.set_ylim(0, self.ymax_UMIs * 1.1)
            ax.axhline(self.ymax_UMIs, linestyle='--', color=color_UMIs)
            ax.vlines(x=self.halfsat_reads_per_cell, ymin=0, ymax=UMI_model.f(
                self.halfsat_reads_per_cell, self.halfsat_reads_per_cell, self.ymax_UMIs), linestyle=':',
                      color=color_UMIs)
            ax.hlines(y=UMI_model.f(self.halfsat_reads_per_cell, self.halfsat_reads_per_cell, self.ymax_UMIs), xmin=0,
                      xmax=self.halfsat_reads_per_cell, linestyle=':', color=color_UMIs)
            ax.text(self.halfsat_reads_per_cell + xmax / 20, 0.65 * UMI_model.f(
                self.halfsat_reads_per_cell, self.halfsat_reads_per_cell, self.ymax_UMIs), 'half-saturation point: \n' +
                    format(self.halfsat_reads_per_cell, ',') + ' reads/cell, ' + '\n' +
                    format(self.halfsat_UMI_counts, ',') + ' UMIs/cell', size=8)
            ax.text(0.1 * xmax, self.ymax_UMIs * 1.01, format(self.ymax_UMIs, ',') + ' UMIs max', size=8)

            if all(v is not None for v in [reads_test, UMIs_test]):
                ax.scatter(reads_test, UMIs_test, label='test')
                ax.legend()
            # label the axes
            ax.set_xlabel('Reads per cell')
            ax.set_ylabel('Unique UMIs detected per cell')

            ax.set_title(self.sample_id)

            plt.subplots_adjust(wspace=0.5)
            plt.tight_layout()
            plt.show()


    def make_UMI_table(self, readsDesired=40000):
        """
        Make table of UMI predictions and other key attributes.
        Args:
            readsDesired (int): number of reads to predict number of UMIs for.

        Returns:
            dataframe of prediction and other UMI_model attributes.

        Raises:
            AttributeError: run fit_UMIs() first to set necessary attributes

        """
        fit_attributes = [self.popt, self.pcov, self.halfsat_reads_per_cell, self.ymax_UMIs,
                          self.halfsat_UMI_counts]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_UMIs() first to set necessary attributes.')

        desiredUniqueUMIs = np.round(UMI_model.f(readsDesired, self.halfsat_reads_per_cell, self.ymax_UMIs)).astype(int)
        UMI_table_data = {'sample_id': [self.sample_id],
                          'reads_per_cell': [int(self.current_reads_per_cell)],
                          'current_UMIs': [int(self.current_UMIs)],
                          'halfsat_reads_per_cell': [self.halfsat_reads_per_cell],
                          'halfsat_UMI_counts': [self.halfsat_UMI_counts],
                          'max_UMIs': [self.ymax_UMIs],
                          'predicted UMIs for ' + str(readsDesired) +
                          ' reads per cell': [desiredUniqueUMIs]
                          }

        # pass column names in the columns parameter
        df = pd.DataFrame.from_dict(UMI_table_data)
        return df


    def predict(self, reads_test):
        """
        Predict UMIs for an array of test read values.
        Args:
            reads_test (array_like): array of reads values to predict UMIs from

        Returns:
            array of predicterd UMIs for given read values

        Raises:
            AttributeError: run fit_UMIs() first to set necessary attributes

        """
        # ensure that model has been fit first
        fit_attributes = [self.popt, self.pcov, self.halfsat_reads_per_cell, self.ymax_UMIs,
                          self.halfsat_UMI_counts]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_UMIs() first to set necessary attributes.')

        UMI_predictions = [np.round(UMI_model.f(read, self.halfsat_reads_per_cell, self.ymax_UMIs)).astype(int)
                           for read in reads_test]
        return UMI_predictions

    def score(self, reads_test, UMIs_test):
        """
        Score UMI predictions by comparing with actual UMI test values using Rsquared.

        Args:
            reads_test (array_like): array of reads values to predict UMIs from
            UMIs_test (array_like): array of UMI values to test predicted UMI values.

        Returns:
            Rsquared value for the goodness of model fit

        Raises:
            AttributeError: run fit_UMIs() first to set necessary attributes
        """
        # ensure that model has been fit first
        fit_attributes = [self.popt, self.pcov, self.halfsat_reads_per_cell, self.ymax_UMIs,
                          self.halfsat_UMI_counts]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_UMIs() first to set necessary attributes.')

        UMI_predictions = [np.round(UMI_model.f(read, self.halfsat_reads_per_cell, self.ymax_UMIs)).astype(int)
                           for read in reads_test]

        RSS = (np.subtract(np.array(UMIs_test), np.array(UMI_predictions)) ** 2).sum()
        ymean = np.array(UMIs_test).mean()
        TSS = (np.subtract(np.array(UMIs_test), np.array(ymean)) ** 2).sum()
        Rsquared = 1 - RSS/TSS
        print('RSS: {}  ymean: {}  TSS: {}  Rsquared: {}'.format(RSS, ymean, TSS, Rsquared))
        return Rsquared


def get_reads_and_UMIs_from_json(jsonPath):
    """
    Scrape information relating to UMI of a sample's metrics summary json file (reads and UMIs only).
    Currently only tested on Cellranger 6.0.0+ metrics_summary_json.json files.

    Args:
        jsonPath (string): path to metrics_summary_json.json file

    Returns:
        array of reads from jsonPath produced from CellRanger count
        array of UMIs from jsonPath produced from CellRanger count

    """
    UMI_dict_abbr, sample_id, current_UMIs, current_reads_per_cell = UMI_model.scrapeUMIinfo(jsonPath)
    reads = np.sort(np.array(list(UMI_dict_abbr.keys())).astype(float))
    UMIs = [UMI_dict_abbr[str(int(key))] for key in reads]
    return reads, UMIs


def aggr_UMI_table(jsonPathList, readsDesired=40000, verbose=False):
    """
    Create table of UMI halfsat values and UMI predictions for multiple samples.

    Args:
        jsonPathList (list): a list of json pathes to samples.
        readsDesired (int): number of reads to estimate the UMIs/cell for.
        showPlot (boolean): whether to print halfsat UMI plots for each sample

    Returns:
        sample_df (pandas DataFrame): DataFrame containing UMI numbers for each sample.

    """

    initial = True
    for i in jsonPathList:
        if initial is True:
            initial_UMI_obj = UMI_model(i)
            initial_UMI_obj.fit_UMIs(verbose)
            sample_df = initial_UMI_obj.make_UMI_table(readsDesired)
            initial = False
        else:
            tmp_UMI_obj = UMI_model(i)
            tmp_UMI_obj.fit_UMIs(verbose)
            tmp_df = tmp_UMI_obj.make_UMI_table(readsDesired)
            sample_df = sample_df.append(tmp_df, ignore_index=True)
    return sample_df
