from bs4 import BeautifulSoup
import re
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os

#__all__ = ['satcurves', 'find_satcurves']

class webSum_model:
    def __init__(self, web_summary_html_file, webSummaryType):
        self.web_summary_html_file = web_summary_html_file
        self.webSummaryType = webSummaryType
        self.reads, self.genes, self.saturations, self.current_sat, self.sampname, self.current_median_genes, \
        self.current_mean_reads = webSum_model.scrape_saturation_stats(web_summary_html_file)
        self.seqSatModel = None
        self.popt_saturation = None
        self.pcov_saturation = None
        self.halfsat_sat_reads_per_cell = None
        self.popt_genes = None
        self.pcov_genes = None
        self.halfsat_genes = None
        self.ymax_genes = None
        self.halfsat_genes_reads_per_cell = None

    @staticmethod
    def mm(x, a, b):
        # Gene saturation and fragment saturation curve.
        # a is the "kd"
        # b is the max genes per cell, or fraction saturation
        # here b=1
        return (b * x / (x + a))

    @staticmethod
    def mm_sat(x, a):
        return (x / (x + a))

    @staticmethod
    def lw(N, X):
        return (1 - (1 - np.exp((-N / X))) * (X / N))

    # reference: https://stackoverflow.com/questions/22277982/how-to-find-50-point-after-curve-fitting-using-numpy
    @staticmethod
    def lw2(N, X, y0=0):
        return (1 - (1 - np.exp((-N / X))) * (X / N) + y0)

    @staticmethod
    def scrape_saturation_stats(web_summary_html_file):
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

        constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n') + 1])

        try:
            seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
            current_sat = [
                name for name in seq_summary_table if 'Sequencing Saturation' in name][0][1]
            sampname = constant_data['summary']['sample']['id']
        except KeyError:
            print('check if GEX web summary html file because attributes not found')
            raise

        reads = np.array(constant_data['summary']['analysis_tab']
                         ['seq_saturation_plot']['plot']['data'][0]['x'])[
                1:]  # slice out zeros to avoid model error
        genes = np.array(constant_data['summary']['analysis_tab']
                         ['median_gene_plot']['plot']['data'][0]['y'])[1:]  # slice out zeros to avoid model error
        saturations = np.array(constant_data['summary']['analysis_tab']
                               ['seq_saturation_plot']['plot']['data'][0]['y'])[
                      1:]  # slice out zeros to avoid model error

        cells_table = constant_data['summary']['summary_tab']['cells']['table']['rows']
        for entry in cells_table:
            if entry[0] == 'Median Genes per Cell':
                current_median_genes = entry[1]
            if entry[0] == 'Mean Reads per Cell':
                current_mean_reads = entry[1]

        return reads, genes, saturations, current_sat, sampname, current_median_genes, current_mean_reads


    def fit_model(self, seqSatModel='lw', verbose=True):
        self.seqSatModel = seqSatModel
        #webSum_model.__scrape_saturation_stats(self)
        print('Model being used for sequencing  saturation is: ', seqSatModel)
        # Fit the curves for sequencing saturation
        if seqSatModel == 'lw':
            self.popt_saturation, self.pcov_saturation = curve_fit(webSum_model.lw, self.reads, self.saturations)
            # Find actual halfsat point for sequencing saturation for for lander-waterman
            from scipy.optimize import brentq
            a = np.min(self.reads)
            # https://stackoverflow.com/questions/53631988/valueerror-fa-and-fb-must-have-different-sign
            # Need a wider range to find the halfsat point
            b = np.max(self.reads) * 100
            y0 = -0.50
            solution = brentq(webSum_model.lw2, a, b, args=(self.popt_saturation[0], y0))
            self.halfsat_sat_reads_per_cell = np.round(solution, 0).astype('int')
        elif seqSatModel == 'mm':
            self.popt_saturation, self.pcov_saturation = curve_fit(webSum_model.mm_sat, self.reads, self.saturations)
            self.halfsat_sat_reads_per_cell = np.round(self.popt_saturation[0], 0).astype('int')
        else:
            raise Exception('Model must be either lander-waterman (lw) or michaelis-menten (mm)')
        ymax_sat = 1

        # Fit the curves for gene saturation
        self.popt_genes, self.pcov_genes = curve_fit(webSum_model.mm, self.reads, self.genes, p0=[22000, 1000])
        self.halfsat_genes_reads_per_cell = np.round(self.popt_genes[0], 0).astype('int')
        self.ymax_genes = np.round(self.popt_genes[1], 0).astype('int')
        self.halfsat_genes = np.round(webSum_model.mm(self.halfsat_genes_reads_per_cell,
                                                      self.halfsat_genes_reads_per_cell, self.ymax_genes)).astype(int)

        if verbose is True:
            print('')
            print('Here is the information pertaining to the sequencing saturation model:')
            print('reads from web_summary path: ', self.reads)
            print('UMIs from web_summary path: ', self.saturations)
            print('popt_saturation: {} and pcov_saturation: {}'.format(self.popt_saturation, self.pcov_saturation))
            print('halfsat_reads_per_cell: {}  ymax_sat: {}  halfsat_saturation: {}'.format(
                self.halfsat_sat_reads_per_cell, ymax_sat, 0.5))
            print('')
            print('Here is the information pertaining to the gene saturation model:')
            print('reads from web_summary path: ', self.reads)
            print('genes from web_summary path: ', self.genes)
            print('popt_genes: {} and pcov_genes: {}'.format(self.popt_genes, self.pcov_genes))
            print('halfsat_reads_per_cell: {}  ymax_genes: {}  halfsat_genes per cell: {}'.format(
                self.halfsat_genes_reads_per_cell, self.ymax_genes, self.halfsat_genes))

    def plot(self, readMax=80000, reads_test=None, saturations_test=None, genes_test=None):
        """
        Plot saturation curves from data scraped from web_summary.html file.
        Only works with Cellranger version 3.x and up.

        Args:
            web_summary_html_file (string): web_summary.html file to scrape.
            webSummaryType (string): either 'GEX' or 'ARC' based on Cellranger pipeline used
            seqSatModel (string): either 'lw' for lander-waterman fit or 'mm' for michaelis menten fit
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
        fit_attributes = [self.popt_saturation, self.pcov_saturation, self.popt_genes, self.pcov_genes]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_model() first to set necessary attributes.')
        else:
            # Plot
            plt.rcParams['figure.dpi'] = 120
            fig, ax = plt.subplots(1, 2, figsize=[6.75, 2.75])

            # constraining to only 3 bins in x-axis for visiblity purposes
            ax[0].locator_params(axis='x', nbins=4)
            ax[1].locator_params(axis='x', nbins=4)

            # Plot
            ax[0].scatter(self.reads, self.saturations, label='train')
            if self.seqSatModel == 'lw':
                ax[0].plot(self.reads, webSum_model.lw(self.reads, *self.popt_saturation), 'r-')
                # extrapolate with the curve fit
                x_more = np.linspace(0, readMax, 100)
                xmax = np.max(x_more)
                ax[0].plot(x_more, webSum_model.lw(x_more, *self.popt_saturation), 'k-')
            elif self.seqSatModel == 'mm':
                ax[0].plot(self.reads, webSum_model.mm_sat(self.reads, *self.popt_saturation), 'r-')
                # extrapolate with the curve fit
                x_more = np.linspace(0, readMax, 100)
                xmax = np.max(x_more)
                ax[0].plot(x_more, webSum_model.mm_sat(x_more, *self.popt_saturation), 'k-')
            else:
                raise ValueError('seqSatModel must be either \'lw\' for lander-waterman or \'mm\' for michaelis-menten')

            # set plot boundaries and add asymptotes and lines
            color_sat = 'purple'
            ymax_sat=1
            ax[0].set_ylim(0, 1.1)
            ax[0].axhline(ymax_sat, linestyle='--', color=color_sat)
            if self.seqSatModel == 'lw':
                ax[0].vlines(x=self.halfsat_sat_reads_per_cell, ymin=0, ymax=webSum_model.lw(
                    self.halfsat_sat_reads_per_cell, self.popt_saturation[0]), linestyle=':', color=color_sat)
                ax[0].hlines(y=webSum_model.lw(self.halfsat_sat_reads_per_cell, self.popt_saturation[0]), xmin=0,
                             xmax=self.halfsat_sat_reads_per_cell, linestyle=':', color=color_sat)
                ax[0].text(self.halfsat_sat_reads_per_cell + xmax / 20, 0.4 * webSum_model.lw(self.halfsat_sat_reads_per_cell,
                                                                                     self.halfsat_sat_reads_per_cell),
                           'half-saturation point: \n' + format(self.halfsat_sat_reads_per_cell, ',') + ' reads/cell', size=8)
            elif self.seqSatModel == 'mm':
                ax[0].vlines(x=self.halfsat_sat_reads_per_cell, ymin=0, ymax=webSum_model.mm_sat(
                    self.halfsat_sat_reads_per_cell, self.popt_saturation[0]), linestyle=':', color=color_sat)
                ax[0].hlines(y=webSum_model.mm_sat(self.halfsat_sat_reads_per_cell, self.popt_saturation[0]), xmin=0,
                             xmax=self.halfsat_sat_reads_per_cell, linestyle=':', color=color_sat)
                ax[0].text(self.halfsat_sat_reads_per_cell + xmax / 20, 0.4 * webSum_model.mm_sat(self.halfsat_sat_reads_per_cell,
                                                                                         self.halfsat_sat_reads_per_cell),
                           'half-saturation point: \n' + format(self.halfsat_sat_reads_per_cell, ',') + ' reads/cell', size=8)
            else:
                raise ValueError('seqSatModel must be either \'lw\' for lander-waterman or \'mm\' for michaelis-menten')
            # pl.show()
            # label the axes
            ax[0].set_xlabel('Reads per cell')
            ax[0].set_ylabel('Sequencing Saturations')


            # Gene saturation curve.
            # a is the "kd"
            # b is the max genes per cell, or fraction saturation
            # here b=1

            # Plot
            ax[1].scatter(self.reads, self.genes, color='red', label='train')
            ax[1].plot(self.reads, webSum_model.mm(self.reads, *self.popt_genes), 'r-')

            # extrapolate with the curve fit
            x_more = np.linspace(0, readMax, 100)
            xmax = np.max(x_more)
            ax[1].plot(x_more, webSum_model.mm(x_more, *self.popt_genes), 'k-')

            # set plot boundaries and add asymptotes and lines
            color_genes = 'r'
            ax[1].set_ylim(0, self.ymax_genes * 1.1)
            ax[1].axhline(self.ymax_genes, linestyle='--', color=color_genes)
            ax[1].vlines(x=self.halfsat_genes_reads_per_cell, ymin=0, ymax=webSum_model.mm(
                self.halfsat_genes_reads_per_cell, self.halfsat_genes_reads_per_cell, self.ymax_genes), linestyle=':',
                         color=color_genes)
            ax[1].hlines(y=webSum_model.mm(self.halfsat_genes_reads_per_cell, self.halfsat_genes_reads_per_cell,
                                           self.ymax_genes), xmin=0, xmax=self.halfsat_genes_reads_per_cell,
                         linestyle=':', color=color_genes)
            ax[1].text(self.halfsat_genes_reads_per_cell + xmax / 20, 0.65 * webSum_model.mm(
                self.halfsat_genes_reads_per_cell, self.halfsat_genes_reads_per_cell, self.ymax_genes),
                       'half-saturation point: \n' +
                       format(self.halfsat_genes_reads_per_cell, ',') + ' reads/cell, ' + '\n' +
                       format(self.halfsat_genes, ',') + ' genes/cell', size=8)
            ax[1].text(0.1 * xmax, self.ymax_genes * 1.01, format(self.ymax_genes, ',') + ' genes max', size=8)
            # label the axes
            ax[1].set_xlabel('Reads per cell')
            ax[1].set_ylabel('Median genes per cell')

            if all(v is not None for v in [reads_test, genes_test, saturations_test]):
                ax[0].scatter(reads_test, saturations_test, label='test')
                ax[0].legend()
                ax[1].scatter(reads_test, genes_test, label='test')
                ax[1].legend()

            fig.tight_layout()

    def predict_seq_saturation(self, reads_test):
        # ensure that model has been fit first
        fit_attributes = [self.popt_saturation, self.pcov_saturation, self.popt_genes, self.pcov_genes]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_model() first to set necessary attributes.')

        if self.seqSatModel == 'lw':
            seq_sat_predictions = [round(webSum_model.lw(read, self.popt_saturation[0]), 4)
                                   for read in reads_test]
            return seq_sat_predictions

        elif self.seqSatModel == 'mm':
            seq_sat_predictions = [round(webSum_model.mm_sat(read, self.popt_saturation[0]), 4)
                                   for read in reads_test]
            return seq_sat_predictions
        else:
            raise ValueError('seqSatModel must be either \'lw\' for lander-waterman or \'mm\' for michaelis-menten')

    def score_seq_saturation(self, reads_test, seq_saturation_test):
        # ensure that model has been fit first
        fit_attributes = [self.popt_saturation, self.pcov_saturation, self.popt_genes, self.pcov_genes]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_model() first to set necessary attributes.')

        if self.seqSatModel == 'lw':
            seq_sat_predictions = [round(webSum_model.lw(read, self.popt_saturation[0]), 4)
                                   for read in reads_test]
        elif self.seqSatModel == 'mm':
            seq_sat_predictions = [round(webSum_model.mm_sat(read, self.popt_saturation[0]), 4)
                                   for read in reads_test]
        else:
            raise ValueError('seqSatModel must be either \'lw\' for lander-waterman or \'mm\' for michaelis-menten')

        RSS = (np.subtract(np.array(seq_saturation_test), np.array(seq_sat_predictions)) ** 2).sum()
        ymean = np.array(seq_saturation_test).mean()
        TSS = (np.subtract(np.array(seq_saturation_test), np.array(ymean)) ** 2).sum()
        Rsquared = 1 - RSS/TSS
        print('RSS: {}  ymean: {}  TSS: {}  Rsquared: {}'.format(RSS, ymean, TSS, Rsquared))
        return Rsquared

    def predict_genes(self, reads_test):
        # ensure that model has been fit first
        fit_attributes = [self.popt_saturation, self.pcov_saturation, self.popt_genes, self.pcov_genes]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_model() first to set necessary attributes.')

        gene_predictions = [np.round(webSum_model.mm(read, self.halfsat_genes_reads_per_cell,
                                                     self.ymax_genes)).astype(int) for read in reads_test]
        return gene_predictions

    def score_genes(self, reads_test, genes_test):
        # ensure that model has been fit first
        fit_attributes = [self.popt_saturation, self.pcov_saturation, self.popt_genes, self.pcov_genes]
        if all(v is None for v in fit_attributes):
            raise AttributeError('Please run fit_model() first to set necessary attributes.')

        gene_predictions = [np.round(webSum_model.mm(read, self.halfsat_genes_reads_per_cell,
                                                     self.ymax_genes)).astype(int) for read in reads_test]

        RSS = (np.subtract(np.array(genes_test), np.array(gene_predictions)) ** 2).sum()
        ymean = np.array(genes_test).mean()
        TSS = (np.subtract(np.array(genes_test), np.array(ymean)) ** 2).sum()
        Rsquared = 1 - RSS / TSS
        print('RSS: {}  ymean: {}  TSS: {}  Rsquared: {}'.format(RSS, ymean, TSS, Rsquared))
        return Rsquared

def get_reads_genes_sats_from_web(web_summary_html):
    reads, genes, saturations, current_sat, sampname, current_median_genes, current_mean_reads = \
        webSum_model.scrape_saturation_stats(web_summary_html)
    return reads, saturations, genes

