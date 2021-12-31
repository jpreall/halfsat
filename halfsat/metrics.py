import re
import json
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

__all__ = ['tableGenerator']


def __attributeScraper(html, webSummaryType='GEX'):
    """
    Generate a dictionary of all attributes and values from a web_summary.html file.

    Scrape attributes from cellranger v3-6 web_summary.html to generate a dictionary of
    attributes and values from all metrics in the json portion of web_summary.html.

    Args:
        html (str): Path to web_summary.html file of desired experiment/sample.
        webSummaryType (string): either 'GEX' or 'ARC' based on Cellranger pipeline used.

    Returns:
        dict: Dictionary of attribute (key) and corresponding value. 2 dictionaries
              if 'ARC' pipeline is used.

    Raises:
        ValueError: when neither 'GEX' nor 'ARC' is supplied.
        KeyError: when webSummaryType does not match the web summary html file.

    """

    f = open(html, encoding="utf8")
    soup = BeautifulSoup(f, "html.parser")
    for line in soup.find('script'):
        if 'const data' in line:
            const_data = line
    f.close()

    # need to include +1 to include end of string (the last })
    constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n')+1])

    if webSummaryType == 'GEX':
        # extract data from JSON tables
        try:
            pipeline_table = constant_data['summary']['summary_tab']['pipeline_info_table']['rows']
            seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
            map_summary_table = constant_data['summary']['summary_tab']['mapping']['table']['rows']
            cells_summary_table = constant_data['summary']['summary_tab']['cells']['table']['rows']

        except KeyError:
            print('check GEX web summary html file because attributes not found')
            raise

        tableList = [pipeline_table, seq_summary_table, map_summary_table, cells_summary_table]

        # create dictionary of all attributes
        attributeDict = {}
        for table in tableList:
            for entry in table:
                if entry[0] not in attributeDict:
                    attributeDict[entry[0]] = entry[1]

        if 'diagnostics' in constant_data['summary']:
            diagnostics_dict = constant_data['summary']['diagnostics']
            attributeDict.update(diagnostics_dict)

        return attributeDict

    elif webSummaryType == 'ARC':
        try:
            joint_pipeline_table = constant_data['joint_pipeline_info_table']['rows']
            joint_metrics_table = constant_data['joint_metrics_table']['rows']

            atac_sequencing_table = constant_data['atac_sequencing_table']['rows']
            atac_cells_table = constant_data['atac_cells_table']['rows']
            atac_mapping_table = constant_data['atac_mapping_table']['rows']
            atac_targeting_table = constant_data['atac_targeting_table']['rows']

            gex_sequencing_table = constant_data['gex_sequencing_table']['rows']
            gex_cells_table = constant_data['gex_cells_table']['rows']
            gex_mapping_table = constant_data['gex_mapping_table']['rows']

        except KeyError:
            print('check ARC web summary html file because attributes not found')
            raise

        arcATAC_tableList = [joint_pipeline_table, joint_metrics_table, atac_sequencing_table,
                             atac_cells_table, atac_mapping_table, atac_targeting_table]

        arcGEX_tableList = [joint_pipeline_table, joint_metrics_table, gex_sequencing_table,
                            gex_sequencing_table, gex_cells_table, gex_mapping_table]

        arcATAC_attributeDict = {}
        for table in arcATAC_tableList:
            for entry in table:
                if entry[0] not in arcATAC_attributeDict:
                    arcATAC_attributeDict[entry[0]] = entry[1]

        arcGEX_attributeDict = {}
        for table in arcGEX_tableList:
            for entry in table:
                if entry[0] not in arcGEX_attributeDict:
                    arcGEX_attributeDict[entry[0]] = entry[1]

        return arcATAC_attributeDict, arcGEX_attributeDict

    else:
        raise ValueError('webSummaryType neither GEX nor ARC')


def tableGenerator(htmlList, webSummaryType='GEX', tableType='full', readsDesired=40000):
    """
    Use list of cellranger (v4+) count web_summary.html () to generate a table.

    Use list of web_summary.html to generate a full table of all cellranger count
    metrics, abbreviated table for data delivery, or a table for repooling purposes.

    Args:
        htmlList (list): list of web_summary.html directories.
        tableType (string): string of one of the following: full, data delivery,
                            repooling.
        readsDesired (int): number of reads per cell desired for repooling.

    Returns:
        pandas DataFrame: DataFrame for table specified.

    Raises:
        AssertionError: htmlList must be of list type.
        ValueError: when neither 'GEX' nor 'ARC' is supplied.
        Exception: for tableType argument, please use one of the following: full,
                   delivery doc, repooling.
        AssertionError: readsDesired must be higher than the mean reads/cell of each
                        sample being repooled.
        Exception: if 'ARC' pipeline is supplied, inform user that 'repooling' and
                   'delivery doc' args have not yet beem implemented.
    """

    assert isinstance(htmlList, list), "Please input a list of html directories"
    ATAC_full_df = pd.DataFrame()
    GEX_full_df = pd.DataFrame()

    if webSummaryType == 'GEX':
        initial = True
        # Create dataframe containing all of the web_summaries in list
        for sample in htmlList:
            sample_dict = __attributeScraper(sample, webSummaryType)
            # https://stackoverflow.com/questions/57631895/dictionary-to-dataframe-error-if-using-all-scalar-values-you-must-pass-an-ind
            sample_df = pd.DataFrame([sample_dict])
            if initial is True:
                full_df = sample_df
                initial = False
            else:
                full_df = full_df.append(sample_df)

    elif webSummaryType == 'ARC':
        initial = True
        # Create dataframe containing all of the web_summaries in list
        for sample in htmlList:
            ATAC_sample_dict, GEX_sample_dict = __attributeScraper(sample, webSummaryType)
            # https://stackoverflow.com/questions/57631895/dictionary-to-dataframe-error-if-using-all-scalar-values-you-must-pass-an-ind
            ATAC_sample_df = pd.DataFrame([ATAC_sample_dict])
            GEX_sample_df = pd.DataFrame([GEX_sample_dict])
            if initial is True:
                ATAC_full_df = ATAC_sample_df
                GEX_full_df = GEX_sample_df
                initial = False
            else:
                ATAC_full_df = ATAC_full_df.append(ATAC_sample_df)
                GEX_full_df = GEX_full_df.append(GEX_sample_df)
    else:
        raise ValueError('webSummaryType neither GEX nor ARC')

    def fullTableMaker(fullDataframe):
        # delete columns like "Sample Description" if empty/full of NaNs
        fullDataframe.replace(r'^\s*$', np.nan, inplace=True, regex=True)
        # reference: https://www.jitsejan.com/find-and-delete-empty-columns-pandas-dataframe.html
        empty_cols = [col for col in fullDataframe.columns if fullDataframe[col].isnull().all()]
        fullDataframe.drop(empty_cols, axis=1, inplace=True)
        fullDataframe = fullDataframe.reset_index(drop=True)
        return fullDataframe

    if tableType == 'full':
        if webSummaryType == 'GEX':
            full_df = fullTableMaker(full_df)
            return full_df
        else:
            # delete columns like "Sample Description" if empty/full of NaNs
            ATAC_full_dframe = fullTableMaker(ATAC_full_df)
            GEX_full_dframe = fullTableMaker(GEX_full_df)
            return ATAC_full_dframe, GEX_full_dframe

    elif tableType == 'delivery doc':
        if webSummaryType == 'GEX':
            deliveryHeaders = ['sample id', 'estimated number of cells', 'mean reads per cell', 'median genes per cell',
                               'number of reads', 'sequencing saturation', 'reads mapped to genome', 'reads mapped confidently to genome',
                               'fraction reads in cells', 'median umi counts per cell', 'tso_frac']
        elif webSummaryType == 'ARC':
            raise Exception(
                "delivery doc not yet implemented for ARC pipelines"
            )
        for col in full_df.columns:
            if col.lower() not in deliveryHeaders:
                del full_df[col]
        full_df = full_df.reset_index(drop=True)
        return full_df

    elif tableType == 'repooling':
        if webSummaryType == 'ARC':
            raise Exception("repooling table not yet implemented for ARC pipelines")

        repoolingHeaders = ['sample id', 'estimated number of cells',
                            'mean reads per cell', 'number of reads']
        for col in full_df.columns:
            if col.lower() not in repoolingHeaders:
                del full_df[col]
        full_df.columns = map(str.lower, full_df.columns)
        full_df['estimated number of cells'] = full_df['estimated number of cells'].str.replace(
            ',', '').astype(int)
        full_df['mean reads per cell'] = full_df['mean reads per cell'].str.replace(
            ',', '').astype(int)
        assert full_df['mean reads per cell'].max(
        ) < readsDesired, "Unable to repool if readsDesired is lower than mean reads/cell of one of the considered samples"
        full_df['number of reads'] = full_df['number of reads'].str.replace(',', '').astype(int)
        full_df['reads needed for ' +
                str(readsDesired) + ' reads per cell'] = readsDesired - full_df['mean reads per cell']
        full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'] = full_df['estimated number of cells'] * \
            full_df['reads needed for ' + str(readsDesired) + ' reads per cell']
        full_df['percent of lane'] = full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'] / sum(
            full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'])
        full_df = full_df.reset_index(drop=True)
        return full_df
    else:
        raise Exception(
            "For tableType argument, please use one of the following: full, delivery doc, repooling")
