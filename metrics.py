import re
import json
from bs4 import BeautifulSoup
import pandas as pd

def attributeScraper(html):
    """
    Generate a dictionary of all attributes and values from a web_summary.html file.

    Scrape attributes from cellranger v3-6 web_summary.html to generate a dictionary of
    attributes and values from all metrics in the json portion of web_summary.html.

    Args:
        html (str): Path to web_summary.html file of desired experiment/sample.

    Returns:
        dict: Dictionary of attribute (key) and corresponding value.

    Raises:
    """

    f = open(html, encoding="utf8")
    soup = BeautifulSoup(f)
    for line in soup.find('script'):
        if 'const data' in line:
            const_data = line
    f.close()

    # need to include +1 to include end of string (the last })
    constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n')+1])

    # extract data from JSON tables
    pipeline_table = constant_data['summary']['summary_tab']['pipeline_info_table']['rows']
    seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
    map_summary_table = constant_data['summary']['summary_tab']['mapping']['table']['rows']
    cells_summary_table = constant_data['summary']['summary_tab']['cells']['table']['rows']

    tableList = [pipeline_table, seq_summary_table, map_summary_table, cells_summary_table]

    # create dictionary of all attributes
    attributeDict = {}
    for table in tableList:
        for entry in table:
            if entry[0] not in attributeDict:
                attributeDict[entry[0]] = entry[1]

    # extract data from JSON diagnostics dictionary if exist
    if 'diagnostics' in constant_data['summary']:
        diagnostics_dict = constant_data['summary']['diagnostics']
        attributeDict.update(diagnostics_dict)

    return attributeDict

def tableGenerator(htmlList, tableType='full', readsDesired=40000):
    """
    Use list of cellranger (v4+) count web_summary.html () to generate a table.

    Use list of web_summary.html to generate a full table of all cellranger count metrics,
    abbreviated table for data delivery, or a table for repooling purposes.

    Args:
        htmlList (list): list of web_summary.html directories.
        tableType (string): string of one of the following: full, data delivery, repooling.
        readsDesired (int): number of reads per cell desired for repooling.

    Returns:
        pandas DataFrame: DataFrame for table specified.

    Raises:
        AssertionError: htmlList must be of list type.
        Exception: For tableType argument, please use one of the following: full, delivery doc, repooling.

    """

    assert isinstance(htmlList, list), "Please input a list of html directories"

    initial = True
    # Create dataframe containing all of the web_summaries in list
    for sample in htmlList:
        sample_dict = attributeScraper(sample)
        # https://stackoverflow.com/questions/57631895/dictionary-to-dataframe-error-if-using-all-scalar-values-you-must-pass-an-ind
        sample_df = pd.DataFrame([sample_dict])
        if initial == True:
            full_df = sample_df
            initial = False
        else:
            full_df = full_df.append(sample_df)

    if tableType == 'full':
        # delete columns like "Sample Description" if empty/full of NaNs
        nan_value = float("NaN")
        full_df.replace("", nan_value, inplace=True)
        # reference: https://www.jitsejan.com/find-and-delete-empty-columns-pandas-dataframe.html
        empty_cols = [col for col in full_df.columns if full_df[col].isnull().all()]
        full_df.drop(empty_cols, axis=1, inplace=True)
        full_df = full_df.reset_index(drop=True)
        return full_df

    elif tableType == 'delivery doc':
        deliveryHeaders = ['sample id', 'estimated number of cells', 'mean reads per cell', 'median genes per cell',
                           'number of reads', 'sequencing saturation', 'reads mapped to genome', 'reads mapped confidently to genome',
                           'fraction reads in cells', 'median umi counts per cell', 'tso_frac']
        for col in full_df.columns:
            if col.lower() not in deliveryHeaders:
                del full_df[col]
        full_df = full_df.reset_index(drop=True)
        return full_df

    elif tableType == 'repooling':
        repoolingHeaders = ['sample id', 'estimated number of cells', 'mean reads per cell', 'number of reads']
        for col in full_df.columns:
            if col.lower() not in repoolingHeaders:
                del full_df[col]
        full_df.columns = map(str.lower, full_df.columns)
        full_df['estimated number of cells'] = full_df['estimated number of cells'].str.replace(',', '').astype(int)
        full_df['mean reads per cell'] = full_df['mean reads per cell'].str.replace(',', '').astype(int)
        full_df['number of reads'] = full_df['number of reads'].str.replace(',', '').astype(int)
        full_df['reads needed for ' + str(readsDesired) + ' reads per cell'] = readsDesired - full_df['mean reads per cell']
        full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'] = full_df['estimated number of cells'] * full_df['reads needed for ' + str(readsDesired) + ' reads per cell']
        full_df['percent of lane'] = full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'] / sum(full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'])
        full_df = full_df.reset_index(drop=True)
        return full_df
    else:
        raise Exception("For tableType argument, please use one of the following: full, delivery doc, repooling")
