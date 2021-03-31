import sys, os
sys.path
sys.path.append('C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat')
import halfSat_webSum


def attributeScraper(sample, repooling = False):
    """
    A short description.

    A bit longer description.

    Args:
        variable (type): description

    Returns:
        type: description

    Raises:
        Exception: description

    """

    import re
    from bs4 import BeautifulSoup
    import numpy as np
    import pandas as pd
    dirname = os.path.dirname(sample)
    if not os.path.exists(sample):
        os.makedirs(sample)
    print(dirname)

    f = open(sample, encoding="utf8")
    soup = BeautifulSoup(f)
    f.close()
    data = soup.find_all(['script'])[0].text

    # Make dictionary
    sampleDict = {}

    # Sample ID
    pattern_sampID = re.compile('(?i)(\"Sample ID\", \"[\w&.\-]+\")')
    hits_sampID = pattern_sampID.findall(data)
    #hits_sampID
    sampID = hits_sampID[0].replace('"','').split(", ")
    #sampID
    sampleDict[sampID[0]] = sampID[1]
    #sampleDict

    # Estimated Number of cells
    pattern_estNumCells = re.compile('(?i)(\"Estimated Number of Cells\", \"\d{1,3}(?:,\d{3})*\")')
    hits_estNumCells = pattern_estNumCells.findall(data)
    #hits_estNumCells
    estNumCells = hits_estNumCells[0].replace('"','').split(", ")
    #estNumCells
    sampleDict[estNumCells[0]] = estNumCells[1]
    #sampleDict

    # Mean Reads per Cell
    pattern_meanReadsCell = re.compile('(?i)(\"Mean Reads per Cell\", \"\d{1,3}(?:,\d{3})*\")')
    hits_meanReadsCell = pattern_meanReadsCell.findall(data)
    #hits_meanReadsCell
    meanReadsCell = hits_meanReadsCell[0].replace('"','').split(", ")
    #meanReadsCell
    sampleDict[meanReadsCell[0]] = meanReadsCell[1]
    #sampleDict

    # Median Genes per Cell
    pattern_medGenesCell = re.compile('(?i)(\"Median Genes per Cell\", \"\d{1,3}(?:,\d{3})*\")')
    hits_medGenesCell = pattern_medGenesCell.findall(data)
    #hits_medGenesCell
    medGenesCell = hits_medGenesCell[0].replace('"','').split(", ")
    #medGenesCell
    if repooling == False:
        sampleDict[medGenesCell[0]] = medGenesCell[1]
        #sampleDict

    # Number of Reads
    pattern_numReads = re.compile('(?i)(\"Number of Reads\", \"\d{1,3}(?:,\d{3})*\")')
    hits_numReads = pattern_numReads.findall(data)
    #hits_numReads
    numReads = hits_numReads[0].replace('"','').split(", ")
    #numReads
    sampleDict[numReads[0]] = numReads[1]
    #sampleDict

    # Sequencing Saturation
    pattern_seqSat = re.compile('(?i)(\"Sequencing Saturation\", \"\d{1,2}(?:.\d{1,2})%\")')
    hits_seqSat = pattern_seqSat.findall(data)
    #hits_seqSat
    seqSat = hits_seqSat[0].replace('"','').split(", ")
    #seqSat
    sampleDict[seqSat[0]] = seqSat[1]
    #sampleDict

    # Reads Mapped to Genome
    pattern_readsMapGenome = re.compile('(?i)(\"Reads Mapped to Genome\", \"\d{1,2}(?:.\d{1,2})%\")')
    hits_readsMapGenome = pattern_readsMapGenome.findall(data)
    #hits_readsMapGenome
    readsMapGenome = hits_readsMapGenome[0].replace('"','').split(", ")
    #readsMapGenome
    if repooling == False:
        sampleDict[readsMapGenome[0]] = readsMapGenome[1]
        #sampleDict

    # Reads Mapped Confidently to Genome
    pattern_readsMapConfGenome = re.compile('(?i)(\"Reads Mapped Confidently to Genome\", \"\d{1,2}(?:.\d{1,2})%\")')
    hits_readsMapConfGenome = pattern_readsMapConfGenome.findall(data)
    #hits_readsMapConfGenome
    readsMapConfGenome = hits_readsMapConfGenome[0].replace('"','').split(", ")
    #readsMapConfGenome
    if repooling == False:
        sampleDict[readsMapConfGenome[0]] = readsMapConfGenome[1]
        #sampleDict

    # Fraction Reads in cells
    pattern_fracReadsCells = re.compile('(?i)(\"Fraction Reads in Cells\", \"\d{1,2}(?:.\d{1,2})%\")')
    hits_fracReadsCells = pattern_fracReadsCells.findall(data)
    #hits_fracReadsCells
    fracReadsCells = hits_fracReadsCells[0].replace('"','').split(", ")
    #fracReadsCells
    if repooling == False:
        sampleDict[fracReadsCells[0]] = fracReadsCells[1]
        #sampleDict

    # Median UMI Counts per Cell
    pattern_medUMI = re.compile('(?i)(\"Median UMI Counts per Cell\", \"\d{1,3}(?:,\d{3})*\")')
    hits_medUMI = pattern_medUMI.findall(data)
    #hits_medUMI
    medUMI = hits_medUMI[0].replace('"','').split(", ")
    #medUMI
    if repooling == False:
        sampleDict[medUMI[0]] = medUMI[1]
        #sampleDict

    return sampleDict

def deliveryDocMaker(listSamp):
    import pandas as pd
    initial = True

    for sample in listSamp:
        sample_dict = attributeScraper(sample, repooling=False)
        # https://stackoverflow.com/questions/57631895/dictionary-to-dataframe-error-if-using-all-scalar-values-you-must-pass-an-ind
        sample_df = pd.DataFrame([sample_dict])
        if initial == True:
            full_df = sample_df
            initial = False
        else:
            full_df = full_df.append(sample_df)

    return full_df
foolist = ["C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/BY01_86T_web.html", "C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/BY01_1247_web.html"]
deliveryDocMaker(foolist)

def repoolingTable(listSamp, readsDesired = 40000):
    import pandas as pd
    initial = True

    for sample in listSamp:
        sample_dict = attributeScraper(sample, repooling=True)
        # https://stackoverflow.com/questions/57631895/dictionary-to-dataframe-error-if-using-all-scalar-values-you-must-pass-an-ind
        sample_df = pd.DataFrame([sample_dict])
        if initial == True:
            full_df = sample_df
            initial = False
        else:
            full_df = full_df.append(sample_df)
    full_df.columns = map(str.lower, full_df.columns)
    full_df['estimated number of cells'] = full_df['estimated number of cells'].str.replace(',', '').astype(int)
    full_df['mean reads per cell'] = full_df['mean reads per cell'].str.replace(',', '').astype(int)
    full_df['number of reads'] = full_df['number of reads'].str.replace(',', '').astype(int)
    full_df['sequencing saturation'] = full_df['sequencing saturation'].str.replace(r'%', r'').astype('float') / 100.0
    full_df['reads needed for ' + str(readsDesired) + ' reads per cell'] = readsDesired - full_df['mean reads per cell']
    full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'] = full_df['estimated number of cells'] * full_df['reads needed for ' + str(readsDesired) + ' reads per cell']
    full_df['percent of next lane'] = full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'] / sum(full_df['total reads needed for ' + str(readsDesired) + ' reads per cell'])
    return full_df

repoolingTable(foolist)

foo = attributeScraper(foolist[0])
foo.keys()
