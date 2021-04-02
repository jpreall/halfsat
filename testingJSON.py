import sys, os
sys.path
sys.path.append('C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat')
import halfSat_webSum

camilaList = ['C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/CC01_1314HRL_web.html',
               'C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/CC01_451NH_web.html',
               'C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/CC01_452HR_web.html',
               'C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/CC01_955H2R_web.html']
import re
import json
from bs4 import BeautifulSoup
import pandas as pd
f = open(camilaList[0], encoding="utf8")
soup = BeautifulSoup(f)
for line in soup.find('script'):
    if 'const data' in line:
        const_data = line
        print(type(const_data))
        print(const_data)
f.close()

constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n')+1]) #need to include +1 to include end of string (the last })
const_data[const_data.find('{'):const_data.find('}\n')]
constant_data

pipeline_table = constant_data['summary']['summary_tab']['pipeline_info_table']['rows']
pipeline_table

seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
seq_summary_table

map_summary_table = constant_data['summary']['summary_tab']['mapping']['table']['rows']
map_summary_table

cells_summary_table = constant_data['summary']['summary_tab']['cells']['table']['rows']
cells_summary_table

attributeDict = {}
def attributeScraper(html):
    import re
    import json
    from bs4 import BeautifulSoup
    import pandas as pd
    f = open(html, encoding="utf8")
    soup = BeautifulSoup(f)
    for line in soup.find('script'):
        if 'const data' in line:
            const_data = line
    f.close()

    constant_data = json.loads(const_data[const_data.find('{'):const_data.find('}\n')+1]) #need to include +1 to include end of string (the last })
    const_data[const_data.find('{'):const_data.find('}\n')]
    constant_data

    pipeline_table = constant_data['summary']['summary_tab']['pipeline_info_table']['rows']
    pipeline_table

    seq_summary_table = constant_data['summary']['summary_tab']['sequencing']['table']['rows']
    seq_summary_table

    map_summary_table = constant_data['summary']['summary_tab']['mapping']['table']['rows']
    map_summary_table

    cells_summary_table = constant_data['summary']['summary_tab']['cells']['table']['rows']
    cells_summary_table

    tableList = [pipeline_table, seq_summary_table, map_summary_table, cells_summary_table]

    attributeDict = {}
    for table in tableList:
        for entry in table:
            if entry[0] not in attributeDict:
                attributeDict[entry[0]] = entry[1]
    return attributeDict

foo = pd.DataFrame([attributeScraper(camilaList[0])])
foo
del foo['Sample Description']
foo.columns[0].lower()
foo['estimated number of cells']
def tableGenerator(htmlList, tableType='full', readsDesired=40000):
    import pandas as pd
    initial = True

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
        return full_df

    elif tableType == 'delivery doc':
        deliveryHeaders = ['sample id', 'estimated number of cells', 'mean reads per cell', 'median genes per cell',
                           'number of reads', 'sequencing saturation', 'reads mapped to genome', 'reads mapped confidently to genome',
                           'fraction reads in cells', 'median umi counts per cell']
        for col in full_df.columns:
            if col.lower() not in deliveryHeaders:
                del full_df[col]
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
        return full_df
tableGenerator(camilaList, tableType='repooling')

for entry in pipeline_table:
    if entry[0].lower() in 'sample id': #check if it is sample id entry
        sampleID = entry
        attributeDict[sampleID[0]] = sampleID[1]


pattern_sampID = re.compile('(?i)(\"Sample ID\", \"[\w&.\-]+\")')
hits_sampID = pattern_sampID.findall(pipeline_table[0])
pipeline_table

if pipeline_table[0][0].lower() in 'sample id':
    pipeline_table[0]
sampleID = pipeline_table[entry].replace('"','').split(", ")

pd.DataFrame(seq_summary_table)





constant_data['summary']['summary_tab']['pipeline_info_table']['rows']
for key in constant_data['summary']['summary_tab']:
    print(key)
