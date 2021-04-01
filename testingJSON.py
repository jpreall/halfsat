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
def attributeScrapper(html):
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
pd.DataFrame([attributeScrapper(camilaList[0])])

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
