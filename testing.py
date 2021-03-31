import sys, os
sys.path
sys.path.append('C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat')
import halfSat_webSum


import re
from bs4 import BeautifulSoup
import numpy as np
#import json
#import pandas as pd

filename = 'C:/Users/bhbri/AppData/Local/Temp/6e4e7afe/bamdev1/grid/preall/home/bhe/pipelineWork/halfsat/BY01_86T_web.html'
dirname = os.path.dirname(filename)
if not os.path.exists(dirname):
    os.makedirs(dirname)
print(dirname)

f = open(filename, encoding="utf8")
soup = BeautifulSoup(f)
print(soup.prettify())
f.close()

data = soup.find_all(['script'])[0].text
data
