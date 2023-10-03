import sys
sys.path.append('/home/user/python-libs')
print (sys.path)

import numpy as np
import pandas as pd
import scanpy as sc

gut = sc.read("gut.h5ad") #epithelial cells only, from https://www.gutcellatlas.org/

gut.obs.to_csv("gut_raw.csv") 
