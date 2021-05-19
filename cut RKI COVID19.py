import numpy as np
import pandas as pd
from projectlib import region_setup

candidate=pd.read_csv('External Data/RKI_COVID19.csv', sep=',', header='infer')
chosen=region_setup(12)[0]
candidate=candidate.loc[candidate['IdLandkreis'].isin(chosen)]
candidate.to_csv("External Data/RKI_COVID19_cut_Region12.csv")