import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib as ppp
import warnings
from matplotlib.dates import DateFormatter, AutoDateLocator

srp_df1 = pd.read_csv('datalog_parsed_first.csv')
srp_df2 = pd.read_csv('datalog_parsed_second.csv')
df_total = pd.concat([srp_df1,srp_df2])

df_total.sort_values(['DOM','time']).to_csv('merged_srp.csv',index=1)

