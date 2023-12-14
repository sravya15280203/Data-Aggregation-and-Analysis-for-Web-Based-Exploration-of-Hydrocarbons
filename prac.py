#import numpy as np
from asyncio.windows_events import NULL
from calendar import c
import pandas as pd


'''
df= pd.read_csv("C:/Users/smile/bulk_csv_upgrade/c10/c10h9_fresh/c10h9_fresh.csv")
df1=np.loadtxt
df.head()
print(df.head(5))
'''

def file_extraction(file_path):
        res_list={}
        with open(file_path,'r') as in_file:
            i=0
            for row in in_file:
                if(row=="QM Software:   orca\n"):
                    i=1
                elif(i==1):
                    row_values = [row_val.strip()
                                  for row_val in row.split(':')]
                    res_list[row_values[0]]=row_values[1]
        return res_list

res_dict=file_extraction("C:/Users/smile/bulk_csv_upgrade/c10/c10h9_fresh/pyar-optimiser.log")

df= pd.read_csv("C:/Users/smile/bulk_csv_upgrade/c10/c10h9_fresh/c10h9_fresh.csv")
df["energy"]=NULL
#print(df["xyz_file"])
#print(df.loc[1,"xyz_file"])
for key in res_dict:
    for x in df.index:
        if(key+".xyz" ==df.loc[x,"xyz_file"]):
            #print("s")
            df.loc[x,"energy"]=res_dict[key]

#print (df.to_string())
#print(df.head())
#df.drop(["sl.no"],axis=0,inplace=True)
#df.reset_index(drop=True,inplace=True)
#df.to_csv("C:/Users/smile/bulk_csv_upgrade/c10/c10h9_fresh/c10h9_fresh.csv",index=False)
#df=pd.read_csv("C:/Users/smile/bulk_csv_upgrade/c10/c10h9_fresh/c10h9_fresh.csv").drop("Unnamed: 0", axis=1,inplace=True)