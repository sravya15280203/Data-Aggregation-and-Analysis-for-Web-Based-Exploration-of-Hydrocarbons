from rdkit import Chem
import pandas as pd
import csv
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

df=pd.read_csv("qm9_Hydrocarbons.csv")
smile_list=[]

inin = 0
for sm in df.smiles:
    inin += 1
    try:
        mol=Chem.MolFromSmiles(sm)
        sms=CalcMolFormula(mol)
        
        smile_list.append(sms)
        #if x<=5:
        #   print(x,sms)
        #x=x+1
    except:
        print('haha')
        print(sm)
        print(inin)
df.insert(1,column="Formula",value=smile_list)
df.to_csv("qm9_Hydrocarbons_2.csv")
print("complete")