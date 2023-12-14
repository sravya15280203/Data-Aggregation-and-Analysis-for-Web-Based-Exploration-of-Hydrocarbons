from rdkit import Chem
import pandas as pd
import csv
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

df=pd.read_csv("Hydrocarbons.csv", )
df.drop("serial No.",inplace=True,axis=1)
'''''
mol1= Chem.MolFromSmiles("c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c(c6ccc7c2ccc8c3ccc9c4ccc%10)ccc5")
print(mol1)
formula1 = CalcMolFormula(mol1)
print(formula1)
'''
mol = Chem.MolFromSmiles("C1=C2C=CC3=CC=C4C=C5C6=C7C(=CC8=C9C%10=C%11C(=CC=C%10C=C8)C=CC%12=C%11C(=C79)C(=C6)C=C%12)C%13=CC%14=C%15C(=C2C3=C4C%15=C5%13)C(=C1)C=C%14")
formula = CalcMolFormula(mol)
print(formula)

smile_list=[]
#print(df.head(5),df.tail(5))
#x=0
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
        None
        '''''
        print('haha')
        print(sm)
        print(inin)
        mol=Chem.MolFromSmiles(sm)
        mol=Chem.MolFromSmiles(mol)
        sms=CalcMolFormula(mol)'''
print(df.loc[[7729]])
print(df.loc[[7730]])
print(df.loc[[7731]])
smile_list.insert(7729,"C48H20")
df.insert(1,column="Formula",value=smile_list)
df.to_csv("HYDROCARBON_2.csv")
print("complete")
#print(df.head(),df.tail(9))

'''
#print("index:0",smile_list[0])
#print("index:last",smile_list[-1],len(smile_list))
#print(df.columns)
#print(df.head(5),df.tail(5))
#df["Formula"]= CalcMolFormula(Chem.MolFromSmiles(smile) for smile in df.smiles)
#print(df)
# e.g. cysteine

columns = ["smiles", "formula"]
with open("Hydrocarbons.csv","r") as in_file:
    with open ("HYDROCARBONS_2.csv","w") as out_file:
        csv_writer = csv.writer(out_file)
        csv_writer.writerow(columns)
        for row in csv.reader(in_file):
            try:
                csv_writer.writerow(row)
            except:
                print("error")
print("done")
'''