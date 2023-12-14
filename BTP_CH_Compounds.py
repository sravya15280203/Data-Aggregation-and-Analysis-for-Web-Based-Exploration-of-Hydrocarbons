from asyncio.windows_events import NULL
import pandas as pd
import csv
import matplotlib.pyplot as plt
import re
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import Chem
import os

class HydroCarbon:
    def MolWeight(self):
        df = pd.read_csv("D:\Bachelor-Thesis-Project\Merging_CH_file.csv")
        #print(df.head(5))
        df["Mol_wt"] = 0
        Final_MolList=[]
        for mol in df.Formula:
            num_list = {"elementC": 0, "elementH": 0}
            i = 0
            j = 0
            n = len(mol)
            carry = ""
            while (i < n or j < n):
                if (j < n and mol[j].isnumeric()):
                    if (j != n):
                        carry += mol[j]
                        j += 1
                elif (j == n or mol[j].isalpha() or mol[j]=='-' or mol[j]=='+'):
                    if (i < j and mol[i].isalpha()):
                        val = 0
                        try:
                            val = int(carry)
                        except Exception as e:
                            val=1
                        if(mol[i].upper()=='C'):
                            num_list["elementC"]+=val
                        elif(mol[i].upper()=='H'):
                            num_list["elementH"]+=val
                        else:
                            pass
                        carry=""
                    i=j
                    j+=1
                else:
                    i=j
                    j+=1
            #Final_MolList.append(num_list["elementC"]*12.02 + num_list["elementH"]*1.007)#---> mol_wt

        df["Mol_wt"]=Final_MolList #----> declaring elements to mol_wt
        #print(df.to_string())    
        #df.drop(["Unnamed: 0"],inplace=True, axis=1)
        df.sort_values(by=["Mol_wt"],inplace=True)
        #df.to_csv("D:\Bachelor-Thesis-Project\Merging_CH_file.csv",index=False) 
        
    def Histogram(self):
        df = pd.read_csv("D:\Bachelor-Thesis-Project\qm9_Hydrocarbons_2.csv")
        Final_MolList=[]
        for mol in df.Formula:
            num_list = {"elementC": 0, "elementH": 0}
            i = 0
            j = 0
            n = len(mol)
            carry = ""
            while (i < n or j < n):
                if (j < n and mol[j].isnumeric()):
                    if (j != n):
                        carry += mol[j]
                        j += 1
                elif (j == n or mol[j].isalpha() or mol[j]=='-' or mol[j]=='+'):
                    if (i < j and mol[i].isalpha()):
                        val = 0
                        try:
                            val = int(carry)
                        except Exception as e:
                            val=1
                        if(mol[i].upper()=='C'):
                            num_list["elementC"]+=val
                        elif(mol[i].upper()=='H'):
                            num_list["elementH"]+=val
                        else:
                            pass
                        carry=""
                    i=j
                    j+=1
                else:
                    i=j
                    j+=1
        
            if num_list["elementC"]<=10 and num_list["elementH"]<=12:#---> condition for histrogram
                Final_MolList.append(num_list["elementC"]*12.02 + num_list["elementH"]*1.007)
        plt.title("PubChem-HYD")
        plt.hist(Final_MolList)
        plt.xlabel("Molecular weight")
        plt.ylabel("No of Molecules")
        #plt.savefig("PyAr.png",transparent=True)
        plt.show()

    def log_file_extraction(file_path):
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
        return(res_list)

    def calling_log_file_extraction(path):
        res_dict=HydroCarbon.log_file_extraction(path)
        df= pd.read_csv("C:/Users/smile/bulk_csv_upgrade/c10/c10h9_fresh/c10h9_fresh.csv")
        df["energy"]=NULL

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
    
    def formula_convo(self):
        df=pd.read_csv("D:\Bachelor-Thesis-Project\Merging_CH_file.csv")
        smile_list=[]

        index = 0
        for smi in df.smiles_format:
            index += 1
            try:
                mol=Chem.MolFromSmiles(smi)
                mol_form=CalcMolFormula(mol)
                
                smile_list.append(mol_form)
                #if x<=5:
                #   print(x,mol_form)
                #x=x+1
            except:
                print(smi)
                print(index)
        #smile_list.insert(7729,"C48H20")--->for hydrocardon2 file one formula goes under exception on index 7729 
        df.insert(4,column="Formula",value=smile_list)
        #print(df.head())
        df.to_csv("D:\Bachelor-Thesis-Project\Merging_CH_file.csv",index=False)
        #print("complete")

    def dat_to_csv(self):
        columns = ["smiles", "serial No."]
        #in_file=open("filepath", 'r')
        #out_file=open("file", 'w')
        with open("C:\\Users\\smile\\Downloads\\chemgiri-pubchem_data\\pubchem.tar\\pubchem\\pubchem.dat", 'r') as in_file:
            with open("pubchem.csv", 'w') as out_file:
                csv_writer = csv.writer(out_file)
                csv_writer.writerow(columns)
                for row in in_file:
                    row_values = [row_val.strip()
                                  for row_val in row.split(',')]
                    csv_writer.writerow(row_values)
        return
    def ch_compounds(self):#---> for 1 file 
        columns = ["smiles", "serial No."]
        with open("pubchem.csv","r") as in_file:
            with open ("Hydrocarbons.csv","w") as out_file:
                csv_writer = csv.writer(out_file)
                csv_writer.writerow(columns)
                for row in csv.reader(in_file):
                    #res = ''.join([i for i in row if not i.isdigit()])
                    #row_line=[val.strip() for val in row.split(',')]
                    #row_values = [row_val.strip()
                                 # for row_val in row.split(',')]
                    try :
                        s1="".join(c for c in row[0] if c.isalpha())
                        #print(s1)
                        if(all(c=="c" or c=="h" or c=="C" or c=="H" for  c in s1)):
                            csv_writer.writerow(row)
                        else:
                            pass   
                    except:
                        pass
                out_file.close()
            in_file.close()         
    def qm9_ch_compounds(self):# for another file
        columns = ["mol_id","smiles","A","B","C","mu","alpha","homo","lumo","gap","r2","zpve","u0","u298","h298","g298","cv","u0_atom","u298_atom","h298_atom","g298_atom"]
        with open("qm9 1.csv","r") as in_file:
            with open ("qm9_Hydrocarbons.csv","w") as out_file:
                csv_writer = csv.writer(out_file)
                csv_writer.writerow(columns)
                for row in csv.reader(in_file):
                    #res = ''.join([i for i in row if not i.isdigit()])
                    #row_line=[val.strip() for val in row.split(',')]
                    #row_values = [row_val.strip()
                                 # for row_val in row.split(',')]
                    try :
                        s1="".join(c for c in row[1] if c.isalpha())
                        #print(s1)
                        if(all(c=="c" or c=="h" or c=="C" or c=="H" for  c in s1)):
                            csv_writer.writerow(row)
                        else:
                            pass   
                    except:
                        pass
                out_file.close()
            in_file.close() 
    
    def merging(self,path):
        
        for i in range(10):
            y=r'D:\bulk_csv_upgrade\combo.csv'
            x="C:\\Users\\smile\\bulk_csv_upgrade\\c{c}\\c{c}h{h}_fresh\\c{c}h{h}_fresh.csv".format(c=10,h=10)# for histogram merge
            z=path #for a,b,c, enthropy...etc merge
        df=pd.concat(map(pd.read_csv,[y,z]),ignore_index=True)
        df.to_csv(y,index=False)
        '''
        df1=pd.read_csv("D:\Bachelor-Thesis-Project\Merging_CH_file.csv")
        df1.drop(["energy1","rem"],inplace=True,axis=1)
        print(df1)
        df1.to_csv("D:\Bachelor-Thesis-Project\Merging_CH_file.csv",index=False)
    
       '''
    def out_file_extraction(self,path):
        global column_name , list_we_need
        out_path=path
        column_name=re.search("fresh.c.*out",out_path)
        column_name=column_name.group().split("\\")
        column_name=column_name[1]
        column_name=column_name.replace(".out",".xyz")

        with open(out_path,"r",encoding = 'utf-8') as input_file:
            content =input_file.readlines()
            list_we_need={}
        
            count=0
            if (" ****ORCA TERMINATED NORMALLY****" in content[-2]):
                for i in content:
                    word="Total Dipole Moment "
                    word1="Rotational spectrum "
                    word2="Total thermal energy"
                    word3="OPTIMIZATION RUN DONE"
                    word4="Total enthalpy"
                    word5="Final entropy term "
                    word6="Final Gibbs free energy"
                    if word in i:    
                        num=re.findall("-?\d+\.\d+",content[count+3])
                    if word1 in i: 
                        num1=re.findall("-?\d+\.\d+",content[count+3])
                    if word2 in i:
                        num2=re.findall("-?\d+\.\d+",content[count])
                    if word3 in i:
                        num3=re.findall("-?\d+\.\d+",content[count-3])
                    if word4 in i:
                        num4=re.findall("-?\d+\.\d+",content[count])
                    if word5 in i:
                        num5= re.findall("-?\d+\.\d+",content[count])
                    if word6 in i:
                        num6=re.findall("-?\d+\.\d+",content[count])    

                    count+=1

                list_we_need["Dipole moment"]=float(num[0])
                list_we_need["A"]=float(num1[0])
                list_we_need["B"]=float(num1[1])
                list_we_need["C"]=float(num1[2])
                list_we_need["Thermal energy"]=float(num2[0])
                list_we_need["Final single point energy"]=float(num3[0])
                list_we_need["Total Enthalpy"]=float(num4[0])
                list_we_need["Entropy"]=float(num5[0])
                list_we_need["Gibbs free energy"]=float(num6[0])
                #print(list_we_need,"\n")
                
            else:
                print("error")
                exit()
    def inserting(self,path1):
        
        df=pd.read_csv(path1)
        #print(column_name)
        
        for i in df.xyz_file:
            if i==column_name:
                df.loc[df.xyz_file==i,["Dipole moment",'A','B','C',"Final single point energy","Thermal energy","Total Enthalpy","Entropy","Gibbs free energy"]]=list_we_need["Dipole moment"],list_we_need["A"],list_we_need["B"],list_we_need["C"],list_we_need["Final single point energy"],list_we_need["Thermal energy"],list_we_need["Total Enthalpy"],list_we_need["Entropy"],list_we_need["Gibbs free energy"]
        #print(df.head())
        df.to_csv(path1,index=False)  
    def path(self):
        global individual_path
        dir_path=[]
        file_path=[]
        out_file_path=[]
        for dirpath,dirnames,filepath in os.walk('D:\\bulk_csv_upgrade'):
            #print(dirnames)
            dir_path.append(dirpath)
            file_path.append(filepath)
            #print(filenames)
    
        for i in file_path:
            for j in i:
                if j.endswith(r'.out'):
                    out_file_path.append(j)
        individual_path=["D:\\bulk_csv_upgrade\c10\c10h9_fresh\{}".format(i) for i in out_file_path]
        
ob1=HydroCarbon()
#ob1.path()
ob1.Histogram()

'''
for i in range(2,11):
    try:
        path="D:\\bulk_csv_upgrade\c10\c10h{}_fresh\\fresh.csv".format(i)
        ob1.merging(path)
    except:
        pass
'''

'''
 #  getting paths through ob1.path for so many file paths

path1=r'D:\bulk_csv_upgrade\c10\c10h9_fresh\fresh.csv'
for particular_path in individual_path:
    try:
        ob1.out_file_extraction(particular_path)
        ob1.inserting(path1)
    except:
        pass
'''
#ob1.merging(path)
#ob1.formula_convo()
#ob1.MolWeight()
