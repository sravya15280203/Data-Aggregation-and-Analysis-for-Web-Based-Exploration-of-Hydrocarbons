
import sqlite3
import csv

class molecules():
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
    def ch_compounds(self):
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
ob1=molecules()
#ob1.dat_to_csv()
ob1.ch_compounds()
print("completed")
