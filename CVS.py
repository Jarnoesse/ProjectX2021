#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import pandas as pd
import sys
import os
import random
#checks for arguments
if len(sys.argv) <= 1:
    print("Please provide FASTA file")
    quit()

#help message
if sys.argv[1] == '-h' or sys.argv[1] == 'help' or sys.argv[1] == '-help':
    print("Usage: CVS.py inputfile.csv outputfile.csv")
    quit()

#INPUT
input = sys.argv[1]
#verifies if the input file exists
if not os.path.exists(input):
    print('Error: file %s does not exist!', input)
    quit()

#verifies if the output was passed in the arguments
output_csv = 'outputCVS.csv'
if len(sys.argv) > 2:
    output = sys.argv[2]

print("--START--")
R =  800 #how many starts
T = 200 #how long the path
#class definition
class Cvs:
    
    #initializer 
    def __init__(self,name,dataframe,R,T,rejects=50,output = output_csv):
        #class name
        self.name      = name
        #class starting dataframe
        self.dataframe = pd.read_csv(dataframe)
        #defining R and T
        self.R = R
        self.T = T
        #reject limit for gradient ascent
        self.rejects_limit = rejects
        #output file
        self.output_csv = output
    
    #this method initializes the indexes and creates the secondary dataframe
    def initialize_indexes(self):
        #randomly selects the first indexes
        self.indexes       = sorted(random.sample(range(len(self.dataframe.columns)),self.number),key = int)
        #create the complementary list of indexes
        self.other_indexes = list(range(0,len(self.dataframe.columns)))
        for element in self.indexes:
            if element in self.other_indexes:
                self.other_indexes.remove(element)
        #creates the secondary dataframe, dropping the NaNs
        self.smaller_dataframe = self.dataframe.iloc[:,self.indexes]
        self.smaller_dataframe = self.smaller_dataframe.dropna()
        self.smaller_dataframe['Combined'] = self.smaller_dataframe.values.tolist()
    
    #this method calculates the parameters
    def calculate_parameters(self):
        self.k = self.smaller_dataframe['Combined'].value_counts()
        self.m = self.k.value_counts()
        self.M = self.k.sum()
    
    #this method calculates Hk
    def calculate_hk(self):
        self.hk = 0
        for element in list(map(list,self.m.items())):
            self.hk += ( ( element[0] * element[1] ) / self.M) * np.log(( element[0] * element[1] ) / self.M)
        self.hk = -self.hk
    
    #this method calculates Hs
    def calculate_hs(self):
        self.hs = 0
        for element in list(map(list,self.m.items())):
            self.hs += ( ( element[0]) / self.M) * np.log(( element[0]) / self.M)
        self.hs = -self.hs
    
    #this method swaps one element from the indexes list with one from the complementary list
    def swap(self):
        condition = False
        while (condition != True):
            i = random.randint(0,len(self.indexes)-1)
            j = random.randint(0,len(self.other_indexes)-1)
            swap = self.indexes[i]
            self.indexes[i] = self.other_indexes[j]
            self.other_indexes[j] = swap
            self.indexes = sorted(self.indexes,key = int)
            
            print(self.indexes)
            #it now controls if the indexes have never been used before
            if len(self.history.iloc[:,0]) <1:
                condition = True
            else:
                condition = True
                for element in self.history.iloc[:,self.number].to_numpy():
                    if element == self.indexes:
                        print(self.history.iloc[:,self.number].to_numpy(),"lista ",element,"uguale")
                        condition = False
        #updates the secondary dataframe
        temporary_df = pd.DataFrame(columns = range(0,self.number))
        temporary_df.loc[len(temporary_df)] = self.indexes
        temporary_df["Combined"] = temporary_df.values.tolist()
        self.history = pd.concat([self.history,temporary_df])
        
        dataframe = self.dataframe.iloc[:,self.indexes].dropna()
        dataframe["Combined"] = dataframe.values.tolist()
        self.smaller_dataframe = dataframe

        
    

        
        
    def task(self,debug = False):

        for i in range(0,self.R):
            self.initialize_indexes()
            self.calculate_parameters()
            self.calculate_hk()
            self.calculate_hs()
            self.hs_ref = self.hs
            self.hk_ref = self.hk
            self.indexes_ref = self.indexes
            self.history = pd.DataFrame(columns = range(0,self.number))
            self.history.loc[len(self.history)] = self.indexes
            self.history["Combined"] = self.history.values.tolist()
            rejects = 0
            runs = 0
            while rejects < self.rejects_limit and runs < self.T:
                self.swap()
                self.calculate_parameters()
                self.calculate_hk()
                self. calculate_hs()
                if (self.hk_ref < self.hk):
                    print("I changed hk=",self.hk_ref," hk* =",self.hk)
                    self.hk_ref = self.hk
                    self.hs_ref = self.hs
                    self.indexes_ref = self.indexes
                    
                else:
                    rejects += 1
                    print("I maintained, rejects =",rejects,"hk=",self.hk_ref," hk* =",self.hk)
                runs += 1
                
                
            self.output.loc[len(self.output)] = [self.indexes,self.hk,self.hs,self.number]
            self.history = pd.DataFrame(columns=self.history.columns)
            
    def run(self):
        self.output = pd.DataFrame(columns = ["Indexes","Hk","Hs","n"])
        self.n = [2,5,10,20,30]
        for self.number in self.n:
            print("Exploring n=",self.number)
            if self.number < self.dataframe.shape[1]:
                self.task()
        print("saving the results...")    
        self.output.to_csv(output_csv,index = None)

CVS = Cvs("CVS",input, R, T)
CVS.run()




