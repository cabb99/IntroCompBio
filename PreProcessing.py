# # Data preprocessing
# It allows to create new variables that are more related to the dependent variables. 
# It is also possible to modify the independent variables to diminish the effect of outliers.
# 
# The tables containing the new data are saved in Qtraining.csv and Qscoring.csv
# 
# ## Making all the variables quantitative
# To make all the variables quantitative we take all the values that are yes/no, pos/neg F/M and convert them to 1 and 0. 
# In the case that there is no data or not conclusive data we use NaN (from numpy). 
# 
# In the special case of cyto.cat that has multiple values we create a column for each unique value. 
# Each column contains 1 or 0 depending depending if the data corresponds to that unique value or not.
# 
# ## Preprocessing data for proteins
# Sometimes a low concentration, as well as a high concentration of the protein can cause illness. 
# If both a high value and a low value can have the same effect on the dependent variables, then there will not be a good correlation. 
# So we are comparing the protein against the mean and squaring the result. (What i mistakenly called normalized)
# Since the mean is expected to be 0 because the data of all the proteins are mean centered, squaring the result should be similar to the "Normalization", what I am calling now "Squared"
# 
# ## Cutoff of Dependent data.
# If a patient lives longer than 2 years or has a remission longer than 2 years is considered as a category.
# Some patients have much longer periods, that may be independent of any variable, because the patient is still alive.
# We don't want this data interfering with the prediction, so we may suppose that all the patients that live longer than 2 years and a half (130 weeks) live just 130 weeks.
# We can apply the same logic for Remission.
# 
# ## To Do
# It is possible to do also a pca analysis (maybe on protenis) and use the pca columns.
# It is also possible to use threeshold values (for example if x=1 if x>0.4 else 0). 
# It is also possible to square the protein data, changing the mean, from 0 to another value.
# This new value may be chossen to maximize the correlation between the dependent and independent variables.


import numpy,pandas
#For the preprocessing we need the data from Dream9.xlsx
Dream9_training=pandas.read_excel('Dream9.xlsx',"trainingData") 
Dream9_scoring=pandas.read_excel('Dream9.xlsx',"scoringData")
Dream9=pandas.concat([Dream9_training,Dream9_scoring])

#Division of types of Variables
All=list(Dream9_training.keys())
Sc=list(Dream9_scoring.keys())

#Dependent variables are present in the training set but not in the scoring set
Dependent=[]
for v in All:
    if v not in Sc:
        Dependent+=[v]

#Categorical variables have discrete values and can't be measured by euclidean distances
Categorical=['SEX', 'PRIOR.MAL', 'PRIOR.CHEMO', 'PRIOR.XRT', 'Infection', 'cyto.cat', 
             'ITD', 'D835', 'Ras.Stat', 'resp.simple', 'Relapse', 'vital.status']

#The last 231 variables are proteins
Protein=All[-231:]



###############
#  Functions  #
###############

def alias(Series,aliases):
    #Changes the values on a Series with aliases as a dict that transform the old values in the new values
    new_Series=(Series=='This creates and new series')
    for key,data in zip(Series.keys(),Series):
        new_Series[key]=False
        for val in aliases:
            try:
                if numpy.isnan(data):
                    new_Series[key]=numpy.nan
            except:
                pass             
            if data==val:
                new_Series[key]=aliases[val]
                break
    return new_Series

def split(Series,All_Data):
    #For Series with multiple values, creates a table with a column for each unique value
    #The value is True for the correct column and False for all the other columns
    D=[]
    for value in All_Data[Series.name].unique():
        q=(Series==value)
        q.name='%s=%s'%(q.name,value)
        D+=[q]
    return pandas.concat(D,axis=1)


def difference_from_mean(Table,All_Data):
    #This function creates a new Table with the values equal to (value-mean)/std
    #Since most of the values are already centered around 0 it woul be better to just take the square?
    D=[]
    for i,var in enumerate(Table.keys()):
        m=All_Data[var].mean()
        std=All_Data[var].std()
        D+=[(Table[var]-m)**2/std]
        D[i].name='%s_Normalized'%var
    return pandas.concat(D,axis=1)

def squared(Table):
    #This function squares all the values on a table
    D=[]
    for i,var in enumerate(Table.keys()):
        D+=[Table[var]**2]
        D[i].name='%s_Squared'%var
    return pandas.concat(D,axis=1)

def cutoff(Series,cutoff):
    #This function makes values above a threeshold equal to the threeshold
    new_Series=Series.copy()
    for key,data in zip(Series.keys(),Series):
        if data>cutoff:
            new_Series[key]=cutoff
        else:
            new_Series[key]=data
    new_Series.name='%s_cut'%Series.name
    return new_Series
    
####################
#  Pre-processing  #
####################    
    
def PreProcess(table,Dream9):
    #Select all variables that are not caregorical
    Tables=[table[[v for v in table.keys() if v not in Categorical]]]
    
    #Convert yes/no to 1/0
    Alias_Dict={'SEX':{'F':1},'PRIOR.MAL':{'YES':1},'PRIOR.CHEMO':{'YES':1},'PRIOR.XRT':{'YES':1},
                'Infection':{'Yes':1},'ITD':{'POS':1,'ND':numpy.nan},'D835':{'POS':1,'ND':numpy.nan},
                'Ras.Stat':{'POS':1,'NotDone':numpy.nan},'resp.simple':{'CR':1},'Relapse':{'Yes':1},
                'vital.status':{'A':1}}
    Aliased=[]
    for key in Alias_Dict:
        if key in table.keys():
            Aliased+=[alias(table[key],Alias_Dict[key])]
    Tables+=[pandas.concat(Aliased,axis=1)]
    
    #Split data that has multiple values
    Tables+=[split(table['cyto.cat'],Dream9)]
    
    #Create new data for protein
    #Tables+=[difference_from_mean(table[Protein],Dream9)]
    
    #Create new data for protein
    Tables+=[squared(table[Protein])]
    
    #Join everything
    
    Cut=[]
    for key in ['Overall_Survival','Remission_Duration']:
        if key in table.keys():
            Cut+=[cutoff(table[key],130)]
    if len(Cut)>0:
        Tables+=[pandas.concat(Cut,axis=1)]

    return pandas.concat(Tables,axis=1)

if __name__=='__main__':
   
    #This part is to test the function

    #Open the data and read in pandas
    Dream9_training=pandas.read_excel('Dream9.xlsx',"trainingData")
    Dream9_scoring=pandas.read_excel('Dream9.xlsx',"scoringData")
    Dream9=pandas.concat([Dream9_training,Dream9_scoring])

    #Create the new tables
    Q_training=PreProcess(Dream9_training,Dream9)
    Q_scoring=PreProcess(Dream9_scoring,Dream9)

    #Save the tables as csv
    Q_training.to_csv('Qtraining.csv')
    Q_scoring.to_csv('Qscoring.csv')

    #Number of columns and rows of new Table
    print Q_training.shape
    print Q_scoring.shape
