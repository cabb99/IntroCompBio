
# coding: utf-8

# In[1]:

import numpy,pandas
from sklearn.decomposition import PCA
import numpy as np

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

def make_pca(Table,All_Data,n,name='PCA'):
    pca = PCA(n_components=n)
    pca.fit(All_Data[Table.keys()])
    trans_PCA=pca.transform(Table)
    return pandas.DataFrame(trans_PCA,columns=['%s_%i'%(name,i+1) for i in range(n)],index=Table.index)

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

def binned(Series,bins=[0,52,104]):
    #This function will bin the results from Remission and Overall Survival as expected    
    bins = numpy.array(bins)
    digitized = list(numpy.digitize(Series, bins))
    for i,v in enumerate(Series):
        if numpy.isnan(v):
            digitized[i]=numpy.nan
    return pandas.Series(digitized,index=Series.index,name='%s_binned'%Series.name)*52-26
####################
#  Pre-processing  #
####################    
    
def PreProcess(table,Dream9):
    #Select all variables that are not Categorical
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
    Tables+=[make_pca(table[Protein],Dream9,200)]
    Tables+=[make_pca(squared(table[Protein]),squared(Dream9[Protein]),200,name='PCA_Sq')]
    
    Cut=[]
    for key in ['Overall_Survival','Remission_Duration']:
        if key in table.keys():
            Cut+=[cutoff(table[key],130)]
    if len(Cut)>0:
        Tables+=[pandas.concat(Cut,axis=1)]
        
    Bin=[]
    for key in ['Overall_Survival','Remission_Duration']:
        if key in table.keys():
            Bin+=[binned(table[key])]
    if len(Bin)>0:
        Tables+=[pandas.concat(Bin,axis=1)]
    
    
    #Join everything
    return pandas.concat(Tables,axis=1)

Q_Dependent=Dependent+['Overall_Survival_cut','Remission_Duration_cut','Overall_Survival_binned','Remission_Duration_binned']
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
    #A=binned(Dream9_training['Remission_Duration'])


# In[2]:

Q_training=PreProcess(Dream9_training,Dream9)
Q_scoring=PreProcess(Dream9_scoring,Dream9)
Q_Dependent=[v for v in Q_training.keys() if v not in Q_scoring.keys()]


# In[3]:

if __name__=='__main__':
    #Correlation for Continuous variables
    Corr=pandas.DataFrame()
    for Variable in Q_Dependent:
        C=Q_training[[t for t in Q_training.keys() if (t not in Q_Dependent)]+[Variable]].corr()[Variable][:-1]
        Corr=Corr.append(C)
    #Write correlation as csv
    Corr.T.to_csv('Correlations.csv')


# In[4]:

if __name__=='__main__':
    #Most important Variables in Correlation
    for Variable in Q_Dependent:
        A=Corr.T[Variable]**2
        A.sort(ascending=False)
        print Corr[A.head().index].T[Variable]


# In[5]:

if __name__=='__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib
    get_ipython().magic(u'matplotlib inline')

    x = -Q_training['Age.at.Dx']
    y = Q_training['HGB']
    z = -Q_training['PCA_Sq_151']
    area = Q_training['cyto.cat=diploid']*25+75 # 0 to 15 point radiuses
    colors = Q_training['resp.simple']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(x.name)
    ax.set_ylabel(y.name)
    ax.set_zlabel(z.name)
    ax.scatter(x, y, z,s=area, c=colors,cmap=matplotlib.cm.coolwarm_r)
    #Blue is 1 and red is 0
    plt.show()


# In[6]:

#Calculate how much information in gained on each column
#Calculate the entropy of the subset
def information_gain(Table,Dependent,Independent):
    Table=Table[Table[Dependent].notnull()]
    freq=[]
    for dval in Table[Dependent].unique():
        freq+=[sum(Table[Dependent]==dval)]
    Freq=[float(f)/sum(freq) for f in freq]
    E=0
    for f in Freq:
        E+=-f*np.log(f)/np.log(2)
    #print 'Subset Entropy:', E
    Vars=[]

    #Calculate the entropy of each variable
    for ind in Independent:     
        if ind in Categorical:
            IG=E
            for ival in Table[ind].unique():
                if np.isnan(ival):
                    continue
                SubTable=Table[Table[ind]==ival]
                #print SubTable
                freq=[]
                for dval in Table[Dependent].unique():
                    freq+=[sum(SubTable[Dependent]==dval)]
                Freq=[float(f)/sum(freq) for f in freq]
                #print Freq
                ES=0
                for f in Freq:
                    ES+=-f*np.log(f)/np.log(2) if f<>0 else 0
                #print ES
                IG-=float(len(SubTable))/len(Table)*ES
            #print 'Information gain from %s: %f'%(ind,IG)
            Vars+=[(IG,ind)]
        else:
            Threeshold=[]
            prev_SubTableA_len=0
            for ival in np.arange(min(Table[ind]),max(Table[ind]),(max(Table[ind])-min(Table[ind]))/500.0):
                IG=E
                SubTableA=Table[Table[ind]<ival]
                SubTableB=Table[Table[ind]>=ival]
                if len(SubTableA)<1 or len(SubTableB)<1:
                    continue
                if len(SubTableA)==prev_SubTableA_len:
                    continue
                else:
                    prev_SubTableA_len=len(SubTableA)
                freq=[]
                for dval in Table[Dependent].unique():
                    freq+=[sum(SubTableA[Dependent]==dval)]
                Freq=[float(f)/sum(freq) for f in freq]
                #print Freq
                ES=0
                for f in Freq:
                    ES+=-f*np.log(f)/np.log(2) if f<>0 else 0
                #print ES
                IG-=float(len(SubTableA))/len(Table)*ES
                #print SubTable
                freq=[]
                for dval in Table[Dependent].unique():
                    freq+=[sum(SubTableB[Dependent]==dval)]
                Freq=[float(f)/sum(freq) for f in freq]
                #print Freq
                ES=0
                for f in Freq:
                    ES+=-f*np.log(f)/np.log(2) if f<>0 else 0
                #print ES
                IG-=float(len(SubTableB))/len(Table)*ES
                Threeshold+=[(IG,ival)]
            Threeshold.sort(reverse=True)
            #print Threeshold
            #break
            #print 'Information gain from %s: %f at theeshold:%f'%(ind,Threeshold[0][0],Threeshold[0][1])
            if len(Threeshold)>0:
                Vars+=[(Threeshold[0][0],ind,Threeshold[0][1])]
            else:
                Vars+=[(0,ind)]
    Information_gain=pandas.Series([v[0] for v in Vars],index=[v[1] for v in Vars],name='Information Gain')
    Threesholds=pandas.Series([v[2] for v in Vars if len(v)>2],index=[v[1] for v in Vars if len(v)>2],name='Threeshold')
    return pandas.concat([Information_gain,Threesholds],axis=1)

if __name__=='__main__':
    Q_Cat=['resp.simple','Relapse','vital.status','Overall_Survival_binned','Remission_Duration_binned']
    Ts=[]
    for Variable in Q_Cat:
        print Variable
        Independent=[v for v in Q_training.keys() if v in Q_scoring.keys()]
        Ts+=[information_gain(Q_training,Variable,Independent)]
    Information_Gain=pandas.concat(Ts,keys=Q_Cat,axis=1)
    Information_Gain.to_csv('InformationGain.csv')


# In[7]:

#Most important variables in Information Gain
if __name__=='__main__':
    for Variable in Q_Cat:
        A=Information_Gain[Variable]
        A=A.sort('Information Gain',ascending=False)
        A.name=Variable
        print Variable
        print A.head()


# In[8]:

if __name__=='__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib
    get_ipython().magic(u'matplotlib inline')

    x = -Q_training['Ras.Stat']
    y = Q_training['NPM1.3542']
    z = -Q_training['FIBRINOGEN']
    area = Q_training['PCA_Sq_38']*25+75 # 0 to 15 point radiuses
    colors = Q_training['resp.simple']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(x.name)
    ax.set_ylabel(y.name)
    ax.set_zlabel(z.name)
    ax.scatter(x, y, z,s=area, c=colors,cmap=matplotlib.cm.coolwarm_r)
    #Blue is 1 and red is 0
    plt.show()

