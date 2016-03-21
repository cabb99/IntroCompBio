
# coding: utf-8

# In[1]:

import PreProcessing
from Prediction import AIPrediction
import pandas,json


# In[2]:

#Define Prediction name
Pname='LinearRegression'

#Select starting variables. If you want to start again, make sure to erase 'Good_variables_%s.json'%Pname 
Good_variables={'pCR': [u'KDR_Squared', u'ATF3', u'RPS6_Squared', u'cyto.cat=Misc', u'GATA3', u'CDKN2A_Squared', u'NF2.pS518', u'CASP9.cl315_Squared', u'IGFBP2', u'SMAD3_Squared', u'PRKAA1_2.pT172_Squared', u'HDAC3_Squared', u'CLPP', u'PRIOR.MAL', u'ATG7', u'cyto.cat=diploid', u'DLX1_Squared', u'MSI2', u'CCNE2', u'NPM1.3542', u'ARC', u'cyto.cat=21', u'ITGAL', u'SMAD2_Squared', u'RPS6.pS240_244', u'MYC', u'LCK_Squared', u'ITGA2', u'GAPDH', u'CCNE1', u'PA2G4.pT70_Squared', u'cyto.cat=-7', u'MTOR.pS2448_Squared', u'CD44', u'PRKCB.II_Squared', u'MAP2K1_2.pS217_221_Squared', u'BAD.pS136_Squared', u'CASP9.cl330', u'GSKA_B.pS21_9', u'CTSG', u'FOXO3_Squared', u'TGM2', u'STAT3.pS727', u'CASP8_Squared', u'PIK3CA', u'RPS6', u'SFN', u'PTK2_Squared', u'ZNF296_Squared', u'PRKCD.pT507', u'Age.at.Dx', u'STMN1_Squared', u'YWHAZ_Squared', u'HSPB1', u'STMN1', u'PDK1.pS241_Squared', u'CDK1', u'MAPK9'],
                'pRelapse': [u'cyto.cat=t9;22', 'IGFBP2_Squared', u'CCND3', u'KIT_Squared', u'PTEN.pS380T382T383', u'BCL2_Squared', u'BAK1_Squared', u'SMAD5.pS463_Squared', 'MDM2', 'ARC', u'PTPN11_Squared', u'H3histon_Squared', u'PA2G4.pS65_Squared', 'HDAC1_Squared', u'EIF2S1.pS51._Squared'], 
                'OS': [u'PRIOR.MAL', u'ARC', u'cyto.cat=diploid', u'H3histon', u'Age.at.Dx', u'PTGS2_Squared', u'SMAD4', u'PA2G4.pS65', u'STMN1', u'EIF2AK2', u'H3K27Me3', u'HSP90AA1_B1'], 
                'Remission': [u'CASP9.cl330', u'ERG', u'ALBUMIN', u'CASP3.cl175', u'TP53', u'RPS6KB1.pT389', u'PLAC1', u'JMJD6', u'SMAD3_Squared', u'ERG_Squared', u'TRIM24', u'Age.at.Dx', u'HSPA1A_L', u'ATG7_Squared', u'ARC_Squared', u'STAT3.pS727', u'CBL_Squared', u'BIRC5_Squared', u'ARC', u'YWHAE', u'SMAD5.pS463', u'BRAF_Squared', u'MTOR.pS2448_Squared']}

#Select Dependent Variables
Variables_for_prediction={'pCR':'resp.simple',
                          'pRelapse':'Relapse',
                          'OS':'Overall_Survival_cut',
                          'Remission':'Remission_Duration_cut'}

#Select scoring measure to use for minimization
Fast_measure='BAC' #Select between BAC PCC or Scr. You can select two Measures, 
                       #the second measure will overwrite the first one. ('BAC','Auroc')
Slow_measure={'pCR':'BAC',
              'pRelapse':'BAC',
              'OS':'BAC',
              'Remission':'BAC'}

#Select max time in hours to run the prediction
max_time_h= 5/60.0#Max number of hours that you can run the simulation


# In[3]:

#Read prediction bias
Var_value=pandas.DataFrame.from_csv('Variable_value.csv')


# In[4]:

if True: #Use Previous Good Variables
    #Open Good variables if it exists
    try:
        with open('Good_variables_%s.json'%Pname) as handle:
            Good_variables=json.loads(handle.read())
    except IOError:
        print 'Did not found previous good variables'
        pass
    
if True:#Use previous prefiction bias
    #Open Variable values if it exists
    try:
        pandas.DataFrame.from_csv('Var_value_%s.csv'%Pname)
    except IOError:
        print 'Did not found previous variable values'
        pass


# In[5]:

#Define prediction class
from sklearn import linear_model
Pred=AIPrediction(pivot=Good_variables,method=linear_model.LinearRegression())
print Pred.result(rep=10)


# In[6]:

#Define Best Prediction
Pred.pivot=Good_variables
groups=Pred.create_groups(50)
Best_scores=Pred.score(groups=groups,Measure=tuple(set(Slow_measure.values())))
Best_result=Pred.accuracy(groups=groups,Measure=Fast_measure)


# In[7]:

#Define Biased Sample
import random,numpy
def biased_sample(variables,k,bias):
    var=list(variables)
    Bias=bias[variables].copy()
    n = len(variables)
    bias=numpy.cumsum(Bias[variables])
    if len(bias)==0:
        bias=[1 for i in range(n)]
    if not 0 <= k <= n:
        #print population
        k=n
        #raise ValueError("sample larger than population")
    r = random.random()
    result = [None] * k # Initialize result
    selected = set()
    max_freq=bias[-1]
    for i in xrange(k):
        #print i
        r = random.random() * max_freq
        #r=max_freq
        for j,b in enumerate(bias):
            #print r,j,b
            if r<b:
                break
        selected.add(var[j])
        result[i] = var[j]
        #print r,b,j,var[j],'\n',bias
        var.pop(j)
        bias=numpy.cumsum(Bias[var])
        if len(bias)>0:
            max_freq=bias[-1]
        else:
            break
    return result


# In[8]:

#Evolution algorithm
import time,json

Stop=time.time()+60*60*max_time_h #Stop time in seconds
i=0 #Iteration number counter
while time.time()<Stop:  
    Update_Best=False #Only update best if conditions apply
    i=i+1
    #print i
    #Select Test variables
    Test_variables={}
    change_variables={}
    change={}
    for key in Good_variables.keys():
        Rnd=random.random()
        if Rnd>0.65: #Sometimes drop a variable
            l=len(Good_variables[key])
            n=min(l-random.randint(1,5),1)
            change.update({key:'remove'})
            new_variables=biased_sample(Good_variables[key],n,1/(Var_value[key]+1))
            change_variables.update({key:[v for v in Good_variables[key] if v not in new_variables]})
            Test_variables.update({key:new_variables})
        elif Rnd>0.1: #Sometimes add a variable
            varis=[v for v in Var_value[key].index if v not in Good_variables[key]]
            l=len(varis)
            n=min(random.randint(1,5),l)
            change.update({key:'add'})
            change_variables.update({key:biased_sample(varis,n,Var_value[key])})
            new_variables=Good_variables[key]+change_variables[key]
            Test_variables.update({key:new_variables})
        else: #Sometimes choose something completely random
            l=min(len(Var_value[key].index),len(Good_variables[key])*2)
            change.update({key:'try'})
            new_variables=biased_sample(varis,random.randint(1,l),Var_value[key])
            change_variables.update({key:new_variables})
            Test_variables.update({key:new_variables})
        #print key,change, change_variables[key]

    #Update Test_variables
    Pred.pivot=Test_variables
    try:    
        Test_result=Pred.accuracy(groups=groups[:10],Measure=Fast_measure) #Fast measure
    except KeyError:
        groups=Pred.create_groups(50)        
        continue
    
    #Test against best prediction
    Next=True
    rew={}
    for key in Good_variables.keys():
        rew.update({key:1.0})
        if Test_result[key]>Best_result[key]:
            Next=False
            if change[key]=='add':
                rew.update({key:1.05})
            elif change[key]=='remove':
                rew.update({key:0.95})
            elif change[key]=='try':
                rew.update({key:1.02})
            else:
                raise 'Undefined change'       
        else:
            if change[key]=='add':
                rew.update({key:0.9})
            elif change[key]=='remove':
                rew.update({key:1.1})
            elif change[key]=='try':
                rew.update({key:0.98})
            else:
                raise 'Undefined change'
    
    #Update bias weights
    for key in Good_variables.keys():
        for var in change_variables[key]:
            if Var_value[key][var]>0.02 and Var_value[key][var]<2.0:
                Var_value[key][var]*=rew[key]
            if Var_value[key][var]<0.02 and Var_value[key][var]>0:
                Var_value[key][var]=0.02
            if Var_value[key][var]>2.0:
                Var_value[key][var]=2.0
   
    #Add another bias to allow unused variables to be used again
    #To Be Done

    #Write to csv
    Var_value.to_csv('Var_value_%s.csv'%Pname)
    
    #Continue with next if it did not improve
    if Next:
        continue
    
    #Test how many times it has a better result
    Test_scores=Pred.score(groups=groups,Measure=tuple(set(Slow_measure.values())))
    for key in Test_scores.keys():
        conf=0.0
        new_BACs=Test_scores[key][Slow_measure[key]]
        old_BACs=Best_scores[key][Slow_measure[key]]
        new_BACs.sort()
        old_BACs.sort()
        pos=0.0
        tot=float(len(new_BACs)**2)
        for new in new_BACs:
            for k,old in enumerate(old_BACs):
                if old>new:
                    pos+=k
                    break
        conf=pos/tot
        #If more than 60% are better, update
        if conf>0.6:
            with open('Log_%s.txt'%Pname,'a+') as handle:
                handle.write('Updated: '+key+'\n')
            Good_variables.update({key:Test_variables[key]})
            Update_Best=True
    
    #Update variables
    if Update_Best:
        Pred.pivot=Good_variables
        groups=Pred.create_groups(50)
        Best_result=Pred.accuracy(groups=groups)
        Best_scores=Pred.score(groups=groups,Measure=tuple(set(Slow_measure.values())))
        with open('Log_%s.txt'%Pname,'a+') as handle:
            handle.write(str(Good_variables)+'\n')
            handle.write(str(Best_result)+'\n')
        with open('Good_variables_%s.json'%Pname,'w+') as handle:
            handle.write(json.dumps(Good_variables))


# In[9]:

print 'iter',i
Pred.pivot=Good_variables
print Pred.result(rep=10)


# In[10]:

Pred.predict()

