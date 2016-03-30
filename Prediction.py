
# coding: utf-8

# In[32]:

import numpy
import random
import pandas
from scipy.stats import norm
import PreProcessing

#Create the new tables
Training_set=PreProcessing.Q_training
Scoring_set=PreProcessing.Q_scoring

Variables_for_prediction={'pCR':'resp.simple','pRelapse':'Relapse','OS':'Overall_Survival_binned','Remission':'Remission_Duration_binned'}

class Prediction:
    ########################
    ##   Initialization   ##
    ########################
    def __init__(self,training=Training_set,PredVar=Variables_for_prediction):
        self.training=training
        self.PredVar=PredVar
    
    ######################
    ##  Group Creation  ##
    ######################
    
    def create_groups(self,rep,n=5):
        #Creates rep*n training and testing groups
        Groups=[]
        for i in range(rep):
            Groups+=self.split(n)
        return Groups
    
    def split(self,n=5):
        #Divides the training groups n times, to have training and testing groups.
        #The testing group will contain 1/n of the data
        #It returns a list of dictionaries.
        #To call a group use for example groups[0]['test']
        #That will return the first testing group
        training_keys=self.training.T.keys().get_values().copy()
        random.shuffle(training_keys)
        sublist=numpy.array_split(training_keys,n)
        groups=[]
        for i in range(n):
            train=[]
            test=[]
            for j in range(n):
                sl=list(sublist[j])
                if j<>i:
                    train+=sl
                else:
                    test+=sl
            train.sort()
            test.sort()
            groups+=[{'train':self.training.loc[train],'test':self.training.loc[test]}]
        return groups
    
    
    ####################
    ##   Prediction   ##
    ####################
    
    def predict(self,scoring=Scoring_set,out='prediction.csv'):
        #Main prediction function.
        #Returns the result of the prediction for the Scoring set as a csv
        #Also prints the results
        with open(out,'w+') as handle:
            handle.write('#Patient_id, pCR, pRelapse, Remission, OS\n')
            print '#Patient_id, pCR, pRelapse, Remission, OS'
            Data=zip(scoring.T.keys(),
                     self.pCR(self.training,scoring),
                     self.pRelapse(self.training,scoring),
                     self._bin(self.Remission(self.training,scoring)),
                     self._bin(self.OS(self.training,scoring)))
            for d in Data:
                handle.write(','.join(str(k) for k in d)+'\n')
                print ','.join(str(k) for k in d)
    
    
    def pCR(self,training,testing): 
        #Returns the probability of Complete Remission ('resp.simple')
        return self.pPred(training,testing,'pCR')
        
    def pRelapse(self,training,testing): 
        #Returns the probability of Relapse ('Relapse')
        return self.pPred(training,testing,'pRelapse')
        
    def Remission(self,training,testing): 
        #Returns the estimated Remission Duration ('Remission_Duration') in weeks (not binned)
        return self.qPred(training,testing,'Remission')
        
    def OS(self,training,testing): 
        #Returns the estimated Overall Survival ('Overall_Survival') in weeks (not binned)
        return self.qPred(training,testing,'OS')
    
    def pPred(self,training,testing,var):
        #Probability prediction, should return a value between 0 and 1.
        #Used for pCR and pRelapse
        avg_p=sum(training[self.PredVar[var]] == True)/float(len(training))
        return pandas.Series([avg_p+random.random()*1E-6 for i in range(len(testing))],index=testing.index)
        
    def qPred(self,training,testing,var):
        #Quantitative prediction, returns the number of weeks.
        #Used for Remission and OS
        #It will be binned after this, so if this is a binned prediction, multiply the result by 52
        count=numpy.bincount(self._bin(training[self.PredVar[var]]))
        count[0]=0 #bin 0 is nan, do not count nan
        val=range(len(count))
        val.sort(key=lambda i: count[i],reverse=True)
        mode=(val[0]-1)*52.0+26
        return pandas.Series([mode+random.random()*1E-2 for i in range(len(testing))],index=testing.index)
    
    def _bin(self,data,bins=[0,52,104]):
        #This function will bin the results from Remission and Overall Survival as expected    
        bins = numpy.array(bins)
        digitized = numpy.digitize(data.convert_objects(convert_numeric=True), bins)
        for i,v in enumerate(data):
            if numpy.isnan(v):
                digitized[i]=0
        return pandas.Series(digitized,index=data.index)

    
    ######################
    ## Scoring methods  ##
    ######################
    
    def accuracy(self,rep=1,groups=[],Measure='BAC'):
        #Returns the Balanced accuracy prediction for all the functions
        groups=self.create_groups(rep) if groups==[] else groups
        S=self.score(Measure=Measure,groups=groups)
        if type(Measure)<>tuple:
            Measure=[Measure]
        Acc={s:numpy.nan for s in S}
        for M in Measure:
            Acc.update({s:numpy.mean(S[s][M]) for s in S if M in S[s].keys()})
        return Acc
    
    def result(self,rep=1,groups=[],alpha=0.5,Measure='All'):    
        #Returns the scoring calculations for the prediction
        #You can select wich scoring to return with Measure
        #alpha allows to select the scoring that is lower than alpha*100% of the scores, supossing a normal distribution.
        #If alpha is 0.5 it will return the mean
        groups=self.create_groups(rep) if groups==[] else groups
        S=self.score(groups,Measure)
        T=pandas.concat([pandas.Series({key:numpy.mean(S[v][key])-norm.ppf(alpha)*numpy.std(S[v][key]) for key in S[v]},name=v) for v in S],axis=1)
        return T.T
        #return  {'pCR':pCR_Error,'pRelapse':pRelapse_Error,'Remission':Remission_Error,'OS':OS_Error}
        #pRelapse_Error
    
    def score(self,groups=[],Measure='All'):
        #Main scoring function
        #Returns a dictionary for all the scores
        #You can select a Measure
        #Will create groups if no group given.

        #Define groups if not defined
        groups=self.split() if groups==[] else groups
        
        #Initialize score dictionary
        lg=range(len(groups))
        S={'pCR':{},'pRelapse':{},'Remission':{},'OS':{}}
        [S[v].update({'BAC':[0 for r in lg],'PCC':[0 for r in lg],'Scr':[0 for r in lg]})for v in S]
        [S[v].update({'AUROC':[0 for r in lg]})for v in ['pCR','pRelapse']]
        [S[v].update({'CI':[0 for r in lg]})for v in ['OS','Remission']]
        #Independent=[v for v in self.training.keys() if v not in Dependent]
        
        if type(Measure)==tuple:
            Measures=Measure
        else:
            Measures=[Measure]
            
        for Measure in Measures:    
            for i,g in enumerate(groups):
                train=g['train']
                test=g['test']#[Independent]

                #Train and evaluate
                pCR_values=self.pCR(train,test)
                pRelapse_values=self.pRelapse(train,test)
                Remission_values=self.Remission(train,test)
                OS_values=self.OS(train,test)

                #Compare the expeced values to the obtained values
                #test=g['test']

                #Calculate the BAC scores
                if Measure in ['All','BAC']:
                    S['pCR']['BAC'][i]=self._pBAC(pCR_values,test['resp.simple'])
                    S['pRelapse']['BAC'][i]=self._pBAC(pRelapse_values,test['Relapse'])
                    S['Remission']['BAC'][i]=self._cBAC(self._bin(Remission_values),self._bin(test['Remission_Duration']))
                    S['OS']['BAC'][i]=self._cBAC(self._bin(OS_values),self._bin(test['Overall_Survival']))
                    if Measure=='BAC':
                        continue

                #Calculate the PCC scores
                if Measure in ['All','PCC']:
                    S['pCR']['PCC'][i]=self._pPCC(pCR_values,test['resp.simple'])
                    S['pRelapse']['PCC'][i]=self._pPCC(pRelapse_values,test['Relapse'])
                    S['Remission']['PCC'][i]=self._pPCC(Remission_values,test['Remission_Duration'])
                    S['OS']['PCC'][i]=self._pPCC(OS_values,test['Overall_Survival'])
                    if Measure=='PCC':
                        continue

                #Calculate the AUROC scores
                if Measure in ['All','AUROC']:
                    S['pCR']['AUROC'][i]=self._AUROC(pCR_values,test['resp.simple'])
                    S['pRelapse']['AUROC'][i]=self._AUROC(pRelapse_values,test['Relapse'])
                    if Measure=='AUROC':
                        continue

                #Calculate the CI scores
                if Measure in ['All','CI']:
                    S['Remission']['CI'][i]=self._RCI(Remission_values,test)
                    S['OS']['CI'][i]=self._OSCI(OS_values,test)
                    if Measure=='CI':
                        continue

                #Expected Score
                if Measure in ['All','Scr']:
                    S['pCR']['Scr'][i]=-self._expectedScore(pCR_values,test['resp.simple'])
                    S['pRelapse']['Scr'][i]=-self._expectedScore(pRelapse_values,test['Relapse'])
                    S['Remission']['Scr'][i]=-self._expectedScore(self._bin(Remission_values),self._bin(test['Remission_Duration']),norm=2)
                    S['OS']['Scr'][i]=-self._expectedScore(self._bin(OS_values),self._bin(test['Overall_Survival']),norm=2)
                    if Measure=='Scr':
                        continue
        
        return S
    
    def _pBAC(self,predicted,expected):
        TP=float(((predicted>=0.5) & (expected>=0.5)).dropna().sum())
        TN=float(((predicted<0.5) & (expected<0.5)).dropna().sum())
        P=max(float((expected>=0.5).dropna().sum()),1E-64)
        N=max(float((expected<0.5).dropna().sum()),1E-64)
        return (TP/P+TN/N)/2

    
    def _cBAC(self,predicted,expected):
        TV=[]
        a=set(expected)
        for val in set(a):
            if val>0:
                TP=float(((predicted==val) & (expected==val)).sum())
                P=max(float((expected==val).sum()),1E-64)
                TV+=[TP/P]
        return numpy.mean(TV)
        
    def _AUROC(self,predicted,expected):
        AUC=0
        TPRl=0
        FPRl=0
        K=list(set(numpy.concatenate((predicted,[1.0,0.0]))))
        K.sort(reverse=True)
        for k in K:
            T=(predicted>=k)
            TP=((expected==1) & (T==1)).sum()
            FP=((expected==0) & (T==1)).sum()
            TPR=TP/max(float(expected.sum()),1E-64)#may be 0 sometimes
            FPR=FP/max(float((expected==0).sum()),1E-64)#may be 0 sometimes   
            AUC+=TPRl*(FPR-FPRl)
            FPRl=FPR
            TPRl=TPR
        return AUC
            
    
    def _OSCI(self,predicted,expected):        
        c=0.0
        H=0.0
        for i,ai,pi,Ai in zip(range(len(predicted)),predicted,expected['Overall_Survival'],expected['vital.status']):
            if numpy.isnan(ai) or numpy.isnan(pi):
                continue
            for j,aj,pj,Aj in zip(range(len(predicted)),predicted,expected['Overall_Survival'],expected['vital.status']):
                if i>=j:
                    continue
                if numpy.isnan(aj) or numpy.isnan(pj):
                    continue
                if ai<=aj and Ai: #The patient i has smaller Survival but is still alive
                    continue
                if aj==ai and Aj: #The patient j has smaller Survival but is still alive
                    continue
                
                if numpy.sign(round(ai-aj,5))==numpy.sign(round(pi-pj,5)):
                    H+=1
                    c+=1
                else:
                    c+=1
        return H/c
                    
    def _RCI(self,predicted,expected):
        c=0.0
        H=0.0
        for i,ai,pi,Ai in zip(range(len(predicted)),predicted,expected['Remission_Duration'],expected['Relapse']):
            if numpy.isnan(ai) or numpy.isnan(pi):
                continue
            for j,aj,pj,Aj in zip(range(len(predicted)),predicted,expected['Remission_Duration'],expected['Relapse']):
                if i>=j:
                    continue
                if numpy.isnan(aj) or numpy.isnan(pj): #no value on Relapse
                    continue
                if ai<=aj and not Ai: #The patient i has smaller Remission but has not Relapsed
                    continue
                if aj==ai and not Aj: #The patient j has smaller Remission but has not Relapsed
                    continue
                
                if numpy.sign(round(ai-aj,5))==numpy.sign(round(pi-pj,5)):
                    H+=1
                    c+=1
                else:
                    c+=1
        return H/c
    
    def _pPCC(self,predicted,expected): #Pearson Correlation Coeficient
        A=pandas.concat([predicted,expected],axis=1,keys=['pred','exp'])
        A=A.dropna() #drop data that has na as a result
        p=A['pred'].mean()
        a=A['exp'].mean()
        sp=max(((A['pred']-p)**2).sum()**0.5,1E-64) #Sometimes this is 0, then S=0 too
        sa=max(((A['exp']-a)**2).sum()**0.5,1E-64) #Sometimes this is 0, then S=0 too
        S=(A['pred']-p)*(A['exp']-a)
        return (S.sum()/sp/sa+1)/2
    
    def _expectedScore(self,predicted,expected,norm=1):
        A=pandas.concat([predicted,expected],axis=1,keys=['pred','exp'])
        len1=len(A)
        A=A.dropna() #drop data that has na as a result
        len2=len(A)
        missing_values=len1-len2
        p=A['pred']
        a=A['exp']
        #print len(p),len(predicted),len(Scoring_set),len(Scoring_set)/float(len(predicted))
        return (sum(((p-a)/norm+missing_values)**2)*len(Scoring_set)/float(len(predicted)))**0.5 #*len(p)/len(predicted)*len(Scoring_set)

if __name__=='__main__':
   
    #This part is to test the function
    Dummy=Prediction()
    print Dummy.result(alpha=0.5,Measure='All',rep=5) 
    print dir(Dummy)


# In[38]:

#Dummy.accuracy(Measure=('AUROC','AUROC'))


# In[39]:

#print Dummy.accuracy(Measure='PCC',rep=5) 


# In[68]:

from sklearn import datasets, linear_model
import numpy
import random
import pandas
from scipy.stats import norm
import PreProcessing 

#Create the new tables
Training_set=PreProcessing.Q_training
Scoring_set=PreProcessing.Q_scoring
Dependent=PreProcessing.Q_Dependent
Independent=[v for v in Training_set.keys() if v not in Dependent]

Variables_for_prediction={'pCR':'resp.simple','pRelapse':'Relapse','OS':'Overall_Survival_binned','Remission':'Remission_Duration_binned'}


#Best set of variables
Good_variables={'pCR': [u'KDR_Squared', u'ATF3', u'RPS6_Squared', u'cyto.cat=Misc', u'GATA3', u'CDKN2A_Squared', u'NF2.pS518', u'CASP9.cl315_Squared', u'IGFBP2', u'SMAD3_Squared', u'PRKAA1_2.pT172_Squared', u'HDAC3_Squared', u'CLPP', u'PRIOR.MAL', u'ATG7', u'cyto.cat=diploid', u'DLX1_Squared', u'MSI2', u'CCNE2', u'NPM1.3542', u'ARC', u'cyto.cat=21', u'ITGAL', u'SMAD2_Squared', u'RPS6.pS240_244', u'MYC', u'LCK_Squared', u'ITGA2', u'GAPDH', u'CCNE1', u'PA2G4.pT70_Squared', u'cyto.cat=-7', u'MTOR.pS2448_Squared', u'CD44', u'PRKCB.II_Squared', u'MAP2K1_2.pS217_221_Squared', u'BAD.pS136_Squared', u'CASP9.cl330', u'GSKA_B.pS21_9', u'CTSG', u'FOXO3_Squared', u'TGM2', u'STAT3.pS727', u'CASP8_Squared', u'PIK3CA', u'RPS6', u'SFN', u'PTK2_Squared', u'ZNF296_Squared', u'PRKCD.pT507', u'Age.at.Dx', u'STMN1_Squared', u'YWHAZ_Squared', u'HSPB1', u'STMN1', u'PDK1.pS241_Squared', u'CDK1', u'MAPK9'],
                'pRelapse': [u'cyto.cat=t9;22', 'IGFBP2_Squared', u'CCND3', u'KIT_Squared', u'PTEN.pS380T382T383', u'BCL2_Squared', u'BAK1_Squared', u'SMAD5.pS463_Squared', 'MDM2', 'ARC', u'PTPN11_Squared', u'H3histon_Squared', u'PA2G4.pS65_Squared', 'HDAC1_Squared', u'EIF2S1.pS51._Squared'], 
                'OS': [u'PRIOR.MAL', u'ARC', u'cyto.cat=diploid', u'H3histon', u'Age.at.Dx', u'PTGS2_Squared', u'SMAD4', u'PA2G4.pS65', u'STMN1', u'EIF2AK2', u'H3K27Me3', u'HSP90AA1_B1'], 
                'Remission': [u'CASP9.cl330', u'ERG', u'ALBUMIN', u'CASP3.cl175', u'TP53', u'RPS6KB1.pT389', u'PLAC1', u'JMJD6', u'SMAD3_Squared', u'ERG_Squared', u'TRIM24', u'Age.at.Dx', u'HSPA1A_L', u'ATG7_Squared', u'ARC_Squared', u'STAT3.pS727', u'CBL_Squared', u'BIRC5_Squared', u'ARC', u'YWHAE', u'SMAD5.pS463', u'BRAF_Squared', u'MTOR.pS2448_Squared']}


#Variables for prediction
Variables_for_prediction={'pCR':'resp.simple','pRelapse':'Relapse','OS':'Overall_Survival_cut','Remission':'Remission_Duration_cut'}

class AIPrediction(Prediction):
    def __init__(self,training=Training_set,pivot=Good_variables,
                 PredVar=Variables_for_prediction,method=linear_model.LinearRegression(),
                 binned=False):
        self.training=training
        self.pivot=pivot
        self.ols=method
        self.PredVar=PredVar
        self.binned=binned
        
    def create_groups(self,rep):
        Groups=[]
        for i in range(rep):
            Groups+=self.split()
        New_Groups=[]
        for g in Groups:
            training=g['train']
            testing=g['test']
            Accept=True
            for Dependent in self.PredVar.values():
                A=pandas.concat([training[Independent],pandas.DataFrame(training[Dependent])],axis=1,keys=['ind','dep'])
                A=A.dropna()
                if len(A['dep',])<=0:
                    Accept=False
            if Accept:
                New_Groups+=[g]
        return New_Groups
    
    def pPred(self,training,testing,dep):
        ind=self.pivot[dep] #Select independent variables
        #ind=pandas.concat([training[ind],testing[ind]]).dropna(axis=1).keys() #drop variables that have na
        A=pandas.concat([training[ind],pandas.DataFrame(training[self.PredVar[dep]])],axis=1,keys=['ind','dep'])
        A=A.dropna() #drop data that has na as a result
        global S0,S1
        S0=A['ind',]
        S1=A['dep',][self.PredVar[dep]]
        self.ols.fit(A['ind',],A['dep',][self.PredVar[dep]]) #train
        test=testing[ind].dropna()
        Results=self.ols.predict(test)#.T[0] #predict
        #if len(Results.shape)==2:
        #    Results=Results.T[0]
        for i,val in enumerate(Results):
            if val>1:
                Results[i]=1
            if val<0:
                Results[i]=0
        R=pandas.Series(Results,index=test.index)
        return pandas.Series(R,index=testing.index)
        
    def qPred(self,training,testing,dep):
        ind=self.pivot[dep] #Select independent variables
        #ind=pandas.concat([training[ind],testing[ind]]).dropna(axis=1).keys() #drop variables that have na
        A=pandas.concat([training[ind],pandas.DataFrame(training[self.PredVar[dep]])],axis=1,keys=['ind','dep'])
        A=A.dropna() #drop data that has na as a result
        global S0,S1
        S0=A['ind',]
        S1=A['dep',][self.PredVar[dep]]
        if self.binned:
            S1=self._bin(S1)
        self.ols.fit(S0,S1) #train
        test=testing[ind].dropna()
        Results=self.ols.predict(test)#.T[0] #predict
        if self.binned:
            Results=Results*52-26
        #if len(Results.shape)==2:
        #    Results=Results.T[0]
        for i,val in enumerate(Results):
            if val<0:
                Results[i]=0
        R=pandas.Series(Results,index=test.index)
        return pandas.Series(R,index=testing.index)
    
    
if __name__=='__main__':
    LR=AIPrediction(pivot=Good_variables)
    print LR.result(rep=10)


# In[38]:

if __name__=='__main__':
    LR=AIPrediction(pivot=Good_variables,method=linear_model.RidgeCV(alphas=[0.1, 1.0, 10.0]))
    print LR.result(rep=10)


# In[39]:

if __name__=='__main__':
    from sklearn.kernel_ridge import KernelRidge
    LR=AIPrediction(pivot=Good_variables,method=KernelRidge(alpha=1.0))
    print LR.result(rep=10)


# In[69]:

if __name__=='__main__':
    from sklearn.svm import SVC
    LR=AIPrediction(pivot=Good_variables,method=SVC(),binned=True)
    print LR.result(rep=10)
    #A=SVC()
    #A.fit(Training_set[Independent[:3]],Training_set[['resp.simple']+['Relapse']])
    #A.predict(Scoring_set[Independent[:3]]).T


# In[70]:

if __name__=='__main__':
    from sklearn.ensemble import RandomForestRegressor
    LR=AIPrediction(pivot=Good_variables,method=RandomForestRegressor())
    print LR.result(rep=10)


# In[74]:

if __name__=='__main__':
    from sklearn.linear_model import SGDClassifier
    LR=AIPrediction(pivot=Good_variables,method=SGDClassifier(),binned=True)
    print LR.result(rep=10)


# In[77]:

if __name__=='__main__':
    from sklearn.linear_model import Perceptron
    LR=AIPrediction(pivot=Good_variables,method=Perceptron(),binned=True)
    print LR.result(rep=10)

