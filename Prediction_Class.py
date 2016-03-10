# # Prediction
# In order to predict I created a Prediction Class. The good thing about Classes is that you can create new classes that can heredate the propierties and methods of other Classes. That means that we have don't have to write this part of the code again.
# 
# ## Group creation
# 
# The class contains a function to divide the 
# 
#     Prediction.split(n)
#     
# Splits the training set in n parts. n-1 parts will go to the training set, while 1 part will be on the testing set. We create n groups, so each of the n parts can be in the testing set. For example if we split the testing sets in 3 parts (n=3) We will be returned 3 groups that will contain:
# 
# |Group|Testing|Training|
# |-----|-------|--------|
# |0    |1      |2,3     |
# |1    |2      |1,3     |
# |2    |3      |1,2     |
# 
# The output is a list of dictionaries. Here are some examples of how to call the groups:
#     
#     groups=Prediction.split(n) #Assign the group to a variable
#         [{'test':<Table>,'train'<Table>},{'test':<Table>,'train'<Table>},...]
#     g=groups[0] #Returns the frist dictionary
#         {'test':<Table>,'train'<Table>}
#     g['train'] #Returns the training set
# 
#     Prediction.create_groups(rep,[n])
#     
# Will create n groups rep times. This method will invoke split rep times creating a bigger list. It is useful to have more testing and training groups to have a more reliable scoring.
#     
# ## Predicts
# 
# ### Main prediction function
#     Prediction.predict(scoring=Scoring_set,out='prediction.csv')
#     
# Will give results about the Scoring set using the full training set with the prediction methods used.
# 
# ### Prediction methods
# 
#     Prediction.pCR(training,testing): 
# Returns a Series containing the probability of Complete Remission ('resp.simple')
#         
#     Prediction.pRelapse(training,testing): 
# Returns a Series containing the probability of Relapse ('Relapse')
#         
#     Prediction.Remission(training,testing): 
# Returns a Series with the estimated Remission Duration ('Remission_Duration') in weeks (not binned)
#         
#     Prediction.OS(training,testing): 
# Returns the estimated Overall Survival ('Overall_Survival') in weeks (not binned)
# 
# ### Reduced prediction methods.
# The prediction algorithms for pCR and pRelapse should be similar, as well as the algorithms gor Remission and Overall Survival. We can have only two prediction algorithms: pPred for the probabilities prediction and qPred for the time Prediction
# 
#     Prediction.pPred(training,testing,var):
# Returns the Probability prediction as a Series. 
# Should return a value between 0 and 1.
# Used for pCR and pRelapse
#         
#     Prediction.qPred(self,training,testing,var):
# Quantitative prediction, returns the number of weeks.
# Used for Remission and OS
# It will be binned after this, so if this is a binned prediction, multiply the result by 52
#  
#     Prediction._bin(data,bins=[0,52,104])
# Bins the results as 0,1 2,3. The 0 values correspond to numpy.nan values.
# 
# ## Scoring functions
#     Prediction.accuracy(self,rep=1,groups=[]):
# Returns the Balanced accuracy prediction for all the functions
#     
#     Prediction.result(self,groups=[],alpha=0.5,Measure='All',rep=1):    
# Returns the scoring calculations for the prediction as a Table
# You can select wich scoring to return with Measure
# alpha allows to select the scoring that is lower than alpha*100% of the scores, supossing a normal distribution.
# If alpha is 0.5 it will return the mean
# 
#     Prediction.score(self,groups=[],Measure='All'):
# Main scoring function
# Returns a dictionary for all the scores
# You can select a Measure
# Will create groups if no group given.
# 
# ### Balanced accuracy (BAC)
# The Balanced accuracy is an improved scoring function over the accuracy. While the accuracy measures the number of correct gueses over the entire number of guesses, the balanced accuracy measures the number of correct guesses in each subgroup of possible outcomes and then obtains a mean. For example in a population where P=99 and N=1, a random guess where we say that all the population is P would give 99% accuracy, while only having 50% balanced accuracy.
# 
# |Total population| Positive (P) | Negative (N) |
# |--|--|--|
# |Predicted positive| True positive (TP)| False positive (FP)| 
# |Predicted necative| False negative(FN)| True negative (TN)|
# 
# #### Probabilities
# For the Balanced accuracy we define a that the predicted positive values are all those values that have a probability greater or equal than 0.5, the rest being predicted negative.
# 
# $$BAC=\frac{TP}{P}+\frac{TN}{N}$$
# 
# #### Categorical
# BAC is not strictly defined for continuous variables, so to estimate this scoring method we used the binned results. We have three variables, so a True negative is not defined. We define TP as correct guesses, where the test bin is the same as the predicted bin.
# 
# $$BAC=\frac{1}{n} \sum_{n} \frac{TP}{P}$$
# 
# ### Area under the receiver operating characteristic curve (AUROC)
# AUROC is a complimentary analysis to BAC. The ROC curve is a curve that compares the tradeoff between sensitivity and specificity. It create a cutoff by supossing a cutoff (k) on the prediction, where values greater or equal than k will be considered positive, and values smaller than k are considered negative. 
# 
# At k=1 all the predictions are negative, and there are no false negatives nor True positives (We are at point (0,0). As k decreases the number of True positives also increases. 
# 
# If there is a k value where all True positives have been found and there are no false positives, we will be at the point (1,0). If this point is reached AUROC will be 1. 
# 
# When k=0 all predictions are positive and we will arrive at point (1,1).
# 
# AUROC measures the area under this curve and is not defined for the Categorical dependent variables, since we don't have a parameter "k" to modify.
# <img src=https://upload.wikimedia.org/wikipedia/commons/8/8c/Receiver_Operating_Characteristic.png>
# 
# ### Pearson correlation coefficient (PCC)
# Measures the correlation between the predicted values and the actual values. It originally return a value between -1 and 1, but is normalized to be between 0 and 1. For Remission duration and Overall Survival it computes with the real results (not binned). If your prediction return binned results you can unbin them by multiplying by 52 and adding 26.
# 
# $$PCC = \frac{\sum_{i=1}^{n}(p_i-\overline{p})(a_i-\overline{a})}{\sqrt{\sum_{i=1}^{n}(p_i-\overline{p})^2}\sqrt{\sum_{i=1}^{n}(a_i-\overline{a})^2}}$$
# 
# $$PCC_{norm} = (PCC + 1)/2$$
# 
# ### Concordance Index (CI)
# For each pair of results, it calculates if the predicted results are in the same order than the expected results. For example: for a pair (51,60),(64,62), where 51 and 64 are the predicted values, and 60 and 62 are the expected values. If $51<64$ and $60<62$ (The same relationsip), then $h_{ij}=0$.
# 
# $$\sum_{i<j}^{N}h(i,j)$$
# 
# The pairs that are counted are not censored pairs. A pair is censored if the Overall survival of i is lower than j, but i is alive. If j is equal to i and i or j is alive it is also censored. Nevertheless, if the Overall Survival of j is lower and i and j is alive it is not censored (...Please tell me if you can understand that). A similar relation is applied to Relapse Status with Remission, so it does not matter whether the patients were alive or dead.

# In[9]:

import numpy
import random
import pandas
from scipy.stats import norm

#Open the data and read in pandas
Dream9_training=pandas.read_excel('Dream9.xlsx',"trainingData")
Dream9_scoring=pandas.read_excel('Dream9.xlsx',"scoringData")
Dream9=pandas.concat([Dream9_training,Dream9_scoring])

#Create the new tables
from PreProcessing import PreProcess
Training_set=PreProcess(Dream9_training,Dream9)
Scoring_set=PreProcess(Dream9_scoring,Dream9)

Variables_for_prediction={'pCR':'resp.simple','pRelapse':'Relapse','OS':'Overall_Survival','Remission':'Remission_Duration'}

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
        digitized = numpy.digitize(data, bins)
        for i,v in enumerate(data):
            if numpy.isnan(v):
                digitized[i]=0
        return digitized

    
    ######################
    ## Scoring methods  ##
    ######################
    
    def accuracy(self,rep=1,groups=[]):
        #Returns the Balanced accuracy prediction for all the functions
        groups=self.create_groups(rep) if groups==[] else groups
        S=self.score(Measure='BAC',groups=groups)
        return {s:numpy.mean(S[s]['BAC']) for s in S}
    
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
        [S[v].update({'BAC':[0 for r in lg],'PCC':[0 for r in lg]})for v in S]
        [S[v].update({'AUROC':[0 for r in lg]})for v in ['pCR','pRelapse']]
        [S[v].update({'CI':[0 for r in lg]})for v in ['OS','Remission']]
        #Independent=[v for v in self.training.keys() if v not in Dependent]
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

if __name__=='__main__':
   
    #This part is to test the function
    Dummy=Prediction()
    print Dummy.result(alpha=0.5,Measure='All',rep=5)
    print Dummy.__doc__    
    print dir(Dummy)



