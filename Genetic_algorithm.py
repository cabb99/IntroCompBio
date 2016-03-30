#!/usr/bin/python
import PreProcessing
from Prediction import AIPrediction
import pandas,json

##############
# Parameters #
##############


#Define Prediction name
Pname='LinearRegressionv2'

#Select Dependent Variables
Variables_for_prediction={'pCR':'resp.simple',
                          'pRelapse':'Relapse',
                          'OS':'Overall_Survival_cut',
                          'Remission':'Remission_Duration_cut'}

#Number of entities in the population
Entities=500

#max_time_h= 5/60.0#Max number of hours that you can run the simulation

#Select scoring measure to use for minimization
Fast_measure='Scr' #Select between BAC PCC or Scr. You can select two Measures, 
                       #the second measure will overwrite the first one. ('BAC','Auroc')
n_groups =3 #The number of groups tested on each iteration. Select a small number 
            #for the genetic algorithm and a big number to choose 
            #the final Good variables for Prediction

#Slow_measure={'pCR':'BAC',
#              'pRelapse':'BAC',
#              'OS':'BAC',
#              'Remission':'BAC'}

#Define prediction class
#from sklearn.kernel_ridge import KernelRidge
#Prediction_Method=KernelRidge
#Prediction_arguments={'alpha':'1.0'}

###Linear Regression
from sklearn import datasets, linear_model
Prediction_Method=linear_model.LinearRegression
Prediction_arguments={}
binned=False

###Random Forest Classifier
#from sklearn.ensemble import RandomForestClassifier
#Prediction_Method=RandomForestClassifier
#Prediction_arguments={}
#binned=True
paralel=True
threads=8

###########
# Program #
###########

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


#Read prediction bias
Initial_Var_value=pandas.DataFrame.from_csv('Variable_value.csv')


# # Reading of previous results
# To make the algorithm easier to use, in each iteration it saves the best prediction.
# When it restarts it opens the results of the precious iteration.
# The saved files are:
# - Var_value_???.csv #Contains information about the bias used to test the variables. If you want some variables to not be tested you can use 0. To use some variables more frequently use 1.
# - Population_???.json #The last population created in the simulation
# - Generation_method_???.json #The methods used to create the last population. It allows to calculate how to generate the new population based on its performance.
# - Best_variables_???.json #Contains the best variables used in the last simulation. If the number of groups is small is not recomended to use this variables, allow at least 5 iterations with a big group.
# If the files are not present it will create new variables using a biased sampling.
# It will also include the best variables if present.

if True:#Use previous prefiction bias
    #Open Variable values if it exists
    try:
        Var_value=pandas.DataFrame.from_csv('Var_value_%s.csv'%Pname)
        print 'Found previous variable values'     
    except IOError:
        Var_value=pandas.DataFrame.from_csv('Variable_value.csv')
        print 'Did not found previous variable values'
        pass

if True: #Use Previous Population
    #Open Population if it exists
    try:
        with open('Population_%s.json'%Pname) as handle:
            Population=json.loads(handle.read())
        print 'Found previous Population'     
    except IOError:
        Population=[{key:biased_sample(Var_value[key].index,random.randint(5,25),Var_value[key]) for key in Var_value.keys()} for i in range(Entities)]
        Generation_method=[{key:'Sample' for key in Var_value.keys()}]*Entities
        print 'Did not found previous Population'
        pass
    
if True: #Use Previous Generation_method
    #Open GGeneration_method if it exists
    try:
        with open('Generation_method_%s.json'%Pname) as handle:
            Generation_method=json.loads(handle.read())
        print 'Found previous Generation methods'    
    except IOError:
        Generation_method=[{key:'Sample' for key in Var_value.keys()}]*Entities
        print 'Did not found previous Generation methods'
        pass
    
if True: #Use Previous Good variables
    #Open Population if it exists
    try:
        with open('Good_variables_%s.json'%Pname) as handle:
            Good_variables=json.loads(handle.read())
            Population+=[Good_variables]
            Generation_method+=[{key:'Elitist' for key in Var_value.keys()}]
        print 'Found previous Good variables'    
    except IOError:
        print 'Did not found previous Good variables'
        pass


# # Genetic algorithm
# 
# Evaluates the prediction a small number of times (10 times) and sorts the variables as a function of the score.
# The variables with the best scores are then transformed in the following ways:
# 1) Addition of variables
# Adds 1 to 5 new variables chosen from a biased sampling.
# 2) Reduction of variables
# Removes 1 to 5 variables chosen by a biased sampling. If the number of variables is 5 it does not remove any more variable.
# 3) Change of variables
# Changes 1 to 3 variables.
# 4) Cross-over
# Selects two sets of variables and selects some variables from the two sets with a biased sampling. The lenght of the variables is at least the number of variables in the smallest set.
# 5) Elitist
# It selects without change the bests sets.
# 
# The probability of selecting one set of variables is given by the function:
# 
# random^3
# 
# We use an adaptative genetic algorithm, and the method to create the new population is calculated by evaluating the performance of the last iteration. Each methods receive a number of points, evaluated by the following function:
# 
# \begin{equation} e^\frac{-9i}{l} \end{equation}
# 
# Where i is the ranking in the sorted population and l the number of entities in the population.
# 
# We also specify that a maximum number of deletions and elitist selection to force the evaluation of new sets.
# There is also a minimum number of new sets created by each method.

# In[ ]:

import time
keys=Variables_for_prediction.keys()
keys.sort()

methods=['Crossover','Mutation','Deletion','Addition','Sample','Elitist']
Method_freq=pandas.concat([pandas.Series([0 for i in range(len(methods))],name=key,index=methods) for key in keys],axis=1)

try:
    iteration=0
    with open('Results_%s.log'%Pname) as Log:
        for line in Log:
            if len(line)>0 and line[0]<>'#':
                iteration=int(line.split()[0])
                
except IOError:
    iteration=0

with open('Results_%s.log'%Pname,'a+') as Log:
    Log.write('#%s\n'%time.ctime())
    Log.write('#%s %s %s %s %s\n'%tuple(['iteration']+keys))

def pool_process(x):
    pivot=x[0]
    groups=x[1]
    #print pivot
    #print len(groups)
    Pred=AIPrediction(pivot=pivot,method=Prediction_Method(*Prediction_arguments),binned=binned,PredVar=Variables_for_prediction)      
    try:
        result=Pred.accuracy(groups=groups,Measure=Fast_measure)
    except KeyError:
        result=numpy.nan
    return result

if paralel:
    from multiprocessing import Pool    
    pool = Pool(threads)
    
while True:

    #Evaluate population
    Pred=AIPrediction(method=Prediction_Method(),binned=binned)
    groups=Pred.create_groups(n_groups)
    #Best_scores=[]
    
    if paralel:
        Results=pool.map(pool_process, zip(Population,[groups]*len(Population)))
    else:
        Results=[]        
        for arg in zip(Population,[groups]*len(Population)):
            Results+=[pool_process(arg)]
    #print Results
    if numpy.nan in Results:
        continue    
    Sorted_population={}
    Sorted_generation_method={}
    Best_result={}
    Average_result={}
    Good_variables={}
    for key in keys:
        key_results=[(r[key],p[key],m[key]) for r,p,m in zip(Results,Population,Generation_method)]
        key_results.sort(reverse=True)
        Sorted_population.update({key:[p for r,p,m in key_results]})
        Sorted_generation_method.update({key:[m for r,p,m in key_results]})
        Best_result.update({key:key_results[0][0]})
        Average_result.update({key:numpy.mean([r[0] for r in key_results])})
        Good_variables.update({key:key_results[0][1]})
    with open('Results_%s.log'%Pname,'a+') as Log:
        Log.write('%i %4.6f %4.6f %4.6f %4.6f\n'%tuple([iteration]+[Best_result[key] for key in keys]))
        print '%i %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f'%tuple([iteration]+[Best_result[key] for key in keys]+[Average_result[key] for key in keys])
    with open('Good_variables_%s.json'%Pname,'w+') as handle:
        handle.write(json.dumps(Good_variables))
    #Adaptative genetic algorithm (the rates of crossover and mutation vary according to their efficacy)
    #print Method_freq
    import numpy
    def Fpoints(x,max_x):
        return numpy.exp(-x/float(max_x)*9.0)

    for dep in keys:
        l=float(len(Sorted_generation_method[dep]))
        m=float(len(Method_freq[dep].keys()))
        max_points=sum(Fpoints(i,l) for i in range(int(l)))
        for method in Method_freq[dep].keys():
            points=0
            c=0
            for i,gm in enumerate(Sorted_generation_method[dep]): 
                if gm==method:
                    points+=Fpoints(i,l)
                    #print l,m
                    c=c+1
            #print points
            #expected_points=max_points/float(c) if c>0 else 0
            if method<>'Elitist' and method<>'Deletion':
                Method_freq[dep][method]=int(round(Entities/20+points/max_points*(l-m*(Entities/20))))
            else:
                Method_freq[dep][method]=min(int(round(Entities/20+points/max_points*(l-m*(Entities/20)))),Entities/10)
        for i in range(Entities+10):
            if sum(Method_freq[dep])<Entities:
                Method_freq[dep][random.sample(Method_freq[dep].keys(),1)[0]]+=1
            elif sum(Method_freq[dep])>Entities:
                meth=random.sample(Method_freq[dep].keys(),1)[0]
                if Method_freq[dep][meth]>1:
                    Method_freq[dep][meth]-=1
            else:
                break
            assert i<Entities+8,'Can not obtain just %i'%meth
    #print method,points
    #Method_freq[dep][method]=1+int(points/max_points/c*(Entities-2*len(Method_freq[dep].keys())))

    #Generate new population
    ps={dep:[] for dep in Method_freq.keys()}
    ms={dep:[] for dep in Method_freq.keys()}
    variables=list(Var_value.index)
    for dep in keys:
        for method in Method_freq[dep].keys():
            for i in range(int(Method_freq[dep][method])):
                ms[dep]+=[method]
                if method=='Elitist':
                    #Select the bests predictions
                    ps[dep]+=[Sorted_population[dep][i]]
                elif method=='Crossover':
                    #Select two variables groups and cross them
                    A=int((random.random()**2)*len(Sorted_population[dep])/2)
                    B=A
                    while B==A:
                        B=int((random.random()**2)*len(Sorted_population[dep])/2)
                    ParentA=Sorted_population[dep][A]
                    ParentB=Sorted_population[dep][B]
                    min_l=min(len(ParentA),len(ParentB))
                    max_l=len(set(ParentA+ParentB))
                    New_vars=random.sample(list(set(ParentA+ParentB)),random.randint(min_l,max_l))
                    ps[dep]+=[New_vars]
                elif method=='Addition':
                    #Select one variable group and add  some new variables
                    A=int((random.random()**3)*len(Sorted_population[dep])/2)
                    ParentA=Sorted_population[dep][A] #Select the variable groups
                    var2add=[v for v in Var_value[dep].index if v not in ParentA] #Variables to add
                    l=len(var2add)
                    n=min(random.randint(1,5),l)
                    new_variables=biased_sample(var2add,n,Var_value[dep])
                    #new_variables=Good_variables[key]+change_variables[key]
                    #Test_variables.update({key:new_variables})
                    ps[dep]+=[ParentA+new_variables]
                elif method=='Deletion':
                    #Select one variable group and add  some new variables
                    A=int((random.random()**3)*len(Sorted_population[dep])/2)
                    ParentA=Sorted_population[dep][A] #Select the variable groups
                    l=len(ParentA)
                    n=max(l-random.randint(1,5),5)
                    new_variables=biased_sample(ParentA,n,1/(Var_value[dep]+1))
                    ps[dep]+=[new_variables]
                elif method=='Mutation':
                    #Select one variable group and add and remove the same number of variables
                    #Select one variable group and add  some new variables
                    A=int((random.random()**3)*len(Sorted_population[dep])/2)
                    ParentA=Sorted_population[dep][A] #Select the variable groups
                    var2add=[v for v in Var_value[dep].index if v not in ParentA]#Variables to add
                    l=len(ParentA)
                    n=max(l-random.randint(1,3),3)
                    new_variables=biased_sample(var2add,l-n,Var_value[dep])
                    New_Parent=biased_sample(ParentA,n,1/(Var_value[dep]+1))
                    ps[dep]+=[New_Parent+new_variables]                    
                elif method=='Sample':
                    ps[dep]+=[biased_sample(variables,random.randint(5,25),Var_value[dep])]

    Population=[dict(zip(ps.keys(),v)) for v in zip(*[ps[key] for key in keys])]
    Generation_method=[dict(zip(ms.keys(),v)) for v in zip(*[ms[key] for key in keys])]
    with open('Population_%s.json'%Pname,'w+') as handle:
        handle.write(json.dumps(Population))
    with open('Generation_method_%s.json'%Pname,'w+') as handle:
        handle.write(json.dumps(Generation_method))
    iteration=iteration+1
    
    


# In[ ]:

print 'iter',i
Pred.pivot=Good_variables
print Pred.result(rep=10)
Pred.predict()


# In[148]:

Method_freq


# In[131]:

for p,m in zip(ps['pRelapse'],ms['pRelapse']):
    print len(p)-len(set(p)),m
len(Population),len(Generation_method)


# In[8]:

Sorted_generation_method['Remission']=[random.sample(['Crossover','Mutation','Deletion','Addition','Sample','Elitist'],1)[0] for i in range(50)]


# In[70]:

for dep in Method_freq.keys():
    l=float(len(Sorted_generation_method[dep]))
    m=float(len(Method_freq[dep].keys()))
    max_points=sum(Fpoints(i,l) for i in range(int(l)))
    for method in Method_freq[dep].keys():
        points=0
        c=0
        for i,gm in enumerate(Sorted_generation_method[dep]): 
            if gm==method:
                points+=Fpoints(i,l)
                print l,m
                c=c+1
        #print points
        #expected_points=max_points/float(c) if c>0 else 0
        Method_freq[dep][method]=int(round(1+points/max_points*(l-m)))
    for i in range(Entities+10):
        if sum(Method_freq[dep])<Entities:
            Method_freq[dep][random.sample(Method_freq[dep].keys(),1)[0]]+=1
        elif sum(Method_freq[dep])>Entities:
            meth=random.sample(Method_freq[dep].keys(),1)[0]
            if Method_freq[dep][meth]>1:
                Method_freq[dep][meth]-=1
        else:
            break
        assert i<Entities+8,'Can not obtain just %i'%meth
Method_freq
        


# In[39]:

T=[]
for i in range(10000):
    T+=[int(random.random()**3*50/2)]
import matplotlib as mpl
get_ipython().magic(u'matplotlib inline')
mpl.pyplot.hist(T,bins=25)


# In[94]:

import time
print time.ctime()
dir(time)


# In[6]:

#Define Best Prediction
Pred.pivot=Good_variables
groups=Pred.create_groups(10)
Best_scores=Pred.score(groups=groups,Measure=tuple(set(Slow_measure.values())))
Best_result=Pred.accuracy(groups=groups,Measure=Fast_measure)


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
        Test_result=Pred.accuracy(groups=groups[:2],Measure=Fast_measure) #Fast measure
    except KeyError:
        groups=Pred.create_groups(10)
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
            if Var_value[key][var]<0.02:
                Var_value[key][var]=0.02
            if Var_value[key][var]>2.0:
                Var_value[key][var]=2.0
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
        groups=Pred.create_groups(10)
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

