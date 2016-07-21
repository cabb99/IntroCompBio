#!/usr/bin/python
import PreProcessing
from Prediction import AIPrediction
import pandas,json

##############
# Parameters #
##############

Predictor='KernelRidge'
Repetition=1
paralel=True
threads=8
max_iterations=1000
Entities=1000  #Number of entities in the population

##################
# Initialization #
##################

#Define Prediction name
Pname='%s_Rep%i'%(Predictor,Repetition)

#Select Dependent Variables
Variables_for_prediction={'pCR':'resp.simple',
                          'pRelapse':'Relapse',
                          'OS':'Overall_Survival',
                          'Remission':'Remission_Duration'}



#max_time_h= 5/60.0#Max number of hours that you can run the simulation

#Select scoring measure to use for minimization
Fast_measure='BAC' #Select between BAC PCC or Scr. You can select two Measures, 
                       #the second measure will overwrite the first one. ('BAC','Auroc')
n_groups =2 #The number of groups tested on each iteration. Select a small number 
            #for the genetic algorithm and a big number to choose 
            #the final Good variables for Prediction

#Slow_measure={'pCR':'BAC',
#              'pRelapse':'BAC',
#              'OS':'BAC',
#              'Remission':'BAC'}

#Define prediction class



###Prediction methodology
#import sklearn

#from sklearn.ensemble import RandomForestRegressor
if Predictor=='RandomForestRegresor':
    from sklearn.ensemble import RandomForestRegressor
    Prediction_Method=RandomForestRegressor
    Prediction_arguments={}
    binned=False
elif Predictor=='RandomForestClassifier':
    from sklearn.ensemble import RandomForestClassifier    
    Prediction_Method=RandomForestClassifier
    Prediction_arguments={}
    binned=True
elif Predictor=='LinearRegression':
    from sklearn.linear_model import LinearRegression      
    Prediction_Method=LinearRegression
    Prediction_arguments={}
    binned=False
elif Predictor=='KernelRidge':
    #Does not work    
    from sklearn.kernel_ridge import KernelRidge 
    Prediction_Method=KernelRidge
    Prediction_arguments={'alpha':'1.0'}
    binned=False




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

#Use previous prefiction bias
try: #Open Variable values if it exists
    Var_value=pandas.DataFrame.from_csv('Var_value_%s.csv'%Pname)
    print 'Found previous variable values'     
except IOError:
    Var_value=pandas.DataFrame.from_csv('Variable_value.csv')
    print 'Did not found previous variable values'
    pass

#Use Previous Population
try: #Open Population if it exists
    with open('Population_%s.json'%Pname) as handle:
        Population=json.loads(handle.read())
    print 'Found previous Population'     
except IOError:
    Population=[{key:biased_sample(Var_value[key].index,random.randint(5,25),Var_value[key]) for key in Var_value.keys()} for i in range(Entities)]
    Generation_method=[{key:'Sample' for key in Var_value.keys()}]*Entities
    print 'Did not found previous Population'
    pass
    
#Use Previous Generation_method

try: #Open Generation_method if it exists
    with open('Generation_method_%s.json'%Pname) as handle:
        Generation_method=json.loads(handle.read())
    print 'Found previous Generation methods'    
except IOError:
    Generation_method=[{key:'Sample' for key in Var_value.keys()}]*Entities
    print 'Did not found previous Generation methods'
    pass
    
#Use Previous Good variables
    
try: #Open Population if it exists
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
    
    if iteration>=max_iterations-1:
        n_groups=20

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
    
    if iteration>=max_iterations:
        break

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
                Method_freq[dep][method]=int(round(Entities/20+points/max_points*(Entities-m*(Entities/20))))
            else:
                Method_freq[dep][method]=min(int(round(Entities/20+points/max_points*(Entities-m*(Entities/20)))),Entities/10)
        for i in range(Entities+10):
            if sum(Method_freq[dep])<Entities:
                Method_freq[dep][random.sample(Method_freq[dep].keys(),1)[0]]+=1
            elif sum(Method_freq[dep])>Entities:
                meth=random.sample(Method_freq[dep][Method_freq[dep]>1].keys(),1)[0]
                if Method_freq[dep][meth]>1:
                    Method_freq[dep][meth]-=1
            else:
                break
            assert i<Entities+8,'Can not obtain just %i'%Entities
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
                    new_variables=biased_sample(ParentA,n,Var_value[dep])
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
                    New_Parent=biased_sample(ParentA,n,Var_value[dep])
                    ps[dep]+=[New_Parent+new_variables]                    
                elif method=='Sample':
                    ps[dep]+=[biased_sample(variables,random.randint(5,25),Var_value[dep])]

    Population=[dict(zip(keys,v)) for v in zip(*[ps[key] for key in keys])]
    Generation_method=[dict(zip(keys,v)) for v in zip(*[ms[key] for key in keys])]
    with open('Population_%s.json'%Pname,'w+') as handle:
        handle.write(json.dumps(Population))
    with open('Generation_method_%s.json'%Pname,'w+') as handle:
        handle.write(json.dumps(Generation_method))
    iteration=iteration+1

Pred=AIPrediction(pivot=Good_variables,method=Prediction_Method(),binned=binned,PredVar=Variables_for_prediction)
Pred.pivot=Good_variables
groups=Pred.create_groups(50)
print Pred.result(groups=groups)
Pred.predict(out='Prediction_%s.csv'%Pname)
