from SALib.sample.saltelli import sample
from SALib.analyze.sobol import analyze
from SALib.test_functions import Ishigami
from SALib.analyze import sobol
import numpy as np
import math
import tale as simple_model
import dcas9 as dcas
import pprint
import timeit
import random
import matplotlib.pyplot as plt
from multiprocessing import Pool, TimeoutError
from multiprocessing import Process,Queue, current_process, freeze_support





#-------SETUP---------------------------------------------------
#[modelname]_[setup]_[samplesize]
#models available: tale, dcas
#setups available: tale->all, ky
#                  dcas->all, main, dcas, sgrna, ky
#samples: int
name = 'dcas_ky_2000'
logspace=False
PROC_NUM=4
c=[5,100,5]
#dont forget to change the parameters that you send to your model
#---------------------------------------------------------------

def main(name,logspace):

    name = name.lower()
    model, par, parameters, bounds, samples_to_take = setups(name)
    problem = {'num_vars' : len(parameters),
               'names' : parameters,
               'bounds' : bounds
                }

    params = sample(problem, samples_to_take, calc_second_order=True)
    print(len(params))
    if logspace:
        params = log_distribute(params)

    # Create queues
    done = Queue()
    paramsdiv=[]
    p=[]
    remaining=0
    for i in range (0,PROC_NUM):
        x=params.shape[0]/PROC_NUM
        remaining+= x - int(x)
        x=int(x)
        y=params[x*i:x*(i+1),:] #filling the y with the x parameter sets

        if i == PROC_NUM-1 and remaining!=0:
            y=np.append(y,params[-int(remaining):,:], axis=0)
        paramsdiv.append(y)
        p.append(Process(target=evaluate,args=(paramsdiv[i],model,par,done,i)))
        p[i].start()
    print(np.asarray(paramsdiv).shape)

    Ydiv=np.zeros(PROC_NUM).tolist() #no of processes
    count=0
    while True:
        temp,no =done.get()
        Ydiv[no]=temp.tolist()
        count+=1
        if count==PROC_NUM:
            break
    Y=[]
    for i in range(0,PROC_NUM):
        Y+=Ydiv[i]

    Y=np.asarray(Y)
    print(Y.shape)
    #print(len(params[:,2]))
    #plt.hist(params[:,2], bins='auto')  # arguments are passed to np.histogram
    #plt.title("Histogram with 'auto' bins")
    #plt.show()

    Si = sobol.analyze(problem, Y, print_to_console=True, calc_second_order=True)
    print(Si)
    write_file(Si,parameters,name,logspace)

'''
model= 'tale'
parameters= ['kon','koff','yR','yG','yM','aR','aGmax','bR','bG']
bounds= [[0.001,10],
         [0.001,10],
         [0.01,1],
         [0.01,1],
         [0.01,1],
         [0.1,5],
        [0.1,5],
        [0.01,5],
        [0.01,5]]
samples_to_take = 10000
#dont forget to change the parameters that you send to your model'''

def log_distribute(params):
    new=[]
    for set in params:
        temp=[]
        for i,x in enumerate(set):
            temp.append(math.log(x,10))
        new.append(temp)
    #mean and variance
    mean=[]
    std=[]
    for i in range(len(params[0])):
        mean.append(np.mean(new[:][i]))
        std.append(np.std(new[:][i]))
    #normal variate
    for i,x in enumerate(new):
        for j,y in enumerate(x):
            Z = random.gauss(0, 1)
            new[i][j]=10**(mean[j] + std[j]*Z)
    return np.asarray(new)


def evaluate(values,model,par,done,no):
    Y = np.zeros([values.shape[0]])
    print("Start evaluation method with " + str(len(values)) + " parameter sets")
    start = timeit.default_timer()
    for i, X in enumerate(values):
        if model == 'tale': Y[i]=simple_model.main(X[0],X[1],X[2])
        if model == 'dcas': Y[i]=dcas.main(c,par,X[0],X[1],X[2],X[3],X[4],X[5])
        if i%100==0 and i != 0 :
            time_100= timeit.default_timer() - start
            remaining_time=(int((len(values)-i)/100)*time_100)/60
            print('progress: ' + str((i/float(len(values)))*100) + ' %')
            print('time for 100 calculations: ' + str(time_100) + ' sec')
            print('remaining time: ' + str(remaining_time) + ' min\n')
            start = timeit.default_timer()
    done.put([Y,no])

def write_file(Si, parameters, name, logspace):
    if logspace==False:filename = name + '.csv'
    else: filename = name + '_logspace.csv'

    columns=("Param1,Param2,S1,S1con,S2,S2con,ST,STcon")
    file = open(filename, 'w')
    file.write(columns + '\n')
    S1=Si['S1']
    S1con=Si['S1_conf']
    ST=Si['ST']
    STcon=Si['ST_conf']
    S2=Si['S2'].tolist()
    S2con=Si['S2_conf'].tolist()
    for i,x in enumerate(S1):
        file.write(parameters[i] + ', - ,' + str(S1[i]) + ',' + str(S1con[i]) + ', - , - ,' + str(ST[i]) + ',' + str(STcon[i]) + '\n')
    for i in range(0,len(parameters)):
        for j in range(0,len(parameters)):
            if math.isnan(S2[i][j]) == False:
                file.write(parameters[i] + ',' + parameters[j] + ', - , - ,' + str(S2[i][j]) + ',' + str(S2con[i][j]) + ', - , - ,' + '\n')
    file.close()


def setups(conf):
    model,par,samp = conf.split('_')
    sample=samp
    if model=='dcas':
        if par == 'dcas':
            parameters= ['kdton','kdtoff','kon2','koff2','Dox','at','bt','bd','dmt','dt','dmd','dd','Tmax']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [1, 50],
                 [0.1, 5],
                 [0.1, 5],
                 [0.1, 5],
                [0.1,0.3],
                [0.001,0.3],
                [0.1,0.3],
                [0.001,0.3],
                [0.1, 5]]
        if par == 'dcasred':
            parameters= ['kdton','kdtoff','kon2','koff2']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10]]
        if par == 'sgrna':
            parameters= ['khron','khroff','kon1','koff1','Rham','ah','bh','dmh','dh','ds','hmax']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [95000, 100000],
                 [0.1, 5],
                 [0.1, 5],
                [0.1,0.3],
                [0.001,0.3],
                [0.001,0.3],
                [0.1, 5]]
        if par == 'ky':
            parameters= ['kdson','kdsoff','kon3','koff3','kon1','Rham'] #['kdton','kdson','kdsoff','kdtoff','kon2','koff2','khron','kon3','koff3','khroff','kon1','koff1','dh','ds','dd']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.00001, 10],
                 [50, 175]]
    elif model=='tale':
        if par == 'all':
            parameters= ['kon','koff','yR','yG','yMR','yMG','aR','aGmax','bR','bG']
            b= [[0.001,10],
                 [0.001,10],
                 [0.001,0.3],
                 [0.001,0.3],
                 [0.1,0.3],
                 [0.1,0.3],
                 [0.1,5],
                [0.1,5],
                [0.01,5],
                [0.01,5]]

        elif par == 'ky':
            parameters = ['kon', 'koff', 'yR']
            b = [[0.001,10],
                 [0.001,10],
                 [0.001,0.3]]
    return model,par,parameters,b,int(sample)
if __name__ == '__main__':
    main(name,logspace)
