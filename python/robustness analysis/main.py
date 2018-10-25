from SALib.sample.latin import sample
import numpy as np
import math
import simple_model
import FinalDcas_doses_without_rate as dcas
import timeit
from multiprocessing import Pool, TimeoutError
from multiprocessing import Process,Queue, current_process, freeze_support
import pprint



name= 'dcas_ky_5000'
c=[5,100,5]
PROC_NUM=4

def main(name,c):
    name = name.lower()
    model, par, parameters, bounds, samples_to_take = setups(name)
    problem = {'num_vars' : len(parameters),
               'names' : parameters,
               'bounds' : bounds
                }
    params = sample(problem, samples_to_take)

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
        p.append(Process(target=evaluate,args=(paramsdiv[i],model,par,c,done,i)))
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


    #Y = evaluate(params, model, par, c)

    write_file(Y, parameters, name, c, params)

def evaluate(values,model,par, c, done, no):
    Y = np.zeros(( values.shape[0], int((c[1]/c[0])) ))
    print("Start evaluation method with " + str(len(values)) + " parameter sets")
    start = timeit.default_timer()
    counter=0
    for i, X in enumerate(values):
        counter+=1
        print(counter)
        if model == 'tale': Y[i]=simple_model.main(c, par, X[0],X[1],X[2])
        if model == 'dcas': Y[i]=dcas.main(c, par,X[0],X[1])
        if i%100==0 and i != 0 :
            time_100= timeit.default_timer() - start
            remaining_time=(int((len(values)-i)/100)*time_100)/60
            print('progress: ' + str((i/float(len(values)))*100) + ' %')
            print('time for 100 calculations: ' + str(time_100) + ' sec')
            print('remaining time: ' + str(remaining_time) + ' min\n')
            start = timeit.default_timer()
    done.put([Y,no])

def write_file(Y, parameters, name, c, params):
    filename = name + '.csv'
    columns=(",".join(parameters))
    for i in range(c[0],c[1]+1,c[2]):
        columns+= ',' + 'G/c=' + str(i)
    columns+= ',E,S'
    file = open(filename, 'w')
    file.write(columns + '\n')

    mean,std=mean_std(Y)
    for i,x in enumerate(Y):
        E,S= Error_Strength(x)
        file.write(",".join(map(str, params[i])))
        file.write(',' + ",".join(map(str, x)) + ',')
        file.write(str(E) + ',' + str(S) + "\n")
    write2(file,mean,Y,parameters)
    write2(file,std,Y,parameters)
    file.close()

def write2(file,x,Y,parameters):
    for i in range(len(parameters)):file.write(',')
    for i in range (Y.shape[1]): file.write(str(x[i]) + ',')
    file.write('\n')

def mean_std(Y):
    mean=[]
    std=[]
    for i in range(Y.shape[1]):
        mean.append(np.mean(Y[:,i]))
        if i == 3:
            counter=0
            for j in range (Y.shape[0]):
                counter+=j
        std.append(np.std(Y[:,i]))
    return mean,std

def Error_Strength(goi):
    E=(abs((goi[-1]-goi[0]))/goi[0])
    return E, E/(E+1)

def setups(conf):
    model,par,samp = conf.split('_')
    sample=samp
    if model=='dcas':
        if par == 'dcasred':
            parameters= ['kdton','kdtoff','kon2','koff2']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10]]
        if par == 'dcas':
            parameters = ['kdton','kdtoff','kon2','koff2','Dox','at','bt','dmt','dt','dmd','dd','Tmax']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [1, 50],
                 [0.1, 5],
                 [0.1, 5],
                [0.001,1],
                [0.001,1],
                [0.001,1],
                [0.001,1],
                [0.1, 5]]
        if par == 'sgrna':
            parameters= ['khron','khroff','kon1','koff1','dh']
            b= [[0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001, 10],
                 [0.001,0.3]]
        if par == 'ky':
            parameters= ['kon1','Rham']
            b= [[0.00001, 0.01],
                [20,200]]
    elif model=='tale':
        if par == 'diffn':
            parameters= ['kon','koff','aGmax']
            b= [[0.001,10],
                [0.001,10],
                [1,5]]#less than one is unrealistic for sfGFP
        if par == 'all':
            parameters= ['kon','koff','yR','yG','yM','aR','bR']
            b= [[0.001,10],
                 [0.001,10],
                 [0.001,0.3],
                 [0.001,0.3],
                 [0.001,0.3],
                 [0.1,5],
                [0.1,5]]

        elif par == 'ky':
            parameters = ['kon', 'koff', 'yR']
            b = [[0.001,10],
                 [0.001,10],
                 [0.001,0.3]]
    return model,par,parameters,b,int(sample)
if __name__ == '__main__':
    main(name,c)
