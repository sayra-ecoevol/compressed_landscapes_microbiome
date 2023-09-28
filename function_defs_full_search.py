#!/usr/bin/env python
# coding: utf-8
'''
@author : sayra_ecoevol
'''
from scipy import *
import numpy as np
import time, math, json, pprint, requests, urllib, io, os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
colors = sns.color_palette()
from matplotlib.backends import backend_pdf as bpd
from community_simulator import *
from community_simulator import Community,essentialtools,usertools,visualization,analysis
from community_simulator.usertools import *
from community_simulator.visualization import *
from itertools import *
from random import *

def validate_simulation_sa743(com_in,N0_init,pro_time,m,step):
    """
    Check accuracy, convergence, and noninvadability of community instance com_in.
    N0 indicates which species were present at the beginning of the simulation.
    """
    com = com_in.copy()
    failures = np.sum(np.isnan(com.N.iloc[0]))
    survive = com.N>0
    com.N[survive] = 1
    if type(com.params) is not list:
        params_list = [com.params for k in range(len(com.N.T))]
    else:
        params_list = com.params
    dlogNdt_survive = pd.DataFrame(np.asarray(list(map(com.dNdt,com.N.T.values,com.R.T.values,params_list))).T,
                                   index=com.N.index,columns=com.N.columns)
    dlogNdt_survive.to_csv(r"data_acc/dNdt_for_step =" +str(step)+"_after_time_"+str(pro_time)+ "for gen = " +str(m)+".csv", sep='\t',float_format='%.4f')
    dlogNdt_survive.to_html(r"html/datah_acc/dNdt_for_step =" +str(step)+"_after_time_"+str(pro_time)+ "for gen = " +str(m)+".html")
    accuracy = np.max(abs(dlogNdt_survive))
    com.N[N0_init>0] = 1
    com.N[survive] = 0
    dlogNdt_extinct = pd.DataFrame(np.asarray(list(map(com.dNdt,com.N.T.values,com.R.T.values,params_list))).T,
                               index=com.N.index,columns=com.N.columns)

    invaders = np.sum(dlogNdt_extinct>0)
    return accuracy.mean(),(invaders>0).sum()

def readassumptions(filepath):
    file = open(r"psoc_setup.txt", "r")
    contents = file.read()
    dictionary0 = eval(contents)
    file.close()
    n_families = dictionary0['n_families'];
    n_species = dictionary0['n_species'];
    n_classes =dictionary0['n_classes'];
    n_resources = dictionary0['n_resources'];
    maintenance = dictionary0['maintenance'];
    tau = dictionary0['tau'];
    nspecies_sampled =dictionary0['nspecies_sampled'];
    scale_set = dictionary0['scale_set'] #this is needed for conversion to actual cell counts from the coded abundances
    iteri = dictionary0['iterations']
    propagation_time = dictionary0['prop_times']
    #assumptions in the model
    file = open(r"psoc_assump.txt", "r")
    contents = file.read()
    a_default = eval(contents)
    file.close()
    global assumptions
    assumptions = a_default
    print(assumptions)
    return assumptions, iteri, maintenance, tau, scale_set

def getdata():
    user='sayra-ecoevol'
    pao='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'
    github_session = requests.Session()
    github_session.auth = (user, pao)
    csv_url = r"https://raw.githubusercontent.com/sayra-ecoevol/CFO2020-/main/consumer-choice/phi_global_truncated.csv"
    download = github_session.get(csv_url).content
    commfunction1 = pd.read_csv(io.StringIO(download.decode('utf-8')), error_bad_lines=False)
    csv_url = r"https://raw.githubusercontent.com/sayra-ecoevol/CFO2020-/main/consumer-choice/resource_stoichiometric_matrix_20_0.3.csv"
    download = github_session.get(csv_url).content
    Dd = pd.read_csv(io.StringIO(download.decode('utf-8')), error_bad_lines=False,header=None)
    csv_url = r"https://raw.githubusercontent.com/sayra-ecoevol/CFO2020-/main/consumer-choice/alpha_global_truncated.csv"
    download = github_session.get(csv_url).content
    commfunction2 = pd.read_csv(io.StringIO(download.decode('utf-8')), error_bad_lines=False)
    return commfunction1, Dd, commfunction2

def mkdirects(okay_path,muuc,siggc,replicaa):
        for folder in ["dataN","dataR","dataS","data_acc","graphs","para", "html","results"]:
            path = str(okay_path)+"/"+str(muuc)+"_"+str(siggc)+"_"+str(replicaa)+"/"+str(folder)
            try:
                os.makedirs(path)
            except OSError:
                print ("yolo");
            else:
                foo=1
            for folder2 in ["resh","datahN","datahR","datahS","datah_acc","parah"]:
                path = str(okay_path)+"/"+str(muuc)+"_"+str(siggc)+"_"+str(replicaa)+"/html/" + str(folder2)
                try:
                    os.makedirs(path)
                except OSError:
                    foo=0
                else:
                    foo=1
def mkdirects2(okay_path,muuc,siggc,replicaa):
    for folder in ["dataN","dataR","dataS","data_acc","graphs","para", "html","results"]:
        path = str(okay_path)+"/type2_"+str(muuc)+"_"+str(siggc)+"_"+str(replicaa)+"/"+str(folder)
        try:
            os.makedirs(path)
        except OSError:
            print ("yolo");
        else:
            foo=1
        for folder2 in ["resh","datahN","datahR","datahS","datah_acc","parah"]:
            path = str(okay_path)+"/type2_"+str(muuc)+"_"+str(siggc)+"_"+str(replicaa)+"/html/" + str(folder2)
            try:
                os.makedirs(path)
            except OSError:
                foo=0
            else:
                foo=1



def mkdirects3(okay_path,muuc,siggc,replicaa):
    for folder in ["dataN","dataR","dataS","data_acc","graphs","para", "html","results"]:
        path = str(okay_path)+"/type3_"+str(muuc)+"_"+str(siggc)+"_"+str(replicaa)+"/"+str(folder)
        try:
            os.makedirs(path)
        except OSError:
            print ("yolo");
        else:
            foo=1
        for folder2 in ["resh","datahN","datahR","datahS","datah_acc","parah"]:
            path = str(okay_path)+"/type3_"+str(muuc)+"_"+str(siggc)+"_"+str(replicaa)+"/html/" + str(folder2)
            try:
                os.makedirs(path)
            except OSError:
                foo=0
            else:
                foo=1

def savepara(okay_path,c_param,Stot,M,para_muc,para_sigc,com):
    com_in = com.copy()
    o = open("para/plateparams.txt",'w')
    print("\n".join("{}\t{}".format(k, v) for k, v in com_in.params[0].items()), file=o)
    o.close()
    plt.figure()
    fig,ax=plt.subplots()
    sns.heatmap(c_param,vmin=0,square=True,linewidths=.5,xticklabels=False,yticklabels=False,cbar=True,ax=ax)
    ax.set_title('consumer matrix'+ str(Stot)+"_"+str(M)+"_"+str(para_muc)+"_"+str(para_sigc)+"_"+"binary-gamma")
    plt.savefig(r"graphs/consumption_map"+ ".png")
    plt.close();

    outp = pd.DataFrame(com_in.N)
    outp.to_csv(r"dataN/N0.csv")
    outp.to_html(r"html/datahN/N0.html")

    outp = pd.DataFrame(com_in.R)
    outp.to_csv(r"dataR/R0.csv")
    outp.to_html(r"html/datahR/R0.html")

    file = open(str(okay_path)+"/psoc_setup.txt", "r")
    contents = file.read()
    dictionary0 = eval(contents)
    file.close()

    f = open(r"para/setup.csv", "w")
    f.write("{\n")
    for k in dictionary0.keys():
        f.write("'{}':'{}'\n".format(k, dictionary0[k]))
    f.write("}")
    f.close()
    return print('done')

def derange(s):
    import random
    d=s[:]
    while np.any([a==b for a,b in zip(d,s)]):random.shuffle(d)
    return d

def step2generator(list1,totspecies):
    list2 = derange(list1)
    rd2list = []
    for i in range(len(list1)):
        rd2list.append((list1[i],list2[i]))
    pairs = set(rd2list)
    finalset = set((a,b) if a<=b else (b,a) for a,b in pairs)
    lenn = len(finalset)
    while lenn !=totspecies:
        list2 = derange(list1)
        rd2list = []
        for i in range(len(list1)):
            rd2list.append((list1[i],list2[i]))
        pairs = set(rd2list)
        finalset = set((a,b) if a<=b else (b,a) for a,b in pairs)
        lenn = len(finalset)
    rd2list = list(finalset)
    return rd2list

from itertools import combinations
def rSubset(arr, r):
    return list(combinations(arr, r))

def step3generator(list1,totspecies):
    import random
    listf = rSubset(list1,3)
    sample = totspecies
    rdlist = []
    rdlist = random.sample(range(0, len(listf)), sample)
    rd3list = []
    for i in range(len(rdlist)):
        rd3list.append(listf[rdlist[i]])

    tval = []
    for i in list1:
        res = i in chain(*rd3list)
        if res == False :
            t=0
            tval.append(t)
        else :
            t=1
            tval.append(t)

    test = np.sum(tval)
    if test!=totspecies:
        #print(rd3list)
        #print(test)
        while test!= totspecies:
            rdlist = []
            rdlist = random.sample(range(0, len(listf)), sample)
            rd3list = []
            for i in range(len(rdlist)):
                rd3list.append(listf[rdlist[i]])
            print(rd3list)
            tval = []
            for i in list1:
                res = i in chain(*rd3list)
                if res == False :
                    t=0
                    print(i)
                    tval.append(t)
                else :
                    t=1
                    tval.append(t)
            test = np.sum(tval)
            print(test)
    return rd3list

def almost_matches(array, xx, rtol=1e-03, atol=1e-10):
    answer = []
    cols = array.shape[1]
    for y in range(cols):
        if abs(array[0][y]-xx) <= (atol + rtol * abs(xx)):
            answer.append(y)
    return answer

def almost_matches2(array, xx, rtol=1e-04, atol=0.0):
    answer = []
    cols = array.shape[1]
    for y in range(cols):
        if abs(array[0][y]-xx) <= (atol + rtol * abs(xx)):
            answer.append(y)
    return answer

def stepsearch(stepsize,inplate,n_species,iterii,noffwells,scale_sett,toll,typee,phiii,alphaaa,globey,list_species):
    commfunc = np.zeros((iterii,1))
    generation = np.zeros((iterii,1))
    failure = np.zeros((iterii,1))
    chosenwell = np.zeros((iterii,1))
    commfunc[0]=0
    neutraldirs = np.zeros((iterii,1))
    neutraldirs[0]=0
    community_in = inplate.copy()
    community_refresh = inplate.copy()
    propagation_time = [100,200]
    fac = 2
    start_time = time.time()
    o = open(r"results/validate.csv",'a')
    print((str(stepsize)+"-step validate simulation"), file=o)
    o.close()
    tol_c = 1e-7
    for k in range(1,iterii):
        #for i in range(1,noffwells):
        #    community_in.N['W'+str(i)]=community_in.N['W0']

        if stepsize == 1:
            print(stepsize)
            for i in range(0,noffwells-1):
                community_in.N.iat[i,i+1]= 1 - community_in.N.iat[i,0]

        if stepsize == 2:
            x = step2generator(list_species,n_species)
            for j in range(0,n_species):
                community_in.N.iat[x[j][0],j+1]= 1-community_in.N.iat[x[j][0],0]
                community_in.N.iat[x[j][1],j+1]= 1-community_in.N.iat[x[j][1],0]

        if stepsize == 3:
            y = step3generator(list_species,n_species)
            for j in range(0,n_species):
                community_in.N.iat[y[j][0],j+1]= 1-community_in.N.iat[y[j][0],0]
                community_in.N.iat[y[j][1],j+1]= 1-community_in.N.iat[y[j][1],0]
                community_in.N.iat[y[j][2],j+1]= 1-community_in.N.iat[y[j][2],0]

        N0 = community_in.N
        community_in.N.to_csv(r"dataN/N0_"+str(stepsize)+"-step_for_gen = "+ str(k)+".csv", sep='\t')
        community_in.N.to_html(r"html/datahN/N0_"+str(stepsize)+"_for_gen = "+ str(k) +".html")
        community_in.scale=scale_sett
        community_in.Propagate(propagation_time[0])
        community_in.N.to_csv(r"dataN/N_"+str(stepsize)+"-step_after_ "+str(propagation_time[0])+"units_for_gen = "+ str(k)+".csv", sep='\t')
        community_in.N.to_html(r"html/datahN/N_"+str(stepsize)+"-step_after_ "+str(propagation_time[0])+"units_for_gen = "+ str(k)+".html")
        community_in.R.to_csv(r"dataR/R_"+str(stepsize)+"-step_after_ "+str(propagation_time[0])+"units_for_gen = "+ str(k)+".csv", sep='\t')
        community_in.R.to_html(r"html/datahR/R_"+str(stepsize)+"-step_after_ "+str(propagation_time[0])+"units_for_gen = "+ str(k)+".html")
        accf = float(validate_simulation_sa743(community_in,N0,propagation_time[0],-1,stepsize)[0])
        check = float(validate_simulation_sa743(community_in,N0,propagation_time[0],-1,stepsize)[1])
        o = open(r"results/validate.csv",'a')
        print(propagation_time[0],accf,check, file=o)
        o.close()
        if check > 0:
            o = open(r"results/validate.csv",'a')
            print("invasion problem after t= "+str(propagation_time[0])+" time steps", file=o)
            o.close()
        #print('k is ',k)
        i=1
        sum = propagation_time[0]
        counttime = 0
        while accf > toll and counttime < 3600:
            start_timelocal = time.time()
            community_in.Propagate(propagation_time[i])
            print(propagation_time[i])
            sum = sum + propagation_time[i]
            print('sum is ',sum)
            community_in.N.to_csv(r"dataN/N_"+str(stepsize)+"-step_after_ "+str(sum)+" steps_for_gen = "+ str(k)+".csv", sep='\t')
            community_in.N.to_html(r"html/datahN/N_"+str(stepsize)+"-step_after_ "+str(sum)+" steps_for_gen = "+ str(k)+".html")
            community_in.R.to_csv(r"dataR/R_"+str(stepsize)+"-step_after_ "+str(sum)+" steps_for_gen = "+ str(k)+".csv", sep='\t')
            community_in.R.to_html(r"html/datahR/R_"+str(stepsize)+"-step_after_ "+str(sum)+" steps_for_gen = "+ str(k)+".html")
            accf = float(validate_simulation_sa743(community_in,N0,sum,k,stepsize)[0])
            print('acc is ',accf)
            check = float(validate_simulation_sa743(community_in,N0,sum,k,stepsize)[1])
            o = open(r"results/validate.csv",'a')
            print(sum,accf,check, file=o)
            o.close()
            if check > 0:
                o = open(r"results/validate.csv",'a')
                print("invasion problem after t= "+str(propagation_time[i])+" time steps", file=o)
                o.close()
            propagation_time.append(propagation_time[i]*fac)
            i=i+1
            time_ellocal = time.time() - start_timelocal
            counttime = counttime + time_ellocal

        community_in.N.to_csv(r"dataS/N_"+str(stepsize)+"-step_steady_for gen="+str(k)+".csv", sep='\t')
        community_in.N.to_html(r"html/datahS/N_"+str(stepsize)+"-step_steady_for_gen = "+ str(k)+".html")
        community_in.R.to_csv(r"dataS/R_"+str(stepsize)+"-step_steady_for_gen = "+ str(k)+".csv", sep='\t')
        community_in.R.to_html(r"html/datahS/R_"+str(stepsize)+"-step_steady_for_gen = "+ str(k)+".html")
        if typee == 1 :
            arr = phiii
            functionarray = arr.dot(community_in.N)
        elif typee == 2:
            arr = alphaaa
            functionarray = arr.dot(community_in.N)
        elif typee == 3:
            function = []
            for i in range(0,noffwells):
                                          function.append(0.031*community_in.N.iat[2,i]*community_in.N.iat[9,i]*community_in.N.iat[12,i]*community_in.N.iat[22,i]*community_in.N.iat[30,i]*community_in.N.iat[40,i])
            functionarray = np.asarray([function])
        outp = pd.DataFrame(functionarray.T)
        order = np.argsort(functionarray)
        print(functionarray)
        print(functionarray.shape)
        highest = order[0,noffwells-1]
        print(functionarray[0,highest])
        bigg = functionarray[0,highest]
        print(bigg)
        if abs(bigg) < 1:
            choices = almost_matches2(functionarray,bigg,rtol=1e-4,atol=0.0)
        else :
            choices = almost_matches(functionarray,bigg,rtol = 1e-3,atol = 1e-10)
        print((choices))
        neutraldirs[k] = len(choices)
        import random
        nwell = random.choice(choices)
        print(nwell)
        outp = pd.DataFrame(functionarray.T)
        outp.to_html(r"html/resh/functionarray_"+str(stepsize)+"step_allwells for gen = "+str(k)+".html")
        outp = pd.DataFrame(order.T)
        outp.to_html(r"html/resh/ordered_functionarray_"+str(stepsize)+"step_allwells for gen = "+str(k)+".html")
        comm = nwell
        chosenwell[k]=comm
        commfunc[k] = functionarray[0,comm]
        N0upd = N0['W'+str(nwell)]
        community_in = community_refresh
        for j in range(n_species):
            community_in.N['W0'] = N0upd

        generation[k]=k+1
        o = open(r'results\validate.csv','a')
        print("generation= "+str(k), file=o)
        o.close()
        diffcomm = commfunc[k]-commfunc[k-1]
        if diffcomm < 0 and abs(diffcomm) > .05*(commfunc[k-1]):
            o = open(r'results\validate.csv','a')
            print("for generation="+str(k)+"the community function decreased!", file=o)
            o.close()
        #if abs(diffcomm) < tol_c and k > 15:
        #    break
        #    print(diffcomm)
    print(propagation_time)
    outp = pd.DataFrame(np.asarray(np.concatenate((generation,commfunc,chosenwell,neutraldirs,failure), axis=1)))
    outp.to_csv(r"results/commfunc_run_step"+str(stepsize)+".csv", sep='\t')
    outp.to_html(r"html/resh/commfunc_run_step"+str(stepsize)+".html")
    time_el = time.time() - start_time
    with open(r"results/commfunc_run_step"+str(stepsize)+".csv", 'a') as file:
        file.write('time elapsed = ' + str(time_el) + 'seconds')

    from matplotlib.ticker import MaxNLocator
    plt.figure()
    plt.plot(generation,commfunc,'o',color='red')
    plt.xlabel('generation')
    plt.ylabel('community function')
    plt.grid(True, color = "black", linewidth = "1.2", linestyle = "--")
    plt.savefig(r"results/run_"+str(stepsize)+ "-step.png")
    return community_in, choices
