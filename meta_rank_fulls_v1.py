
'''
@author : sayra_ecoevol
'''
from scipy import *
import numpy as np
import time, math, json, pprint, requests, urllib, io, os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import random as rd
from itertools import *
colors = sns.color_palette()
from matplotlib.backends import backend_pdf as bpd
from community_simulator import *
from community_simulator import Community,essentialtools,usertools,visualization,analysis
from community_simulator.usertools import *
from community_simulator.visualization import *
from function_defs_full_search import *

assumpt = readassumptions(r"psoc_assump.txt")[0]
nofwells = assumpt['n_wells']
iteri = readassumptions(r"psoc_setup.txt")[1]
scale_set = readassumptions(r"psoc_setup.txt")[4]
assumpt['n_wells'] = 2**int(assumpt['SA'][0])
#assumpt['n_wells'] = 10
#print(int(assumpt['SA'][0]))
nofwells = assumpt['n_wells']

#set tolerance for the ODE integrator
tolg = 5e-2
#supply all resources
#make initial state
init_state_global = MakeInitialState(assumpt)
#make params
params = MakeParams(assumpt)
#make matrices c and D
c,D = MakeMatrices(assumpt)
#save initial state of the wells
outp = pd.DataFrame(init_state_global[0])
#output as html
outp.to_html(r"globalN0.html")
#extract Stot and M from the initial state
N0=init_state_global[0]
R0=init_state_global[1]
Stot=len(N0)
M = len(R0)

#get global parameters from github
user='sayra-ecoevol'
pao='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'
github_session = requests.Session()
github_session.auth = (user, pao)

#get phis
commfunc1 = getdata()[0]
#get D
Dsim = getdata()[1]
#get alphas
commfunc2 = getdata()[2]
#assign phis to a new array
phii = np.zeros((1,Stot))
for j in range(Stot):
    phii[0,j] = commfunc1.iat[j,0]

#assign alphas to a new array
alphaa = np.zeros((1,Stot))
for j in range(Stot):
    alphaa[0,j] = commfunc2.iat[j,0]



#this is needed for 2,3-step greedy walk, not for full search.
list1 = []
for i in range(Stot):
    list1.append(i)

x1 = step2generator(list1,Stot)
y1 = step3generator(list1,Stot)

new_path = os.getcwd()
print ("The current working directory is %s" % new_path)
#bunch of experiments begin
assumpt['R0_food'] = 100.0
for para in [[15,2.3]]:
    para_muc = para[0]
    para_sigc = para[1]
    for replica in range(10):
        init_state_same = MakeInitialState(assumpt)
        init_state_same[1].iloc[5] = 100.0
        init_state_same[1].iloc[0] = 0.0
        R0 = init_state_same[1]
        mkdirects(new_path,para_muc,para_sigc,replica)
        init_state = init_state_global
        def dNdt(N,R,params):
                return MakeConsumerDynamics(assumpt)(N,R,params)
        def dRdt(N,R,params):
            return MakeResourceDynamics(assumpt)(N,R,params)
        dynamics = [dNdt,dRdt]
        user='sayra-ecoevol'
        pao='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'
        github_session = requests.Session()
        github_session.auth = (user, pao)
        #import consumer preference matrix from github
        csv_url = r"https://raw.githubusercontent.com/sayra-ecoevol/CFO2020-/main/consumer-choice/consumption_matrix_globe_"+str(Stot)+"_"+str(M)+"_"+str(para_muc)+"_"+str(para_sigc)+"_replica_" + str(replica)+".csv"
        print(csv_url)
        download = github_session.get(csv_url).content
        c_param = pd.read_csv(io.StringIO(download.decode('utf-8')), error_bad_lines=False, header = None)
        print(c_param)
        params_thisexp=[{'c':c_param,
                    'm':assumpt['m'],
                    'w':1.0,
                    'D':Dsim,
                    'g':1.0,
                    'l':0.6,
                    'R0':R0.values[:,k],
                    'tau':assumpt['tau']
                    } for k in range(nofwells)]
        plate1 = Community(init_state_same,dynamics,params_thisexp,parallel=False)
        change_path = str(new_path)+"/"+str(para_muc)+"_"+str(para_sigc)+"_"+str(replica)
        os.chdir(change_path)
        savepara(new_path,c,Stot,M,para_muc,para_sigc,plate1)
        n_species = Stot
        tot = 2**n_species
        form = '0'+str(n_species)+'b'
        for j in range(tot):
            y=format(j, form)
            for i in range(n_species):
                plate1.N.iat[i,j]=y[i]
        print(plate1.N)

        iterii=2
        typef = 1
        stepsearch(0,plate1,Stot,iterii,nofwells,scale_set,tolg,typef,phii,alphaa,y1,list1)
