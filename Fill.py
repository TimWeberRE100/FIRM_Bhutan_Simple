# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:17:37 2020
@author: cheng + tim
"""

from Input import *
from Simulation import Reliability
import numpy as np
import datetime as dt

def fill_deficit(deficit,india_imports,india_limit,india_annual,impflag,eff,step):
    idx = np.where(deficit > 0)[0]
    for idd, i in np.ndenumerate(idx):
        d = deficit[i]
        t = i
        count = 0
        #print("--------------------")
        #print(idd, " of ", len(idx))
        try:
            while d > 0 and t >= 0 and count < step:
                #print("t = ",t)
                year = t // 8760
                start = year * 8760
                end = (year+1) * 8760
                if t == i - 1:
                    d = d / eff
                if impflag:
                    remaining = india_annual - sum(india_imports[start:end])
                    assert remaining >= 0
                    hydro_c = min(india_imports[t] + d, india_limit, india_imports[t] + remaining)
                    d = d - (hydro_c - india_imports[t])
                    india_imports[t] = hydro_c
                    if remaining == 0:
                        print("Year", year, " annual limit met")
                        t = start - 1
                    else:
                        idxx = np.where(india_imports < india_limit)[0]
                        t = sorted([i for i in idxx if i < t])[-1]
                count += 1
        except:
            continue
    return india_imports

def save(imp,suffix):
    np.savetxt('Results/Dispatch_IndiaImports' + suffix, imp, fmt='%f', delimiter=',', newline='\n', header='India Imports')
    
def maxx(x):
    return np.reshape(x, (-1, 8760)).sum(axis=-1).max()/1e6

def mean(x):
    return x.sum()/years/1e6

def Analysis(optimisation_x,suffix):
    starttime = dt.datetime.now()
    print('Deficit fill starts at', starttime)

    S = Solution(optimisation_x)
    
    Deficit_energy1, Deficit_power1, Deficit1, DischargePH1, DischargePond1, Spillage1 = Reliability(S, baseload=baseload, india_imports=np.zeros(intervals), daily_pondage=daily_pondage, pond_hours=pond_hours)
    Max_deficit1 = np.reshape(Deficit1, (-1, 8760)).sum(axis=-1) # MWh per year
    PIndia = Deficit1.max() * pow(10, -3) # GW

    GIndia = resolution * (Max_deficit1).max() / efficiencyPH    

    print("GIndia_max:", GIndia/1e6)
    
    GIndia2 = Deficit1.sum() / years / 0.8
    
    print("GIndia_mean:", GIndia2/1e6)

    Indiamax = energy*pow(10,9)
    
    print("IMPORTS")
    print("------------------------------")
    india_imports = np.zeros(intervals)
    Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=india_imports, daily_pondage=daily_pondage, pond_hours=pond_hours)
    imp = fill_deficit(Deficit,india_imports,sum(S.CInter)*1e3,Indiamax,True,0.8,168)
    Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=imp, daily_pondage=daily_pondage, pond_hours=pond_hours)
    print("India generation:", maxx(imp))
    print("Remaining deficit:", Deficit.sum()/1e6)
    step = 1
    while Deficit.sum() > allowance*years and step < 50:
        imp = fill_deficit(Deficit,imp,sum(S.CInter)*1e3,Indiamax,True,0.8,168)
        Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=imp, daily_pondage=daily_pondage, pond_hours=pond_hours)
        step += 1
    print("India generation max:", maxx(imp))
    print("India generation mean:", mean(imp))
    print("Remaining deficit final:", Deficit.sum()/1e6)

    save(imp,suffix)

    endtime = dt.datetime.now()
    print('Deficit fill took', endtime - starttime)

    from Statistics import Information
    Information(optimisation_x,imp)

    return True

if __name__=='__main__':
    suffix = "_Super_existing_20_True.csv"
    optimisation_x = np.genfromtxt('Results/Optimisation_resultx{}'.format(suffix), delimiter=',')
    Analysis(optimisation_x,suffix)