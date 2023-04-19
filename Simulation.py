# To simulate energy supply-demand balance based on long-term, high-resolution chronological data
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

import numpy as np

def Reliability(solution, hydro, start=None, end=None):
    """Deficit = Simulation.Reliability(S, hydro=...)"""

    ###### CALCULATE NETLOAD FOR EACH INTERVAL ######
    Netload = (solution.MLoad.sum(axis=1) - solution.GPV.sum(axis=1) - solution.GInter.sum(axis=1) - solution.exports)[start:end] \
        - hydro[start:end] # - solution.GWind.sum(axis=1); Sj-ENLoad(j, t), MW
    length = len(Netload)
    
    solution.hydro = hydro # MW

    ###### CREATE STORAGE SYSTEM VARIABLES ######
    Pcapacity_PH = sum(solution.CPHP) * pow(10, 3) # S-CPHP(j), GW to MW
    Scapacity_PH = solution.CPHS * pow(10, 3) # S-CPHS(j), GWh to MWh
    efficiencyPH, resolution = (solution.efficiencyPH, solution.resolution)

    DischargePH, ChargePH, StoragePH = map(np.zeros, [length] * 3)
    Deficit_energy, Deficit_power = map(np.zeros, [length] * 2)

    for t in range(length):
        ###### INITIALISE INTERVAL ######
        Netloadt = Netload[t]
        Storage_PH_t1 = StoragePH[t-1] if t>0 else 0.5 * Scapacity_PH

        ##### UPDATE STORAGE SYSTEMS ######
        Discharge_PH_t = min(max(0, Netloadt), Pcapacity_PH, Storage_PH_t1 / resolution)
        Charge_PH_t = min(-1 * min(0, Netloadt), Pcapacity_PH, (Scapacity_PH - Storage_PH_t1) / efficiencyPH / resolution)
        Storage_PH_t = Storage_PH_t1 - Discharge_PH_t * resolution + Charge_PH_t * resolution * efficiencyPH

        DischargePH[t] = Discharge_PH_t
        ChargePH[t] = Charge_PH_t
        StoragePH[t] = Storage_PH_t

        diff1 = Netloadt - Discharge_PH_t + Charge_PH_t
        
        ###### DETERMINE DEFICITS ######
        if diff1 <= 0:
            Deficit_energy[t] = 0
            Deficit_power[t] = 0
        elif (Discharge_PH_t == Pcapacity_PH):
            Deficit_energy[t] = 0
            Deficit_power[t] = diff1
        elif (Discharge_PH_t == Storage_PH_t1 / resolution):
            Deficit_energy[t] = diff1
            Deficit_power[t] = 0    

    Deficit = Deficit_energy + Deficit_power
    Spillage = -1 * np.minimum(Netload + ChargePH - DischargePH, 0)

    ###### ERROR CHECKING ######
    assert 0 <= int(np.amax(StoragePH)) <= Scapacity_PH, 'Storage below zero or exceeds max storage capacity'
    assert np.amin(Deficit) > -0.1, 'DeficitD below zero'
    assert np.amin(Spillage) >= 0, 'Spillage below zero'

    ###### UPDATE SOLUTION OBJECT ######
    solution.DischargePH, solution.ChargePH, solution.StoragePH = (DischargePH, ChargePH, StoragePH)
    solution.Deficit_energy, solution.Deficit_power, solution.Deficit, solution.Spillage = (Deficit_energy, Deficit_power, Deficit, Spillage)

    return Deficit_energy, Deficit_power, Deficit, DischargePH