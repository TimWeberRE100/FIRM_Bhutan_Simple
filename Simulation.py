# To simulate energy supply-demand balance based on long-term, high-resolution chronological data
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

import numpy as np

def Reliability(solution, baseload, india_imports, daily_pondage, start=None, end=None):
    """Deficit = Simulation.Reliability(S, hydro=...)"""

    ###### CALCULATE NETLOAD FOR EACH INTERVAL ######
    Netload = (solution.MLoad.sum(axis=1) - solution.GPV.sum(axis=1) - baseload.sum(axis=1))[start:end] \
                - india_imports[start:end] # - solution.GWind.sum(axis=1); Sj-ENLoad(j, t), MW
    length = len(Netload)
    
    solution.india_imports = india_imports # MW

    ###### CREATE STORAGE SYSTEM VARIABLES ######
    Pcapacity_PH = sum(solution.CPHP) * pow(10, 3) # S-CPHP(j), GW to MW
    Scapacity_PH = solution.CPHS * pow(10, 3) # S-CPHS(j), GWh to MWh
    Pcapacity_Pond = sum(solution.CHydro_Pond)
    #Scapacity_Pond = [4*x for x in solution.CHydro_Pond]
    efficiencyPH, resolution = (solution.efficiencyPH, solution.resolution)

    DischargePH, ChargePH, StoragePH, DischargePond, StoragePond = map(np.zeros, [length] * 5)
    Deficit_energy, Deficit_power = map(np.zeros, [length] * 2)
    
    Storage_Pond_t1 = [0] * len(solution.CHydro_Pond)

    daily_pondage_divided = daily_pondage.sum(axis=1) / 24

    for t in range(length):
        """ if t%1000 == 0:
            print(t) """
        ###### INITIALISE INTERVAL ######
        Netloadt = Netload[t]
        Storage_PH_t1 = StoragePH[t-1] if t>0 else 0.5 * Scapacity_PH

        if t % 24 == 0:
            Storage_Pond_t1 = daily_pondage_divided[t]
        else:
            Storage_Pond_t1 += daily_pondage_divided[t]

        # Calculate pond discharge
        Discharge_Pond_t = np.minimum(np.maximum(0, Netloadt), Pcapacity_Pond)
        Discharge_Pond_t = np.minimum(Discharge_Pond_t, Storage_Pond_t1/resolution)
        Storage_Pond_t = Storage_Pond_t1 - Discharge_Pond_t * resolution

        DischargePond[t] = Discharge_Pond_t
        StoragePond[t] = Storage_Pond_t

        # Update storage systems
        Netloadt = Netloadt - Discharge_Pond_t
        
        """ for i in range(0,len(solution.CHydro_Pond)):
            if t % 24 == 0:
                Storage_Pond_t1[i] = daily_pondage[t,i] / 16
            else:
                Storage_Pond_t1[i] += daily_pondage[t,i] / 16
            
            ###### CALCULATE POND DISCHARGE #####
            Discharge_Pond_t_i = min(max(0,Netloadt), Pcapacity_Pond[i], Storage_Pond_t1[i] / resolution)
            Storage_Pond_t_i = Storage_Pond_t1[i] - Discharge_Pond_t_i*resolution

            DischargePond[t,i] = Discharge_Pond_t_i
            StoragePond[t,i] = Storage_Pond_t_i

        
        #print(t, DischargePond.sum(axis=1))
        #print(t, StoragePond)

        ##### UPDATE STORAGE SYSTEMS ######
        Netloadt = Netloadt - DischargePond.sum(axis=1)[t] #### CHECK THIS IS THE CORRECT AXIS """
        #print(t, Netload)
        #print(t, Netloadt)
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
    solution.DischargePH, solution.ChargePH, solution.StoragePH, solution.DischargePond, solution.StoragePond = (DischargePH, ChargePH, StoragePH, DischargePond, StoragePond)
    solution.Deficit_energy, solution.Deficit_power, solution.Deficit, solution.Spillage = (Deficit_energy, Deficit_power, Deficit, Spillage)

    return Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage