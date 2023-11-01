# To simulate energy supply-demand balance based on long-term, high-resolution chronological data
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

import numpy as np

def Reliability(solution, baseload, india_imports, daily_pondage, start=None, end=None):
    """Deficit = Simulation.Reliability(S, hydro=...)"""

    ###### CALCULATE NETLOAD FOR EACH INTERVAL ######
    Netload = (solution.MLoad.sum(axis=1) - solution.GPV.sum(axis=1) - baseload.sum(axis=1))[start:end] \
                - india_imports # - solution.GWind.sum(axis=1); Sj-ENLoad(j, t), MW
    length = len(Netload)
    
    solution.india_imports = india_imports # MW

    ###### CREATE STORAGE SYSTEM VARIABLES ######
    Pcapacity_PH = sum(solution.CPHP) * pow(10, 3) # S-CPHP(j), GW to MW
    Scapacity_PH = solution.CPHS * pow(10, 3) # S-CPHS(j), GWh to MWh
    Pcapacity_Pond = sum(solution.CHydro_Pond) * pow(10, 3)
    #Scapacity_Pond = [4*x for x in solution.CHydro_Pond]
    efficiencyPH, resolution = (solution.efficiencyPH, solution.resolution)

    DischargePH, ChargePH, StoragePH, DischargePond, StoragePond = map(np.zeros, [length] * 5)
    Deficit_energy, Deficit_power = map(np.zeros, [length] * 2)
    
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
            Storage_Pond_t1 = StoragePond[t-1] + daily_pondage_divided[t]

        # Calculate pond discharge
        Discharge_Pond_t = np.minimum(np.maximum(0, Netloadt), Pcapacity_Pond)
        Discharge_Pond_t = np.minimum(Discharge_Pond_t, Storage_Pond_t1/resolution)
        Storage_Pond_t = Storage_Pond_t1 - Discharge_Pond_t * resolution

        DischargePond[t] = Discharge_Pond_t
        StoragePond[t] = Storage_Pond_t

        if (t - 1 % 24 == 0):
            pond_daily_exports = StoragePond[t]
            counter = 0
            while ((pond_daily_exports > 0.1) & (counter < 24)):
                discharge_pond_original = DischargePond[t - counter]
                DischargePond[t - counter] = np.minimum(Pcapacity_Pond, discharge_pond_original + StoragePond[t]/resolution)
                exported_pond_t = DischargePond[t - counter] - discharge_pond_original
                StoragePond[t - counter] = StoragePond[t - counter] - exported_pond_t * resolution

                pond_daily_exports -= exported_pond_t
                counter+=1

        """ if t==0:
            print(Netloadt, Pcapacity_Pond, Discharge_Pond_t) 
 """
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
    Netload = Netload - DischargePond
    Spillage = -1 * np.minimum(Netload + ChargePH - DischargePH, 0)

    #print(DischargePond)

    ###### ERROR CHECKING ######
    assert 0 <= int(np.amax(StoragePH)) <= Scapacity_PH, 'Storage below zero or exceeds max storage capacity'
    assert np.amin(Deficit) > -0.1, 'DeficitD below zero'
    assert np.amin(Spillage) >= 0, 'Spillage below zero'

    ###### UPDATE SOLUTION OBJECT ######
    solution.DischargePH, solution.ChargePH, solution.StoragePH, solution.DischargePond, solution.StoragePond = (DischargePH, ChargePH, StoragePH, DischargePond, StoragePond)
    solution.Deficit_energy, solution.Deficit_power, solution.Deficit, solution.Spillage = (Deficit_energy, Deficit_power, Deficit, Spillage)

    return Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage

if __name__ == '__main__':
    from Input import *
    from Network import Transmission 

    suffix = "_Super_existing_20_True.csv"
    Optimisation_x = np.genfromtxt('Results/Optimisation_resultx{}'.format(suffix), delimiter=',')
    
    # Initialise the optimisation
    S = Solution(Optimisation_x)

    CIndia = np.nan_to_num(np.array(S.CInter))

    # Simulation with only baseload
    Deficit_energy1, Deficit_power1, Deficit1, DischargePH1, DischargePond1, Spillage1 = Reliability(S, baseload=baseload, india_imports=np.zeros(intervals), daily_pondage=daily_pondage)
    Max_deficit1 = np.reshape(Deficit1, (-1, 8760)).sum(axis=-1) # MWh per year
    PIndia = Deficit1.max() * pow(10, -3) # GW

    GIndia = resolution * (Max_deficit1).max() / efficiencyPH

    PenEnergy = 0
    PenPower = 0
    
    # Simulation with baseload, all existing capacity, and all hydrogen
    Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=np.ones(intervals) * CIndia.sum() * pow(10,3), daily_pondage=daily_pondage)

    # Deficit penalty function
    PenDeficit = max(0, Deficit.sum() * resolution - S.allowance)

    # India import profile
    india_imports = np.clip(Deficit1, 0, CIndia.sum() * pow(10,3))

    # Simulation using the existing capacity generation profiles - required for storage average annual discharge
    Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=india_imports, daily_pondage=daily_pondage)

    # Discharged energy from storage systems
    GPHES = DischargePH.sum() * resolution / years * pow(10,-6) # TWh per year

    # Transmission capacity calculations
    TDC = Transmission(S) if 'Super' in node else np.zeros((intervals, len(TLoss))) # TDC: TDC(t, k), MW
    CDC = np.amax(abs(TDC), axis=0) * pow(10, -3) # CDC(k), MW to GWeakeak

    # Transmission penalty function
    PenDC = 0

    # Average annual electricity generated by existing capacity
    GHydro = resolution * (baseload.sum() + DischargePond.sum()) / efficiencyPH / years # MWh
    
    # Average annual electricity imported through external interconnections
    GIndia = resolution * india_imports.sum() / years / efficiencyPH # MWh

    # Assume all spillage, curtailment, and Hydro_CH2 generation is exported to india
    export_annual = (Spillage.sum() + indiaExportProfiles.sum()) * resolution / years
    Ghydro_CH2 = indiaExportProfiles.sum() * resolution / years

    # Levelised cost of electricity calculation
    cost = factor * np.array([sum(S.CPV), 0, GIndia * pow(10,-6), sum(S.CPHP), S.CPHS, GPHES] + list(CDC) + [sum(S.CPV), 0, (GHydro + Ghydro_CH2) * pow(10, -6), 0, 0]) # $b p.a.
    cost = cost.sum()
    loss = np.sum(abs(TDC), axis=0) * TLoss
    loss = loss.sum() * pow(10, -9) * resolution / years # PWh p.a.
    LCOE = cost / abs(energy - loss) 

    print(LCOE, energy, loss, GIndia * pow(10,-6) * factor[2] / (GIndia*pow(10,-9)), (GHydro + Ghydro_CH2) * pow(10, -6) * factor[18] / ((GHydro + Ghydro_CH2)*pow(10,-9)), len(CDC))