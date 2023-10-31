# Load profiles and generation mix data (LPGM) & energy generation, storage and transmission information (GGTA)
# based on x/capacities from Optimisation and flexible from Dispatch
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

from Input import *
from Simulation import Reliability
from Network import Transmission

import numpy as np
import datetime as dt

def Debug(solution):
    """Debugging"""

    Load, PV, India = (solution.MLoad.sum(axis=1), solution.GPV.sum(axis=1), solution.MIndia.sum(axis=1))
#    Wind = solution.GWind.sum(axis=1)
    Baseload = solution.MBaseload.sum(axis=1)
    Pond = solution.MPond.sum(axis=1)

    DischargePH, ChargePH, StoragePH = (solution.DischargePH, solution.ChargePH, solution.StoragePH)
    Deficit_energy, Deficit_power, Deficit, Spillage = (solution.Deficit_energy, solution.Deficit_power, solution.Deficit, solution.Spillage)

    PHS = solution.CPHS * pow(10, 3) # GWh to MWh
    efficiencyPH = solution.efficiencyPH

    for i in range(intervals):
        # Energy supply-demand balance
#        assert abs(Load[i] + ChargePH[i] + ChargeB[i] + Spillage[i]
#                   - PV[i] - Inter[i] - Wind[i] - Baseload[i] - Peak[i] - DischargePH[i] + DischargeB[i] - Deficit[i]) <= 1
        assert abs(Load[i] + ChargePH[i] + Spillage[i] \
                   - PV[i] - India[i] - Baseload[i] - Pond[i] - DischargePH[i] - Deficit[i]) <= 1

        # Discharge, Charge and Storage
        if i==0:
            assert abs(StoragePH[i] - 0.5 * PHS + DischargePH[i] * resolution - ChargePH[i] * resolution * efficiencyPH) <= 1
        else:
            assert abs(StoragePH[i] - StoragePH[i - 1] + DischargePH[i] * resolution - ChargePH[i] * resolution * efficiencyPH) <= 1

        # Capacity: PV, wind, Discharge, Charge and Storage
        try:
            assert np.amax(PV) - sum(solution.CPV) * pow(10, 3) <= 0.1, print("PV",np.amax(PV) - sum(solution.CPV) * pow(10, 3))
#            assert np.amax(Wind) <= sum(solution.CWind) * pow(10, 3), print(np.amax(Wind) - sum(solution.CWind) * pow(10, 3))
            assert np.amax(India) - sum(solution.CInter) * pow(10,3) <= 0.1

            assert np.amax(DischargePH) - sum(solution.CPHP) * pow(10, 3) <= 0.1, print("DischargePH",np.amax(DischargePH) - sum(solution.CPHP) * pow(10, 3))
            assert np.amax(ChargePH) - sum(solution.CPHP) * pow(10, 3) <= 0.1, print("ChargePH",np.amax(ChargePH) - sum(solution.CPHP) * pow(10, 3))
            assert np.amax(StoragePH) - solution.CPHS * pow(10, 3) <= 0.1, print("StoragePH",np.amax(StoragePH) - solution.CPHS * pow(10, 3))
        except AssertionError:
            pass

    print('Debugging: everything is ok')

    return True

def LPGM(solution):
    """Load profiles and generation mix data"""

    Debug(solution)

    C = np.stack([(solution.MLoad).sum(axis=1),
                  solution.MBaseload.sum(axis=1), solution.MPond.sum(axis=1), solution.MIndia.sum(axis=1), solution.GPV.sum(axis=1), #solution.GWind.sum(axis=1),
                  solution.DischargePH, solution.Deficit, -1 * (solution.Spillage + indiaExportProfiles), -1 * solution.ChargePH,
                  solution.StoragePH,
                  solution.CHTH, solution.THTS, solution.TSSA, solution.SAZH, solution.ZHPE, solution.PEMO, solution.IN1CH, solution.IN2TS, solution.IN3SA, solution.IN4PE])

    C = np.around(C.transpose())

    datentime = np.array([(dt.datetime(firstyear, 1, 1, 0, 0) + x * dt.timedelta(minutes=60 * resolution)).strftime('%a %-d %b %Y %H:%M') for x in range(intervals)])
    C = np.insert(C.astype('str'), 0, datentime, axis=1)

    header = 'Date & time,Operational demand,' \
             'RoR Hydropower (MW),Pond Hydropower (MW),India Imports (MW),Solar photovoltaics (MW),PHES-Discharge (MW),Energy deficit (MW),India Exports (MW),PHES-Charge (MW),' \
             'PHES-Storage (MWh),' \
             'CHTH,THTS,TSSA,SAZH,ZHPE,PEMO,IN1CH,IN2TS,IN3SA,IN4PE'

    np.savetxt('Results/LPGM_{}_{}_{}_{}_Network.csv'.format(node,scenario,percapita,export_flag), C, fmt='%s', delimiter=',', header=header, comments='')

    if 'Super' in node:
        header = 'Date & time,Operational demand,' \
                 'RoR Hydropower (MW),Pond Hydropower (MW),India Imports (MW),Solar photovoltaics (MW),PHES-Discharge (MW),Energy deficit (MW),India Exports (MW),'\
                 'Transmission,PHES-Charge (MW),' \
                 'PHES-Storage'

        Topology = solution.Topology[np.where(np.in1d(Nodel, coverage) == True)[0]]

        for j in range(nodes):
            MCH2_exports = np.zeros(intervals) if j != 1 else indiaExportProfiles

            C = np.stack([(solution.MLoad)[:, j],
                          solution.MBaseload[:, j], solution.MPond[:, j],solution.MIndia[:, j], solution.MPV[:, j], #solution.MWind[:, j],
                          solution.MDischargePH[:, j], solution.MDeficit[:, j], -1 * (solution.MSpillage[:, j] + MCH2_exports), Topology[j], 
                          -1 * solution.MChargePH[:, j],
                          solution.MStoragePH[:, j]])
            C = np.around(C.transpose())

            C = np.insert(C.astype('str'), 0, datentime, axis=1)
            np.savetxt('Results/LPGM_{}_{}_{}_{}_{}.csv'.format(node,scenario,percapita, export_flag,solution.Nodel[j]), C, fmt='%s', delimiter=',', header=header, comments='')

    print('Load profiles and generation mix is produced.')

    return True

def GGTA(solution):
    """GW, GWh, TWh p.a. and A$/MWh information"""
    # Import cost factors
    factor = np.genfromtxt('Data/factor.csv', dtype=None, delimiter=',', encoding=None)
        
    factor = dict(factor)

    # Import capacities [GW, GWh] from the least-cost solution
    CPV, CInter, CPHP, CPHS = (sum(solution.CPV), sum(solution.CInter), sum(solution.CPHP), solution.CPHS) # GW, GWh
#    CWind = sum(solution.CWind)
    CapHydro = CHydro_max.sum() # GW

    # Import generation energy [GWh] from the least-cost solution
    GPV, GHydro, GIndia = map(lambda x: x * pow(10, -6) * resolution / years, (solution.GPV.sum(), solution.MBaseload.sum() + solution.MPond.sum(), solution.MIndia.sum())) # TWh p.a.
    DischargePH = solution.DischargePH.sum()
    Ghydro_CH2 = indiaExportProfiles.sum() * resolution / years
#    GWind = solution.GWind.sum()
    CFPV = GPV / CPV / 8.76 if CPV != 0 else 0
#    CFWind = GWind / CWind / 8.76
    
    # Calculate the annual costs for each technology
    CostPV = factor['PV'] * CPV # A$b p.a.
#    CostWind = factor['Wind'] * CWind # A$b p.a.
    CostHydro = factor['Hydro'] * (GHydro + Ghydro_CH2) # A$b p.a.
    CostPH = factor['PHP'] * CPHP + factor['PHS'] * CPHS + factor['PHES-VOM'] * DischargePH * pow(10, -6) * resolution / years # A$b p.a.
    CostIndia = factor['India'] * GIndia # A$b p.a.
#    if scenario>=21:
#        CostPH -= factor['LegPH']

    CostT = np.array([factor['CHTH'], factor['THTS'], factor['TSSA'], factor['SAZH'], factor['ZHPE'], factor['PEMO'], factor['IN1CH'], factor['IN2TS'], factor['IN3SA'], factor['IN4PE']])
    CostDC, CostAC, CDC, CAC = [],[],[],[]

    for i in range(0,len(CostT)):
        CostDC.append(CostT[i]) if dc_flags[i] else CostAC.append(CostT[i])
        CDC.append(solution.CDC[i]) if dc_flags[i] else CAC.append(solution.CDC[i])
    CostDC, CostAC, CDC, CAC = [np.array(x) for x in [CostDC, CostAC, CDC, CAC]]
    
    CostDC = (CostDC * CDC).sum() if len(CDC) > 0 else 0 # A$b p.a.
    CostAC = (CostAC * CAC).sum() if len(CAC) > 0 else 0 # A$b p.a.

#    if scenario>=21:
#        CostDC -= factor['LegINTC']

    CostAC += factor['ACPV'] * CPV # + factor['ACWind'] * CWind # A$b p.a.
    
    # Calculate the average annual energy demand
    Energy = (MLoad).sum() * pow(10, -9) * resolution / years # PWh p.a.
    Exports = indiaExportProfiles.sum() + solution.MSpillage.sum()
    Loss = np.sum(abs(solution.TDC), axis=0) * TLoss
    Loss = Loss.sum() * pow(10, -9) * resolution / years # PWh p.a.

    # Calculate the levelised cost of elcetricity at a network level
    LCOE = (CostPV + CostIndia + CostHydro + CostPH + CostDC + CostAC) / (Energy - Loss) # + CostWind / (Energy - Loss)
    LCOEPV = CostPV / (Energy - Loss)
#    LCOEWind = CostWind / (Energy - Loss)
    LCOEIndia = CostIndia / (Energy - Loss)
    LCOEHydro = CostHydro / (Energy - Loss)
    LCOEPH = CostPH / (Energy - Loss)
    LCOEDC = CostDC / (Energy - Loss)
    LCOEAC = CostAC / (Energy - Loss)
    
    # Calculate the levelised cost of generation
#    LCOG = (CostPV + CostWind + CostHydro + CostBio) * pow(10, 3) / (GPV + GWind + GHydro + GBio)
    LCOG = (CostPV + CostHydro + CostIndia) * pow(10, 3) / (GPV + GHydro + GIndia)
    LCOGP = CostPV * pow(10, 3) / GPV if GPV!=0 else 0
#    LCOGW = CostWind * pow(10, 3) / GWind if GWind!=0 else 0
    LCOGH = CostHydro * pow(10, 3) / GHydro if GHydro!=0 else 0
    LCOGI = CostIndia * pow(10, 3) / GIndia if GIndia != 0 else 0

    # Calculate the levelised cost of balancing
    LCOB = LCOE - LCOG
    LCOBS_P = CostPH / (Energy - Loss)
    LCOBT = (CostDC + CostAC) / (Energy - Loss)
    LCOBL = LCOB - LCOBS_P - LCOBT

    print('Levelised costs of electricity:')
    print('\u2022 LCOE:', LCOE)
    print('\u2022 LCOG:', LCOG)
    print('\u2022 LCOB:', LCOB)
    print('\u2022 LCOG-PV:', LCOGP, '(%s)' % CFPV)
#    print('\u2022 LCOG-Wind:', LCOGW, '(%s)' % CFWind)
    print('\u2022 LCOG-Hydro:', LCOGH)
    print('\u2022 LCOG-External_Imports:', LCOGI)
    print('\u2022 LCOB-PHES_Storage:', LCOBS_P)
    print('\u2022 LCOB-Transmission:', LCOBT)
    print('\u2022 LCOB-Spillage & loss:', LCOBL)

    size = 20 + len(list(solution.CDC))
    D = np.zeros((1, size))
    header = 'Annual demand (PWh),Annual Energy Losses (PWh),' \
             'PV Capacity (GW),PV Avg Annual Gen (GWh),Hydro Capacity (GW),Hydro Avg Annual Gen (GWh),Inter Capacity (GW),India Avg Annual Imports (GWh),India Avg Annual Exports (GWh),' \
             'PHES-PowerCap (GW),PHES-EnergyCap (GWh),' \
             'CHTH,THTS,TSSA,SAZH,ZHPE,PEMO,IN1CH,IN2TS,IN3SA,IN4PE,' \
             'LCOE,LCOG,LCOB,LCOG_PV, LCOG_Hydro,LCOG_IndiaImports,LCOBS_PHES,LCOBT,LCOB_LossesExports'
    D[0, :] = [Energy * pow(10, 3), Loss * pow(10, 3), CPV, GPV, CapHydro, GHydro, CInter, GIndia,Exports] \
              + [CPHP, CPHS] \
              + list(solution.CDC) \
              + [LCOE, LCOG, LCOB, LCOGP, LCOGH, LCOGI, LCOBS_P, LCOBT, LCOBL] # + [CWind, GWind, LCOGW]

    np.savetxt('Results/GGTA_{}_{}_{}_{}.csv'.format(node,scenario,percapita,export_flag), D, fmt='%f', delimiter=',',header=header)
    print('Energy generation, storage and transmission information is produced.')

    return True

def Information(x, flexible):
    """Dispatch: Statistics.Information(x, Flex)"""

    start = dt.datetime.now()
    print("Statistics start at", start)

    S = Solution(x)
    Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=flexible, daily_pondage=daily_pondage)
    
    try:
        assert Deficit.sum() * resolution < 0.1, 'Energy generation and demand are not balanced.'
    except AssertionError:
        pass
    
    assert np.reshape(baseload.sum(axis=1) + S.DischargePond, (-1, 8760)).sum(axis=-1).max() <= Hydromax, "Hydro generation exceeds requirement"

    #S.TDC = Transmission(S, output=True) if 'APG' in node else np.zeros((intervals, len(TLoss))) # TDC(t, k), MW
    S.TDC = Transmission(S, output=True)
    S.CDC = np.amax(abs(S.TDC), axis=0) * pow(10, -3) # CDC(k), MW to GW
    
    S.CHTH, S.THTS, S.TSSA, S.SAZH, S.ZHPE, S.PEMO, S.IN1CH, S.IN2TS, S.IN3SA, S.IN4PE = map(lambda k: S.TDC[:, k], range(S.TDC.shape[1]))

    if 'Super' not in node:
        S.MPV = S.GPV
    #    S.MWind = S.GWind if S.GWind.shape[1]>0 else np.zeros((intervals, 1))
        S.MIndia = S.GIndia
        S.MDischargePH = np.tile(S.DischargePH, (nodes, 1)).transpose()
        S.MDeficit = np.tile(S.Deficit, (nodes, 1)).transpose()
        S.MChargePH = np.tile(S.ChargePH, (nodes, 1)).transpose()
        S.MStoragePH = np.tile(S.StoragePH, (nodes, 1)).transpose()
        S.MSpillage = np.tile(S.Spillage, (nodes, 1)).transpose()

    #S.MBaseload = np.clip(S.MHydro, None, CHydro * pow(10, 3)) # GHydro(t, j), GW to MW

    S.MPHS = S.CPHS * np.array(S.CPHP) * pow(10, 3) / sum(S.CPHP) # GW to MW

    # S.CHTH, S.THTS, S.TSSA, S.SAZH, S.ZHPE, S.PEMO, S.IN1CH, S.IN2TS, S.IN3SA, S.IN4PE
    S.Topology = np.array([S.IN1CH + S.CHTH,                    # CH
                  -1 * (S.CHTH + S.THTS),                       # TH
                  S.THTS + S.TSSA + S.IN2TS,                    # TS
                  -1 * (S.TSSA + S.IN3SA + S.SAZH),             # SA
                  S.SAZH + S.ZHPE,                              # ZH
                  -1 * (S.ZHPE + S.IN4PE + S.PEMO),             # PE
                  S.PEMO,                                       # MO
                  -1 * S.IN1CH,                                 # IN1
                  -1 * S.IN2TS,                                 # IN2
                  S.IN3SA,                                      # IN3
                  S.IN4PE])                                     # IN4

    LPGM(S)
    GGTA(S)

    end = dt.datetime.now()
    print("Statistics took", end - start)

    return True

if __name__ == '__main__':
    suffix="_Super_existing_20_True.csv"
    Optimisation_x = np.genfromtxt('Results/Optimisation_resultx{}'.format(suffix), delimiter=',')
    flexible = np.genfromtxt('Results/Dispatch_Flexible{}'.format(suffix), delimiter=',', skip_header=1)
    Information(Optimisation_x, flexible)