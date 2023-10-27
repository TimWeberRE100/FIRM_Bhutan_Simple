# A transmission network model to calculate inter-regional power flows
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

import numpy as np

def Transmission(solution, output=False):
    """TDC = Network.Transmission(S)"""

    Nodel, PVl, Interl, Hydrol = (solution.Nodel, solution.PVl, solution.Interl, solution.Hydrol)
#    Windl = solution.Windl
    intervals, nodes = (solution.intervals, solution.nodes)
    
    CHydro_Pond = solution.CHydro_Pond
    pondfactor = np.tile(CHydro_Pond, (intervals, 1)) / sum(CHydro_Pond) if sum(CHydro_Pond) != 0 else 0
    MPond_long = np.tile(solution.DischargePond, (len(CHydro_Pond), 1)).transpose() * pondfactor 

    MPV, MBaseload, MPond = map(np.zeros, [(nodes, intervals)] * 3)
#    MWind = map(np.zeros, [(nodes, intervals)] * 1)
    for i, j in enumerate(Nodel):
        MPV[i, :] = solution.GPV[:, np.where(PVl==j)[0]].sum(axis=1)
#        MWind[i, :] = solution.GWind[:, np.where(Windl==j)[0]].sum(axis=1)
        MBaseload[i, :] = solution.baseload[:, np.where(Hydrol==j)[0]].sum(axis=1)
        MPond[i, :] = MPond_long[:, np.where(Hydrol==j)[0]].sum(axis=1)
        """ if solution.node=='Super':
            MInter[i, :] = solution.GInter[:, np.where(Interl==j)[0]].sum(axis=1) """
    MPV, MBaseload, MPond = (MPV.transpose(), MBaseload.transpose(), MPond.transpose()) # Sij-GPV(t, i), Sij-GWind(t, i), MW
#    MWind = MWind.transpose()
    
    MLoad = solution.MLoad # EOLoad(t, j), MW

    defactor = MLoad / MLoad.sum(axis=1)[:, None]
    MDeficit = np.tile(solution.Deficit, (nodes, 1)).transpose() * defactor # MDeficit: EDE(j, t)

    M_minFactors = np.full((intervals, nodes), pow(10,-9)) # Matrix of 10^(-9) required to distribute spillage between nodes when no solar generation
    MPW = MPV + M_minFactors # + MWind
    spfactor = np.divide(MPW, MPW.sum(axis=1)[:, None], where=MPW.sum(axis=1)[:, None]!=0)
    MSpillage = np.tile(solution.Spillage, (nodes, 1)).transpose() * spfactor # MSpillage: ESP(j, t)

    CPHP = solution.CPHP
    pcfactor = np.tile(CPHP, (intervals, 1)) / sum(CPHP) if sum(CPHP) != 0 else 0
    MDischargePH = np.tile(solution.DischargePH, (nodes, 1)).transpose() * pcfactor # MDischarge: DPH(j, t)
    MChargePH = np.tile(solution.ChargePH, (nodes, 1)).transpose() * pcfactor # MCharge: CHPH(j, t)

    CIndia = np.append(np.array([0]*(nodes-len(solution.Interl))), np.nan_to_num(np.array(solution.CInter))) # GW
    india_imports = solution.india_imports # MW
    if CIndia.sum() == 0:
        ifactor = np.tile(CIndia, (intervals, 1))
    else:
        ifactor = np.tile(CIndia, (intervals, 1)) / CIndia.sum()
    MIndia = np.tile(india_imports, (nodes, 1)).transpose() * ifactor

    """ exportNodes = np.array([0,0,0,0,0,0,0,CHydro[0]+CHydro[1], CHydro[2], CHydro[3]+CHydro[4], CHydro[5]+CHydro[6]]) # Electricity exported to each node proportional to the capacity of hydro nearest that interconnection
    efactor = np.tile(exportNodes, (intervals,1)) / sum(exportNodes) if sum(exportNodes) != 0 else 0
    MExport = np.tile(solution.exports, (nodes, 1)).transpose() * efactor """

    print(MLoad.shape,MChargePH.shape,MSpillage.shape,MPV.shape,MIndia.shape,MBaseload.shape,MPond.shape,MDischargePH.shape,MDeficit.shape)

    MImport = MLoad + MChargePH + MSpillage \
              - MPV - MIndia - MBaseload - MPond - MDischargePH - MDeficit # - MWind; EIM(t, j), MW
    
    coverage = solution.coverage
    if len(coverage) > 1:
        # Imorts into external nodes
        IN1CH = -1 * MImport[:, np.where(Nodel=='IN1')[0][0]] if 'IN1' in coverage else np.zeros(intervals)
        IN2TS = -1 * MImport[:, np.where(Nodel=='IN2')[0][0]] if 'IN2' in coverage else np.zeros(intervals)
        IN3SA = MImport[:, np.where(Nodel=='IN3')[0][0]] if 'IN3' in coverage else np.zeros(intervals)
        IN4PE = MImport[:, np.where(Nodel=='IN4')[0][0]] if 'IN4' in coverage else np.zeros(intervals)
        
        # Imports into outer internal nodes
        PEMO = MImport[:, np.where(Nodel=='MO')[0][0]] if 'MO' in coverage else np.zeros(intervals)

        # Imports into inner internal nodes
        CHTH = MImport[:, np.where(Nodel=='CH')[0][0]] - IN1CH if 'CH' in coverage else np.zeros(intervals)
        ZHPE = -1 * MImport[:, np.where(Nodel=='PE')[0][0]] - IN4PE - PEMO if 'PE' in coverage else np.zeros(intervals)
        
        SAZH = MImport[:, np.where(Nodel=='ZH')[0][0]] - ZHPE if 'ZH' in coverage else np.zeros(intervals)
        TSSA = -1 * MImport[:, np.where(Nodel=='SA')[0][0]] - SAZH - IN3SA if 'SA' in coverage else np.zeros(intervals)
        THTS = MImport[:, np.where(Nodel=='TS')[0][0]] - TSSA - IN2TS if 'TS' in coverage else np.zeros(intervals)

        # Check the final node
        THTS1 = -1 * MImport[:, np.where(Nodel=='TH')[0][0]] - CHTH if 'TH' in coverage else np.zeros(intervals)

        assert abs(THTS - THTS1).max() <= 0.1, print('THTS Error', abs(THTS - THTS1).max())

        TDC = np.array([CHTH, THTS, TSSA, SAZH, ZHPE, PEMO, IN1CH, IN2TS, IN3SA, IN4PE]).transpose() # TDC(t, k), MW   
    
    else:
        TDC = np.zeros((intervals, len(solution.TLoss)))

    if output:
        MStoragePH = np.tile(solution.StoragePH, (nodes, 1)).transpose() * pcfactor # SPH(t, j), MWh
        solution.MPV, solution.MInter, solution.MHydro = (MPV, MInter, MHydro)
#        solution.MWind = MWind        
        solution.MDischargePH, solution.MChargePH, solution.MStoragePH = (MDischargePH, MChargePH, MStoragePH)
        solution.MDeficit, solution.MSpillage = (MDeficit, MSpillage)

    return TDC