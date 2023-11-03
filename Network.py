# A transmission network model to calculate inter-regional power flows
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

import numpy as np

def Transmission(solution, domestic_only=False, export_only=False, output=False):
    """TDC = Network.Transmission(S)"""

    Nodel, PVl, Interl, Hydrol, Windl, expl = (solution.Nodel, solution.PVl, solution.Interl, solution.Hydrol, solution.Windl, solution.expl)
    intervals, nodes, inters = (solution.intervals, solution.nodes, len(Interl))
    
    CHydro_Pond = solution.CHydro_Pond
    pondfactor = np.tile(CHydro_Pond, (intervals, 1)) / sum(CHydro_Pond) if sum(CHydro_Pond) != 0 else 0
    MPond_long = np.tile(solution.DischargePond, (len(CHydro_Pond), 1)).transpose() * pondfactor 

    MPV, MBaseload, MPond, MWind = map(np.zeros, [(nodes, intervals)] * 4)
    for i, j in enumerate(Nodel):
        MPV[i, :] = solution.GPV[:, np.where(PVl==j)[0]].sum(axis=1)
        MWind[i, :] = solution.GWind[:, np.where(Windl==j)[0]].sum(axis=1)
        MBaseload[i, :] = solution.baseload[:, np.where(Hydrol==j)[0]].sum(axis=1)
        MPond[i, :] = MPond_long[:, np.where(Hydrol==j)[0]].sum(axis=1)
    MPV, MBaseload, MPond, MWind = (MPV.transpose(), MBaseload.transpose(), MPond.transpose(), MWind.transpose()) # Sij-GPV(t, i), Sij-GWind(t, i), MW
    
    MLoad = solution.MLoad # EOLoad(t, j), MW

    defactor = MLoad / MLoad.sum(axis=1)[:, None]
    MDeficit = np.tile(solution.Deficit, (nodes, 1)).transpose() * defactor # MDeficit: EDE(j, t)

    CIndia = np.append(np.array([0]*(nodes-len(solution.Interl))), np.nan_to_num(np.array(solution.CInter))) # GW

    CPHP = solution.CPHP
    pcfactor = np.tile(CPHP, (intervals, 1)) / sum(CPHP) if sum(CPHP) != 0 else 0
    MDischargePH = np.tile(solution.DischargePH, (nodes, 1)).transpose() * pcfactor # MDischarge: DPH(j, t)
    MChargePH = np.tile(solution.ChargePH, (nodes, 1)).transpose() * pcfactor # MCharge: CHPH(j, t)

    india_imports = solution.india_imports # MW
    if CIndia.sum() == 0:
        ifactor = np.tile(CIndia, (intervals, 1))
    else:
        ifactor = np.tile(CIndia, (intervals, 1)) / CIndia.sum()
    MIndia = np.tile(india_imports, (nodes, 1)).transpose() * ifactor

    efactor = np.array([0,0,0,0,0,0,0,1,0,0,0])
    ch2factor = np.array([0,1,0,0,0,0,0,0,0,0,0])
    MExport = np.tile(solution.indiaExportProfiles, (nodes, 1)).transpose() * efactor
    MHydro_CH2 = np.tile(solution.indiaExportProfiles, (nodes, 1)).transpose() * ch2factor

    if solution.export_flag:
        CHydro_nodes = np.zeros(nodes)
        for j in range(0,len(Nodel)):
            CHydro_nodes[j] = solution.CHydro_max[expl==Nodel[j]].sum()
        expfactor = np.tile(CHydro_nodes, (intervals, 1)) / sum(CHydro_nodes) if sum(CHydro_nodes) != 0 else 0
        MSpillage_exp = np.tile(solution.Spillage, (nodes, 1)).transpose() * expfactor # MSpillage: ESP(j, t)
        MSpillage = np.zeros((nodes, intervals)).transpose()

        BaseloadDomestic = MChargePH.sum(axis=1) + MLoad.sum(axis=1) - MIndia.sum(axis=1) - MDischargePH.sum(axis=1) - MDeficit.sum(axis=1)
        BaseloadDomestic[BaseloadDomestic < 0] = 0
        BaseloadExports = MBaseload.sum(axis=1) - BaseloadDomestic
        BaseloadExports[BaseloadExports < 0] = 0
        b1factor = np.divide(MBaseload, MBaseload.sum(axis=1)[:, None], where=MBaseload.sum(axis=1)[:, None]!=0)
        MBaseload_exp  = np.tile(BaseloadExports, (nodes, 1)).transpose() * b1factor

        SolarDomestic = MChargePH.sum(axis=1) + MLoad.sum(axis=1) - MIndia.sum(axis=1) - MBaseload.sum(axis=1) - MDischargePH.sum(axis=1) - MDeficit.sum(axis=1)
        SolarDomestic[SolarDomestic < 0] = 0
        SolarExports = MPV.sum(axis=1) - SolarDomestic        
        SolarExports[SolarExports < 0] = 0
        s1factor = np.divide(MPV, MPV.sum(axis=1)[:, None], where=MPV.sum(axis=1)[:, None]!=0)
        MPV_exp  = np.tile(SolarExports, (nodes, 1)).transpose() * s1factor
        
        WindDomestic = MChargePH.sum(axis=1) + MLoad.sum(axis=1) - MIndia.sum(axis=1) - MBaseload.sum(axis=1) - MPV.sum(axis=1) - MDischargePH.sum(axis=1) - MDeficit.sum(axis=1)
        WindDomestic[WindDomestic < 0] = 0
        WindExports = MWind.sum(axis=1) - WindDomestic      
        WindExports[WindExports < 0] = 0  
        w1factor = np.divide(MWind, MWind.sum(axis=1)[:, None], where=MWind.sum(axis=1)[:, None]!=0)
        MWind_exp  = np.tile(WindExports, (nodes, 1)).transpose() * w1factor

        PondDomestic = MChargePH.sum(axis=1) + MLoad.sum(axis=1) - MPV.sum(axis=1) - MWind.sum(axis=1) - MIndia.sum(axis=1) - MBaseload.sum(axis=1) - MDischargePH.sum(axis=1) - MDeficit.sum(axis=1)
        PondDomestic[PondDomestic < 0] = 0
        PondExports = MPond.sum(axis=1) - PondDomestic
        PondExports[PondExports < 0]
        p1factor = np.divide(MPond, MPond.sum(axis=1)[:, None], where=MPond.sum(axis=1)[:, None]!=0)
        MPond_exp  = np.tile(PondExports, (nodes, 1)).transpose() * p1factor
        
    else:
        M_minFactors = np.full((intervals, nodes), pow(10,-9)) # Matrix of 10^(-9) required to distribute spillage between nodes when no solar generation
        MPW = MPV + M_minFactors + MWind + MPond
        spfactor = np.divide(MPW, MPW.sum(axis=1)[:, None], where=MPW.sum(axis=1)[:, None]!=0)
        MSpillage = np.tile(solution.Spillage, (nodes, 1)).transpose() * spfactor # MSpillage: ESP(j, t)
        MSpillage_exp = np.zeros((nodes, intervals)).transpose()

    if domestic_only:
        MImport = MLoad + MChargePH + MSpillage + \
                MBaseload_exp + MPV_exp + MWind_exp + MPond_exp \
                - MPV - MWind - MIndia - MBaseload - MPond - MDischargePH - MDeficit # EIM(t, j), MW
    elif export_only:
        MImport = MSpillage_exp + MExport - MHydro_CH2 - MBaseload_exp - MPV_exp - MWind_exp - MPond_exp
    else:
        MImport = MLoad + MChargePH + MSpillage_exp + MSpillage + MExport \
              - MPV - MWind - MIndia - MBaseload - MPond - MDischargePH - MDeficit - MHydro_CH2 # EIM(t, j), MW
    
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
        solution.MPV, solution.MIndia, solution.MBaseload, solution.MPond, solution.MWind = (MPV, MIndia, MBaseload, MPond, MWind)
        solution.MDischargePH, solution.MChargePH, solution.MStoragePH = (MDischargePH, MChargePH, MStoragePH)
        solution.MDeficit, solution.MSpillage, solution.MSpillage_exp, solution.MExport = (MDeficit, MSpillage, MSpillage_exp, MExport)
        solution.MBaseload_exp, solution.MPV_exp, solution.MWind_exp, solution.MPond_exp = (MBaseload_exp, MPV_exp, MWind_exp, MPond_exp)

    return TDC