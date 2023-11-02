# Modelling input and assumptions
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

import numpy as np
from Optimisation import scenario, node, percapita, export_flag, import_flag

######### DEBUG ##########
scenario = 'existing'
node = 'Super'
percapita = 3
export_flag = True
import_flag = True
#########################

###### NODAL LISTS ######
Nodel = np.array(['CH', 'TH', 'TS', 'SA', 'ZH', 'PE', 'MO', 'IN1', 'IN2', 'IN3', 'IN4'])
PVl =   np.array(['CH']*1 + ['TH']*1 + ['TS']*1 + ['SA']*1 + ['ZH']*1 + ['PE']*1 + ['MO']*1)
#Windl = np.array(['ME']*1 + ['SB']*1 + ['TE']*1 + ['PA']*1 + ['SE']*1 + ['PE']*1 + ['JO']*1 + ['KT']*1 + ['KD']*1 + ['SW']*1)
Interl = np.array(['IN1']*1 + ['IN2']*1 + ['IN3']*1 + ['IN4']*1) if ((node=='Super') & import_flag) else np.array([]) # Add external interconnections if ASEAN Power Grid scenario
resolution = 1

###### DATA IMPORTS ######
MLoad = np.genfromtxt('Data/electricity{}.csv'.format(percapita), delimiter=',', skip_header=1, usecols=range(4, 4+len(Nodel))) # EOLoad(t, j), MW
TSPV = np.genfromtxt('Data/pv.csv', delimiter=',', skip_header=1, usecols=range(4, 4+len(PVl))) # TSPV(t, i), MW
#TSWind = np.genfromtxt('Data/wind.csv', delimiter=',', skip_header=1, usecols=range(4, 4+len(Windl))) # TSWind(t, i), MW

if 'uni'in scenario:
    assets = np.genfromtxt('Data/assets_uni100.csv'.format(scenario), dtype=None, delimiter=',', encoding=None)[1:, 3:].astype(float)
    constraints = np.genfromtxt('Data/constraints_existing.csv'.format(scenario), dtype=None, delimiter=',', encoding=None)[1:, 3:].astype(float)
    if scenario == 'uni200':
        assets = 2*assets
    elif scenario == 'uni50':
        assets = 0.5*assets
    elif scenario == 'uni0':
        assets = 0*assets
else:
    assets = np.genfromtxt('Data/assets_{}.csv'.format(scenario), dtype=None, delimiter=',', encoding=None)[1:, 3:].astype(float)
    constraints = np.genfromtxt('Data/constraints_{}.csv'.format(scenario), dtype=None, delimiter=',', encoding=None)[1:, 3:].astype(float)

if scenario == 'existing':
    hydrol = np.array(['CH']*2+['MO']*1+['TH']*1+['TS']*1+['ZH']*1)
elif scenario == 'construction':
    hydrol = np.array(['CH']*3+['MO']*3+['TH']*3+['TS']*1+['ZH']*3)
elif scenario == 'all':
    hydrol = np.array(['CH']*6+['MO']*5+['TH']*3+['TS']*2+['ZH']*4+['PE']*1)

""" print(assets)
print(constraints)
print(hydrol) """

CHydro_max, CHydro_RoR, CHydro_Pond = [assets[:, x] * pow(10, -3) for x in range(assets.shape[1])] # CHydro(j), MW to GW
EHydro = constraints[:, 0] # GWh per year
hydroProfiles = np.genfromtxt('Data/RoR_{}.csv'.format(scenario), delimiter=',', skip_header=1, usecols=range(4,4+len(hydrol)), encoding=None).astype(float)

""" print(CHydro_max)
print(CHydro_RoR)
print(CHydro_Pond) """

indiaExportProfiles = hydroProfiles[:,1] # Tala power station is full export to India

#print(indiaExportProfiles)

for i in range(0,len(hydroProfiles[0])):
    hydroProfiles[i,1] = 0
CHydro_Pond[1] = 0
#print(hydroProfiles)

baseload = np.ones((MLoad.shape[0], len(CHydro_RoR)))
#print(baseload)
for i in range(0,MLoad.shape[0]):
    for j in range(0,len(CHydro_RoR)):
        baseload[i,j] = min(hydroProfiles[i,j],CHydro_RoR[j]*pow(10,3)) if CHydro_Pond[j] != 0 else hydroProfiles[i,j]

daily_pondage = np.zeros((MLoad.shape[0], len(CHydro_RoR)))
for i in range(0,MLoad.shape[0]):
    for j in range(0,len(CHydro_RoR)):
        daily_pondage[i,j] = sum(hydroProfiles[i:i+23,j] - baseload[i:i+23,j]) if i % 24 == 0 else daily_pondage[i-1,j]


#print(baseload.shape)
#print(baseload[0:23,0])

""" if export_flag == 'export':
    exports = MLoad.sum(axis=1) - baseload # Export excess hydro to India, assume no export of solar PV
    exports[exports > 0] = 0
elif export_flag == 'no_export':
    exports = np.zeros(MLoad.shape[0])
else:
    print("Export flag error") """

###### CONSTRAINTS ######
# Energy constraints
Hydromax = EHydro.sum() * pow(10,3) # GWh to MWh per year

# Transmission constraints
""" externalImports = 0.05 if node=='Super' else 0
CDC7max, CDC8max, CDC9max, CDC10max = 4 * [externalImports * MLoad.sum() / MLoad.shape[0] / 1000] # 5%: External interconnections: THKD, INSE, PHSB, MW to GW """

###### TRANSMISSION LOSSES ######

# HVDC backbone scenario
dc_flags = np.array([True,True,True,True,True,True,True,True,True,True])
    
""" # HVAC backbone scenario
dc_flags = np.array([False,False,False,False,False,False,False,False,False,False]) """
    
TLoss = []
TDistances = [75, 75, 40, 35, 75, 50, 80, 90, 50, 55] # ['CHTH', 'THTS', 'TSSA', 'SAZH', 'ZHPE', 'PEMO', 'IN1CH', 'IN2TS', 'IN3SA', 'IN4PE']
for i in range(0,len(dc_flags)):
    TLoss.append(TDistances[i]*0.03) if dc_flags[i] else TLoss.append(TDistances[i]*0.07)
TLoss = np.array(TLoss)* pow(10, -3)

###### STORAGE SYSTEM CONSTANTS ######
efficiencyPH = 0.8

###### COST FACTORS ######
factor = np.genfromtxt('Data/factor.csv', delimiter=',', usecols=1)

###### SIMULATION PERIOD ######
firstyear, finalyear, timestep = (2013, 2022, 1)

###### SCENARIO ADJUSTMENTS #######
# Node values
if 'Super' == node:
    coverage = Nodel

""" else:
    coverage = np.array([node])

    MLoad = MLoad[:, np.where(np.in1d(Nodel, coverage)==True)[0]]
    TSPV = TSPV[:, np.where(np.in1d(Nodel, coverage)==True)[0]]
    #TSWind = TSWind[:, np.where(np.in1d(Nodel, coverage)==True)[0]]

    CBaseload = CBaseload[np.where(np.in1d(Nodel, coverage)==True)[0]]
    CHydro = CHydro[np.where(np.in1d(Nodel, coverage)==True)[0]]
    CPeak = CHydro - CBaseload # GW

    EHydro = EHydro[np.where(np.in1d(Nodel, coverage)==True)[0]]

    Hydromax = EHydro.sum() * pow(10,3) # GWh to MWh per year

    baseload = np.ones(MLoad.shape[0]) * CBaseload.sum() * 1000 # GW to MW

    Nodel, PVl, Interl = [x[np.where(np.in1d(x, coverage)==True)[0]] for x in (Nodel, PVl, Interl)]
#    Nodel, PVl, Windl, Interl = [x[np.where(np.in1d(x, coverage)==True)[0]] for x in (Nodel, PVl, Windl, Interl)]
 """
###### DECISION VARIABLE LIST INDEXES ######
intervals, nodes = MLoad.shape
years = int(resolution * intervals / 8760)
pzones = TSPV.shape[1] # Solar PV and wind sites
# wzones = TSWind.shape[1]
# pidx, widx, phidx, bidx = (pzones, pzones + wzones, pzones + wzones + nodes, pzones + wzones + 2*nodes) # Index of solar PV (sites), wind (sites), pumped hydro power (service areas), and battery power (service areas)
pidx, phidx = (pzones, pzones + nodes) # Index of solar PV (sites), wind (sites), pumped hydro power (service areas)
inters = len(Interl) # Number of external interconnections
iidx = phidx + 1 + inters # Index of external interconnections, noting pumped hydro energy (network) and battery energy (network) decision variables after the index of battery power

###### NETWORK CONSTRAINTS ######
energy = (MLoad).sum() * pow(10, -9) * resolution / years # PWh p.a.
#export_annual = exports.sum() * pow(10,-9) * resolution / years #PWh p.a.
contingency_ph = list(0.25 * (MLoad).max(axis=0) * pow(10, -3))[:(nodes)] # MW to GW

#manage = 0 # weeks
#allowance = MLoad.sum(axis=1).max() * 0.05 * manage * 168 * efficiencyPH # MWh
allowance = min(0.00002*np.reshape(MLoad.sum(axis=1), (-1, 8760)).sum(axis=-1)) # Allowable annual deficit of 0.002%, MWh

###### DECISION VARIABLE UPPER BOUNDS ######
pv_ub = [20.] * pzones
phes_ub = [20.] * nodes
phes_s_ub = [200.]
inters_ub = [20.] * inters if import_flag else []

class Solution:
    """A candidate solution of decision variables CPV(i), CWind(i), CPHP(j), S-CPHS(j)"""

    def __init__(self, x):
        self.x = x
        self.MLoad = MLoad
        self.intervals, self.nodes = (intervals, nodes)
        self.resolution = resolution
        self.baseload = baseload
        self.indiaExportProfiles = indiaExportProfiles
        self.daily_pondage = daily_pondage

        self.CPV = list(x[: pidx]) # CPV(i), GW
#        self.CWind = list(x[pidx: widx]) # CWind(i), GW
        self.GPV = TSPV * np.tile(self.CPV, (intervals, 1)) * pow(10, 3) # GPV(i, t), GW to MW
#        self.GWind = TSWind * np.tile(self.CWind, (intervals, 1)) * pow(10, 3) # GWind(i, t), GW to MW

#        self.CPHP = list(x[widx: phidx]) # CPHP(j), GW
        self.CPHP = list(x[pidx: phidx]) # CPHP(j), GW
        self.CPHS = x[phidx] # S-CPHS(j), GWh
        self.efficiencyPH = efficiencyPH

        self.CInter = list(x[phidx+1: ]) if node == 'Super' else len(Interl)*[0] #CInter(j), GW
        self.GIndia = np.tile(self.CInter, (intervals, 1)) * pow(10,3) # GInter(j, t), GW to MW

        self.Nodel, self.PVl, self.Hydrol = (Nodel, PVl, hydrol)
        self.Interl = Interl
#        self.Windl = Windl
        self.node = node
        self.scenario = scenario
        self.export_flag = export_flag
        self.import_flag = import_flag
        self.allowance = allowance
        self.coverage = coverage
        self.TLoss = TLoss

        """ self.CBaseload, self.CPeak = (CBaseload, CPeak)
        self.CHydro = CHydro # GW """
        self.CHydro_RoR = CHydro_RoR
        self.CHydro_Pond = CHydro_Pond
        self.CHydro_max = CHydro_max
        

    def __repr__(self):
        """S = Solution(list(np.ones(64))) >> print(S)"""
        return 'Solution({})'.format(self.x)