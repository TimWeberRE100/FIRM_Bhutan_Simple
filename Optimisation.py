# To optimise the configurations of energy generation, storage and transmission assets
# Copyright (c) 2019, 2020 Bin Lu, The Australian National University
# Licensed under the MIT Licence
# Correspondence: bin.lu@anu.edu.au

from scipy.optimize import differential_evolution
from argparse import ArgumentParser
import datetime as dt
import csv

parser = ArgumentParser()
parser.add_argument('-i', default=400, type=int, required=False, help='maxiter=4000, 400')
parser.add_argument('-p', default=10, type=int, required=False, help='popsize=2, 10')
parser.add_argument('-m', default=0.5, type=float, required=False, help='mutation=0.5')
parser.add_argument('-r', default=0.3, type=float, required=False, help='recombination=0.3')
parser.add_argument('-e', default=20, type=int, required=False, help='per-capita electricity = 3, 6, 20 MWh/year; prefix 2 for flat curves; 3100, 3200 for universal flat')
parser.add_argument('-n', default='Super', type=str, required=False, help='Super, CH, TH...')
parser.add_argument('-s', default='existing', type=str, required=False, help='all, construction, existing, uni200, uni100, uni50, uni0')
parser.add_argument('-z', default='export', type=str, required=False, help='export, no_export')
args = parser.parse_args()

scenario = args.s
node = args.n
percapita = args.e
export_flag = (args.z == 'export')

from Input import *
from Simulation import Reliability
from Network import Transmission

def F(x):
    """This is the objective function."""

    # Initialise the optimisation
    S = Solution(x)

    CIndia = np.nan_to_num(np.array(S.CInter))
    #print("India: ", CIndia," GW")

    # Simulation with only baseload
    Deficit_energy1, Deficit_power1, Deficit1, DischargePH1, DischargePond1, Spillage1 = Reliability(S, baseload=baseload, india_imports=np.zeros(intervals), daily_pondage=daily_pondage)
    Max_deficit1 = np.reshape(Deficit1, (-1, 8760)).sum(axis=-1) # MWh per year
    PIndia = Deficit1.max() * pow(10, -3) # GW

    GIndia = resolution * (Max_deficit1).max() / efficiencyPH

    PenPower = max(0, PIndia - CIndia.sum()) * pow(10,3)
    PenEnergy = 0
    
    # Simulation with baseload, all existing capacity, and all hydrogen
    Deficit_energy, Deficit_power, Deficit, DischargePH, DischargePond, Spillage = Reliability(S, baseload=baseload, india_imports=np.ones(intervals) * CIndia.sum() * pow(10,3), daily_pondage=daily_pondage)

    # Deficit penalty function
    PenDeficit = max(0, Deficit.sum() * resolution - S.allowance)

    # India import profile
    india_imports = np.clip(Deficit1, 0, CIndia.sum() * pow(10,3)) # MW

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
    GHydro = resolution * (baseload.sum() + DischargePond.sum()) / efficiencyPH / years
    
    # Average annual electricity imported through external interconnections
    GIndia = resolution * india_imports.sum() / years / efficiencyPH

    # Assume all spillage, curtailment, and Hydro_CH2 generation is exported to india
    export_annual = (Spillage.sum() + indiaExportProfiles.sum()) * resolution / years
    Ghydro_CH2 = indiaExportProfiles.sum() * resolution / years

    # Levelised cost of electricity calculation
    cost = factor * np.array([sum(S.CPV), 0, GIndia * pow(10,-6), sum(S.CPHP), S.CPHS, GPHES] + list(CDC) + [sum(S.CPV), 0, (GHydro + Ghydro_CH2) * pow(10, -6), 0, 0]) # $b p.a.
    cost = cost.sum()
    loss = np.sum(abs(TDC), axis=0) * TLoss
    loss = loss.sum() * pow(10, -9) * resolution / years # PWh p.a.
    LCOE = cost / abs(energy - loss) 

    #print("Costs: ",energy, loss, cost, LCOE)
    
    ########### INCLUSION OF EXPORT ENERGY IN LCOE CALCULATION?###############
    ### IF NO - REMOVE GHYDRO_CH2 FROM COSTS CALCULATION, REMOVE SPILLAGE FROM ENERGY CALCULATION
    
    with open('Results/record_{}_{}_{}_{}.csv'.format(node,scenario,percapita, export_flag), 'a', newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(np.append(x,[PenDeficit+PenEnergy+PenPower+PenDC,PenDeficit,PenEnergy,PenPower,LCOE]))

    Func = LCOE + PenDeficit + PenEnergy + PenPower + PenDC
    
    return Func

if __name__=='__main__':
    starttime = dt.datetime.now()
    print("Optimisation starts at", starttime)

#    lb = [0.]       * pzones + [0.]     * wzones + contingency_ph   + contingency_b     + [0.]      + [0.]     + [0.]    * inters + [0.] * nodes
#    ub = [10000.]   * pzones + [300]    * wzones + [10000.] * nodes + [10000.] * nodes  + [100000.] + [100000] + [1000.] * inters + [50.] * nodes

    lb = [0.] * pzones + [0.] * nodes + [0.] + [0.] * inters
    ub = pv_ub + phes_ub + phes_s_ub + inters_ub

    # start = np.genfromtxt('Results/init.csv', delimiter=',')

    result = differential_evolution(func=F, bounds=list(zip(lb, ub)), tol=0, # init=start,
                                    maxiter=args.i, popsize=args.p, mutation=args.m, recombination=args.r,
                                    disp=True, polish=False, updating='deferred', workers=-1) ###### CHANGE WORKERS BACK TO -1

    with open('Results/Optimisation_resultx_{}_{}_{}_{}.csv'.format(node,scenario,percapita,export_flag), 'w', newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(result.x)

    endtime = dt.datetime.now()
    print("Optimisation took", endtime - starttime)

    from Dispatch import Analysis
    Analysis(result.x,'_{}_{}_{}_{}.csv'.format(node,scenario,percapita,export_flag))
