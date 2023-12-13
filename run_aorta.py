import vessel
import pickle
import os
import time
import numpy as np

def saveVessel(vess):
    with open('vessel.pickle', 'wb') as file:
        pickle.dump(vess,file)
    return

def loadVessel():
    with open('vessel.pickle', 'rb') as file:
        vess = pickle.load(file)
    return vess

os.system("mpiexec python3 utils_init_vessel.py")
startTime = time.time()

if os.path.exists('vessel.pickle'):
    simulation_vessel = loadVessel()
    simulation_vessel.startTime = simulation_vessel.currTime
else:
    simulation_vessel = vessel.Vessel(radius=0.0647, thickness=0.004, length=1.5, numLen=48, numCirc=24, numRad=4)
    simulation_vessel.outletPressure = 1333.22*105
    simulation_vessel.inletFlow = -97.0
    simulation_vessel.gnr_step_size = 1
    simulation_vessel.gnr_max_days = 720
    simulation_vessel.damping = 1e7
    simulation_vessel.penalty = 1e9
    simulation_vessel.tolerance = 1e-3
    simulation_vessel.simulationExecutable = "~/svFSI-build/svFSI-build/mysvfsi"
    simulation_vessel.setInputFileValues()
    simulation_vessel.numProcessorsSolid = 72
    simulation_vessel.smoothAttributesValue = 0.1
    simulation_vessel.aneurysm = 1
    simulation_vessel.zod = 0.3
    simulation_vessel.zapex = 0
    simulation_vessel.thetaod = np.pi/2
    simulation_vessel.thetaapex = np.pi
    os.system('mkdir -p ' + simulation_vessel.outputDir)
    os.system('mkdir -p ' + 'meshIterations')
    os.system('mkdir -p ' + 'meshResults')
    os.system('mkdir -p ' + 'simulationResults')
    os.system('mkdir -p ' + 'materialResults')
    simulation_vessel.initializeVessel()
    simulation_vessel.estimateFluidIteration()
    simulation_vessel.runSolidIteration()
    simulation_vessel.estimateFluidIteration()
    simulation_vessel.runMaterialIteration()
    simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
    simulation_vessel.writeStatus(simulation_vessel.currTime, "SG")
    simulation_vessel.incrementIteration()
    saveVessel(simulation_vessel)

while simulation_vessel.timeStep < simulation_vessel.total_time_steps:
    while simulation_vessel.residual > simulation_vessel.tolerance or simulation_vessel.timeIter < 5:
        simulation_vessel.runSolidIteration()
        simulation_vessel.estimateFluidIteration()
        simulation_vessel.runMaterialIteration()
        simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
        simulation_vessel.writeStatus(simulation_vessel.currTime)
        simulation_vessel.incrementIteration()
        saveVessel(simulation_vessel)
    simulation_vessel.runSystoleDiastole()
    simulation_vessel.timeIter = 0
    simulation_vessel.residual = simulation_vessel.tolerance*10.0
    simulation_vessel.incrementTimestep()
    saveVessel(simulation_vessel)
