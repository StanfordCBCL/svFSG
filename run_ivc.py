import vessel
import pickle
import os
import time

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
    simulation_vessel = vessel.Vessel(radius=0.8573, thickness=0.0743, length=0.8573, numLen=4, numCirc=24, numRad=4)
    simulation_vessel.gnr_step_size = 4.0
    #simulation_vessel.inletFlow = 1.2*simulation_vessel.inletFlow
    #simulation_vessel.outletPressure = (1/1.2)*simulation_vessel.outletPressure
    #simulation_vessel.adaptiveMesh = True
    #simulation_vessel.tevg = 1
    simulation_vessel.damping = 1e4
    simulation_vessel.penalty = 1e8
    simulation_vessel.tolerance = 1e-4
    simulation_vessel.numProcessorsSolid = 4
    simulation_vessel.numProcessorsFluid = 4
    simulation_vessel.simulationExecutable = "~/svFSI-build/svFSI-build/mysvfsi"
    simulation_vessel.setInputFileValues()
    os.system('mkdir -p ' + simulation_vessel.outputDir)
    os.system('mkdir -p ' + 'meshIterations')
    os.system('mkdir -p ' + 'meshResults')
    os.system('mkdir -p ' + 'simulationResults')
    os.system('mkdir -p ' + 'materialResults')
    simulation_vessel.initializeVessel()

while simulation_vessel.timeStep < simulation_vessel.total_time_steps:
    while simulation_vessel.residual > simulation_vessel.tolerance or simulation_vessel.timeIter < 3:
        simulation_vessel.runSolidIteration()
        simulation_vessel.runFluidIteration()
        simulation_vessel.runMaterialIteration()
        simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
        simulation_vessel.writeStatus(simulation_vessel.currTime)
        simulation_vessel.incrementIteration()
        saveVessel(simulation_vessel)
    simulation_vessel.timeIter = 0
    simulation_vessel.residual = simulation_vessel.tolerance*10.0
    simulation_vessel.incrementTimestep()
    saveVessel(simulation_vessel)
