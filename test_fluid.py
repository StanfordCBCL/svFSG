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

if os.path.exists('vessel.pickle'):
    simulation_vessel = loadVessel()
    simulation_vessel.startTime = simulation_vessel.currTime
else:
    simulation_vessel = vessel.Vessel(radius=0.857, thickness=0.07, length=0.857*2, numLen=8, numCirc=12)
    simulation_vessel.simulationExecutable = "-np 24 ~/svFSI-build/svFSI-build/mysvfsi"
    simulation_vessel.initializeVessel()
    simulation_vessel.runFluidIteration()
    os.system('mkdir -p ' + simulation_vessel.outputDir)
    os.system('mkdir -p ' + 'meshIterations')



startTime = time.time()

while simulation_vessel.timeStep < simulation_vessel.max_days:
    simulation_vessel.timeIter = 0
    simulation_vessel.residual = simulation_vessel.tolerance*10.0
    while simulation_vessel.residual > simulation_vessel.tolerance or simulation_vessel.timeIter < 3:
        simulation_vessel.runSolidIteration()
        simulation_vessel.runFluidIteration()
        simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
        simulation_vessel.writeStatus(simulation_vessel.currTime)
        simulation_vessel.incrementIteration()
        saveVessel(simulation_vessel)
    simulation_vessel.incrementTimestep()
    saveVessel(simulation_vessel)
