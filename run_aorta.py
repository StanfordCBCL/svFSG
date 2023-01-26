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
    simulation_vessel = vessel.Vessel(gnrStepSize=20.0, vesselType="segmentation", segmentationName="aorta_fsg", numLen=36, numCirc=24)
    simulation_vessel.estimateOuterSegmentation = True
    simulation_vessel.thicknessRatio = 0.0795
    simulation_vessel.outletPressure = 1333.3*100
    simulation_vessel.inletFlow = -97.0
    simulation_vessel.damping = 8e5
    simulation_vessel.penalty = 1e8
    simulation_vessel.simulationExecutable = "-np 24 ~/svFSI-build/svFSI-build/mysvfsi"
    simulation_vessel.setInputFileValues()
    os.system('mkdir -p ' + simulation_vessel.outputDir)
    os.system('mkdir -p ' + 'meshIterations')
    os.system('mkdir -p ' + 'meshResults')
    simulation_vessel.initializeVessel()

startTime = time.time()

while simulation_vessel.timeStep < simulation_vessel.max_days:
    while simulation_vessel.residual > simulation_vessel.tolerance or simulation_vessel.timeIter < 3:
        simulation_vessel.runSolidIteration()
        simulation_vessel.estimateFluidIteration()
        simulation_vessel.runMaterialIteration()
        simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
        simulation_vessel.writeStatus(simulation_vessel.currTime)
        simulation_vessel.incrementIteration()
        saveVessel(simulation_vessel)
        simulation_vessel.damping = 1e5
        simulation_vessel.setInputFileValues()

    simulation_vessel.timeIter = 0
    simulation_vessel.residual = simulation_vessel.tolerance*10.0
    simulation_vessel.incrementTimestep()
    saveVessel(simulation_vessel)
