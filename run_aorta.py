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
    simulation_vessel = vessel.Vessel(vesselType="segmentation", segmentationName="aorta_fsg", numLen=48, numCirc=24)
    simulation_vessel.estimateOuterSegmentation = True
    simulation_vessel.thicknessRatio = 0.0795
    simulation_vessel.outletPressure = 1333.3*100*1.4
    simulation_vessel.inletFlow = -97.0
    simulation_vessel.gnr_step_size = 20
    simulation_vessel.gnr_max_days = 7200
    simulation_vessel.numProcessorsSolid = 48
    simulation_vessel.numProcessorsFluid = 128
    simulation_vessel.damping = 1e4
    simulation_vessel.penalty = 5e8
    simulation_vessel.simulationExecutable = "~/svFSI-build/svFSI-build/mysvfsi"
    simulation_vessel.flipInlet = True
    simulation_vessel.setInputFileValues()
    simulation_vessel.smoothAttributesValue = 0.1
    os.system('mkdir -p ' + simulation_vessel.outputDir)
    os.system('mkdir -p ' + 'meshIterations')
    os.system('mkdir -p ' + 'meshResults')
    os.system('mkdir -p ' + 'simulationResults')
    simulation_vessel.initializeVessel()
    simulation_vessel.runFluidIteration()
    simulation_vessel.runSolidIteration()
    simulation_vessel.runMaterialIteration()
    simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
    simulation_vessel.writeStatus(simulation_vessel.currTime, "FSG")
    simulation_vessel.incrementIteration()
    saveVessel(simulation_vessel)

while simulation_vessel.timeStep < simulation_vessel.total_time_steps:
    while simulation_vessel.residual > simulation_vessel.tolerance or simulation_vessel.timeIter < 3:
        while simulation_vessel.residual > simulation_vessel.tolerance:
            simulation_vessel.runSolidIteration()
            simulation_vessel.skipFluidIteration()
            simulation_vessel.runMaterialIteration()
            simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
            simulation_vessel.writeStatus(simulation_vessel.currTime, "SG")
            simulation_vessel.incrementIteration()
            saveVessel(simulation_vessel)
        simulation_vessel.runFluidIteration()
        simulation_vessel.runSolidIteration()
        simulation_vessel.runMaterialIteration()
        simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
        simulation_vessel.writeStatus(simulation_vessel.currTime, "FSG")
        simulation_vessel.incrementIteration()
        saveVessel(simulation_vessel)
    simulation_vessel.timeIter = 0
    simulation_vessel.residual = simulation_vessel.tolerance*10.0
    simulation_vessel.incrementTimestep()
    saveVessel(simulation_vessel)
