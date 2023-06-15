import vessel
import pickle
import os
import time
import numpy as np
import pyvista as pv

def saveVessel(vess):
    with open('vessel.temp', 'wb') as file:
        pickle.dump(vess,file)
    os.system('mv vessel.temp vessel.pickle')
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
    simulation_vessel = vessel.Vessel(vesselType="segmentation", segmentationName="line", numLen=42, numCirc=32, numRad=6, rotateSegmentation=True)
    simulation_vessel.radius = 1.1
    simulation_vessel.gnr_step_size = 2.0
    simulation_vessel.zcenter = 0.895
    simulation_vessel.estimateOuterSegmentation = True
    simulation_vessel.constantThickness = True
    simulation_vessel.thickness = 0.0743
    #simulation_vessel.outletPressure = 13000
    simulation_vessel.rotation_matrix = np.loadtxt('ivcRotationMatrix.dat')
    simulation_vessel.tevg = 1
    simulation_vessel.damping = 1e4
    simulation_vessel.penalty = 1e9
    simulation_vessel.simulationExecutable = "~/svFSI-build/svFSI-build/mysvfsi"
    simulation_vessel.numProcessorsSolid = 48
    simulation_vessel.numProcessorsFluid = 96
    simulation_vessel.flipContours = True
    simulation_vessel.flipInlet = True
    simulation_vessel.smoothAttributesValue = 0.1
    simulation_vessel.setInputFileValues()
    os.system('mkdir -p ' + simulation_vessel.outputDir)
    os.system('mkdir -p ' + 'meshIterations')
    os.system('mkdir -p ' + 'meshResults')
    os.system('mkdir -p ' + 'simulationResults')
    os.system('mkdir -p ' + 'materialResults')
    simulation_vessel.initializeVessel()
    simulation_vessel.runFluidIteration()
    simulation_vessel.runSolidIteration()
    simulation_vessel.runMaterialIteration()
    simulation_vessel.currTime = time.time() - startTime + simulation_vessel.startTime
    simulation_vessel.writeStatus(simulation_vessel.currTime, "FSG")
    simulation_vessel.incrementIteration()
    saveVessel(simulation_vessel)

#simulation_vessel.vesselReference = pv.read('mesh_4_0.vtu')
#simulation_vessel.timeIter = 1

while simulation_vessel.timeStep < simulation_vessel.total_time_steps:
    while simulation_vessel.residual > simulation_vessel.tolerance or simulation_vessel.timeIter < 3 or simulation_vessel.skippedFluid == True:
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
