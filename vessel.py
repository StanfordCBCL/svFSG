from utils import *
import time
from cvessel import cvessel
import sympy

# This class provides functionality to running an FSG simulation
class Vessel():

    def __init__(self, **kwargs):

        self.vesselType = kwargs.get("vesselType", "cylinder")
        self.radius = kwargs.get("radius", 1.0)
        self.thickness = kwargs.get("thickness", 0.1)
        self.length = kwargs.get("length", 10.0)
        self.torusRadius = kwargs.get("torusRadius", 10.0)
        self.torusFraction = kwargs.get("torusFraction", 0.5)
        self.vesselSegmentations = kwargs.get("vesselSegmentations", "")
        self.vesselReference = None
        self.vesselSolid = None
        self.vesselFluid = None
        self.solidResult = None
        self.fluidResult = None
        self.sigma_h = None
        self.tau_h = None
        self.zod = 1
        self.zapex = 0
        self.thetaod = 1
        self.thetaapex = 0

        self.nG = 8
        self.total_time_steps = 0
        self.gnr_max_days = kwargs.get("gnr_max_days", 720)
        self.gnr_step_size = kwargs.get("gnr_step_size", 1.0)
        self.prefix = kwargs.get("prefix", "./")
        self.numRad = kwargs.get("numRad", 4)
        self.numCirc = kwargs.get("numCirc", 12)
        self.numLen = kwargs.get("numLen", 12)
        self.aneurysm = kwargs.get("aneurysm", 0)
        self.tevg = kwargs.get("tevg", 0)
        self.timeStep = 0
        self.timeIter = 0
        self.omega = 0.5
        self.inletFlow = kwargs.get("inletFlow", -20)
        self.outletPressure = kwargs.get("outletPressure", 6150)
        self.residual = 1.0
        self.residualType = None
        self.viscosity = kwargs.get("viscosity", 0.04)
        self.simulationInputDirectory = kwargs.get("simulationInputDirectory", "FolderSimulationInputFiles")
        self.simulationExecutable = kwargs.get("simulationExecutable","~/svFSI-build/svFSI-build/mysvfsi")
        self.vesselName = kwargs.get("vesselName", "vessel")
        self.resultNum = kwargs.get("resultNum",1000)
        self.resultDir = kwargs.get("resultDir","results")
        self.fluidDir = kwargs.get("fluidDir","fluid-results")
        self.outputDir = kwargs.get("outputDir","Outputs")
        self.estimateWSS = kwargs.get("estimateWSS", False)
        self.predictMethod = kwargs.get("predictMethod", "aitken")
        self.damping = kwargs.get("damping", 1e4)
        self.penalty = kwargs.get("penalty", 1e8)
        self.startTime = 0.0
        self.currTime = 0.0
        self.tolerance = kwargs.get("tolerance", 1e-2)
        self.segmentationName = kwargs.get("segmentationName", None)
        self.estimateOuterSegmentation = kwargs.get("estimateOuterSegmentation",False)
        self.thicknessRatio = kwargs.get("thicknessRatio",0.1)
        self.constantThickness = kwargs.get("constantThickness",False)
        self.adaptiveMesh = kwargs.get("adaptiveMesh",False)
        self.rotation_matrix = np.eye(3)
        self.scaling_value = 1.0
        self.numProcessorsSolid = None
        self.numProcessorsFluid = None
        self.numProcessorsFluidSolid = None
        self.zcenter = 0.0
        self.nq = 8
        self.iq_eps = 1e-12
        self.mat_W = []
        self.mat_V = []
        self.mat_D = []
        self.cvessels = []
        self.flipContours = False
        self.flipInlet = False
        self.averageStress = True
        self.averageVolume = True
        self.solidLinearSolverType = "GMRES"
        self.smoothAttributesValue = 0.0
        self.anastomosisGraftRatio = 1.0
        self.anastomosisConnectionRatio = 1.0
        self.anastomosisLengthRatio = 0.5
        self.anastomosisTransitionLength = 0.25
        self.skippedFluid = False

    def propogateIQNILS(self):
        self.mat_W = []
        self.mat_V = []
        self.mat_D = []

        pointLocatorSolid = vtk.vtkPointLocator()
        pointLocatorSolid.SetDataSet(self.solidResult)
        pointLocatorSolid.BuildLocator()

        numPts = self.vesselReference.GetNumberOfPoints()
        numCells = self.vesselReference.GetNumberOfCells()

        dcurr = []
        dprev = []

        for q in range(numPts):
            originalCoordinate = np.array(self.vesselReference.GetPoint(q))
            displacement_prev = np.array(self.vesselReference.GetPointData().GetArray("displacements").GetTuple3(q))
            currentCoordinate = originalCoordinate + displacement_prev

            pointIdSolid = pointLocatorSolid.FindClosestPoint(currentCoordinate)
            solidCoordinate = np.array(self.solidResult.GetPoint(pointIdSolid))
            solidDispacement = np.array(self.solidResult.GetPointData().GetArray("Displacement").GetTuple3(pointIdSolid))

            displacement_curr = solidCoordinate + solidDispacement - originalCoordinate
            residual_curr = displacement_curr - displacement_prev
            self.vesselReference.GetPointData().GetArray("residual_curr").SetTuple(q, residual_curr)

            dcurr.append(displacement_curr)
            dprev.append(displacement_prev)


        rcurr = np.array(self.vesselReference.GetPointData().GetArray("residual_curr")).flatten()
        rprev = np.array(self.vesselReference.GetPointData().GetArray("residual_prev")).flatten()

        dcurr = np.array(dcurr).flatten()
        dprev = np.array(dprev).flatten()

        self.mat_W = []
        self.mat_V = []
        self.mat_D = []
        self.mat_D.append(dcurr)
        vnew =  dprev + 0.5*rcurr
        dnew = vnew.reshape((-1,3))

        self.timeIter = 1


        return


    def writeStatus(self, currTime, extra=""):
        with open('svDriverIterations','a') as f:
            print("%d %d %.3e %s %5.4f %5.2f %s" %(self.timeStep,self.timeIter,self.residual, self.residualType, self.omega, currTime, extra), file=f)

    def setTime(self, timeVal):
        self.time = timeVal
        return

    def incrementTimestep(self):
        self.timeStep+=1
        return

    def incrementIteration(self):
        self.timeIter+=1
        return

    def runFluidIteration(self):
        self.updateSolid()
        self.saveSolid()
        self.updateFluid()
        self.saveFluid()
        self.runFluid()
        self.updateFluidResults()
        self.appendIterfaceResult()
        self.skippedFluid = False
        return

    def estimateFluidIteration(self):
        self.estimateIterfaceResult()
        return

    def skipFluidIteration(self):
        numCells = self.vesselReference.GetNumberOfCells()
        for q in range(numCells):
            wss_prev = self.vesselReference.GetCellData().GetArray('wss_curr').GetTuple1(q)
            self.vesselReference.GetCellData().GetArray('wss_prev').SetTuple1(q,wss_prev)
        self.skippedFluid = True
        return

    def runFluidSolidIteration(self):
        self.updateSolid()
        self.saveSolid()
        self.updateFluid()
        self.saveFluid()
        self.runFluidSolid()
        self.updateFluidSolidResults()
        self.appendSolidResult()
        self.appendIterfaceResult()
        self.skippedFluid = False
        return

    def runSolidIteration(self):
        self.updateSolid()
        self.saveSolid()
        self.runSolid()

        with open(self.resultDir+'/histor.dat') as f:
            if 'NaN' in f.read():
                print("Simulation has NaN! Reducing omega.", file=sys.stderr)
                self.appendReducedResult()
                return

        self.updateSolidResults()
        self.appendSolidResult()
        return

    def runMaterialIteration(self):
        self.updateMaterial()
        self.updateReference()
        self.saveReference()
        self.checkResidual()
        return

    def checkResidual(self):
        tolTypes = ["disp", "sigma_inv", "wss", "jac"]
        tol1 = 100.0*np.max(np.linalg.norm(self.vesselReference.get_array('residual_curr'),axis=1))
        tol2 = np.max(abs((self.vesselReference.get_array('inv_prev')-self.vesselReference.get_array('inv_curr'))/self.sigma_h))
        tol3 = np.max(abs((self.vesselReference.get_array('wss_prev')-self.vesselReference.get_array('wss_curr'))/self.tau_h))
        tol4 = np.max(abs(self.vesselReference.get_array('varWallProps')[:,37::45] - 1.0))
        tolVals = np.array([tol1,tol2,tol3,tol4])
        self.residual = np.max(tolVals)
        self.residualType = tolTypes[tolVals.argmax()]
        return

    def runSystoleDiastole(self):
        os.system('rm -rf '+self.resultDir +'/*')
        if self.timeIter == 0 and self.timeStep == 0:
            if self.numProcessorsSolid is not None:
                os.system("srun -n " + str(self.numProcessorsSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_systole_mm.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_systole_mm.mfs")
            os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/solid_systole_' + str(self.timeStep) + '.vtu')
            if self.numProcessorsSolid is not None:
                os.system("srun -n " + str(self.numProcessorsSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_diastole_mm.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_diastole_mm.mfs")
            os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/solid_diastole_' + str(self.timeStep) + '.vtu')
        else:
            if self.numProcessorsSolid is not None:
                os.system("srun -n " + str(self.numProcessorsSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_systole.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_systole.mfs")
            os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/solid_systole_' + str(self.timeStep) + '.vtu')
            if self.numProcessorsSolid is not None:
                os.system("srun -n " + str(self.numProcessorsSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_diastole.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_diastole.mfs")
            os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/solid_diastole_' + str(self.timeStep) + '.vtu')
        print("Solid simulation finished.")
        return


    def runSolid(self):
        os.system('rm -rf '+self.resultDir +'/*')
        if self.timeIter == 0 and self.timeStep == 0:
            if self.numProcessorsSolid is not None:
                os.system("srun -n " + str(self.numProcessorsSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_mm.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_mm.mfs")
        else:
            if self.numProcessorsSolid is not None:
                os.system("srun -n " + str(self.numProcessorsSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_aniso.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_aniso.mfs")
        print("Solid simulation finished.")
        os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/solid_' + str(self.timeStep) + '.vtu')
        return

    def runFluid(self):
        os.system('rm -rf '+self.resultDir +'/*')
        if self.numProcessorsFluid is not None:
            os.system("srun -n " + str(self.numProcessorsFluid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_fluid.mfs")
        else:
            os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_fluid.mfs")
        with open(self.resultDir+'/histor.dat') as f:
            if 'NaN' in f.read():
                raise RuntimeError("Simulation has NaN!")
        print("Fluid simulation finished.")
        os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/fluid_' + str(self.timeStep) + '.vtu')
        return

    def runFluidSolid(self):
        os.system('rm -rf '+self.resultDir +'/*')
        if self.timeIter == 0 and self.timeStep == 0:
            if self.numProcessorsFluidSolid is not None:
                os.system("srun -n " + str(self.numProcessorsFluidSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_mm.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_mm.mfs")
        else:
            if self.numProcessorsFluidSolid is not None:
                os.system("srun -n " + str(self.numProcessorsFluidSolid) + " " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_aniso.mfs")
            else:
                os.system("srun " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_aniso.mfs")
        with open(self.resultDir+'/histor.dat') as f:
            if 'NaN' in f.read():
                raise RuntimeError("Simulation has NaN!")
        print("FSI simulation finished.")
        os.system('cp ' + self.resultDir + '/result_'+str(self.resultNum)+'.vtu simulationResults/fluidsolid_' + str(self.timeStep) + '.vtu')

        return

        
    def initializeVessel(self):
        if self.vesselType == "cylinder":
            self.vesselReference = self.initializeCylinder()
            self.estimateIterfaceResult()
        elif self.vesselType == "segmentation":
            self.vesselReference = self.initializeSegmentation()
            self.estimateIterfaceResult()
        elif self.vesselType == "anastomosis":
            self.vesselReference = self.initializeAnastomosis()
            self.estimateIterfaceResult()
        elif self.vesselType == "torus":
            self.vesselReference = self.initializeTorus()
            self.estimateIterfaceResult()
        else:
            print("TODO: Implement new initializer")
        self.saveReference()

        self.total_time_steps = self.gnr_max_days//self.gnr_step_size
        self.cvessels = [cvessel() for i in range(self.vesselReference.GetNumberOfCells() * self.nG)]

        return

    def updateMaterial(self):
        numCells = self.vesselReference.GetNumberOfCells()
        input_array = []

        for q in range(numCells):
            sigma_inv = self.vesselReference.GetCellData().GetArray('inv_curr').GetTuple1(q)
            tauw_wss = self.vesselReference.GetCellData().GetArray('wss_curr').GetTuple1(q)
            
            simulate = int(self.vesselReference.GetCellData().GetArray('Simulate').GetTuple1(q))
            aneurysmValue = self.vesselReference.GetCellData().GetArray('aneurysmValue').GetTuple1(q)
            tevgValue = self.vesselReference.GetCellData().GetArray('tevgValue').GetTuple1(q)

            #Rotate into GnR membrane frame
            e_r = self.vesselReference.GetCellData().GetArray('e_r').GetTuple(q)
            e_t = self.vesselReference.GetCellData().GetArray('e_t').GetTuple(q)
            e_z = self.vesselReference.GetCellData().GetArray('e_z').GetTuple(q)
            Q = np.array((e_r,e_t,e_z))

            defGrad_mem = []
            #Gauss point values
            for p in range(self.nG):

                defGrad = self.vesselReference.GetCellData().GetArray('defGrad').GetTuple(q)
                defGrad_g = defGrad[p*9:(p+1)*9]
                defGrad_g = np.reshape(defGrad_g, (3,3))
                defGrad_mem_g = np.matmul(Q,np.matmul(defGrad_g, np.transpose(Q)))
                defGrad_mem_g = np.reshape(defGrad_mem_g,9)

                #Queue process

                self.cvessels[q*self.nG + p].prefix = self.outputDir
                self.cvessels[q*self.nG + p].name = 'python_'+str(q)+'_'+str(p)
                self.cvessels[q*self.nG + p].number = q*self.nG + p
                self.cvessels[q*self.nG + p].restart = self.timeStep
                self.cvessels[q*self.nG + p].iteration = self.timeIter
                self.cvessels[q*self.nG + p].simulate = simulate
                self.cvessels[q*self.nG + p].num_days = self.gnr_max_days
                self.cvessels[q*self.nG + p].step_size = self.gnr_step_size
                self.cvessels[q*self.nG + p].sigma_inv = sigma_inv
                self.cvessels[q*self.nG + p].tauw_wss = tauw_wss
                self.cvessels[q*self.nG + p].tevg = tevgValue
                self.cvessels[q*self.nG + p].aneurysm = aneurysmValue

                for i in range(9):
                    self.cvessels[q*self.nG + p].F[i] = defGrad_mem_g[i]

                defGrad_mem = np.hstack((defGrad_mem,defGrad_mem_g))

            self.vesselReference.GetCellData().GetArray('defGrad_mem').SetTuple(q,defGrad_mem)



        numrank = int(np.loadtxt('numrank'))
        cvessels_split = np.array_split(self.cvessels,numrank)
        for i in range(numrank):
            cvessel_file = open("materialResults/cvessel_array_out_"+str(i)+".dat","wb")
            pickle.dump(cvessels_split[i], cvessel_file)
            cvessel_file.close()

        print("Running points...")
        time1 = time.time()
        os.system("mpiexec python3 utils_run_vessel.py")
        time2 = time.time()
        print("Time to run_vessel: " + str(time2 - time1))

        for i in range(numrank):
            cvessel_file = open("materialResults/cvessel_array_in_"+str(i)+".dat","rb")
            cvessels_temp = pickle.load(cvessel_file)
            for vess_temp in cvessels_temp:
                self.cvessels[vess_temp.number] = vess_temp
            cvessel_file.close()
        
        print("Parsing points...")

        for q in range(numCells):
            stiffness_mem_g = []
            sigma_mem_g = []
            J_target_g = []
            J_curr_g = []
            for p in range(self.nG):

                output_data = self.cvessels[q*self.nG + p].out_array

                J_curr =  np.linalg.det(np.reshape(output_data[48:57], (3,3)))
                J_target = output_data[46]/output_data[47]

                stiffness = output_data[1:37]
                sigma = output_data[37:46]

                # Dont forget units 
                stiffness = np.array(stiffness)*10.0
                sigma = np.array(sigma)*10.0

                sigma_inv = output_data[57]
                wss = output_data[58]

                stiffness_mem_g = stiffness_mem_g + stiffness.tolist()
                sigma_mem_g = sigma_mem_g + sigma.tolist()
                J_target_g = J_target_g + [J_target]
                J_curr_g = J_curr_g + [J_curr]


            self.vesselReference.GetCellData().GetArray('stiffness_mem').SetTuple(q,stiffness_mem_g)
            self.vesselReference.GetCellData().GetArray('sigma_mem').SetTuple(q,sigma_mem_g)
            self.vesselReference.GetCellData().GetArray('J_target').SetTuple(q,J_target_g)
            self.vesselReference.GetCellData().GetArray('J_curr').SetTuple(q,J_curr_g)

        return

    def updateSolid(self):
        self.vesselSolid = self.vesselReference.warp_by_vector("displacements")
        arrayNames = self.vesselReference.array_names
        for name in arrayNames:
            if name not in ["GlobalNodeID", "varWallProps", "GlobalElementID", "InnerRegionID", "OuterRegionID", "DistalRegionID", "ProximalRegionID", "StructureID", "Pressure", "Coordinate"]:
                if name in self.vesselSolid.point_data:
                    self.vesselSolid.point_data.remove(name)
                if name in self.vesselSolid.cell_data:
                    self.vesselSolid.cell_data.remove(name)
        return

    def updateFluid(self):
        surf = self.vesselSolid.extract_surface()
        inner = thresholdModel(surf, 'InnerRegionID',0.5,1.5)
        if self.fluidResult is not None:
            tempFluid = self.generateFluidMesh(inner)
            self.vesselFluid = interpolateSolution(self.fluidResult, tempFluid)
        else:
            self.vesselFluid = self.generateFluidMesh(inner)


    def saveSolid(self):
        vol = self.vesselSolid

        os.makedirs(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces', exist_ok=True)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-complete.mesh.vtu',vol)
        surf = vol.extract_surface()
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-complete.mesh.vtp',surf)
        outer = thresholdModel(surf, 'OuterRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_outer.vtp',outer)
        distal = thresholdModel(surf, 'DistalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_inlet.vtp',distal)
        proximal = thresholdModel(surf, 'ProximalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_outlet.vtp',proximal)
        inner = thresholdModel(surf, 'InnerRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_inner.vtp',inner)


        ring_distal = thresholdModel(inner, 'DistalRegionID',0.5,1.5, allScalars=False)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/ring_inlet.vtp',ring_distal)

        ring_proximal = thresholdModel(inner, 'ProximalRegionID',0.5,1.5, allScalars=False)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/ring_outlet.vtp',ring_proximal)

        return

    def saveFluid(self):
        fluidVol = self.vesselFluid

        os.makedirs(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces', exist_ok=True)
        # Add global node id to inner volume
        numPtsVol = fluidVol.GetNumberOfPoints()
        numCellsVol = fluidVol.GetNumberOfCells()

        # Need to order ids to match globalNodeID order
        fluidVol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPtsVol,numPtsVol).astype(int),name="GlobalNodeID"))
        fluidVol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCellsVol,numCellsVol).astype(int),name="GlobalElementID"))

        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-complete.mesh.vtu',fluidVol)
        innerVolSurf = pv.wrap(fluidVol).extract_surface()
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-complete.mesh.vtp',innerVolSurf)

        innerProx = thresholdModel(innerVolSurf, 'ProximalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces/lumen_inlet.vtp',pv.wrap(innerProx))
        innerDist = thresholdModel(innerVolSurf, 'DistalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces/lumen_outlet.vtp',innerDist)
        innerWall = thresholdModel(innerVolSurf, 'OuterRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces/lumen_wall.vtp',innerWall)

        return

    def saveReference(self):
        save_data(self.prefix + 'meshIterations/mesh_' + str(self.timeStep) + '_' + str(self.timeIter) + '.vtu', self.vesselReference)
        save_data(self.prefix + 'meshResults/mesh_' + str(self.timeStep) + '.vtu', self.vesselReference)
        return


    def updateFluidSolidResults(self):
        resultdata = read_data(self.resultDir+'/result_'+str(self.resultNum)+'.vtu', file_format="vtu")
        self.solidResult = thresholdModel(resultdata,'Domain_ID', 1.5, 2.5, extract=False, cell=True)
        self.fluidResult = thresholdModel(resultdata,'Domain_ID', 0.5, 1.5, extract=False, cell=True)

        return

    def updateSolidResults(self):
        self.solidResult = read_data(self.resultDir+'/result_'+str(self.resultNum)+'.vtu', file_format="vtu")
        return

    def updateFluidResults(self):
        self.fluidResult = read_data(self.resultDir+'/result_'+str(self.resultNum)+'.vtu', file_format="vtu")
        return


    def appendIterfaceResult(self):
        """
        Add results of fluid simulation
        """
        numCells = self.vesselReference.GetNumberOfCells()

        fluidSurface = self.fluidResult.extract_surface()
        pointLocatorFluid = vtk.vtkPointLocator()
        pointLocatorFluid.SetDataSet(fluidSurface)
        pointLocatorFluid.BuildLocator()

        innerReference = thresholdModel(self.vesselReference.extract_surface(), 'InnerRegionID',0.5,1.5)
        numPointsInner = innerReference.GetNumberOfPoints()

        arrayNames = innerReference.array_names
        for name in arrayNames:
            if name not in ["displacements", "fluidStressQueryID", "WSS", "Pressure"]:
                if name in innerReference.point_data:
                    innerReference.point_data.remove(name)
                if name in self.vesselSolid.cell_data:
                    innerReference.cell_data.remove(name)

        for q in range(numPointsInner):
            fluidCoordinate = np.array(innerReference.GetPoint(q)) + np.array(innerReference.GetPointData().GetArray("displacements").GetTuple3(q))
            fluidId = int(pointLocatorFluid.FindClosestPoint(fluidCoordinate))
            innerReference.GetPointData().GetArray('WSS').SetTuple1(q,np.linalg.norm(fluidSurface.GetPointData().GetArray('WSS').GetTuple3(fluidId)))
            innerReference.GetPointData().GetArray('Pressure').SetTuple1(q,fluidSurface.GetPointData().GetArray('Pressure').GetTuple1(fluidId))

        if self.smoothAttributesValue:
            innerReference = smoothAttributes(innerReference, self.smoothAttributesValue, 100)

        pointLocatorInner = vtk.vtkPointLocator()
        pointLocatorInner.SetDataSet(innerReference)
        pointLocatorInner.BuildLocator()

        for q in range(numCells):
            fluidStressId = int(self.vesselReference.GetCellData().GetArray('fluidStressQueryID').GetTuple1(q))
            cell = self.vesselReference.GetCell(fluidStressId)
            cellPts = cell.GetPointIds()

            # For each cell, use trapezoidal integration to compute WSS
            cellWSS = 0.0
            numberOfPoints = 0
            for p in range(cellPts.GetNumberOfIds()):
                ptId = cellPts.GetId(p)
                if self.vesselReference.GetPointData().GetArray('InnerRegionID').GetTuple1(ptId) == 1:
                    innerCoordinate = np.array(self.vesselReference.GetPoint(ptId))
                    innerId = int(pointLocatorInner.FindClosestPoint(innerCoordinate))
                    self.vesselReference.GetPointData().GetArray('Pressure').SetTuple1(ptId,innerReference.GetPointData().GetArray('Pressure').GetTuple1(innerId))
                    cellWSS += innerReference.GetPointData().GetArray('WSS').GetTuple1(innerId)
                    numberOfPoints += 1
            #Average WSS in cell
            cellWSS *= 1/float(numberOfPoints)

            wss_prev = self.vesselReference.GetCellData().GetArray('wss_curr').GetTuple1(q)
            self.vesselReference.GetCellData().GetArray('wss_prev').SetTuple1(q,wss_prev)
            self.vesselReference.GetCellData().GetArray('wss_curr').SetTuple1(q,cellWSS)

        return


    def estimateIterfaceResult(self):
        """
        Estimates the wall shear stress and assigns pressure without
        fluid simulation result
        """

        # Calculate slice centers
        centers = np.zeros((self.numLen+1,3))
        inner = thresholdModel(self.vesselReference.extract_surface(), 'InnerRegionID',0.5,1.5)
        numInnerPts = inner.GetNumberOfPoints()
        for q in range(numInnerPts):
            pointCoordinate = np.array(inner.GetPoint(q)) + np.array(inner.GetPointData().GetArray("displacements").GetTuple3(q))
            pointSlice = int(inner.GetPointData().GetArray('sliceIds').GetTuple1(q))
            centers[pointSlice] = centers[pointSlice] + pointCoordinate/(self.numCirc)


        numCells = self.vesselReference.GetNumberOfCells()

        for q in range(numCells):
            fluidStressId = int(self.vesselReference.GetCellData().GetArray('fluidStressQueryID').GetTuple1(q))
            cell = self.vesselReference.GetCell(fluidStressId)
            cellPts = cell.GetPointIds()

            # For each cell, use trapezoidal integration to compute WSS
            radius = 0.0
            numberOfPoints = 0
            for p in range(cellPts.GetNumberOfIds()):
                ptId = cellPts.GetId(p)
                if self.vesselReference.GetPointData().GetArray('InnerRegionID').GetTuple1(ptId) == 1:
                    fluidCenter = centers[int(self.vesselReference.GetPointData().GetArray('sliceIds').GetTuple1(ptId))]
                    fluidCoordinate = np.array(self.vesselReference.GetPoint(ptId)) + np.array(self.vesselReference.GetPointData().GetArray("displacements").GetTuple3(ptId))
                    pointRadius = np.linalg.norm(fluidCoordinate - fluidCenter)
                    self.vesselReference.GetPointData().GetArray('Pressure').SetTuple1(ptId,self.outletPressure)
                    radius += pointRadius
                    numberOfPoints += 1

            #Average WSS in cell
            radius *= 1/float(numberOfPoints)
            cellWSS = 4.0*self.viscosity*abs(self.inletFlow)/(np.pi*(radius**3.0))


            wss_prev = self.vesselReference.GetCellData().GetArray('wss_curr').GetTuple1(q)
            self.vesselReference.GetCellData().GetArray('wss_prev').SetTuple1(q,wss_prev)
            self.vesselReference.GetCellData().GetArray('wss_curr').SetTuple1(q,cellWSS)

        return

    def appendReducedResult(self):

        numPts = self.vesselReference.GetNumberOfPoints()
        numCells = self.vesselReference.GetNumberOfCells()

        self.omega = self.omega/2.0

        # Calculate cauchy green tensor
        for q in range(numPts):
            rcurr = np.array(self.vesselReference.GetPointData().GetArray("residual_curr").GetTuple3(q))
            dcurr = np.array(self.vesselReference.GetPointData().GetArray("displacements").GetTuple3(q))
            displacement = dcurr - self.omega*rcurr

            self.vesselReference.GetPointData().GetArray("displacements").SetTuple(q, displacement)

        self.vesselReference = computeGaussValues(self.vesselReference,"displacements")

        return

    def appendSolidResult(self):

        if self.smoothAttributesValue:
            self.solidResult = smoothAttributes(self.solidResult,self.smoothAttributesValue,100)

        pointLocatorSolid = vtk.vtkPointLocator()
        pointLocatorSolid.SetDataSet(self.solidResult)
        pointLocatorSolid.BuildLocator()

        numPts = self.vesselReference.GetNumberOfPoints()
        numCells = self.vesselReference.GetNumberOfCells()

        time1 = time.time()

        dcurr = []
        dprev = []

        for q in range(numPts):
            originalCoordinate = np.array(self.vesselReference.GetPoint(q))
            displacement_prev = np.array(self.vesselReference.GetPointData().GetArray("displacements").GetTuple3(q))
            currentCoordinate = originalCoordinate + displacement_prev

            pointIdSolid = pointLocatorSolid.FindClosestPoint(currentCoordinate)
            solidCoordinate = np.array(self.solidResult.GetPoint(pointIdSolid))
            solidDispacement = np.array(self.solidResult.GetPointData().GetArray("Displacement").GetTuple3(pointIdSolid))

            displacement_curr = solidCoordinate + solidDispacement - originalCoordinate
            residual_curr = displacement_curr - displacement_prev
            self.vesselReference.GetPointData().GetArray("residual_curr").SetTuple(q, residual_curr)

            dcurr.append(displacement_curr)
            dprev.append(displacement_prev)


        rcurr = np.array(self.vesselReference.GetPointData().GetArray("residual_curr")).flatten()
        rprev = np.array(self.vesselReference.GetPointData().GetArray("residual_prev")).flatten()

        dcurr = np.array(dcurr).flatten()
        dprev = np.array(dprev).flatten()



        if self.predictMethod == "none":
            dnew = dcurr.reshape((-1, 3))


        elif self.predictMethod == "aitken":
            if self.timeIter > 1:
                diff = rcurr - rprev
                self.omega = -self.omega*np.dot(rprev,diff)/np.dot(diff,diff)
                #if self.skippedFluid == False:
                #    self.omega = 0.25
            elif self.timeStep == 0 and self.timeIter == 0:
                self.omega = 1.0
            else:
                self.omega = 0.5

            if self.omega > 2.0:
                self.omega = 2.0
            elif self.omega < 0.1:
                self.omega = 0.1

            vnew =  dprev + self.omega*rcurr
            dnew = vnew.reshape((-1,3))


        elif self.predictMethod == "iqnils":

            if self.timeIter == 0:
                self.mat_W = []
                self.mat_V = []
                self.mat_D = []
                self.mat_D.append(dcurr)
                vnew =  dprev + 0.5*rcurr
                dnew = vnew.reshape((-1,3))
            elif self.timeIter == 1:
                self.mat_D.append(dcurr)
                vnew =  dprev + 0.5*rcurr
                dnew = vnew.reshape((-1,3))
            else:
                self.mat_D.append(dcurr)
                self.mat_W.append(self.mat_D[-1] - self.mat_D[-2])
                self.mat_V.append(rcurr - rprev)

                # trim to max number of considered vectors
                self.mat_V = self.mat_V[-self.nq:]
                self.mat_W = self.mat_W[-self.nq:]

                # remove linearly dependent vectors
                while True:
                    # QR decomposition
                    qq, rr = np.linalg.qr(np.array(self.mat_V[:self.nq]).T)

                    # tolerance for redundant vectors
                    i_eps = np.where(
                        np.abs(np.diag(rr)) < self.iq_eps
                    )[0]
                    if not np.any(i_eps):
                        break

                    print("Filtering " + str(len(i_eps)) + " time steps")
                    for i in reversed(i_eps):
                        self.mat_V.pop(i)
                        self.mat_W.pop(i)

                # solve for coefficients
                bb = np.linalg.solve(rr.T, -np.dot(np.array(self.mat_V), rcurr))
                cc = np.linalg.solve(rr, bb)

                # update
                vnew = dcurr + np.dot(np.array(self.mat_W).T, cc)
                dnew = vnew.reshape((-1, 3))


        # Calculate cauchy green tensor
        for q in range(numPts):
            rcurr = np.array(self.vesselReference.GetPointData().GetArray("residual_curr").GetTuple3(q))
            self.vesselReference.GetPointData().GetArray("residual_prev").SetTuple(q, rcurr)
            self.vesselReference.GetPointData().GetArray("displacements").SetTuple(q, dnew[q])

        self.vesselReference = computeGaussValues(self.vesselReference,"displacements")

        # Get stress invariant
        for q in range(numCells):
            cell = self.solidResult.GetCell(q)
            cellPts = cell.GetPointIds()

            # Elementwise values
            fluidStressId = int(self.vesselReference.GetCellData().GetArray('fluidStressQueryID').GetTuple1(q))
            solidStressId = int(self.vesselReference.GetCellData().GetArray('solidStressQueryID').GetTuple1(q))

            sigma_inv = 0.0

            # For each cell, use trapezoidal integration to compute sigma. Choose to average radially or not.
            if self.averageStress:
                stressRange = range(fluidStressId, solidStressId)
            else:
                stressRange = range(q, q+1)

            for p in stressRange:
                cell = self.solidResult.GetCell(p)
                cellPts = cell.GetPointIds()
                cellSigma = 0.0
                for r in range(cellPts.GetNumberOfIds()):
                    ptId = cellPts.GetId(r)
                    pointSigma = self.solidResult.GetPointData().GetArray('Cauchy_stress').GetTuple6(ptId)           
                    cellSigma += (pointSigma[0]+pointSigma[1]+pointSigma[2])
                cellSigma *= 1/float(cellPts.GetNumberOfIds())
                sigma_inv = sigma_inv + cellSigma

            sigma_inv = 0.1*sigma_inv
            if self.averageStress:
                sigma_inv = sigma_inv/(solidStressId-fluidStressId)

            inv_prev = self.vesselReference.GetCellData().GetArray("inv_curr").GetTuple1(q)
            self.vesselReference.GetCellData().GetArray("inv_prev").SetTuple1(q, inv_prev)
            self.vesselReference.GetCellData().GetArray("inv_curr").SetTuple1(q, sigma_inv)

        return


    def appendFluidResult(self):

        numPoints = self.vesselFluid.GetNumberOfPoints()

        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(self.fluidResult)
        pointLocator.BuildLocator()

        velocity_array = np.zeros((numPoints,3))
        pressure_array = np.zeros((numPoints,1))

        for q in range(numPoints):
            coordinate = self.vesselFluid.GetPoint(q)
            pointIdSource = pointLocator.FindClosestPoint(coordinate)

            velocity_array[q,:] = self.fluidResult.GetPointData().GetArray('Velocity').GetTuple(pointIdSource)
            pressure_array[q,:] = self.fluidResult.GetPointData().GetArray('Pressure').GetTuple1(pointIdSource)

        self.vesselFluid.GetPointData().AddArray(pv.convert_array(np.array(velocity_array),name="Velocity"))
        self.vesselFluid.GetPointData().AddArray(pv.convert_array(np.array(self.outletPressure),name="Pressure"))

        return


    def updateReference(self):

        numCells = self.vesselReference.GetNumberOfCells()
        numPoints = self.vesselReference.GetNumberOfPoints()

        array_names = self.vesselReference.array_names

        dv_max = 0

        if "temp_array" in array_names:
            self.vesselReference.rename_array("varWallProps","material_global")
            self.vesselReference.rename_array("temp_array","varWallProps")

        for q in range(numCells):
            e_r = self.vesselReference.GetCellData().GetArray('e_r').GetTuple(q)
            e_t = self.vesselReference.GetCellData().GetArray('e_t').GetTuple(q)
            e_z = self.vesselReference.GetCellData().GetArray('e_z').GetTuple(q)
            Q = np.array((e_r,e_t,e_z))

            J_curr = self.vesselReference.GetCellData().GetArray('J_curr').GetTuple(q)
            inv_curr = self.vesselReference.GetCellData().GetArray('inv_curr').GetTuple1(q)

            sigma_gnr_g = []
            stiffness_g = []
            varWallProps_g = []
            p_est_g = []

            J_c = 0.0
            p_est_c = 0.0
            for p in range(self.nG):
                J_c += J_curr[p]/self.nG

                sigma_gnr_mem = self.vesselReference.GetCellData().GetArray('sigma_mem').GetTuple(q)[p*9:(p+1)*9]
                sigma_gnr = rotateStress(sigma_gnr_mem,Q)
                sigma_gnr_g = sigma_gnr_g + sigma_gnr.tolist()

                p_est = (10.0*inv_curr-(sigma_gnr[0]+sigma_gnr[1]+sigma_gnr[2]))/3.0
                p_est_g = p_est_g + [p_est]
                p_est_c += p_est/self.nG


            for p in range(self.nG):
                J_target = self.vesselReference.GetCellData().GetArray('J_target').GetTuple(q)[p]
                stiffness_mem = self.vesselReference.GetCellData().GetArray('stiffness_mem').GetTuple(q)[p*36:(p+1)*36]
                stiffness = rotateStiffness(stiffness_mem,Q)
                stiffness_g = stiffness_g + stiffness.tolist()

                if self.averageVolume:
                    J_final = J_target/J_c
                else:
                    J_final = J_target/J_curr[p]

                varWallProps_g = np.hstack((varWallProps_g,np.hstack((stiffness,p_est_g[p],J_final,sigma_gnr_g[p*6:(p+1)*6],p))))

                if np.abs(J_final - 1) > dv_max:
                    dv_max = np.abs(J_final - 1)

            self.vesselReference.GetCellData().GetArray('varWallProps').SetTuple(q,varWallProps_g)
            self.vesselReference.GetCellData().GetArray('p_est').SetTuple(q,p_est_g)
            self.vesselReference.GetCellData().GetArray('sigma_gnr').SetTuple(q,sigma_gnr_g)

            self.vesselReference.GetCellData().GetArray('J_curr').SetTuple(q,J_curr)
            self.vesselReference.GetCellData().GetArray('J_c').SetTuple1(q,J_c)
            self.vesselReference.GetCellData().GetArray('p_est_c').SetTuple1(q,p_est_c)

        """
        print("Editing data...")

        if dv_max > 0.05:
            print("Scaling dV_max: " + str(dv_max))
            for q in range(numCells):
                varWallProp = np.array(self.vesselReference.GetCellData().GetArray('varWallProps').GetTuple(q))
                for i in range(37,360,45):
                    varWallProp[i] = 1 + (varWallProp[i]-1)/(dv_max/0.05)
                self.vesselReference.GetCellData().GetArray('varWallProps').SetTuple(q,varWallProp)


        print("Finished editing data...")
        """

        return



    def initializeCylinder(self):
        print("Initializing cylindrical vessel...")

        with open('FolderVesselConfigurationFiles/Native_in_handshake_') as f:
            lines = f.readlines()[1:]
            nativeIn = []
            for line in lines:
                nativeIn.append([float(x) for x in line.split()])

        self.sigma_h = nativeIn[13][0]
        self.tau_h = nativeIn[13][1]

        #Build material array (constant for now)
        materialArray = [ nativeIn[7][0],nativeIn[4][0],nativeIn[4][1],nativeIn[4][2],nativeIn[2][0]*10.0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[5][2],nativeIn[2][4]*10.0,nativeIn[2][5],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][3],nativeIn[5][3],nativeIn[2][6]*10.0,nativeIn[2][7],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][4],nativeIn[5][4],nativeIn[2][8]*10.0,nativeIn[2][9],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][5],nativeIn[5][5],nativeIn[2][10]*10.0,nativeIn[2][11],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][1],nativeIn[5][1],nativeIn[2][2]*10.0,nativeIn[2][3],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[16][0],nativeIn[14][2],nativeIn[14][1],0,0,0,0,0,0,0,0,0,0]

        #********************
        ang1 = nativeIn[3][4]
        ang2 = nativeIn[3][5]

        points = np.empty([(self.numCirc+1)*(self.numLen+1)*(self.numRad+1),3])

        structureIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        innerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        outerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        proximalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        distalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        sliceIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))

        proximalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
        distalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1
        innerIds[0:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1
        outerIds[self.numRad:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1

        fluidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)
        solidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)

        tevgValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        aneurysmValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        simulateIds = np.ones(self.numCirc*self.numLen*self.numRad)

        e_r = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_t = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_z = np.zeros((self.numCirc*self.numLen*self.numRad,3))

        e_ma = np.zeros((self.numCirc*self.numLen*self.numRad,98))

        aFac = 1.5
        
        for i in range(self.numLen+1):
            for j in range(self.numCirc+1):
                for k in range(self.numRad+1):
                    xPt = (self.radius + self.thickness*k/self.numRad)*np.cos(2.0*j*np.pi/self.numCirc)
                    yPt = (self.radius + self.thickness*k/self.numRad)*np.sin(2.0*j*np.pi/self.numCirc)
                    if self.adaptiveMesh:
                        if self.length*i/self.numLen - self.length/2.0 == 0:
                            zPt = 0.0
                        else:
                            zPt = -(1 - 2.0*i/(self.numLen))/abs(1 - 2.0*i/(self.numLen))*self.length*0.5*(np.exp(aFac*abs(1 - 2.0*i/(self.numLen)))-1.0)/(np.exp(aFac)-1.0)
                    else:
                        zPt = self.length*i/self.numLen - self.length/2.0
                    structureIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] = i*(self.numCirc)*(self.numRad+1) + (j%self.numCirc)*(self.numRad+1) + k + 1
                    points[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k,:] = [xPt, yPt, zPt]
                    sliceIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] =  i


        for i in range(self.numLen):
            for j in range(self.numCirc):
                for k in range(self.numRad):

                    theta = 2.0*(j+0.5)*np.pi/self.numCirc

                    xPt = np.cos(theta)
                    yPt = np.sin(theta)
                    if self.adaptiveMesh:
                        if self.length*(i+0.5)/self.numLen - self.length/2.0 == 0:
                            zPt = 0.0
                        else:
                            zPt = -(1 - 2.0*(i+0.5)/(self.numLen))/abs(1 - 2.0*(i+0.5)/(self.numLen))*self.length*0.5*(np.exp(aFac*abs(1 - 2.0*(i+0.5)/(self.numLen)))-1.0)/(np.exp(aFac)-1.0)
                    else:
                        zPt = self.length*(i+0.5)/self.numLen - self.length/2.0

                    if self.adaptiveMesh:
                        if self.length*(i+0.5)/self.numLen - self.length/2.0 == 0:
                            zPt = 0.0
                        else:
                            zPt = -(1 - 2.0*(i+0.5)/(self.numLen))/abs(1 - 2.0*(i+0.5)/(self.numLen))*self.length*0.5*(np.exp(aFac*abs(1 - 2.0*(i+0.5)/(self.numLen)))-1.0)/(np.exp(aFac)-1.0)
                    else:
                        zPt = self.length*(i+0.5)/self.numLen - self.length/2.0

                    v_r = np.array([xPt,yPt,0])
                    v_t = np.array([-yPt,xPt,0])
                    v_z = np.array([0,0,1])

                    e_r[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_r
                    e_t[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_t
                    e_z[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_z

                    A = np.array([v_r,v_t,v_z])
                    # Is this a rotation matrix?
                    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                                          np.eye(3, dtype=float)) > 1e-6):
                        print(A)
                        raise RuntimeError('Matrix *A* does not describe a rotation.')

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,:] = materialArray

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,5:8]   = v_r
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,8:11]  = v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,11:14] = v_z

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,18:21] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,32:35] =  np.cos(0.0*np.pi/180.0)*v_z + np.sin(0.0*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,46:49] =  np.cos(ang1*np.pi/180.0)*v_z + np.sin(ang1*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,60:63] =  np.cos(ang2*np.pi/180.0)*v_z + np.sin(ang2*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,74:77] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,88:91] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t

                    fluidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc
                    solidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc + self.numRad

                    if self.tevg:
                        tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getTEVGValue([xPt,yPt,zPt], self.radius)
                        if tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] > 0:
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,4] = 200000000
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,1] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,2] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,3] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,14] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,28] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,42] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,56] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,70] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,84] = 0
                    if self.aneurysm:
                        aneurysmValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getAneurysmValue(zPt, self.zod, self.zapex, theta, self.thetaod, self.thetaapex)

        vol = pv.StructuredGrid()
        vol.points = points
        vol.dimensions = [self.numRad+1, self.numCirc+1, self.numLen+1]

        numCells = vol.GetNumberOfCells()
        vol.GetPointData().AddArray(pv.convert_array((proximalIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((distalIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((innerIds).astype(int),name="InnerRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outerIds).astype(int),name="OuterRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((structureIds).astype(int),name="StructureID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="GlobalElementID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="CellStructureID"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="varWallProps"))
        vol.GetPointData().AddArray(pv.convert_array((sliceIds).astype(int),name="sliceIds"))

        #Exchange these to gauss point quantities
        #Matrices/arrays
        vol.GetCellData().AddArray(pv.convert_array(np.zeros((self.numCirc*self.numLen*self.numRad,45*8)).astype(float),name="temp_array"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad_mem"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(36),8),(numCells,1)).astype(float),name="stiffness_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(9),8),(numCells,1)).astype(float),name="sigma_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(6),8),(numCells,1)).astype(float),name="sigma_gnr"))
        #Scalar values
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_curr"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_target"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([0],8),(numCells,1)).astype(float),name="p_est"))

        vol.GetCellData().AddArray(pv.convert_array(np.ones(numCells).astype(float),name="J_c"))
        vol.GetCellData().AddArray(pv.convert_array(np.zeros(numCells).astype(float),name="p_est_c"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_prev"))

        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_prev"))


        vol.GetCellData().AddArray(pv.convert_array(e_r.astype(float),name="e_r"))
        vol.GetCellData().AddArray(pv.convert_array(e_t.astype(float),name="e_t"))
        vol.GetCellData().AddArray(pv.convert_array(e_z.astype(float),name="e_z"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="material_global"))
        vol.GetCellData().AddArray(pv.convert_array(fluidStressQueryIds.astype(int),name="fluidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(solidStressQueryIds.astype(int),name="solidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(simulateIds.astype(int),name="Simulate"))
        vol.GetCellData().AddArray(pv.convert_array(aneurysmValue.astype(float),name="aneurysmValue"))
        vol.GetCellData().AddArray(pv.convert_array(tevgValue.astype(float),name="tevgValue"))
        vol = vol.cast_to_unstructured_grid()
        vol = clean(vol, tolerance=1e-10)
        numPts = vol.GetNumberOfPoints()
        vol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPts,numPts).astype(int),name="GlobalNodeID"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="displacements"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_curr"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_prev"))
        vol.GetPointData().AddArray(pv.convert_array(np.zeros(numPts).astype(float),name="Pressure"))
        vol.GetPointData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numPts).astype(float),name="WSS"))
        vol.GetPointData().AddArray(pv.convert_array(vol.points.astype(float),name="Coordinate"))

        return vol


    def initializeTorus(self):
        print("Initializing toric vessel...")

        with open('FolderVesselConfigurationFiles/Native_in_handshake_') as f:
            lines = f.readlines()[1:]
            nativeIn = []
            for line in lines:
                nativeIn.append([float(x) for x in line.split()])

        self.sigma_h = nativeIn[13][0]
        self.tau_h = nativeIn[13][1]

        #Build material array (constant for now)
        materialArray = [ nativeIn[7][0],nativeIn[4][0],nativeIn[4][1],nativeIn[4][2],nativeIn[2][0]*10.0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[5][2],nativeIn[2][4]*10.0,nativeIn[2][5],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][3],nativeIn[5][3],nativeIn[2][6]*10.0,nativeIn[2][7],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][4],nativeIn[5][4],nativeIn[2][8]*10.0,nativeIn[2][9],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][5],nativeIn[5][5],nativeIn[2][10]*10.0,nativeIn[2][11],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][1],nativeIn[5][1],nativeIn[2][2]*10.0,nativeIn[2][3],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[16][0],nativeIn[14][2],nativeIn[14][1],0,0,0,0,0,0,0,0,0,0]

        #********************

        passive_stress = np.array([0, 0, 0])
        for i in range(5):
            angle = nativeIn[3][i+1]
            lamda_prestretch = nativeIn[5][i+1]
            F_alpha_ntau_s = lamda_prestretch * np.array([0, np.sin(angle*np.pi/180.0), np.cos(angle*np.pi/180.0)])
            sigma_alpha_ntau_s = hat_S_alpha(lamda_prestretch, nativeIn[2][2 + i*2]*10.0, nativeIn[2][3 + i*2])*np.power(F_alpha_ntau_s,2)
            passive_stress = passive_stress + nativeIn[7][i+1]*sigma_alpha_ntau_s
        g_e = nativeIn[2][0]*10.0
        eps_e = nativeIn[7][0]

        points = np.empty([(self.numCirc+1)*(self.numLen+1)*(self.numRad+1),3])

        structureIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        innerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        outerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        proximalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        distalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        sliceIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        
        if self.flipInlet:
            distalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
            proximalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1
        else:
            proximalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
            distalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1

        innerIds[0:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1
        outerIds[self.numRad:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1

        fluidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)
        solidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)

        tevgValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        aneurysmValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        simulateIds = np.ones(self.numCirc*self.numLen*self.numRad)

        e_r = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_t = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_z = np.zeros((self.numCirc*self.numLen*self.numRad,3))

        e_ma = np.zeros((self.numCirc*self.numLen*self.numRad,98))


        for i in range(self.numLen+1):
            theta = self.torusFraction * 2 * np.pi * i / self.numLen  # Azimuthal angle (circumference)
            for j in range(self.numCirc+1):
                phi = -2 * np.pi * j / self.numCirc
                for k in range(self.numRad+1):

                    xPt = self.radius * np.sin(phi) + (np.sin(phi)*(self.thickness*k/self.numRad))
                    yPt = (self.torusRadius + self.radius * np.cos(phi)) * np.cos(theta) + ((np.cos(phi) * np.cos(theta))*(self.thickness*k/self.numRad))
                    zPt = (self.torusRadius + self.radius * np.cos(phi)) * np.sin(theta) + ((np.cos(phi) * np.sin(theta))*(self.thickness*k/self.numRad))

                    structureIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] = i*(self.numCirc)*(self.numRad+1) + (j%self.numCirc)*(self.numRad+1) + k + 1
                    points[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k,:] = [xPt, yPt, zPt]
                    sliceIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] =  i

        #lambda_r = nativeIn[4][0]
        #lambda_t = nativeIn[4][1]
        #lambda_z = nativeIn[4][2]
        #outletPressure = self.outletPressure
        #sigma_z = self.outletPressure*self.radius/(2*self.thickness)
        for i in range(self.numLen):
            theta = self.torusFraction * 2 * np.pi * (i + 0.5) / self.numLen  # Azimuthal angle (circumference)
            for j in range(self.numCirc):
                phi = -2 * np.pi * (j + 0.5) / self.numCirc
                #sigma_t_inner = (self.outletPressure*self.radius/self.thickness)*((self.torusRadius - self.radius*np.sin(phi + np.pi/2)/2)/(self.torusRadius - self.radius*np.sin(phi + np.pi/2)))
                #sigma_t_outer = (self.outletPressure*self.radius/self.thickness)*((self.torusRadius + self.radius*np.sin(phi + np.pi/2)/2)/(self.torusRadius + self.radius*np.sin(phi + np.pi/2)))
                for k in range(self.numRad):

                    #frac = (k+0.5)/self.numRad

                    xPt = self.radius * np.sin(phi) + (np.sin(phi)*(self.radius + self.thickness*(k + 0.5)/self.numRad))
                    yPt = (self.torusRadius + self.radius * np.cos(phi)) * np.cos(theta) + ((np.cos(phi) * np.cos(theta))*(self.radius + self.thickness*(k + 0.5)/self.numRad))
                    zPt = (self.torusRadius + self.radius * np.cos(phi)) * np.sin(theta) + ((np.cos(phi) * np.sin(theta))*(self.radius + self.thickness*(k + 0.5)/self.numRad))

                    v_r = np.array([np.sin(phi),np.cos(phi) * np.cos(theta),np.cos(phi) * np.sin(theta)])
                    v_z = np.array([0,-np.sin(theta),np.cos(theta)])
                    v_t = np.cross(v_r, v_z)

                    e_r[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_r
                    e_t[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_t
                    e_z[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_z

                    A = np.array([v_r,v_t,v_z])
                    # Is this a rotation matrix?
                    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                                          np.eye(3, dtype=float)) > 1e-6):
                        print(A)
                        raise RuntimeError('Matrix *A* does not describe a rotation.')

                    #sigma_t = (1-frac)*sigma_t_inner + frac*(sigma_t_outer)
                    #sigma_r = -(1-frac)*self.outletPressure

                    #target = [sigma_r, sigma_t, sigma_z]

                    #optimal_lambda_values = minimize_difference(passive_stress, g_e, eps_e, frac, self.outletPressure, target, [nativeIn[4][1],nativeIn[4][2]])

                    #print(optimal_lambda_values)

                    #materialArray[1:4] = [1/(optimal_lambda_values[0]*optimal_lambda_values[1]), optimal_lambda_values[0], optimal_lambda_values[1]]

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,:] = materialArray

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,5:8]   = v_r
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,8:11]  = v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,11:14] = v_z

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,18:21] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,32:35] =  np.cos(0.0*np.pi/180.0)*v_z + np.sin(0.0*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,46:49] =  np.cos(41.94*np.pi/180.0)*v_z + np.sin(41.94*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,60:63] =  np.cos(318.06*np.pi/180.0)*v_z + np.sin(318.06*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,74:77] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,88:91] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t

                    fluidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc
                    solidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc + self.numRad

                    if self.tevg:
                        tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getTEVGValue([xPt,yPt,zPt], self.radius)
                        if tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] > 0:
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,4] = 200000000
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,1] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,2] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,3] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,14] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,28] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,42] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,56] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,70] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,84] = 0
                    if self.aneurysm:
                        aneurysmValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getAneurysmValue([xPt,yPt,zPt], self.radius)

        vol = pv.StructuredGrid()
        vol.points = points
        vol.dimensions = [self.numRad+1, self.numCirc+1, self.numLen+1]

        numCells = vol.GetNumberOfCells()
        vol.GetPointData().AddArray(pv.convert_array((proximalIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((distalIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((innerIds).astype(int),name="InnerRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outerIds).astype(int),name="OuterRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((structureIds).astype(int),name="StructureID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="GlobalElementID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="CellStructureID"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="varWallProps"))
        vol.GetPointData().AddArray(pv.convert_array((sliceIds).astype(int),name="sliceIds"))

        #Exchange these to gauss point quantities
        #Matrices/arrays
        vol.GetCellData().AddArray(pv.convert_array(np.zeros((self.numCirc*self.numLen*self.numRad,45*8)).astype(float),name="temp_array"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad_mem"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(36),8),(numCells,1)).astype(float),name="stiffness_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(9),8),(numCells,1)).astype(float),name="sigma_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(6),8),(numCells,1)).astype(float),name="sigma_gnr"))
        #Scalar values
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_curr"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_target"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([0],8),(numCells,1)).astype(float),name="p_est"))

        vol.GetCellData().AddArray(pv.convert_array(np.ones(numCells).astype(float),name="J_c"))
        vol.GetCellData().AddArray(pv.convert_array(np.zeros(numCells).astype(float),name="p_est_c"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_prev"))

        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_prev"))


        vol.GetCellData().AddArray(pv.convert_array(e_r.astype(float),name="e_r"))
        vol.GetCellData().AddArray(pv.convert_array(e_t.astype(float),name="e_t"))
        vol.GetCellData().AddArray(pv.convert_array(e_z.astype(float),name="e_z"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="material_global"))
        vol.GetCellData().AddArray(pv.convert_array(fluidStressQueryIds.astype(int),name="fluidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(solidStressQueryIds.astype(int),name="solidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(simulateIds.astype(int),name="Simulate"))
        vol.GetCellData().AddArray(pv.convert_array(aneurysmValue.astype(float),name="aneurysmValue"))
        vol.GetCellData().AddArray(pv.convert_array(tevgValue.astype(float),name="tevgValue"))
        vol = vol.cast_to_unstructured_grid()
        vol = clean(vol, tolerance=1e-10)
        numPts = vol.GetNumberOfPoints()
        vol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPts,numPts).astype(int),name="GlobalNodeID"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="displacements"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_curr"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_prev"))
        vol.GetPointData().AddArray(pv.convert_array(np.zeros(numPts).astype(float),name="Pressure"))
        vol.GetPointData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numPts).astype(float),name="WSS"))
        vol.GetPointData().AddArray(pv.convert_array(vol.points.astype(float),name="Coordinate"))

        return vol



    def initializeAnastomosis(self):
        print("Initializing anastomosis vessel...")

        with open('FolderVesselConfigurationFiles/Native_in_handshake_') as f:
            lines = f.readlines()[1:]
            nativeIn = []
            for line in lines:
                nativeIn.append([float(x) for x in line.split()])

        self.sigma_h = nativeIn[13][0]
        self.tau_h = nativeIn[13][1]

        #Build material array (constant for now)
        materialArray = [ nativeIn[7][0],nativeIn[4][0],nativeIn[4][1],nativeIn[4][2],nativeIn[2][0]*10.0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[5][2],nativeIn[2][4]*10.0,nativeIn[2][5],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][3],nativeIn[5][3],nativeIn[2][6]*10.0,nativeIn[2][7],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][4],nativeIn[5][4],nativeIn[2][8]*10.0,nativeIn[2][9],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][5],nativeIn[5][5],nativeIn[2][10]*10.0,nativeIn[2][11],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][1],nativeIn[5][1],nativeIn[2][2]*10.0,nativeIn[2][3],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[16][0],nativeIn[14][2],nativeIn[14][1],0,0,0,0,0,0,0,0,0,0]

        #********************
        ang1 = nativeIn[3][4]
        ang2 = nativeIn[3][5]

        points = np.empty([(self.numCirc+1)*(self.numLen+1)*(self.numRad+1),3])

        structureIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        innerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        outerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        proximalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        distalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        sliceIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))

        proximalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
        distalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1
        innerIds[0:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1
        outerIds[self.numRad:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1

        fluidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)
        solidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)

        tevgValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        aneurysmValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        simulateIds = np.ones(self.numCirc*self.numLen*self.numRad)

        e_r = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_t = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_z = np.zeros((self.numCirc*self.numLen*self.numRad,3))

        e_ma = np.zeros((self.numCirc*self.numLen*self.numRad,98))

        aFac = 1.5
        
        l1 = self.length*(1-self.anastomosisLengthRatio)/2 - self.anastomosisTransitionLength
        l9 = l1

        l2 = 5*self.anastomosisTransitionLength/12
        l3 = 1*self.anastomosisTransitionLength/6
        l4 = 5*self.anastomosisTransitionLength/12

        l5 = self.length*self.anastomosisLengthRatio

        l6 = 5*self.anastomosisTransitionLength/12
        l7 = 1*self.anastomosisTransitionLength/6
        l8 = 5*self.anastomosisTransitionLength/12

        wl = self.radius*self.anastomosisConnectionRatio
        wr = self.radius*self.anastomosisConnectionRatio

        for i in range(self.numLen+1):
            for j in range(self.numCirc+1):
                for k in range(self.numRad+1):

                    if self.adaptiveMesh:
                        if self.length*i/self.numLen - self.length/2.0 == 0:
                            axi = 0.0
                        else:
                            axi = -(1 - 2.0*i/(self.numLen))/abs(1 - 2.0*i/(self.numLen))*self.length*0.5*(np.exp(aFac*abs(1 - 2.0*i/(self.numLen)))-1.0)/(np.exp(aFac)-1.0)
                    else:
                        axi = self.length*i/self.numLen - self.length/2.0

                    zs, rs, ths = sympy.symbols('zs rs ths')

                    # ***
                    if axi < -self.length/ 2 + l1:
                        r_inner = self.radius 
                        fp = rs

                    elif axi < -self.length/ 2 + l1 + l2:
                        r_inner = 0.5 * (self.radius - 0.5 * wl*2.0) *    np.sin((np.pi * (axi - (0.5*(l2 + l3)))) / l2) + 0.5 * (self.radius + 0.5 * wl*2.0)
                        fp = rs - 0.5 * (self.radius - 0.5 * wl*2.0) * sympy.sin((np.pi * (zs  - (0.5*(l2 + l3)))) / l2) + 0.5 * (self.radius + 0.5 * wl*2.0)

                    elif axi < -self.length/ 2 + l1 + l2 + l3:
                        r_inner = wl
                        fp = rs

                    elif axi < -self.length/ 2 + l1 + l2 + l3 + l4:
                        r_inner = 0.5 * (self.radius*self.anastomosisGraftRatio - 0.5 * wl*2.0) *    np.sin((np.pi * ( axi + (0.5*(l4 + l5)))) / l4) + 0.5 * (self.radius*self.anastomosisGraftRatio + 0.5 * wl*2.0)
                        fp = rs - 0.5 * (self.radius*self.anastomosisGraftRatio - 0.5 * wl*2.0) * sympy.sin((np.pi * (zs + (0.5 * (l4 + l5)))) / l4) + 0.5 * (self.radius*self.anastomosisGraftRatio + 0.5 * wl*2.0)

                    elif axi < -self.length/ 2 + l1 + l2 + l3 + l4 + l5:
                        r_inner = self.radius*self.anastomosisGraftRatio
                        fp = rs

                    elif axi < -self.length/ 2 + l1 + l2 + l3 + l4 + l5 + l6:
                        r_inner = 0.5 * (self.radius*self.anastomosisGraftRatio - 0.5 * wr*2.0) * np.sin((np.pi * -(axi - 0.5 * (l5 + l6))) / l6) + 0.5 * (self.radius*self.anastomosisGraftRatio + 0.5 * wr*2.0)
                        fp = rs - 0.5 * (self.radius*self.anastomosisGraftRatio - 0.5 * wr*2.0) * sympy.sin((np.pi * -(zs - 0.5 * (l5 + l6))) / l6) + 0.5 * (self.radius*self.anastomosisGraftRatio + 0.5 * wr*2.0)

                    elif axi < -self.length/ 2 + l1 + l2 + l3 + l4 + l5 + l6 + l7:
                        r_inner = wr
                        fp = rs

                    elif axi < -self.length/ 2 + l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8:
                        r_inner = 0.5 * (self.radius  - 0.5 * wr*2.0) * np.sin(   (np.pi * -(axi + 0.5 * (l7 + l8))) / l8) + 0.5 * (self.radius + 0.5 * wr*2.0)
                        fp = rs - 0.5 * (self.radius  - 0.5 * wr*2.0) * sympy.sin((np.pi * -(zs  + 0.5 * (l7 + l8))) / l8) + 0.5 * (self.radius + 0.5 * wr*2.0)

                    else:
                        r_inner = self.radius 
                        fp = rs


                    dfdz = sympy.diff(fp, zs)
                    dfdr = sympy.diff(fp, rs)

                    # cylindrical coordinate system
                    cir = 2.0*j*np.pi/self.numCirc
                    radial_coord = float(dfdr.subs([(rs, self.radius), (zs, axi)]))
                    # angular_coord = cir
                    z_coord = float(dfdz.subs([(rs, self.radius), (zs, axi)]))
                    mag = np.sqrt(radial_coord**2 + z_coord**2)
                    # vector_norm = (1/mag)*[radial_coord, angular_coord, z_coord]
                    # import pdb; pdb.set_trace()
                    rad = (
                        r_inner + self.thickness * k / self.numRad * mag/sympy.cos(sympy.atan(z_coord/radial_coord))
                    )
                    # print(mag/sympy.cos(sympy.atan(z_coord/radial_coord)))
                    # NOTE: can change axi here and around line 254 to shift coordinate system along z-axis

                    rad_normal = np.array([radial_coord * np.cos(cir), radial_coord * np.sin(cir), z_coord])
                    rad_normal = rad_normal/np.linalg.norm(rad_normal)


                    points[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k,:] = rad_normal*self.thickness*k/self.numRad + np.array([r_inner * np.cos(cir), r_inner * np.sin(cir), axi])

                    structureIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] = i*(self.numCirc)*(self.numRad+1) + (j%self.numCirc)*(self.numRad+1) + k + 1
                    sliceIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] =  i


        coords=[[0, 1, 0],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 0],
                [0, 0, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 0, 0]]


        for i in range(self.numLen):
            for j in range(self.numCirc):
                for k in range(self.numRad):

                    cellPts = []

                    for coord in coords:
                        cellPts.append(points[(i+coord[0])*(self.numCirc+1)*(self.numRad+1) + (j+coord[1])*(self.numRad+1) + (k+coord[2]),:])

                    cellPts = np.array(cellPts)

                    cellCenter = np.mean(cellPts,axis=0)

                    xPt = cellCenter[0]
                    yPt = cellCenter[1]
                    zPt = cellCenter[2]

                    v_r = np.mean(cellPts[[1,2,5,6]],axis=0) - np.mean(cellPts[[0,3,4,7]],axis=0)
                    v_r = v_r/np.linalg.norm(v_r)
                    v_t = np.mean(cellPts[[4,5,6,7]],axis=0) - np.mean(cellPts[[0,1,2,3]],axis=0)
                    v_t = v_t/np.linalg.norm(v_t)
                    v_z = np.mean(cellPts[[2,3,6,7]],axis=0) - np.mean(cellPts[[0,1,4,5]],axis=0)
                    v_z = v_z/np.linalg.norm(v_z)

                    v_r = np.cross(v_z,v_t)
                    v_r = v_r/np.linalg.norm(v_r)
                    #v_t = np.cross(v_z,v_r)
                    #v_t = v_t/np.linalg.norm(v_t)
                    v_z = np.cross(v_t,v_r)
                    v_z = v_z/np.linalg.norm(v_z)

                    e_r[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_r
                    e_t[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_t
                    e_z[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_z

                    A = np.array([v_r,v_t,v_z])
                    # Is this a rotation matrix?
                    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                                          np.eye(3, dtype=float)) > 1e-6):
                        print(A)
                        raise RuntimeError('Matrix *A* does not describe a rotation.')

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,:] = materialArray

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,5:8]   = v_r
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,8:11]  = v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,11:14] = v_z

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,18:21] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,32:35] =  np.cos(0.0*np.pi/180.0)*v_z + np.sin(0.0*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,46:49] =  np.cos(ang1*np.pi/180.0)*v_z + np.sin(ang1*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,60:63] =  np.cos(ang2*np.pi/180.0)*v_z + np.sin(ang2*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,74:77] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,88:91] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t

                    fluidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc
                    solidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc + self.numRad

                    if self.tevg:
                        tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getTEVGValue([xPt,yPt,zPt],self.length*self.anastomosisLengthRatio/2.0+self.anastomosisTransitionLength/2.0, self.zcenter)
                        if tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] > 0:
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,4] = 200000000
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,1] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,2] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,3] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,14] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,28] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,42] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,56] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,70] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,84] = 0
                    if self.aneurysm:
                        aneurysmValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getAneurysmValue([xPt,yPt,zPt], self.radius)


        vol = pv.StructuredGrid()
        vol.points = points
        vol.dimensions = [self.numRad+1, self.numCirc+1, self.numLen+1]

        numCells = vol.GetNumberOfCells()
        vol.GetPointData().AddArray(pv.convert_array((proximalIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((distalIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((innerIds).astype(int),name="InnerRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outerIds).astype(int),name="OuterRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((structureIds).astype(int),name="StructureID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="GlobalElementID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="CellStructureID"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="varWallProps"))
        vol.GetPointData().AddArray(pv.convert_array((sliceIds).astype(int),name="sliceIds"))

        #Exchange these to gauss point quantities
        #Matrices/arrays
        vol.GetCellData().AddArray(pv.convert_array(np.zeros((self.numCirc*self.numLen*self.numRad,45*8)).astype(float),name="temp_array"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad_mem"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(36),8),(numCells,1)).astype(float),name="stiffness_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(9),8),(numCells,1)).astype(float),name="sigma_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(6),8),(numCells,1)).astype(float),name="sigma_gnr"))
        #Scalar values
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_curr"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_target"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([0],8),(numCells,1)).astype(float),name="p_est"))

        vol.GetCellData().AddArray(pv.convert_array(np.ones(numCells).astype(float),name="J_c"))
        vol.GetCellData().AddArray(pv.convert_array(np.zeros(numCells).astype(float),name="p_est_c"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_prev"))

        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_prev"))


        vol.GetCellData().AddArray(pv.convert_array(e_r.astype(float),name="e_r"))
        vol.GetCellData().AddArray(pv.convert_array(e_t.astype(float),name="e_t"))
        vol.GetCellData().AddArray(pv.convert_array(e_z.astype(float),name="e_z"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="material_global"))
        vol.GetCellData().AddArray(pv.convert_array(fluidStressQueryIds.astype(int),name="fluidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(solidStressQueryIds.astype(int),name="solidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(simulateIds.astype(int),name="Simulate"))
        vol.GetCellData().AddArray(pv.convert_array(aneurysmValue.astype(float),name="aneurysmValue"))
        vol.GetCellData().AddArray(pv.convert_array(tevgValue.astype(float),name="tevgValue"))
        vol = vol.cast_to_unstructured_grid()
        vol = clean(vol, tolerance=1e-10)
        numPts = vol.GetNumberOfPoints()
        vol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPts,numPts).astype(int),name="GlobalNodeID"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="displacements"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_curr"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_prev"))
        vol.GetPointData().AddArray(pv.convert_array(np.zeros(numPts).astype(float),name="Pressure"))
        vol.GetPointData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numPts).astype(float),name="WSS"))
        vol.GetPointData().AddArray(pv.convert_array(vol.points.astype(float),name="Coordinate"))

        return vol



    def initializeSegmentation(self):
        """
        Initialize points and coordinate system for segmentation
        """
        print("Initializing segmented vessel...")

        #Store ctgr slice points
        ctgrPtsInner = parseCTGR(self.segmentationName+'_inner.ctgr', addFinal = True)
        numSegsInner = np.shape(ctgrPtsInner)[0]

        for i in range(np.shape(ctgrPtsInner)[0]):
            for j in range(np.shape(ctgrPtsInner)[1]):
                    ctgrPtsInner[i,j] = self.scaling_value*np.matmul(self.rotation_matrix,ctgrPtsInner[i,j])

        if not self.estimateOuterSegmentation:
            ctgrPtsOuter = parseCTGR(self.segmentationName+'_outer.ctgr')
            numSegsOuter = np.shape(ctgrPtsOuter)[0]
            if numSegsInner != numSegsOuter:
                raise RuntimeError('Number of inner segmenetations does not equal number of outer segmentations!')
            for i in range(np.shape(ctgrPtsOuter)[0]):
                for j in range(np.shape(ctgrPtsOuter)[1]):
                        ctgrPtsOuter[i,j] = self.scaling_value*np.matmul(self.rotation_matrix,ctgrPtsOuter[i,j])

        self.numLen = int(np.max([1,round((self.numLen) / (numSegsInner-1))]) * (numSegsInner-1))

        #Align inner data to each other
        #Flip first contour to prevent negative jacobian (seems necessary on only come ctgrs?)
        if self.flipContours:
            ctgrPtsInner[0] = flipContour(ctgrPtsInner[0])
        for i in range(np.shape(ctgrPtsInner)[0]-1):
            ctgrPtsInner[i+1] = alignContours(ctgrPtsInner[i],ctgrPtsInner[i+1])

        if not self.estimateOuterSegmentation:
            #Align outer data to inner data
            for i in range(np.shape(ctgrPtsInner)[0]):
                ctgrPtsOuter[i] = alignContours(ctgrPtsInner[i],ctgrPtsOuter[i])


        # Now reparameterize to desired number of points
        ctgrPtsInner = interpolateSplineArray(ctgrPtsInner,periodic=True,numPts=self.numCirc)
        ctgrPtsInner = interpolateSplineArray(np.moveaxis(ctgrPtsInner,0,1),numPts=self.numLen+1,redistribute=False)
        ctgrPtsInner = np.moveaxis(ctgrPtsInner,0,1)

        if not self.estimateOuterSegmentation:
            ctgrPtsOuter = interpolateSplineArray(ctgrPtsOuter,periodic=True,numPts=self.numCirc)
            ctgrPtsOuter = interpolateSplineArray(np.moveaxis(ctgrPtsOuter,0,1),numPts=self.numLen+1,redistribute=False)
            ctgrPtsOuter = np.moveaxis(ctgrPtsOuter,0,1)


        with open('FolderVesselConfigurationFiles/Native_in_handshake_') as f:
            lines = f.readlines()[1:]
            nativeIn = []
            for line in lines:
                nativeIn.append([float(x) for x in line.split()])


        self.sigma_h = nativeIn[13][0]
        self.tau_h = nativeIn[13][1]

        #Build material array (constant for now)

        materialArray = [ nativeIn[7][0],nativeIn[4][0],nativeIn[4][1],nativeIn[4][2],nativeIn[2][0]*10.0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[5][2],nativeIn[2][4]*10.0,nativeIn[2][5],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][3],nativeIn[5][3],nativeIn[2][6]*10.0,nativeIn[2][7],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][4],nativeIn[5][4],nativeIn[2][8]*10.0,nativeIn[2][9],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][5],nativeIn[5][5],nativeIn[2][10]*10.0,nativeIn[2][11],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][1],nativeIn[5][1],nativeIn[2][2]*10.0,nativeIn[2][3],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[16][0],nativeIn[14][2],nativeIn[14][1],0,0,0,0,0,0,0,0,0,0]

        #********************
        ang1 = nativeIn[3][4]
        ang2 = nativeIn[3][5]

        points = np.empty([(self.numCirc+1)*(self.numLen+1)*(self.numRad+1),3])

        structureIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        innerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        outerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        proximalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        distalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        sliceIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))

        if self.flipInlet:
            distalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
            proximalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1
        else:
            proximalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
            distalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1


        innerIds[0:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1
        outerIds[self.numRad:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1

        fluidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)
        solidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)

        tevgValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        aneurysmValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        simulateIds = np.ones(self.numCirc*self.numLen*self.numRad)

        e_r = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_t = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_z = np.zeros((self.numCirc*self.numLen*self.numRad,3))

        e_ma = np.zeros((self.numCirc*self.numLen*self.numRad,98))

        for i in range(self.numLen+1):
            for j in range(self.numCirc+1):
                for k in range(self.numRad+1):
                    center = np.mean(ctgrPtsInner[i],axis=0)
                    sliceIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] = i
                    if not self.estimateOuterSegmentation:
                        xPt = (k/self.numRad)*(ctgrPtsOuter[i,j%self.numCirc,0] - ctgrPtsInner[i,j%self.numCirc,0]) + ctgrPtsInner[i,j%self.numCirc,0]
                        yPt = (k/self.numRad)*(ctgrPtsOuter[i,j%self.numCirc,1] - ctgrPtsInner[i,j%self.numCirc,1]) + ctgrPtsInner[i,j%self.numCirc,1]
                        zPt = (k/self.numRad)*(ctgrPtsOuter[i,j%self.numCirc,2] - ctgrPtsInner[i,j%self.numCirc,2]) + ctgrPtsInner[i,j%self.numCirc,2]
                    else:
                        p_radius = np.linalg.norm(ctgrPtsInner[i,j%self.numCirc,:] - center)
                        if i == 0 or i == self.numLen:
                            p_v_r = (ctgrPtsInner[i,j%self.numCirc,:] - center)/p_radius
                        else:
                            p_v_t = ctgrPtsInner[i,(j+1)%self.numCirc,:] - ctgrPtsInner[i,j%self.numCirc - 1,:]
                            p_v_z = ctgrPtsInner[i + 1,j%self.numCirc,:] - ctgrPtsInner[i - 1,j%self.numCirc,:]
                            p_v_r = np.cross(p_v_t, p_v_z)
                            p_v_r = p_v_r/np.linalg.norm(p_v_r)

                        if self.constantThickness:
                            p_thickness = self.thickness
                        else:
                            p_thickness = p_radius*self.thicknessRatio
                        xPt = (k/self.numRad)*(p_thickness*p_v_r[0]) + ctgrPtsInner[i,j%self.numCirc,0]
                        yPt = (k/self.numRad)*(p_thickness*p_v_r[1]) + ctgrPtsInner[i,j%self.numCirc,1]
                        zPt = (k/self.numRad)*(p_thickness*p_v_r[2]) + ctgrPtsInner[i,j%self.numCirc,2]
                    structureIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] = i*(self.numCirc)*(self.numRad+1) + (j%self.numCirc)*(self.numRad+1) + k + 1
                    points[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k,:] =  np.array([xPt, yPt, zPt])


        coords=[[0, 1, 0],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 0],
                [0, 0, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 0, 0]]


        for i in range(self.numLen):
            for j in range(self.numCirc):
                for k in range(self.numRad):

                    cellPts = []

                    for coord in coords:
                        cellPts.append(points[(i+coord[0])*(self.numCirc+1)*(self.numRad+1) + (j+coord[1])*(self.numRad+1) + (k+coord[2]),:])

                    cellPts = np.array(cellPts)

                    cellCenter = np.mean(cellPts,axis=0)

                    xPt = cellCenter[0]
                    yPt = cellCenter[1]
                    zPt = cellCenter[2]

                    v_r = np.mean(cellPts[[1,2,5,6]],axis=0) - np.mean(cellPts[[0,3,4,7]],axis=0)
                    v_r = v_r/np.linalg.norm(v_r)
                    v_t = np.mean(cellPts[[4,5,6,7]],axis=0) - np.mean(cellPts[[0,1,2,3]],axis=0)
                    v_t = v_t/np.linalg.norm(v_t)
                    v_z = np.mean(cellPts[[2,3,6,7]],axis=0) - np.mean(cellPts[[0,1,4,5]],axis=0)
                    v_z = v_z/np.linalg.norm(v_z)

                    v_r = np.cross(v_z,v_t)
                    v_r = -v_r/np.linalg.norm(v_r)
                    #v_t = np.cross(v_z,v_r)
                    #v_t = v_t/np.linalg.norm(v_t)
                    v_t = np.cross(v_z,v_r)
                    v_t = -v_t/np.linalg.norm(v_t)

                    e_r[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_r
                    e_t[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_t
                    e_z[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_z

                    A = np.array([v_r,v_t,v_z])
                    # Is this a rotation matrix?
                    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                                          np.eye(3, dtype=float)) > 1e-6):
                        print(A)
                        raise RuntimeError('Matrix *A* does not describe a rotation.')

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,:] = materialArray

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,5:8]   = v_r
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,8:11]  = v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,11:14] = v_z

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,18:21] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,32:35] =  np.cos(0.0*np.pi/180.0)*v_z + np.sin(0.0*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,46:49] =  np.cos(ang1*np.pi/180.0)*v_z + np.sin(ang1*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,60:63] =  np.cos(ang2*np.pi/180.0)*v_z + np.sin(ang2*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,74:77] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,88:91] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t

                    fluidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc
                    solidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc + self.numRad

                    if self.tevg:
                        tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getTEVGValue([xPt,yPt,zPt], self.radius, self.zcenter)
                        if tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] > 0:
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,4] = 200000000
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,1] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,2] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,3] = 1
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,14] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,28] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,42] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,56] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,70] = 0
                            e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,84] = 0
                    if self.aneurysm:
                        aneurysmValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getAneurysmValue([xPt,yPt,zPt], self.radius)

        vol = pv.StructuredGrid()
        vol.points = points
        vol.dimensions = [self.numRad+1, self.numCirc+1, self.numLen+1]

        numCells = vol.GetNumberOfCells()
        vol.GetPointData().AddArray(pv.convert_array((proximalIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((distalIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((innerIds).astype(int),name="InnerRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outerIds).astype(int),name="OuterRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((structureIds).astype(int),name="StructureID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="GlobalElementID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="CellStructureID"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="varWallProps"))
        vol.GetPointData().AddArray(pv.convert_array((sliceIds).astype(int),name="sliceIds"))

        #Exchange these to gauss point quantities
        #Matrices/arrays
        vol.GetCellData().AddArray(pv.convert_array(np.zeros((self.numCirc*self.numLen*self.numRad,45*8)).astype(float),name="temp_array"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad_mem"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(36),8),(numCells,1)).astype(float),name="stiffness_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(9),8),(numCells,1)).astype(float),name="sigma_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(6),8),(numCells,1)).astype(float),name="sigma_gnr"))
        #Scalar values
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_curr"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_target"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([0],8),(numCells,1)).astype(float),name="p_est"))

        vol.GetCellData().AddArray(pv.convert_array(np.ones(numCells).astype(float),name="J_c"))
        vol.GetCellData().AddArray(pv.convert_array(np.zeros(numCells).astype(float),name="p_est_c"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_prev"))

        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_prev"))


        vol.GetCellData().AddArray(pv.convert_array(e_r.astype(float),name="e_r"))
        vol.GetCellData().AddArray(pv.convert_array(e_t.astype(float),name="e_t"))
        vol.GetCellData().AddArray(pv.convert_array(e_z.astype(float),name="e_z"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="material_global"))
        vol.GetCellData().AddArray(pv.convert_array(fluidStressQueryIds.astype(int),name="fluidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(solidStressQueryIds.astype(int),name="solidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(simulateIds.astype(int),name="Simulate"))
        vol.GetCellData().AddArray(pv.convert_array(aneurysmValue.astype(float),name="aneurysmValue"))
        vol.GetCellData().AddArray(pv.convert_array(tevgValue.astype(float),name="tevgValue"))
        vol = vol.cast_to_unstructured_grid()
        vol = clean(vol, tolerance=1e-10)
        numPts = vol.GetNumberOfPoints()
        vol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPts,numPts).astype(int),name="GlobalNodeID"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="displacements"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_curr"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_prev"))
        vol.GetPointData().AddArray(pv.convert_array(np.zeros(numPts).astype(float),name="Pressure"))
        vol.GetPointData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numPts).astype(float),name="WSS"))
        vol.GetPointData().AddArray(pv.convert_array(vol.points.astype(float),name="Coordinate"))

        return vol



    def generateFluidMesh(self, surf):
        print("Making inner volume...")
        pid = 0
        numQuad = self.numCirc // 4 + 1
        numPerim = self.numCirc // 4
        quadFrac = 0.1
        numTrans = 40
        points = []

        inIds = np.zeros((numQuad**2 + (numTrans-1)*self.numCirc)*(self.numLen+1))
        outIds = np.zeros((numQuad**2 + (numTrans-1)*self.numCirc)*(self.numLen+1))
        wallIds = np.zeros((numQuad**2 + (numTrans-1)*self.numCirc)*(self.numLen+1))

        if self.flipInlet:
            inIds[0:numQuad**2] = 1
            inIds[((self.numLen+1)*numQuad**2):((self.numLen+1)*numQuad**2)+ (numTrans-1)*self.numCirc] = 1

            outIds[((self.numLen+1)-1)*numQuad**2:(self.numLen+1)*numQuad**2] = 1
            outIds[-(numTrans-1)*self.numCirc:] = 1
        else:
            outIds[0:numQuad**2] = 1
            outIds[((self.numLen+1)*numQuad**2):((self.numLen+1)*numQuad**2)+ (numTrans-1)*self.numCirc] = 1

            inIds[((self.numLen+1)-1)*numQuad**2:(self.numLen+1)*numQuad**2] = 1
            inIds[-(numTrans-1)*self.numCirc:] = 1

        for ia in range(self.numLen+1):
            wallIds[((self.numLen+1)*numQuad**2) + (numTrans-2)*self.numCirc + ia*(numTrans-1)*(self.numCirc):((self.numLen+1)*numQuad**2) + (numTrans-1)*self.numCirc + ia*(numTrans-1)*(self.numCirc)] = 1


        numPoints = surf.GetNumberOfPoints()
        origIdsDict = dict(zip((surf.get_array('StructureID')-1).tolist(),range(numPoints)))

        # generate quadratic mesh
        for ia in range((self.numLen+1)):
            c1 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*0*(self.numCirc // 4)]))
            c2 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*1*(self.numCirc // 4)]))
            c3 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*2*(self.numCirc // 4)]))
            c4 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*3*(self.numCirc // 4)]))
            center = (c1 + c2 + c3 + c4)/4.0

            M = np.linalg.lstsq(np.array([[-0.5,-0.5,0,1],[0.5,-0.5,0,1],[0.5,0.5,0,1],[-0.5,0.5,0,1]]),np.vstack((c1,c2,c3,c4)))[0]
            for iy in range(numQuad):
                for ix in range(numQuad):

                    rady = quadFrac * ( iy / (numQuad - 1) - 0.5)
                    radx = quadFrac * ( ix / (numQuad - 1) - 0.5)

                    points.append(np.matmul([radx, rady, 0,1],M))
                    pid += 1

        # generate transition mesh
        for ia in range(self.numLen+1):
            for ir in range(1, numTrans):
                for ic in range(self.numCirc):
                    if ic <= numPerim:
                        quadPoint = points[ic + ia*(numQuad**2)]
                    elif numPerim < ic <= numPerim*2:
                        quadPoint = points[(ic-numPerim + 1)*numQuad - 1 + ia*(numQuad**2)]
                    elif numPerim*2 < ic <= numPerim*3:
                        quadPoint = points[numQuad**2 - 1 - (ic-numPerim*2) + ia*(numQuad**2)]
                    elif numPerim*3 < ic <= numPerim*4:
                        quadPoint = points[(numPerim-(ic-numPerim*3))*numQuad  + ia*(numQuad**2)]

                    outerPoint = surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + ic*(self.numRad+1)])

                    points.append(((outerPoint - quadPoint) * ir/(numTrans - 1)) + quadPoint)
                    pid += 1


        coords=[[0, 1, 0],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 0],
                [0, 0, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 0, 0]]

        cid = 0
        cells = []

        # generate quadratic mesh
        for ia in range(self.numLen):
            for iy in range(numQuad - 1):
                for ix in range(numQuad - 1):
                    ids = []
                    for c in coords:
                        ids += [(iy + c[0]) * numQuad + ix + c[1] + (ia + c[2]) * numQuad ** 2]
                    cells.append(ids)
                    cid += 1
     
        # generate transition mesh
        for ia in range(self.numLen):
            for ic in range(self.numCirc):
                ids = []
                for c in coords:
                        if c[1] == 1:
                            # circular side
                            ids += [(ic + c[0])%self.numCirc + (self.numLen+1) * numQuad ** 2 + (ia + c[2]) * (numTrans - 1) * self.numCirc]
                        else:
                            if (ic + c[0]) <= numPerim:
                                ids += [(ic + c[0]) + (ia + c[2])*(numQuad**2)]
                            elif numPerim < (ic + c[0]) <= numPerim*2:
                                ids += [((ic + c[0])-numPerim + 1)*numQuad - 1 + (ia + c[2])*(numQuad**2)]
                            elif numPerim*2 < (ic + c[0]) <= numPerim*3:
                                ids += [numQuad**2 - 1 - ((ic + c[0])-numPerim*2) + (ia + c[2])*(numQuad**2)]
                            elif numPerim*3 < (ic + c[0]) <= numPerim*4:
                                ids += [(numPerim-((ic + c[0])-numPerim*3))*numQuad  + (ia + c[2])*(numQuad**2)]
                cells.append(ids)
                cid += 1

        # generate circular mesh
        for ia in range(self.numLen):
            for ir in range(numTrans-2):
                for ic in range(self.numCirc):
                    ids = []
                    for c in coords:
                            ids += [(ic + c[0])%self.numCirc + (self.numLen+1) * numQuad ** 2 + (ia + c[2]) * (numTrans - 1) * self.numCirc + (ir + c[1]) * self.numCirc]
                    cells.append(ids)
                    cid += 1


        # each cell is a VTK_HEXAHEDRON
        cellArray = [ [8] + cell for cell in cells]
        cellTypes = np.empty(np.shape(cells)[0], dtype=np.uint8)
        cellTypes[:] = vtk.VTK_HEXAHEDRON

        vol = pv.UnstructuredGrid(np.array(cellArray), np.array(cellTypes), np.array(points))

        vol.GetPointData().AddArray(pv.convert_array((inIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((wallIds).astype(int),name="OuterRegionID"))

        numPts = vol.GetNumberOfPoints()
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="Velocity"))
        vol.GetPointData().AddArray(pv.convert_array(np.zeros(numPts).astype(float),name="Pressure"))
        return vol

    def setInputFileValues(self):

        with open(self.simulationInputDirectory + '/input_aniso.mfs', 'r') as file:
            data = file.readlines()
        data[86] = "      Penalty parameter: " + str(self.penalty) + "\n"
        data[88] = "      Mass damping: " + str(self.damping) + "\n"
        data[118] = "      Value: " + str(self.inletFlow) + "\n"
        data[130] = "      Value: " + str(- self.outletPressure / self.inletFlow) + "\n"
        with open(self.simulationInputDirectory + '/input_aniso.mfs', 'w') as file:
            file.writelines(data)

        with open(self.simulationInputDirectory + '/input_mm.mfs', 'r') as file:
            data = file.readlines()
        data[89] = "      Penalty parameter: " + str(self.penalty) + "\n"
        data[91] = "      Mass damping: " + str(self.damping) + "\n"
        data[120] = "      Value: " + str(self.inletFlow) + "\n"
        data[132] = "      Value: " + str(- self.outletPressure / self.inletFlow) + "\n"
        with open(self.simulationInputDirectory + '/input_mm.mfs', 'w') as file:
            file.writelines(data)

        with open(self.simulationInputDirectory + '/input_fluid.mfs', 'r') as file:
            data = file.readlines()
        data[73] = "      Value: " + str(self.inletFlow) + "\n"
        data[81] = "      Value: " + str(- self.outletPressure / self.inletFlow) + "\n"
        with open(self.simulationInputDirectory + '/input_fluid.mfs', 'w') as file:
            file.writelines(data)

        with open(self.simulationInputDirectory + '/solid_aniso.mfs', 'r') as file:
            data = file.readlines()
        data[59] = "   Penalty parameter: " + str(self.penalty) + "\n"
        data[61] = "   Mass damping: " + str(self.damping) + "\n"
        data[63] = "   LS type: " + str(self.solidLinearSolverType) + " {\n"
        if self.solidLinearSolverType == "GMRES":
            data[65] = "      Max iterations: 10\n"
            data[66] = "      Krylov space dimension: 50\n"
        else:
            data[65] = "      Max iterations: 1000\n"
            data[66] = "      Krylov space dimension: 250\n"
        with open(self.simulationInputDirectory + '/solid_aniso.mfs', 'w') as file:
            file.writelines(data)

        with open(self.simulationInputDirectory + '/solid_mm.mfs', 'r') as file:
            data = file.readlines()
        data[60] = "   Penalty parameter: " + str(self.penalty) + "\n"
        data[62] = "   Mass damping: " + str(self.damping) + "\n"
        #data[58] = "   LS type: " + str(self.solidLinearSolverType) + " {\n"
        # if self.solidLinearSolverType == "GMRES":
        #     data[60] = "      Max iterations: 10\n"
        #     data[61] = "      Krylov space dimension: 50\n"
        # else:
        #     data[60] = "      Max iterations: 1000\n"
        #     data[61] = "      Krylov space dimension: 250\n"
        with open(self.simulationInputDirectory + '/solid_mm.mfs', 'w') as file:
            file.writelines(data)



    """
    def setPressure()

    def setTimestep()

    def setDampening()

    def setPenalty()

    """
