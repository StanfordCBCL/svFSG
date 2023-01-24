import numpy as np
import scipy as sp
import pyvista as pv
import vtk
import os
import pickle
from tqdm import tqdm
import multiprocessing
from multiprocessing import Pool, Manager
from cvessel import cvessel
import time
os.system("taskset -p 0xff %d" % os.getpid())


def read_data(file_name, file_format="vtp", datatype=None):
    """
    Read surface geometry from a file.

    Args:
        file_name (str): Path to input file.
        file_format (str): File format (.vtp, .stl, etc.). 
        datatype (str): Additional parameter for vtkIdList objects.

    Returns:
        polyData (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Output data.
    """

    # Check if file exists
    if not os.path.exists(file_name):
        raise RuntimeError("Could not find file: %s" % file_name)

    # Get reader
    if file_format == 'vtk':
        reader = vtk.vtkPolyDataReader()
    elif file_format == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif file_format == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        raise RuntimeError('Unknown file type %s' % file_format)

    # Read surface geometry.
    reader.SetFileName(file_name)
    reader.Update()
    polydata = reader.GetOutput()

    return pv.wrap(polydata)


def save_data(file_name, data):

    """ Write the given VTK object to a file.

    Args:
        file_name (str): The name of the file to write. 
        data (vtkDataObject): Data to write.
        datatype (str): Additional parameter for vtkIdList objects.
    """

    # Check filename format.
    file_ext = file_name.split(".")[-1]

    data = pv.wrap(data)
    data.points_to_double()

    if file_ext == '':
        raise RuntimeError('The file does not have an extension')

    # Get writer.
    if file_ext == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif file_ext == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        raise RuntimeError('Unknown file type %s' % file_ext)

    # Set file name and the data to write.
    writer.SetFileName(file_name)
    writer.SetInputData(data)
    writer.Update()

    # Write the data.
    writer.Write()

def thresholdModel(data,dataName,low,high,extract=True,cell=False):
    if cell:
        data.GetCellData().SetActiveScalars(dataName)
    else:
        data.GetPointData().SetActiveScalars(dataName)
    t = vtk.vtkThreshold()
    t.SetInputData(data)
    t.ThresholdBetween(low,high)
    if cell:
        t.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, vtk.vtkDataSetAttributes.SCALARS)
    else:
        t.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, vtk.vtkDataSetAttributes.SCALARS)
    t.Update()
    if extract:
        t_surf = vtk.vtkDataSetSurfaceFilter()
        t_surf.SetInputData(t.GetOutput())
        t_surf.Update()
        return pv.wrap(t_surf.GetOutput())
    else:
        return pv.wrap(t.GetOutput())

                
def initializePoolWorker(pool_input_array):
    global output_array
    global cvessels
    global rank
    global numData

    rank = multiprocessing.current_process()._identity[0]-1
    print("Initializing data on rank " + str(rank))
    input_array = pool_input_array[rank]
    numData = len(input_array)

    cvessels = [cvessel] * numData

    for i in range(numData):
        input_data = input_array[i]
        outputDir,num_id,max_days,gnr_step_size, aneurysmValue, tevgValue = input_data
        cvessels[i] = cvessel()
        cvessels[i].initializeVesselHandshake(outputDir,num_id,max_days,gnr_step_size, aneurysmValue, tevgValue)
    
    output_array = [[0]*59]*(numData)

    return

def runPoolWorker(pool_input_array):

    print("Running data on rank " + str(rank))
    input_array = pool_input_array[rank]

    for i in range(numData):
        input_data = input_array[i]
        timeStep, timeIter,sigma_inv,wss,defGrad_mem_g = input_data
        output_array[i] = cvessels[i].updateVesselHandshake(timeStep, timeIter,sigma_inv,wss,defGrad_mem_g)

    return output_array


def clean(data, tolerance=None):
    """Temporary reimplimentation of the algorithm to be faster.
    Has seen 1.1-1.5 times speedup so far
    """
    output = pv.UnstructuredGrid()
    output.Allocate(data.n_cells)

    output.GetPointData().CopyAllocate(data.GetPointData())

    out_cell_data = output.GetCellData()
    out_cell_data.CopyGlobalIdsOn()
    out_cell_data.CopyAllocate(data.GetCellData())

    if tolerance is None:
        tolerance = data.length * 1e-6
    all_nodes = data.points.copy()
    all_nodes = np.around(all_nodes / tolerance)
    unique_nodes, index, ind_nodes = np.unique(all_nodes, return_index=True, return_inverse=True, axis=0)
    unique_nodes *= tolerance
    new_points = pv.vtk_points(unique_nodes)
    output.SetPoints(new_points)
    for name, arr in data.point_arrays.items():
        output[name] = arr[index]

    cell_set = set()
    cell_points = vtk.vtkIdList()
    for i in range(data.n_cells):
        # special handling for polyhedron cells
        cell_type = data.GetCellType(i)
        data.GetCellPoints(i, cell_points)
        if (cell_type == vtk.VTK_POLYHEDRON):
            vtk.vtkUnstructuredGrid.SafeDownCast(data).GetFaceStream(i, cell_points)
            vtk.vtkUnstructuredGrid.ConvertFaceStreamPointIds(cell_points, ind_nodes)
            output.InsertNextCell(data.GetCellType(i), cell_points)
        else:
            if (cell_type == vtk.VTK_POLY_VERTEX or cell_type == vtk.VTK_TRIANGLE_STRIP):
                for j in range(cell_points.GetNumberOfIds()):
                    cell_pt_id = cell_points.GetId(j)
                    new_id = ind_nodes[cell_pt_id]
                    cell_points.SetId(j, new_id)
                new_cell_id = output.InsertNextCell(cell_type, cell_points)
                out_cell_data.CopyData(data.GetCellData(), i, new_cell_id)
                continue

            nn = set()
            for j in range(cell_points.GetNumberOfIds()):
                cell_pt_id = cell_points.GetId(j)
                new_id = ind_nodes[cell_pt_id]
                cell_points.SetId(j, new_id)
                cell_pt_id = cell_points.GetId(j)
                nn.add(cell_pt_id)

                is_unique = not (frozenset(nn) in cell_set)

                # only copy a cell to the output if it is neither degenerate nor duplicate
                if len(nn) == cell_points.GetNumberOfIds() and is_unique:
                    new_cell_id = output.InsertNextCell(data.GetCellType(i), cell_points)
                    out_cell_data.CopyData(data.GetCellData(), i, new_cell_id)
                    cell_set.add(frozenset(nn))

    return output


def getAneurysmValue(point,radius):
    """                                                                                                                                                                                                      
    Get the value of vessel behavior based on point location and radius                                                                                                                                      
    """
    zPt = point[2]
    vesselValue = 0.65*np.exp(-abs(zPt/(1.5/4.0))**2)
    return vesselValue


def getTEVGValue(point, radius):
    """                                                                                                                                                                                                      
    Get the value of vessel behavior based on point location and radius                                                                                                                                      
    """
    value = 0.0
    if abs(point[2]) <= radius:
        value = 1.0*(0.5*np.cos((np.pi/radius)*point[2])+0.5) + 0.25 
    return value

def rotate_elastic_constants(C, A, tol=1e-4):
    """
    Return rotated elastic moduli for a general crystal given the elastic 
    constant in Voigt notation.
    Parameters
    ----------
    C : array_like
        6x6 matrix of elastic constants (Voigt notation).
    A : array_like
        3x3 rotation matrix.
    Returns
    -------
    C : array
        6x6 matrix of rotated elastic constants (Voigt notation).
    """

    A = np.asarray(A)

    # Is this a rotation matrix?
    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                          np.eye(3, dtype=float)) > tol):
        print(A)
        raise RuntimeError('Matrix *A* does not describe a rotation.')

    # Rotate
    return full_3x3x3x3_to_Voigt_6x6(np.einsum('ia,jb,kc,ld,abcd->ijkl',
                                               A, A, A, A,
                                               Voigt_6x6_to_full_3x3x3x3(C)))

def Voigt_6x6_to_full_3x3x3x3(C):
    Voigt_notation = np.array([[0,3,5],[3,1,4],[5,4,2]])
    """
    Convert from the Voigt representation of the stiffness matrix to the full
    3x3x3x3 representation.
    Parameters
    ----------
    C : array_like
        6x6 stiffness matrix (Voigt notation).
    
    Returns
    -------
    C : array_like
        3x3x3x3 stiffness matrix.
    """
    
    C = np.asarray(C)
    C_out = np.zeros((3,3,3,3), dtype=float)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Voigt_i = Voigt_notation[i, j]
                    Voigt_j = Voigt_notation[k, l]
                    C_out[i, j, k, l] = C[Voigt_i, Voigt_j]
    return C_out

def full_3x3x3x3_to_Voigt_6x6(C, tol=1e-4, check_symmetry=True):
    """
    Convert from the full 3x3x3x3 representation of the stiffness matrix
    to the representation in Voigt notation. Checks symmetry in that process.
    """
    # The indices of the full stiffness matrix of (orthorhombic) interest
    Voigt_notation = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]

    C = np.asarray(C)
    Voigt = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            Voigt[i,j] = C[k,l,m,n]
            if check_symmetry:
                # MAJOR SYMMETRY
                #assert abs(Voigt[i,j]-C[m,n,k,l]) < tol, \
                #    '1 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                #    .format(i, j, Voigt[i,j], m, n, k, l, C[m,n,k,l])
                assert abs(Voigt[i,j]-C[l,k,m,n]) < tol, \
                    '2 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                    .format(i, j, Voigt[i,j], l, k, m, n, C[l,k,m,n])
                assert abs(Voigt[i,j]-C[k,l,n,m]) < tol, \
                    '3 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                    .format(i, j, Voigt[i,j], k, l, n, m, C[k,l,n,m])
                #assert abs(Voigt[i,j]-C[m,n,l,k]) < tol, \
                #    '4 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                #    .format(i, j, Voigt[i,j], m, n, l, k, C[m,n,l,k])
                #assert abs(Voigt[i,j]-C[n,m,k,l]) < tol, \
                #    '5 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                #    .format(i, j, Voigt[i,j], n, m, k, l, C[n,m,k,l])
                assert abs(Voigt[i,j]-C[l,k,n,m]) < tol, \
                    '6 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                    .format(i, j, Voigt[i,j], l, k, n, m, C[l,k,n,m])
                #assert abs(Voigt[i,j]-C[n,m,l,k]) < tol, \
                #    '7 Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                #    .format(i, j, Voigt[i,j], n, m, l, k, C[n,m,l,k])

    return Voigt

# Implementation from  https://github.com/libAtoms/matscipy/blob/master/matscipy/elasticity.py

def rotateStress(sigma,Q):
    Q = np.transpose(Q)
    Qinv = np.transpose(Q)

    if np.sometrue(np.abs(np.dot(np.array(Q), np.transpose(np.array(Q))) - 
                          np.eye(3, dtype=float)) > 1e-6):
        print(Q)
        raise RuntimeError('Matrix *Q* does not describe a rotation.')

    sigma = np.reshape(sigma, (3,3))
    sigma = np.matmul(Q,np.matmul(sigma,Qinv))
    # Switch Voigt notation to match svFSI
    if not is_symmetric(sigma):
        print(sigma)
        raise ValueError("Sigma is not symmetric!")
    sigma = np.array([sigma[0,0],sigma[1,1],sigma[2,2],sigma[0,1],sigma[1,2],sigma[2,0]])
    return sigma

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def is_symmetric(a, tol=1e-2):
    return np.all(np.abs(a-a.T) < tol)

def rotateStiffness(DD,Q):
    Q = np.transpose(Q)
    Qinv = np.transpose(Q)

    #Definitely should speed this up
    CC = np.reshape(DD, (6,6))

    CC = rotate_elastic_constants(CC,Q)

    if not is_symmetric(CC):
        print(CC)
        raise ValueError("CC is not symmetric!")
    """
    if not is_pos_def(CC):
        print(CC)
        raise ValueError("CC is not positive definite!")
    """

    CC = np.reshape(CC, (1,36))[0]

    return CC


def getNN(xi):

    N = np.zeros(8)
    Nxi = np.zeros((3,8))

    lx = 1.0 - xi[0];
    ly = 1.0 - xi[1];
    lz = 1.0 - xi[2];
    ux = 1.0 + xi[0];
    uy = 1.0 + xi[1];
    uz = 1.0 + xi[2];

    N[0] = lx*ly*lz/8.0;
    N[1] = ux*ly*lz/8.0;
    N[2] = ux*uy*lz/8.0;
    N[3] = lx*uy*lz/8.0;
    N[4] = lx*ly*uz/8.0;
    N[5] = ux*ly*uz/8.0;
    N[6] = ux*uy*uz/8.0;
    N[7] = lx*uy*uz/8.0;

    Nxi[0,0] = -ly*lz/8.0;
    Nxi[1,0] = -lx*lz/8.0;
    Nxi[2,0] = -lx*ly/8.0;
    Nxi[0,1] =  ly*lz/8.0;
    Nxi[1,1] = -ux*lz/8.0;
    Nxi[2,1] = -ux*ly/8.0;
    Nxi[0,2] =  uy*lz/8.0;
    Nxi[1,2] =  ux*lz/8.0;
    Nxi[2,2] = -ux*uy/8.0;
    Nxi[0,3] = -uy*lz/8.0;
    Nxi[1,3] =  lx*lz/8.0;
    Nxi[2,3] = -lx*uy/8.0;
    Nxi[0,4] = -ly*uz/8.0;
    Nxi[1,4] = -lx*uz/8.0;
    Nxi[2,4] =  lx*ly/8.0;
    Nxi[0,5] =  ly*uz/8.0;
    Nxi[1,5] = -ux*uz/8.0;
    Nxi[2,5] =  ux*ly/8.0;
    Nxi[0,6] =  uy*uz/8.0;
    Nxi[1,6] =  ux*uz/8.0;
    Nxi[2,6] =  ux*uy/8.0;
    Nxi[0,7] = -uy*uz/8.0;
    Nxi[1,7] =  lx*uz/8.0;
    Nxi[2,7] =  lx*uy/8.0;

    return N, Nxi


def computeGaussValues(mesh,name):

    numCells = mesh.GetNumberOfCells()


    nG = 8

    s =  1/np.sqrt(3);
    t = -1/np.sqrt(3);

    xi = np.zeros((3,8))


    xi[0,0] = t; xi[1,0] = t; xi[2,0] = t;
    xi[0,1] = s; xi[1,1] = t; xi[2,1] = t;
    xi[0,2] = s; xi[1,2] = s; xi[2,2] = t;
    xi[0,3] = t; xi[1,3] = s; xi[2,3] = t;
    xi[0,4] = t; xi[1,4] = t; xi[2,4] = s;
    xi[0,5] = s; xi[1,5] = t; xi[2,5] = s;
    xi[0,6] = s; xi[1,6] = s; xi[2,6] = s;
    xi[0,7] = t; xi[1,7] = s; xi[2,7] = s;

    IdM = np.eye(3)
    gaussN = np.zeros(3)
    gaussNx = np.zeros(9)
    allGaussN = np.zeros((numCells,3*nG))
    allGaussNx = np.zeros((numCells,9*nG))

    for q in range(numCells):


        cell = mesh.GetCell(q)
        cellPts = cell.GetPointIds()

        coordinates = []
        field = []
        for p in range(cellPts.GetNumberOfIds()):
            ptId = cellPts.GetId(p)
            field.append(mesh.GetPointData().GetArray(name).GetTuple3(ptId))
            coordinates.append(mesh.GetPoint(ptId))

        field = np.array(field)
        coordinates = np.array(coordinates)

        for g in range(nG):
            N, Nxi = getNN(xi[:,g])
            J = np.matmul(Nxi,coordinates);
            Nx = np.matmul(np.linalg.inv(J),Nxi);
            gaussN = np.matmul(N,field)
            gaussNx = np.ravel(np.matmul(Nx,field) + IdM)
            allGaussN[q, g*3:(g+1)*3] = gaussN
            allGaussNx[q, g*9:(g+1)*9] = gaussNx

    mesh.GetCellData().AddArray(pv.convert_array(allGaussNx,name="defGrad"))

    return mesh

def interpolateSolution(source, target):

    numPoints = target.GetNumberOfPoints()

    pointLocator = vtk.vtkPointLocator()
    pointLocator.SetDataSet(source)
    pointLocator.BuildLocator()

    velocity_array = np.zeros((numPoints,3))
    pressure_array = np.zeros((numPoints,1))

    for q in range(numPoints):
        coordinate = target.GetPoint(q)
        pointIdSource = pointLocator.FindClosestPoint(coordinate)

        target.GetPointData().GetArray('Pressure').SetTuple1(q,source.GetPointData().GetArray('Pressure').GetTuple1(pointIdSource))
        target.GetPointData().GetArray('Velocity').SetTuple(q,source.GetPointData().GetArray('Velocity').GetTuple(pointIdSource))

    return target
