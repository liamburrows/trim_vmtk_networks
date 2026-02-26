import csv
import vtk
#from vmtk import vmtkscripts
import slicer
import numpy as np
import nibabel as nib

def visualise_slicer(surf_obj, node_name = 'object'):
    network_node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
    network_node.SetAndObservePolyData(surf_obj)
    network_node.SetName(node_name)
    network_display = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelDisplayNode')
    network_node.SetAndObserveDisplayNodeID(network_display.GetID())
    network_display.SetColor(.2, .5, .5)
    network_display.SetOpacity(0.5)
    slicer.app.processEvents()

class ExtractCenterlineLogic():

    def __init__(self):
        self.blankingArrayName = 'Blanking'
        self.RadiusArrayName = 'Radius'  # maximum inscribed sphere radius
        #self.RadiusArrayName = 'MaximumInscribedSphereRadius'  # https://github.com/vmtk/vmtk/blob/master/vmtkScripts/vmtkcenterlines.py#L432
        self.groupIdsArrayName = 'GroupIds'
        self.centerlineIdsArrayName = 'CenterlineIds'
        self.tractIdsArrayName = 'TractIds'
        self.topologyArrayName = 'Topology'
        self.marksArrayName = 'Marks'
        self.lengthArrayName = 'Length'
        self.curvatureArrayName = 'Curvature'
        self.torsionArrayName = 'Torsion'
        self.tortuosityArrayName = 'Tortuosity'
        self.frenetTangentArrayName = 'FrenetTangent'
        self.frenetNormalArrayName = 'FrenetNormal'
        self.frenetBinormalArrayName = 'FrenetBinormal'

        self.StopFastMarchingOnReachingTarget = 0

    def read_fcsv(self, fcsv_file, single_list=False):
        """
        Reads end points from an .fcsv file and returns them as a list of (x, y, z) coordinates.
        """
        source_endpoints = []
        target_endpoints = []
        count = 0
        with open(fcsv_file, newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                # Skip lines starting with '#'
                if row[0].startswith('#'):
                    continue
                # Extract x, y, z coordinates (assumed to be columns 1, 2, 3)
                
                x = -float(row[1])
                y = -float(row[2])
                z = float(row[3])
                #endpoints.append([x, y, z])
                label = int(row[9]) 
            
                if label == 0:
                    source_endpoints.append([x, y, z])
                    if len(source_endpoints) > 1:
                        raise ValueError(f'There is more than one end point labelled as the source in fcsv file: {fcsv_file}')
                else:
                    target_endpoints.append([x, y, z])
                count = count+1

        
        if single_list:
            source_flatten = [x for xs in source_endpoints for x in xs]
            target_flatten = [x for xs in target_endpoints for x in xs]
            return source_flatten, target_flatten
        else:
            return source_endpoints, target_endpoints

    def extractCenterline_Ollie(self, surfacePolyData, source_point, end_points):

        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc
        
        surfaceCleaner = vtk.vtkCleanPolyData()
        surfaceCleaner.SetInputData(surfacePolyData)
        surfaceCleaner.Update()

        surfaceTriangulator = vtk.vtkTriangleFilter()
        surfaceTriangulator.SetInputConnection(surfaceCleaner.GetOutputPort())
        surfaceTriangulator.PassLinesOff()
        surfaceTriangulator.PassVertsOff()
        surfaceTriangulator.Update()

        centerlineInputSurface = surfaceTriangulator.GetOutput()

        sourceIdList = vtk.vtkIdList()
        targetIdList = vtk.vtkIdList()
        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(centerlineInputSurface)
        pointLocator.BuildLocator()
        
        sourceId = pointLocator.FindClosestPoint(source_point)
        sourceIdList.InsertNextId(sourceId)
        for target_points in end_points:
            pointId = pointLocator.FindClosestPoint(target_points)
            targetIdList.InsertNextId(pointId)


        centerlineFilter = vtkvmtkComputationalGeometry.vtkvmtkPolyDataCenterlines()
        centerlineFilter.SetInputData(centerlineInputSurface)
        centerlineFilter.SetSourceSeedIds(sourceIdList)
        centerlineFilter.SetTargetSeedIds(targetIdList)
        centerlineFilter.SetRadiusArrayName(self.RadiusArrayName)
        centerlineFilter.SetCostFunction('1/R')
        centerlineFilter.SetFlipNormals(False)
        centerlineFilter.SetAppendEndPointsToCenterlines(0)
        
        enableVoronoiSmoothing = 1
        centerlineFilter.SetSimplifyVoronoi(enableVoronoiSmoothing)

        if self.StopFastMarchingOnReachingTarget == True:
            if targetIdList.GetNumberOfIds() != 1:
                self.PrintError('Parameter Conflict: cannot enable "StopFastMarchingOnReachingTarget" when there is more then one target seed set.')
            else:
                centerlineFilter.SetStopFastMarchingOnReachingTarget(self.StopFastMarchingOnReachingTarget)
        
        centerlineFilter.SetCenterlineResampling(0)
        centerlineFilter.SetResamplingStepLength(1.0)
        centerlineFilter.Update()

        #self.Centerlines = centerlineFilter.GetOutput()
        #self.VoronoiDiagram = centerlineFilter.GetVoronoiDiagram()
        #self.PoleIds = centerlineFilter.GetPoleIds()
        centerlinePolyData = vtk.vtkPolyData()
        centerlinePolyData.DeepCopy(centerlineFilter.GetOutput())

        voronoiDiagramPolyData = vtk.vtkPolyData()
        voronoiDiagramPolyData.DeepCopy(centerlineFilter.GetVoronoiDiagram())
        #self.EikonalSolutionArrayName = centerlineFilter.GetEikonalSolutionArrayName()
        #self.EdgeArrayName = centerlineFilter.GetEdgeArrayName()
        #self.EdgePCoordArrayName = centerlineFilter.GetEdgePCoordArrayName()
        #self.CostFunctionArrayName = centerlineFilter.GetCostFunctionArrayName()
        return centerlinePolyData, voronoiDiagramPolyData
        
    
    def extractCenterline(self, surfacePolyData, source_point, end_points, curveSamplingDistance=1.0):
        """Compute centerline.
        This is more robust and accurate but takes longer than the network extraction.
        :param surfacePolyData:
        :param source_point (x,y,z):
        :param end_points:  A list of tuples: [ (x1,y1,z1), (x2,y2,z2), ... ]
        :return:
        """

        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc

        # Cap all the holes that are in the mesh that are not marked as endpoints
        # Maybe this is not needed.
        capDisplacement = 0.0
        surfaceCapper = vtkvmtkComputationalGeometry.vtkvmtkCapPolyData()
        surfaceCapper.SetInputData(surfacePolyData)
        surfaceCapper.SetDisplacement(capDisplacement)
        surfaceCapper.SetInPlaneDisplacement(capDisplacement)
        surfaceCapper.Update()

        tubePolyData = surfaceCapper.GetOutput()


        sourceIdList = vtk.vtkIdList()
        targetIdList = vtk.vtkIdList()
        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(tubePolyData)
        pointLocator.BuildLocator()
        
        sourceId = pointLocator.FindClosestPoint(source_point)
        sourceIdList.InsertNextId(sourceId)
        for target_points in end_points:
            pointId = pointLocator.FindClosestPoint(target_points)
            targetIdList.InsertNextId(pointId)

        print('checking')
        #slicer.tubePolyData = tubePolyData

        centerlineFilter = vtkvmtkComputationalGeometry.vtkvmtkPolyDataCenterlines()
        centerlineFilter.SetInputData(tubePolyData)
        centerlineFilter.SetSourceSeedIds(sourceIdList)
        centerlineFilter.SetTargetSeedIds(targetIdList)

        centerlineFilter.SetRadiusArrayName(self.RadiusArrayName)
        centerlineFilter.SetCostFunction('1/R')  # this makes path search prefer go through points with large radius
        centerlineFilter.SetFlipNormals(False)
        centerlineFilter.SetAppendEndPointsToCenterlines(0)

        # Voronoi smoothing slightly improves connectivity
        # Unfortunately, Voronoi smoothing is broken if VMTK is used with VTK9, therefore
        # disable this feature for now (https://github.com/vmtk/SlicerExtension-VMTK/issues/34)
        #enableVoronoiSmoothing = (slicer.app.majorVersion * 100 + slicer.app.minorVersion < 413)
        enableVoronoiSmoothing = 1
        centerlineFilter.SetSimplifyVoronoi(enableVoronoiSmoothing)

        if self.StopFastMarchingOnReachingTarget == True:
            if targetIdList.GetNumberOfIds() != 1:
                self.PrintError('Parameter Conflict: cannot enable "StopFastMarchingOnReachingTarget" when there is more then one target seed set.')
            else:
                centerlineFilter.SetStopFastMarchingOnReachingTarget(self.StopFastMarchingOnReachingTarget)


        centerlineFilter.SetCenterlineResampling(0)
        centerlineFilter.SetResamplingStepLength(curveSamplingDistance)
        centerlineFilter.Update()

        centerlinePolyData = vtk.vtkPolyData()
        centerlinePolyData.DeepCopy(centerlineFilter.GetOutput())

        voronoiDiagramPolyData = vtk.vtkPolyData()
        voronoiDiagramPolyData.DeepCopy(centerlineFilter.GetVoronoiDiagram())

        #logging.debug("End of Centerline Computation..")
        return centerlinePolyData, voronoiDiagramPolyData

    def extractCenterline_fcsv(self, surfacePolyData, fcsv_file, curveSamplingDistance=1.0):
        """Compute centerline.
        This is more robust and accurate but takes longer than the network extraction.
        :param surfacePolyData:
        :param endPointsMarkupsNode:
        :return:
        """

        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc

        # Cap all the holes that are in the mesh that are not marked as endpoints
        # Maybe this is not needed.
        capDisplacement = 0.0
        surfaceCapper = vtkvmtkComputationalGeometry.vtkvmtkCapPolyData()
        surfaceCapper.SetInputData(surfacePolyData)
        surfaceCapper.SetDisplacement(capDisplacement)
        surfaceCapper.SetInPlaneDisplacement(capDisplacement)
        surfaceCapper.Update()

        tubePolyData = surfaceCapper.GetOutput()

        seed_method = 1

        if seed_method == 1:
            sourceIdList = vtk.vtkIdList()
            targetIdList = vtk.vtkIdList()
            pointLocator = vtk.vtkPointLocator()
            pointLocator.SetDataSet(tubePolyData)
            pointLocator.BuildLocator()
            
            source_endpoints, target_endpoints = self.read_fcsv(fcsv_file)
            sourceId = pointLocator.FindClosestPoint(source_endpoints[0])
            sourceIdList.InsertNextId(sourceId)
            for target_points in target_endpoints:
                pointId = pointLocator.FindClosestPoint(target_points)
                targetIdList.InsertNextId(pointId)
        else:
            SeedSelector = vmtkscripts.vmtkPointListSeedSelector()
            source_endpoints, target_endpoints = self.read_fcsv(fcsv_file, single_list=True)

            SeedSelector.SourcePoints = source_endpoints
            SeedSelector.TargetPoints = target_endpoints
            SeedSelector.SetSurface(tubePolyData)
            SeedSelector.Execute()
            sourceIdList = SeedSelector.GetSourceSeedIds()
            targetIdList = SeedSelector.GetTargetSeedIds()
            
        print('checking')
        #slicer.tubePolyData = tubePolyData

        centerlineFilter = vtkvmtkComputationalGeometry.vtkvmtkPolyDataCenterlines()
        centerlineFilter.SetInputData(tubePolyData)
        centerlineFilter.SetSourceSeedIds(sourceIdList)
        centerlineFilter.SetTargetSeedIds(targetIdList)

        centerlineFilter.SetRadiusArrayName(self.RadiusArrayName)
        centerlineFilter.SetCostFunction('1/R')  # this makes path search prefer go through points with large radius
        centerlineFilter.SetFlipNormals(False)
        centerlineFilter.SetAppendEndPointsToCenterlines(0)

        # Voronoi smoothing slightly improves connectivity
        # Unfortunately, Voronoi smoothing is broken if VMTK is used with VTK9, therefore
        # disable this feature for now (https://github.com/vmtk/SlicerExtension-VMTK/issues/34)
        #enableVoronoiSmoothing = (slicer.app.majorVersion * 100 + slicer.app.minorVersion < 413)
        enableVoronoiSmoothing = 1
        centerlineFilter.SetSimplifyVoronoi(enableVoronoiSmoothing)

        if self.StopFastMarchingOnReachingTarget == True:
            if targetIdList.GetNumberOfIds() != 1:
                self.PrintError('Parameter Conflict: cannot enable "StopFastMarchingOnReachingTarget" when there is more then one target seed set.')
            else:
                centerlineFilter.SetStopFastMarchingOnReachingTarget(self.StopFastMarchingOnReachingTarget)


        centerlineFilter.SetCenterlineResampling(0)
        centerlineFilter.SetResamplingStepLength(curveSamplingDistance)
        centerlineFilter.Update()

        centerlinePolyData = vtk.vtkPolyData()
        centerlinePolyData.DeepCopy(centerlineFilter.GetOutput())

        voronoiDiagramPolyData = vtk.vtkPolyData()
        voronoiDiagramPolyData.DeepCopy(centerlineFilter.GetVoronoiDiagram())

        #logging.debug("End of Centerline Computation..")
        return centerlinePolyData, voronoiDiagramPolyData

    def preprocess(self, surfacePolyData, targetNumberOfPoints, decimationAggressiveness, subdivide):
        # import the vmtk libraries
        try:
            import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
            import vtkvmtkMiscPython as vtkvmtkMisc
        except ImportError:
            raise ImportError("VMTK library is not found")

        numberOfInputPoints = surfacePolyData.GetNumberOfPoints()
        if numberOfInputPoints == 0:
            raise("Input surface model is empty")
        reductionFactor = (numberOfInputPoints-targetNumberOfPoints) / numberOfInputPoints
        if reductionFactor > 0.0:
            parameters = {}
            inputSurfaceModelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "tempInputSurfaceModel")
            inputSurfaceModelNode.SetAndObserveMesh(surfacePolyData)
            parameters["inputModel"] = inputSurfaceModelNode
            outputSurfaceModelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "tempDecimatedSurfaceModel")
            parameters["outputModel"] = outputSurfaceModelNode
            parameters["reductionFactor"] = reductionFactor
            parameters["method"] = "FastQuadric"
            parameters["aggressiveness"] = decimationAggressiveness
            decimation = slicer.modules.decimation
            cliNode = slicer.cli.runSync(decimation, None, parameters)
            surfacePolyData = outputSurfaceModelNode.GetPolyData()
            slicer.mrmlScene.RemoveNode(inputSurfaceModelNode)
            slicer.mrmlScene.RemoveNode(outputSurfaceModelNode)
            slicer.mrmlScene.RemoveNode(cliNode)

        surfaceCleaner = vtk.vtkCleanPolyData()
        surfaceCleaner.SetInputData(surfacePolyData)
        surfaceCleaner.Update()

        surfaceTriangulator = vtk.vtkTriangleFilter()
        surfaceTriangulator.SetInputData(surfaceCleaner.GetOutput())
        surfaceTriangulator.PassLinesOff()
        surfaceTriangulator.PassVertsOff()
        surfaceTriangulator.Update()

        # new steps for preparation to avoid problems because of slim models (f.e. at stenosis)
        if subdivide:
            subdiv = vtk.vtkLinearSubdivisionFilter()
            subdiv.SetInputData(surfaceTriangulator.GetOutput())
            subdiv.SetNumberOfSubdivisions(1)
            subdiv.Update()
            if subdiv.GetOutput().GetNumberOfPoints() == 0:
                logging.warning("Mesh subdivision failed. Skip subdivision step.")
                subdivide = False

        normals = vtk.vtkPolyDataNormals()
        if subdivide:
            normals.SetInputData(subdiv.GetOutput())
        else:
            normals.SetInputData(surfaceTriangulator.GetOutput())
        normals.SetAutoOrientNormals(1)
        normals.SetFlipNormals(0)
        normals.SetConsistency(1)
        normals.SplittingOff()
        normals.Update()

        return normals.GetOutput()
    
    def extractNetwork(self, surfacePolyData, computeGeometry=False):
        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc

        # Decimate
        # It seems that decimation at this stage is not necessary (decimation in preprocessing is enough).
        # By not decimating here, we can keep th network and centerline extraction results more similar.
        # If network extraction is too slow then one can experiment with this flag.
        decimate = False
        if decimate:
            decimationFilter = vtk.vtkDecimatePro()
            decimationFilter.SetInputData(surfacePolyData)
            decimationFilter.SetTargetReduction(0.99)
            decimationFilter.SetBoundaryVertexDeletion(0)
            decimationFilter.PreserveTopologyOn()
            decimationFilter.Update()

        # Clean and triangulate
        cleaner = vtk.vtkCleanPolyData()
        if decimate:
            cleaner.SetInputData(decimationFilter.GetOutput())
        else:
            cleaner.SetInputData(surfacePolyData)
        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(cleaner.GetOutputPort())
        triangleFilter.Update()
        simplifiedPolyData = triangleFilter.GetOutput()

        # Cut hole at start position
        #if endPointsMarkupsNode and endPointsMarkupsNode.GetNumberOfControlPoints() > 0:
        #    startPosition = [0, 0, 0]
        #    endPointsMarkupsNode.GetNthControlPointPosition(
        #        self.startPointIndexFromEndPointsMarkupsNode(endPointsMarkupsNode), startPosition)
        #else:
        #    # If no endpoints are specific then use the closest point to a corner
        #    bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #    simplifiedPolyData.GetBounds(bounds)
        #    startPosition = [bounds[0], bounds[2], bounds[4]]
        # Cut hole at start position

        # If no endpoints are specific then use the closest point to a corner
        bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        simplifiedPolyData.GetBounds(bounds)
        startPosition = [bounds[0], bounds[2], bounds[4]]

        self.openSurfaceAtPoint(simplifiedPolyData, startPosition)

        # Extract network
        networkExtraction = vtkvmtkMisc.vtkvmtkPolyDataNetworkExtraction()
        networkExtraction.SetInputData(simplifiedPolyData)
        networkExtraction.SetAdvancementRatio(1.05)
        networkExtraction.SetRadiusArrayName(self.RadiusArrayName)
        networkExtraction.SetTopologyArrayName(self.topologyArrayName)
        networkExtraction.SetMarksArrayName(self.marksArrayName)
        networkExtraction.Update()

        if computeGeometry:
            centerlineGeometry = vtkvmtkComputationalGeometry.vtkvmtkCenterlineGeometry()
            centerlineGeometry.SetInputData(networkExtraction.GetOutput())
            centerlineGeometry.SetLengthArrayName(self.lengthArrayName)
            centerlineGeometry.SetCurvatureArrayName(self.curvatureArrayName)
            centerlineGeometry.SetTorsionArrayName(self.torsionArrayName)
            centerlineGeometry.SetTortuosityArrayName(self.tortuosityArrayName)
            centerlineGeometry.SetFrenetTangentArrayName(self.frenetTangentArrayName)
            centerlineGeometry.SetFrenetNormalArrayName(self.frenetNormalArrayName)
            centerlineGeometry.SetFrenetBinormalArrayName(self.frenetBinormalArrayName)
            # centerlineGeometry.SetLineSmoothing(0)
            # centerlineGeometry.SetOutputSmoothedLines(0)
            # centerlineGeometry.SetNumberOfSmoothingIterations(100)
            # centerlineGeometry.SetSmoothingFactor(0.1)
            centerlineGeometry.Update()
            return centerlineGeometry.GetOutput()
        else:
            return networkExtraction.GetOutput()

    def extract_network_skeletons(self, surfacePolyData, nii_path, computeGeometry = False):
        #surfacePolyData - preprocessed mesh!!!!
        
        
        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc
    
        # Decimate
        # It seems that decimation at this stage is not necessary (decimation in preprocessing is enough).
        # By not decimating here, we can keep th network and centerline extraction results more similar.
        # If network extraction is too slow then one can experiment with this flag.
        decimate = False
        if decimate:
            decimationFilter = vtk.vtkDecimatePro()
            decimationFilter.SetInputData(surfacePolyData)
            decimationFilter.SetTargetReduction(0.99)
            decimationFilter.SetBoundaryVertexDeletion(0)
            decimationFilter.PreserveTopologyOn()
            decimationFilter.Update()
    
        # Clean and triangulate
        cleaner = vtk.vtkCleanPolyData()
        if decimate:
            cleaner.SetInputData(decimationFilter.GetOutput())
        else:
            cleaner.SetInputData(surfacePolyData)
        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(cleaner.GetOutputPort())
        triangleFilter.Update()
        simplifiedPolyData = triangleFilter.GetOutput()
    
        # If no endpoints are specific then use the closest point to a corner
        #bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #simplifiedPolyData.GetBounds(bounds)
        #startPosition = [bounds[0], bounds[2], bounds[4]]
    
        # Get start position using max radius of skeletons
        from scipy.ndimage import distance_transform_edt
        from skimage.morphology import skeletonize
        
        referenceVolume = slicer.util.loadVolume(nii_path)
        lbl = nib.load(nii_path)
        label = lbl.get_fdata()
        skeleton = skeletonize(label)
        edt = distance_transform_edt(label)
        radius_map = np.zeros_like(label, dtype=float)
        radius_map[skeleton] = edt[skeleton]
        flat_index = np.argmax(radius_map)
        max_index_tuple = np.unravel_index(flat_index, radius_map.shape)
        max_value = radius_map[max_index_tuple]
        startPosition = convert_voxel_to_RAS(max_index_tuple, referenceVolume)
        
    
        self.openSurfaceAtPoint(simplifiedPolyData, startPosition)
    
        # Extract network
        networkExtraction = vtkvmtkMisc.vtkvmtkPolyDataNetworkExtraction()
        networkExtraction.SetInputData(simplifiedPolyData)
        networkExtraction.SetAdvancementRatio(1.05)
        networkExtraction.SetRadiusArrayName(self.RadiusArrayName)
        networkExtraction.SetTopologyArrayName(self.topologyArrayName)
        networkExtraction.SetMarksArrayName(self.marksArrayName)
        networkExtraction.Update()
    
        if computeGeometry:
            centerlineGeometry = vtkvmtkComputationalGeometry.vtkvmtkCenterlineGeometry()
            centerlineGeometry.SetInputData(networkExtraction.GetOutput())
            centerlineGeometry.SetLengthArrayName(self.lengthArrayName)
            centerlineGeometry.SetCurvatureArrayName(self.curvatureArrayName)
            centerlineGeometry.SetTorsionArrayName(self.torsionArrayName)
            centerlineGeometry.SetTortuosityArrayName(self.tortuosityArrayName)
            centerlineGeometry.SetFrenetTangentArrayName(self.frenetTangentArrayName)
            centerlineGeometry.SetFrenetNormalArrayName(self.frenetNormalArrayName)
            centerlineGeometry.SetFrenetBinormalArrayName(self.frenetBinormalArrayName)
            # centerlineGeometry.SetLineSmoothing(0)
            # centerlineGeometry.SetOutputSmoothedLines(0)
            # centerlineGeometry.SetNumberOfSmoothingIterations(100)
            # centerlineGeometry.SetSmoothingFactor(0.1)
            centerlineGeometry.Update()
            return centerlineGeometry.GetOutput()
        else:
            return networkExtraction.GetOutput()
        
    def openSurfaceAtPoint(self, polyData, holePosition=None, holePointIndex=None):
        '''
        Modifies the polyData by cutting a hole at the given position.
        '''

        if holePointIndex is None:
            pointLocator = vtk.vtkPointLocator()
            pointLocator.SetDataSet(polyData)
            pointLocator.BuildLocator()
            # find the closest point to the desired hole position
            holePointIndex = pointLocator.FindClosestPoint(holePosition)

        if holePointIndex < 0:
            # Calling GetPoint(-1) would crash the application
            raise ValueError("openSurfaceAtPoint failed: empty input polydata")

        # Tell the polydata to build 'upward' links from points to cells
        polyData.BuildLinks()
        # Mark cells as deleted
        cellIds = vtk.vtkIdList()
        polyData.GetPointCells(holePointIndex, cellIds)
        removeFirstCell = True
        if removeFirstCell:
            # remove first cell only (smaller hole)
            if cellIds.GetNumberOfIds() > 0:
                polyData.DeleteCell(cellIds.GetId(0))
                polyData.RemoveDeletedCells()
        else:
            # remove all cells
            for cellIdIndex in range(cellIds.GetNumberOfIds()):
                polyData.DeleteCell(cellIds.GetId(cellIdIndex))
            polyData.RemoveDeletedCells()

    def branchExtractor(self, polyData):
        pass

    def getEndPoints(self, inputNetworkPolyData, startPointPosition):
        '''
        Clips the surfacePolyData on the endpoints identified using the networkPolyData.
        If startPointPosition is specified then start point will be the closest point to that position.
        Returns list of endpoint positions. Largest radius point is be the first in the list.
        '''

        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc


        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(inputNetworkPolyData)
        cleaner.Update()
        network = cleaner.GetOutput()
        network.BuildCells()
        network.BuildLinks(0)

        networkPoints = network.GetPoints()
        radiusArray = network.GetPointData().GetArray(self.RadiusArrayName)

        startPointId = -1
        maxRadius = 0
        minDistance2 = 0

        endpointIds = vtk.vtkIdList()
        for i in range(network.GetNumberOfCells()):
            numberOfCellPoints = network.GetCell(i).GetNumberOfPoints()
            if numberOfCellPoints < 2:
                continue

            for pointIndex in [0, numberOfCellPoints - 1]:
                pointId = network.GetCell(i).GetPointId(pointIndex)
                pointCells = vtk.vtkIdList()
                network.GetPointCells(pointId, pointCells)
                if pointCells.GetNumberOfIds() == 1:
                    endpointIds.InsertUniqueId(pointId)
                    if startPointPosition is not None:
                        # find start point based on position
                        position = networkPoints.GetPoint(pointId)
                        distance2 = vtk.vtkMath.Distance2BetweenPoints(position, startPointPosition)
                        if startPointId < 0 or distance2 < minDistance2:
                            minDistance2 = distance2
                            startPointId = pointId
                    else:
                        # find start point based on radius
                        radius = radiusArray.GetValue(pointId)
                        if startPointId < 0 or radius > maxRadius:
                            maxRadius = radius
                            startPointId = pointId

        endpointPositions = []
        numberOfEndpointIds = endpointIds.GetNumberOfIds()
        if numberOfEndpointIds == 0:
            return endpointPositions
        # add the largest radius point first
        endpointPositions.append(networkPoints.GetPoint(startPointId))
        # add all the other points
        for pointIdIndex in range(numberOfEndpointIds):
            pointId = endpointIds.GetId(pointIdIndex)
            if pointId == startPointId:
                # already added
                continue
            endpointPositions.append(networkPoints.GetPoint(pointId))

        return endpointPositions



