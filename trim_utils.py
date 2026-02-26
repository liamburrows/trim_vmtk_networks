import csv
import vtk
import slicer
import numpy as np
import nibabel as nib


def convert_voxel_to_RAS(voxel_point, referenceVolume):
    # Input: voxel_point = (i,j,k) tuple
    # Output: RAS (world) coordinates as (x,y,z)
    
    ijkToRAS = vtk.vtkMatrix4x4()
    referenceVolume.GetIJKToRASMatrix(ijkToRAS)
    
    voxel = list(voxel_point) + [1.0]
    ras = [0.0, 0.0, 0.0, 1.0]
    ijkToRAS.MultiplyPoint(voxel, ras)
    
    return ras[:3]

def trim_network(network, n, source_endpoint):
    
    pointLocator = vtk.vtkPointLocator()
    pointLocator.SetDataSet(network)
    pointLocator.BuildLocator()
    #print(source_endpoint)
    sourceId = pointLocator.FindClosestPoint(source_endpoint)
    print(sourceId)
    points = network.GetPoints()
    source_coord = points.GetPoint(sourceId)

    trimmed_net = reconstruct_network_with_point_data(network, sourceId, n)

    return trimmed_net, source_coord

def find_current_cell(polydata, point_id):
    cell_ids = vtk.vtkIdList()
    polydata.GetPointCells(point_id, cell_ids)
    if cell_ids.GetNumberOfIds() > 1:
        raise ValueError('There shouldnt be multiple cells attached to a single point id in vmtk networks?')
    else:
        cell_id = cell_ids.GetId(0)
        return cell_id
    

def find_connected_cells_via_coordinates(polydata, point_id):
    """Finds all cells connected to a given point by checking the point coordinates."""
    connected_cells = []
    source_point = polydata.GetPoint(point_id)  # Get the coordinates of the source point

    for i in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(i)
        point_ids = cell.GetPointIds()

        # Loop through each point in the cell and check if the coordinates match the source point
        for j in range(point_ids.GetNumberOfIds()):
            current_point_id = point_ids.GetId(j)
            current_point = polydata.GetPoint(current_point_id)
            
            # Compare the coordinates of the current point with the source point
            if current_point == source_point:
                connected_cells.append(i)
                break  # Stop once a match is found for this cell

    return connected_cells

def copy_cell_with_point_data(polydata, cell_id, output_polydata):
    """Copy a cell (branch) and its point data arrays from one polydata to another."""
    cell = polydata.GetCell(cell_id)
    point_ids = cell.GetPointIds()

    # Add points and create the line for the cell in output_polydata
    new_line = vtk.vtkPolyLine()
    new_line.GetPointIds().SetNumberOfIds(point_ids.GetNumberOfIds())
    
    for i in range(point_ids.GetNumberOfIds()):
        pt_id = point_ids.GetId(i)
        point = polydata.GetPoint(pt_id)
        new_pt_id = output_polydata.GetPoints().InsertNextPoint(point)
        new_line.GetPointIds().SetId(i, new_pt_id)
        

        # Now copy all point data arrays (like radius) to the new polydata
        num_point_arrays = polydata.GetPointData().GetNumberOfArrays()
        for array_idx in range(num_point_arrays):
            point_array = polydata.GetPointData().GetArray(array_idx)
            value = point_array.GetTuple(pt_id)
            output_polydata.GetPointData().GetArray(array_idx).InsertNextTuple(value)

    # Add the new cell to the output polydata
    output_polydata.InsertNextCell(new_line.GetCellType(), new_line.GetPointIds())

    # Now copy all cell data arrays (associated with this specific cell) to the new polydata
    num_cell_arrays = polydata.GetCellData().GetNumberOfArrays()
    for array_idx in range(num_cell_arrays):
        cell_array = polydata.GetCellData().GetArray(array_idx)
        value = cell_array.GetTuple(cell_id)  # Get the cell data for this cell
        output_polydata.GetCellData().GetArray(array_idx).InsertNextTuple(value)

def reconstruct_network_with_point_data(polydata, source_point_id, max_bifurcations):
    """Reconstruct the network up to the nth bifurcation, including point data arrays."""
    # Create a new vtkPolyData to store the result
    output_polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    output_polydata.SetPoints(points)
    lines = vtk.vtkCellArray()
    output_polydata.SetLines(lines)

    # Copy point data arrays from original polydata to the output polydata
    point_data = polydata.GetPointData()
    num_point_arrays = point_data.GetNumberOfArrays()

    # Initialise arrays in output_polydata with the same structure as in the input
    for i in range(num_point_arrays):
        array = point_data.GetArray(i)
        new_array = array.NewInstance()  # Create a new instance of the array
        new_array.SetName(array.GetName())  # Set the same name
        new_array.SetNumberOfComponents(array.GetNumberOfComponents())  # Same number of components
        output_polydata.GetPointData().AddArray(new_array)  # Add it to the output polydata

    cell_data = polydata.GetCellData()
    num_cell_arrays = cell_data.GetNumberOfArrays()
    for i in range(num_cell_arrays):
        array = cell_data.GetArray(i)
        new_array = array.NewInstance()
        new_array.SetName(array.GetName())
        new_array.SetNumberOfComponents(array.GetNumberOfComponents())
        output_polydata.GetCellData().AddArray(new_array)
    
    # BFS-like approach to traverse the network up to n bifurcations
    to_visit = [(source_point_id, 0)]  # (point_id, bifurcation_depth)
    visited_cells = set()

    source_cell = find_current_cell(polydata, source_point_id)
    visited_cells.add(source_cell)
    source_dict = { 'cell' : source_cell, 'gen' : 0, 'parent' : np.nan, 'children' : [], 'side' : 'MPA' }
    dicts = [source_dict]

    copy_cell_with_point_data(polydata, source_cell, output_polydata)
    cell = polydata.GetCell(source_cell)
    point_ids = cell.GetPointIds()
    
    # Look for points at the end of the cell that aren't the current point
    for i in range(point_ids.GetNumberOfIds()):
        next_point_id = point_ids.GetId(i)
        if next_point_id != source_point_id:  # Move to the other end of the cell
            to_visit.append((next_point_id, 0))

    while to_visit:
        point_id, depth = to_visit.pop(0)

        current_cell_id = find_current_cell(polydata, point_id)
        for dict in dicts:
            if dict['cell'] == current_cell_id:
                current_dict = dict
                break
        

        # Find all connected cells (branches) by matching coordinates
        connected_cells = find_connected_cells_via_coordinates(polydata, point_id)

        for cell_id in connected_cells:
            if cell_id not in visited_cells and depth < max_bifurcations:
                current_dict['children'].append(cell_id)
                new_dict = { 'cell' : cell_id, 'gen' : depth + 1, 'parent' : current_dict['cell'], 'children' : [] }
                
                dicts.append(new_dict)
                
                visited_cells.add(cell_id)

                # Copy the current cell and its point data to the new polydata
                copy_cell_with_point_data(polydata, cell_id, output_polydata)

                # Find the points at the other end of the cell and continue traversal
                cell = polydata.GetCell(cell_id)
                point_ids = cell.GetPointIds()
                
                # Look for points at the end of the cell that aren't the current point
                for i in range(point_ids.GetNumberOfIds()):
                    next_point_id = point_ids.GetId(i)
                    if next_point_id != point_id:  # Move to the other end of the cell
                        to_visit.append((next_point_id, depth + 1))


    return output_polydata, dicts

def calc_vol_surf(polydata, cell_id):
    point_data = polydata.GetPointData()
    radius_array = point_data.GetArray('Radius')
    points = polydata.GetPoints()
    
    cell_data = polydata.GetCellData()
    
    cell = polydata.GetCell(cell_id)
    point_ids = cell.GetPointIds()

    vol = 0.0 
    surf = 0.0
    for i in range(point_ids.GetNumberOfIds()-1):
        point_id1, point_id2 = point_ids.GetId(i), point_ids.GetId(i+1)
        r1, r2 = radius_array.GetValue(point_id1), radius_array.GetValue(point_id2)
        coord1, coord2 = np.array(points.GetPoint(point_id1)), np.array(points.GetPoint(point_id2))


            
        h = np.linalg.norm(coord1 - coord2)

        vol_segment = (1/3)*np.pi*h*(r1**2 + r1*r2 + r2**2)
        # (lateral) surface area
        s = np.sqrt(h**2 + (r1-r2)**2)
        a_segment = np.pi * (r1 + r2) * s

        if np.isnan(vol_segment):
            vol_segment = 0
        if np.isnan(a_segment):
            a_segment = 0

        vol += vol_segment
        surf += a_segment


    return vol, surf

def extract_single_cell_as_polydata(networkPolyData, cell_id):
    # Create a cell array and insert the desired cell
    selected_cell = networkPolyData.GetCell(cell_id)
    point_ids = selected_cell.GetPointIds()

    # Create a new set of points and a map from old to new point IDs
    new_points = vtk.vtkPoints()
    id_map = {}
    for i in range(point_ids.GetNumberOfIds()):
        old_id = point_ids.GetId(i)
        coord = networkPolyData.GetPoint(old_id)
        new_id = new_points.InsertNextPoint(coord)
        id_map[old_id] = new_id

    # Create a new polydata object
    new_polydata = vtk.vtkPolyData()
    new_polydata.SetPoints(new_points)

    # Rebuild the cell using the new point IDs
    new_cell_array = vtk.vtkCellArray()
    new_cell = vtk.vtkPolyLine()
    new_cell.GetPointIds().SetNumberOfIds(point_ids.GetNumberOfIds())

    for i in range(point_ids.GetNumberOfIds()):
        old_id = point_ids.GetId(i)
        new_id = id_map[old_id]
        new_cell.GetPointIds().SetId(i, new_id)

    new_cell_array.InsertNextCell(new_cell)
    new_polydata.SetLines(new_cell_array)

    return new_polydata

def find_root_cellID_from_network(networkPolyData):
# THIS FUNCTION WORKS BY CHECKING IF VOLUME IS BIGGEST //AND// IF THE CELL IS AN ENDPOINT - ONLY SUITABLE FOR PA TREE
    
    cell_data = networkPolyData.GetCellData()
    length_array = cell_data.GetArray('Length')
    tort_array = cell_data.GetArray('Tortuosity')
    
    point_data = networkPolyData.GetPointData()
    curv_array = point_data.GetArray('Curvature')
    torsion_array = point_data.GetArray('Torsion')
    radius_array = point_data.GetArray('Radius')
    
    points = networkPolyData.GetPoints()
    
    number_of_cells = networkPolyData.GetNumberOfCells()
    # Generate a list of cell IDs
    cells = list(range(number_of_cells))
    #print(cells)

    max_vol = -1
    for cell_id in cells:
        cell = networkPolyData.GetCell(cell_id)
        point_ids = cell.GetPointIds()

        first = point_ids.GetId(0)
        last = point_ids.GetId(point_ids.GetNumberOfIds() - 1)
        first_connected = find_connected_cells_via_coordinates(networkPolyData, first)
        last_connected = find_connected_cells_via_coordinates(networkPolyData, last)
        is_end_cell = (len(first_connected) == 1 or len(last_connected) == 1)

        if is_end_cell:
            vol, surf = calc_vol_surf(networkPolyData, cell_id)
            if vol > max_vol:
                #print(f'vol = {vol}, cell_id = {cell_id}')
                max_vol = vol
                max_cell_id = cell_id
    
    cell_id = max_cell_id
    return max_cell_id

def reconstruct_network_with_cell_data(polydata, source_cell, max_bifurcations):
    """Reconstruct the network up to the nth bifurcation, including point data arrays."""
    # Create a new vtkPolyData to store the result
    output_polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    output_polydata.SetPoints(points)
    lines = vtk.vtkCellArray()
    output_polydata.SetLines(lines)

    # Copy point data arrays from original polydata to the output polydata
    point_data = polydata.GetPointData()
    num_point_arrays = point_data.GetNumberOfArrays()

    # Initialize arrays in output_polydata with the same structure as in the input
    for i in range(num_point_arrays):
        array = point_data.GetArray(i)
        new_array = array.NewInstance()  # Create a new instance of the array
        new_array.SetName(array.GetName())  # Set the same name
        new_array.SetNumberOfComponents(array.GetNumberOfComponents())  # Same number of components
        output_polydata.GetPointData().AddArray(new_array)  # Add it to the output polydata

    cell_data = polydata.GetCellData()
    num_cell_arrays = cell_data.GetNumberOfArrays()
    for i in range(num_cell_arrays):
        array = cell_data.GetArray(i)
        new_array = array.NewInstance()
        new_array.SetName(array.GetName())
        new_array.SetNumberOfComponents(array.GetNumberOfComponents())
        output_polydata.GetCellData().AddArray(new_array)
    
    # BFS-like approach to traverse the network up to n bifurcations
    #to_visit = [(source_point_id, 0)]  # (point_id, bifurcation_depth)
    to_visit = []
    visited_cells = set()

    #source_cell = find_current_cell(polydata, source_point_id)
    visited_cells.add(source_cell)
    source_dict = { 'cell' : source_cell, 'gen' : 0, 'parent' : np.nan, 'children' : [], 'side' : 'MPA' }
    dicts = [source_dict]

    copy_cell_with_point_data(polydata, source_cell, output_polydata)
    cell = polydata.GetCell(source_cell)
    point_ids = cell.GetPointIds()
    
    # Look for points at the end of the cell that aren't the current point
    for i in range(point_ids.GetNumberOfIds()):
        next_point_id = point_ids.GetId(i)
       # if next_point_id != source_point_id:  # Move to the other end of the cell
            #to_visit.append((next_point_id, 0))
        to_visit.append((next_point_id, 0))

    while to_visit:
        point_id, depth = to_visit.pop(0)

        current_cell_id = find_current_cell(polydata, point_id)
        for dict in dicts:
            if dict['cell'] == current_cell_id:
                current_dict = dict
                break
        
        #if depth > max_bifurcations:
        #    break  # Stop once we reach the desired bifurcation depth

        # Find all connected cells (branches) by matching coordinates
        connected_cells = find_connected_cells_via_coordinates(polydata, point_id)

        for cell_id in connected_cells:
            if cell_id not in visited_cells and depth < max_bifurcations:
                current_dict['children'].append(cell_id)
                new_dict = { 'cell' : cell_id, 'gen' : depth + 1, 'parent' : current_dict['cell'], 'children' : [] }
                
                dicts.append(new_dict)
                
                visited_cells.add(cell_id)

                # Copy the current cell and its point data to the new polydata
                copy_cell_with_point_data(polydata, cell_id, output_polydata)

                # Find the points at the other end of the cell and continue traversal
                cell = polydata.GetCell(cell_id)
                point_ids = cell.GetPointIds()
                
                # Look for points at the end of the cell that aren't the current point
                for i in range(point_ids.GetNumberOfIds()):
                    next_point_id = point_ids.GetId(i)
                    if next_point_id != point_id:  # Move to the other end of the cell
                        to_visit.append((next_point_id, depth + 1))


    return output_polydata, dicts
