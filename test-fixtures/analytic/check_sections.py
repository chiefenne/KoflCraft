"""Read a generated sweep bundle with ParaView's pvpython and save two polygonal sections."""

import json
import sys
from pathlib import Path

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet, vtkPlane
from vtkmodules.vtkFiltersCore import vtkAppendFilter, vtkCutter
from vtkmodules.vtkIOXML import vtkXMLMultiBlockDataReader, vtkXMLUnstructuredGridWriter


def cells(grid):
    return tuple(vtk_to_numpy(a) for a in (
        grid.GetCells().GetOffsetsArray(), grid.GetCells().GetConnectivityArray()
    ))


def read_sizes(grid, axis):
    # Reconstruct each volume from its oriented base and axial displacement.
    # Check connectivity at the remaining stations against the two ends of every prism.
    points = vtk_to_numpy(grid.GetPoints().GetData())
    offsets, ids = cells(grid)
    transverse = [(axis+1) % 3, (axis+2) % 3]
    stations = np.unique(points[:, axis])
    station_count = len(stations) - 1
    count = grid.GetNumberOfCells() // station_count
    nodes = len(points) // len(stations)
    assert grid.GetNumberOfCells() == station_count * count
    assert np.allclose(points.reshape(len(stations), nodes, 3)[:, :, axis], stations[:, None])
    assert np.allclose(points.reshape(len(stations), nodes, 3)[:, :, transverse], points[:nodes, transverse])
    strides = np.diff(offsets).reshape(station_count, count)
    assert np.all(strides == strides[0])
    per_station = int(offsets[count])
    base_ids = ids[:per_station]
    for k, station_ids in enumerate(ids.reshape(station_count, per_station)):
        assert np.array_equal(station_ids - k*nodes, base_ids)
    shared = {}
    pairs = []
    areas = []
    for i in range(count):
        vertices = base_ids[offsets[i]:offsets[i+1]]
        bottom, top = np.split(vertices, 2)
        assert np.array_equal(top-nodes, bottom)
        xy = points[bottom][:, transverse] - points[bottom[0], transverse]
        areas.append(0.5 * np.sum(xy[:, 0]*np.roll(xy[:, 1], -1) - xy[:, 1]*np.roll(xy[:, 0], -1)))
        for a, b in zip(bottom, np.roll(bottom, -1)):
            opposite = shared.pop((b, a), None)
            if opposite is None:
                assert (a, b) not in shared
                shared[(a, b)] = i
            else:
                pairs.append((i, opposite))
    volumes = np.diff(stations)[:, None] * np.array(areas)
    assert np.all(volumes > 0)
    h = np.cbrt(volumes)
    ratios = np.ones_like(h)
    for i, j in pairs:
        r = np.maximum(h[:, i]/h[:, j], h[:, j]/h[:, i])
        ratios[:, i] = np.maximum(ratios[:, i], r)
        ratios[:, j] = np.maximum(ratios[:, j], r)
    r = np.maximum(h[1:]/h[:-1], h[:-1]/h[1:])
    ratios[1:] = np.maximum(ratios[1:], r)
    ratios[:-1] = np.maximum(ratios[:-1], r)
    for name, values in [('CellSizeMm', h), ('MaxNeighbourSizeRatio', ratios)]:
        stored = vtk_to_numpy(grid.GetCellData().GetArray(name))
        assert np.allclose(stored, values.ravel(), rtol=1e-10, atol=0)
    return volumes, ratios, stations


def check_polyhedron_faces(grid, expected):
    # Integrate tetrahedra from each VTK_POLYHEDRON's actual outward face stream.
    refs = grid.GetPolyhedronFaceLocations()
    faces = grid.GetPolyhedronFaces()
    cell_faces = vtk_to_numpy(refs.GetOffsetsArray())
    face_id = vtk_to_numpy(refs.GetConnectivityArray())
    face_owner = np.empty(len(face_id), dtype=np.int64)
    face_owner[face_id] = np.repeat(np.arange(len(expected)), np.diff(cell_faces))
    offsets = vtk_to_numpy(faces.GetOffsetsArray())
    ids = vtk_to_numpy(faces.GetConnectivityArray())
    coords = vtk_to_numpy(grid.GetPoints().GetData())
    starts, nodes = cells(grid)
    origin = coords[nodes[starts[:-1]]]
    result = np.zeros(len(expected))
    for n in np.unique(np.diff(offsets)):
        matches = np.flatnonzero(np.diff(offsets) == n)
        for begin in range(0, len(matches), 30000):
            faces_here = matches[begin:begin+30000]
            owner = face_owner[faces_here]
            p = coords[ids[offsets[faces_here, None] + np.arange(n)]]
            for j in range(1, n-1):
                area = np.cross(p[:, j]-p[:, 0], p[:, j+1]-p[:, 0])
                part = np.einsum('ij,ij->i', p[:, 0]-origin[owner], area) / 6
                np.add.at(result, owner, part)
    poly = vtk_to_numpy(grid.GetCellTypes()) == 42
    assert np.allclose(result[poly], expected[poly], atol=1e-12, rtol=1e-9)


def write_section(grid, axis, fraction, dest):
    bounds = grid.GetBounds()
    origin = [(bounds[2*i]+bounds[2*i+1])/2 for i in range(3)]
    origin[axis] = bounds[2*axis] + fraction*(bounds[2*axis+1]-bounds[2*axis])
    plane = vtkPlane()
    plane.SetOrigin(origin)
    plane.SetNormal(np.eye(3)[axis])
    cut = vtkCutter()
    cut.SetCutFunction(plane)
    cut.SetInputData(grid)
    cut.GenerateTrianglesOff()
    convert = vtkAppendFilter()
    convert.AddInputConnection(cut.GetOutputPort())
    convert.Update()
    result = convert.GetOutput()
    assert result.GetNumberOfCells() and result.GetCellData().GetArray("LayerIndex")
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetInputData(result)
    writer.SetFileName(str(dest))
    if not writer.Write():
        raise OSError(dest)


def inspect(case):
    report = json.loads((case/'report.json').read_text())
    config = report['config']
    shape = config['geometry']
    blind = shape['kind'] == 'block'
    axis = 'xyz'.index(shape['inlet'][0]) if blind else 0
    reader = vtkXMLMultiBlockDataReader()
    reader.SetFileName(str(case/'mesh.vtm'))
    reader.Update()
    dataset = reader.GetOutput()
    blocks = {dataset.GetMetaData(i).Get(vtkCompositeDataSet.NAME()): dataset.GetBlock(i)
              for i in range(dataset.GetNumberOfBlocks())}
    assert set(blocks) == ({'volume','inlet','walls'} if blind else {'volume','inlet','outlet','walls'})
    grid = blocks['volume']
    volumes, ratios, stations = read_sizes(grid, axis)
    check_polyhedron_faces(grid, volumes.ravel())
    assert np.isclose(volumes.sum(), report['volume_mm3'], rtol=1e-10)
    assert ratios.max() <= config['max_neighbour_size_ratio'] + 1e-10
    assert np.isclose(ratios.max(), report['neighbour_size_ratio']['max'], rtol=1e-10)
    layers = vtk_to_numpy(grid.GetCellData().GetArray('LayerIndex'))
    assert np.count_nonzero(layers == 0) == report['core_cells']
    assert set(vtk_to_numpy(grid.GetCellTypes())[layers == 0]) == {42}
    for name, surfaces in blocks.items():
        if name == 'volume':
            continue
        owner = vtk_to_numpy(surfaces.GetCellData().GetArray('OwnerCellId'))
        surface_layers = vtk_to_numpy(surfaces.GetCellData().GetArray('LayerIndex'))
        assert np.array_equal(layers[owner], surface_layers)
        assert np.all(vtk_to_numpy(surfaces.GetCellData().GetArray('PatchId')) == {'inlet':1,'outlet':2,'walls':3}[name])
        xyz = vtk_to_numpy(surfaces.GetPoints().GetData())
        if name == 'walls':
            assert set(surface_layers) == {1}
            if blind:
                offsets, ids = cells(surfaces)
                end = np.isclose(xyz[:, axis], stations[-1])
                cap_faces = np.logical_and.reduceat(end[ids], offsets[:-1])
                assert cap_faces.sum() == volumes.shape[1]
        else:
            plane = stations[0] if name == 'inlet' else stations[-1]
            assert np.allclose(xyz[:, axis], plane)
            assert set(surface_layers) == set(range(config['wall_layers'] + 1))
            offsets, ids = cells(surfaces)
            areas = np.zeros(len(owner))
            for n in np.unique(np.diff(offsets)):
                use = np.where(np.diff(offsets) == n)[0]
                poly = xyz[ids[offsets[use, None] + np.arange(n)]]
                poly -= poly[:, :1].copy()
                area_vector = 0.5*np.cross(poly, np.roll(poly, -1, axis=1)).sum(axis=1)
                areas[use] = area_vector[:, axis] * (-1 if name == 'inlet' else 1)
            assert np.all(areas > 0)
            assert np.allclose(np.bincount(surface_layers, weights=areas), report[name+'_area_by_layer_mm2'])
    write_section(grid, axis, 0.2713, case/'cross-section.vtu')
    write_section(grid, (axis+2) % 3, 0.5031, case/'longitudinal-section.vtu')
    print(case, 'cells', grid.GetNumberOfCells(), 'max size ratio', ratios.max(), 'volume mm3', volumes.sum())


if __name__ == '__main__':
    for arg in sys.argv[1:]:
        inspect(Path(arg))
