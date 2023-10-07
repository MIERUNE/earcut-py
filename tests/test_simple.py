from earcut import deviation, earcut, flatten


def test_flatten():
    data, holeIndices, dim = flatten([[[0, 0], [100, 0], [100, 100], [0, 100]]])
    triangles = earcut(data, holeIndices, dim)
    assert len(triangles) == 2 * 3
    deviation(data, holeIndices, dim, triangles)


def test_empty():
    data = []
    holeIndices = []
    triangles = earcut(data, holeIndices, 2)
    assert len(triangles) == 0
    deviation(data, holeIndices, 2, triangles)


def test_invalid_point():
    data = [100, 200]
    holeIndices = []
    triangles = earcut(data, holeIndices, 2)
    assert len(triangles) == 0
    deviation(data, holeIndices, 2, triangles)


def test_invalid_line():
    data = [0, 0, 100, 200]
    holeIndices = []
    triangles = earcut(data, holeIndices, 2)
    assert len(triangles) == 0
    deviation(data, holeIndices, 2, triangles)


def test_invalid_empty_hole():
    data = [0, 0, 100, 0, 100, 100]
    holeIndices = [3, 3]
    triangles = earcut(data, holeIndices, 2)
    assert len(triangles) == 1 * 3
    # deviation(data, holeIndices, 2, triangles)


def test_invalid_point_hole():
    data = [0, 0, 100, 0, 100, 100, 50, 30]
    holeIndices = [3]
    triangles = earcut(data, holeIndices, 2)
    assert len(triangles) == 3 * 3
    deviation(data, holeIndices, 2, triangles)


def test_invalid_line_hole():
    data = [0, 0, 100, 0, 100, 100, 50, 30, 60, 30]
    holeIndices = [3]
    triangles = earcut(data, holeIndices, 2)
    assert len(triangles) == 5 * 3
    deviation(data, holeIndices, 2, triangles)
