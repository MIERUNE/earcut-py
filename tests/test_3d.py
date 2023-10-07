import numpy as np
import pytest

from earcut import earcut
from earcut.utils_3d import project3d_to_2d


@pytest.fixture
def polygons() -> list[list[np.ndarray]]:
    polygons = [
        [
            np.array(
                [
                    [0, 0, 0],
                    [0, 4, 0],
                    [4, 4, 0],
                    [4, 0, 0],
                ]
            ),
            np.array(
                [
                    [1, 1, 0],
                    [1, 2, 0],
                    [2, 2, 0],
                    [2, 1, 0],
                ]
            ),
            np.array(
                [
                    [2, 2, 0],
                    [2, 3, 0],
                    [3, 3, 0],
                    [3, 2, 0],
                ]
            ),
        ],
        [
            np.array(
                [
                    [0, 0, 0],
                    [0, 0, 4],
                    [0, 4, 4],
                    [0, 4, 0],
                ]
            ),
            np.array(
                [
                    [0, 1, 1],
                    [0, 1, 3],
                    [0, 3, 3],
                    [0, 3, 1],
                ]
            ),
        ],
        [
            np.array(
                [
                    [0, 4, 0],
                    [0, 4, 4],
                    [4, 4, 4],
                    [4, 4, 0],
                ]
            ),
            np.array(
                [
                    [3, 4, 2],
                    [2, 4, 3],
                    [1, 4, 2],
                    [2, 4, 1],
                ]
            ),
        ],
        [
            np.array(
                [
                    [0, 4, 0],
                    [0, 4, 4],
                    [4, 4, 4],
                    [4, 4, 0],
                ]
            ),
            np.array(
                [
                    [3, 4, 2],
                    [2, 4, 3],
                    [1, 4, 2],
                    [2, 4, 1],
                ]
            ),
        ],
        [
            np.array(
                [
                    [4, 0, 0],
                    [4, 4, 0],
                    [7, 4, -3],
                    [7, 0, -3],
                ],
            ),
            np.array(
                [
                    [5, 1, -1],
                    [5, 3, -1],
                    [6, 3, -2],
                    [6, 1, -2],
                ],
            ),
        ],
    ]

    ix = np.linspace(0, 2 * np.pi, 150)
    u = 4 + 0.5 * np.sin(ix * 17)
    x = 4 + u * np.cos(ix)
    z = u * np.sin(ix)
    y = 6 * np.ones(150)
    ix = np.linspace(0, 2 * np.pi, 13)
    in_x = 4 + 2 * np.cos(ix)
    in_z = 2 * np.sin(ix)
    in_y = 6 * np.ones(13)

    polygons.append(
        [
            np.vstack((x, y, z)).T,
            np.vstack((in_x, in_y, in_z)).T,
        ]
    )

    return polygons


def test_triangulate(polygons):
    # Flatten
    for polygon in polygons:
        holeIndices = []
        if len(polygon) > 1:
            hi = polygon[0].shape[0]
            for ring in polygon[1:]:
                holeIndices.append(hi)
                hi += ring.shape[0]
        vertices = np.vstack(polygon)
        flatten_vertices = vertices.flatten()

        # Earcut
        flatten_vertices = project3d_to_2d(flatten_vertices, len(polygon[0]))
        assert flatten_vertices is not None
        if triangles := earcut(flatten_vertices, holeIndices, dim=2):
            triangles = np.asarray(triangles).reshape(-1, 3)
            assert vertices[triangles].shape == (len(triangles), 3, 3)
