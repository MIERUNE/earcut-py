import math
from typing import Optional

import numpy as np


def _normal(a: np.ndarray) -> Optional[np.ndarray]:
    b = np.roll(a, -1, axis=0)
    n = np.average(np.cross(a - b, a + b), axis=0)
    d = np.linalg.norm(n)
    if d < 1e-30:
        return None
    else:
        return n / d


def _project_to_2d(normal: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    nx, ny, nz = normal
    dd = (nx**2 + ny**2) ** 0.5
    if dd < 1e-30:
        if normal[2] > 0:
            return vertices[:, :2]
        else:
            return vertices[:, (1, 0)]
    ax: float = -ny / dd
    ay: float = nx / dd
    theta = math.acos(nz)
    sint = math.sin(theta)
    cost = math.cos(theta)
    s = ax * ay * (1 - cost)
    t = ay * sint
    u = ax * sint
    R = np.asarray(
        [
            ax * ax * (1 - cost) + cost,
            s,
            s,
            ay * ay * (1 - cost) + cost,
            -t,
            u,
        ]
    ).reshape(3, 2)
    return vertices @ R


def project3d_to_2d(data, num_outer: int) -> Optional[np.ndarray]:
    d = np.asarray(data).reshape(-1, 3)
    norm = _normal(d[:num_outer])
    if norm is None:
        return None
    return _project_to_2d(norm, d).flatten()
