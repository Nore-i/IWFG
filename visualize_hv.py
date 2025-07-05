import numpy as np
import matplotlib.pyplot as plt

# ---- 1. same occupancy labelling as before -----------------------
points = np.array([[0.2, 0.8, 0.6],
                   [0.5, 0.4, 0.9],
                   [0.7, 0.2, 0.3],
                   [0.9, 0.6, 0.1],
                   [0.3, 0.5, 0.7],
                   [0.5, 0.5, 0.5]])

ind = np.argsort(points[:, -1])
points = points[ind]

res = 10
lin = np.linspace(0, 1, res)
X, Y, Z = np.meshgrid(lin, lin, lin, indexing='ij')
grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

owner = np.full(grid.shape[0], -1, dtype=int)
for i, p in enumerate(points):
    new = np.all(grid <= p, axis=1) & (owner == -1)
    owner[new] = i

# ---- 2. build a 4-D colour array for ax.voxels -------------------
# occupancy mask
filled = owner.reshape((res, res, res)) >= 0

# colour every cell by its owner; normalise to [0,1] for RGB
col_idx = owner.reshape((res, res, res))
cmap = plt.cm.get_cmap('viridis', len(points))
colors = cmap(col_idx / (len(points)-1))
colors[..., 3] = 0.8                 # global α-channel (80 % opaque)

# ---- 3. render ----------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.voxels(filled, facecolors=colors, edgecolor='k')
ax.set(xlabel='f1', ylabel='f2', zlabel='f3',
       title='Exclusive ΔHV slices rendered as voxels')
plt.tight_layout()
plt.show()
