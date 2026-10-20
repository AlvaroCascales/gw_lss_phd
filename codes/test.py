import hdf5plugin
import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


filename = (
    "/cosmos_storage/simulations/BaccoSims/reorg_preview/"
    "rings/standard/N768_L256/Nenya/0.00/groups_047/"
    "fof_subhalo_history_tab_orph_wweight_047.hdf5"
)


with h5py.File(filename, "r") as f:

    pos = f["Group/GroupPos"][:]

    boxsize = f["Header"].attrs["BoxSize"]
    redshift = f["Header"].attrs["Redshift"]


print("Número de halos:", len(pos))
print("Shape:", pos.shape)
print("Primeras posiciones:")
print(pos[:5])

print("BoxSize:", boxsize)
print("Redshift:", redshift)

print("Mínimos:", pos.min(axis=0))
print("Máximos:", pos.max(axis=0))


# Seleccionamos una lámina fina en la coordenada z
zmin = 0
zmax = 10

mask = (pos[:, 2] >= zmin) & (pos[:, 2] < zmax)
pos_slice = pos[mask]

print("Halos en la lámina:", len(pos_slice))


# Figura x-y
plt.figure(figsize=(8, 8))

plt.scatter(
    pos_slice[:, 0],
    pos_slice[:, 1],
    s=1,
    alpha=0.5
)

plt.xlim(0, boxsize)
plt.ylim(0, boxsize)

plt.xlabel("x")
plt.ylabel("y")

plt.title(
    f"Nenya N768_L256 — halos FoF\n"
    f"z = {redshift:.3f}, {zmin} < Z < {zmax}"
)

plt.gca().set_aspect("equal")

plt.tight_layout()
plt.savefig("halo_distribution.png", dpi=200)

print("Figura guardada en halo_distribution.png")