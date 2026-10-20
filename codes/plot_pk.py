import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


filename = (
    "/cosmos_storage/simulations/BaccoSims/reorg_preview/"
    "rings/standard/N768_L256/Nenya/0.00/powerspec/"
    "powerspec_047_r.txt"
)

# El archivo tiene:
#   4 líneas de cabecera
#   2000 filas primer bloque
#   4 líneas de cabecera
#   2000 filas segundo bloque
#
# Nos quedamos de momento con el segundo bloque.

data = np.loadtxt(
    filename,
    skiprows=2008,
    max_rows=2000
)

print("Shape:", data.shape)
print("Primera fila:")
print(data[0])

k = data[:, 0]

# Potencia adimensional (interpretación preliminar)
delta2 = data[:, 1]

# Número de modos del bin
nmodes = data[:, 4]

# Nos quedamos solo con bins que contienen modos
mask = (
    (nmodes > 0)
    & np.isfinite(delta2)
    & (delta2 > 0)
    & (k <= 5)
)

k = k[mask]
delta2 = delta2[mask]

# Delta^2(k) -> P(k)
pk = 2 * np.pi**2 * delta2 / k**3


plt.figure(figsize=(8, 6))

plt.loglog(k, pk)

plt.xlabel(r"$k\ [h\,\mathrm{Mpc}^{-1}]$")
plt.ylabel(r"$P(k)\ [(h^{-1}\mathrm{Mpc})^3]$")

plt.title("Nenya N768_L256 — nonlinear matter power spectrum")

plt.tight_layout()
plt.savefig("powerspectrum_nenya.png", dpi=200)

print("Número de bins representados:", len(k))
print("k mínimo:", k.min())
print("k máximo:", k.max())
print("Figura guardada en powerspectrum_nenya.png")