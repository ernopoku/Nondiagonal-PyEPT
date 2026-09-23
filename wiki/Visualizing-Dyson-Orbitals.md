# Visualizing Dyson orbitals

This guide takes you from a converged propagator calculation to a three-dimensional orbital plot. It applies to every implemented method, including BD-T1, and to both IP and EA calculations.

## 1. Prepare your environment

Install Nondiagonal PyEPT as described in [[Installation]]. Cube export uses `pyscf.tools.cubegen`, which is included with PySCF. Install [Avogadro](https://avogadro.cc/) separately to view the exported files.

A cube file stores the orbital amplitude on a three-dimensional grid together with the molecular geometry. Export **`pole.dyson_ao`**, the coefficients in the AO basis of the same PySCF `mol` used for the calculation. Do not pass `pole.dyson_mo` to the cube writer.

## 2. Run a calculation and export the orbitals

Save this complete example as `plot_dyson.py`, then run `python plot_dyson.py` in your installed environment.

```python
from pathlib import Path
from pyscf import gto, scf
from pyscf.tools import cubegen
from nondiagonal_ept import EPT

mol = gto.M(
    atom="""
    F  0.0000  0.0000  0.0000
    H  0.0000  0.0000  0.9168
    """,
    unit="Angstrom",
    basis="cc-pvtz",
    verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)
calculation = EPT(mf, "NRL3", frozen=1, sector="ip")
poles = calculation.kernel(targets=[4, 2], tol=1e-9)

output = Path("dyson_orbitals")
output.mkdir(exist_ok=True)

for pole in poles:
    filename = output / (
        f"{calculation.method}_{pole.sector}_MO{pole.target}.cube"
    )
    cubegen.orbital(
        mol, str(filename), pole.dyson_ao,
        nx=100, ny=100, nz=100, margin=5.0,
    )
    print(f"Saved {filename} — PS = {pole.strength:.6f}")

print("Output directory:", output.resolve())
```

The example creates `NRL3_ip_MO4.cube` and `NRL3_ip_MO2.cube` in `dyson_orbitals`, relative to the directory from which you run Python. The labels use zero-based target MO indices. These are correlated Dyson orbitals, not simply the corresponding RHF orbitals.

If you already have `poles = calculation.kernel(...)`, reuse those results and run only the export loop. You do not need to solve the propagator again.

## 3. Open the cube file in Avogadro

1. Open Avogadro and use **File → Open**, or drag a `.cube` file into its window.
2. Recent versions display a surface automatically. Use **Create Surfaces** to adjust it; the menu location depends on your Avogadro version.
3. Select the loaded cube data and display both positive and negative isosurfaces in different colors. A useful starting point is **+0.03 and −0.03** in atomic units of orbital amplitude. Some versions ask for one magnitude and generate both signs.
4. Rotate the molecule and inspect its nodal planes. Reduce the isovalue to display more diffuse regions; increase it to focus on larger amplitudes.
5. Use the viewer's image export or screenshot function to save a figure. Record the isovalue and normalization in the caption.

The two colors represent opposite signs (phases) of the orbital amplitude, **not positive and negative charge**. Reversing every sign leaves the physical orbital unchanged.

[Avogadro cube-file guidance](https://discuss.avogadro.cc/t/how-to-visualize-a-cube-file-in-avogadro/2409/16)

## 4. Export results from several methods

Use the same `mol`, `mf`, and `output` setup above. If you already have a `results` dictionary from `run_methods`, skip the calculation call and reuse it.

```python
from nondiagonal_ept import run_methods

results = run_methods(
    mf, ["NRL3", "NRQ3", "NRP3"],
    frozen=1, sector="ip", targets=[4, 2], tol=1e-9,
)

for method, poles in results.items():
    for pole in poles:
        filename = output / f"{method}_{pole.sector}_MO{pole.target}.cube"
        cubegen.orbital(
            mol, str(filename), pole.dyson_ao,
            nx=100, ny=100, nz=100, margin=5.0,
        )
        print(f"Saved {filename} — PS = {pole.strength:.6f}")
```

Method and sector names in filenames distinguish the results. Repeating an export with the same filename overwrites that file; use separate directories for different geometries or basis sets.

## 5. Choose the normalization

The default Dyson orbital retains its pole strength:

$$
\int |\phi_D(\mathbf r)|^2\,d\mathbf r = \mathrm{PS}.
$$

`pole.strength` is the Python attribute for PS. Keep the original coefficients when you want to retain this amplitude information.

To compare **shapes independently of PS**, divide by the square root of PS before exporting. This changes the norm to one and does not change the stored result:

```python
import numpy as np

for pole in poles:
    if pole.strength <= 1e-12:
        print(f"Skipping MO {pole.target}: PS is too small to normalize reliably")
        continue
    shape_coeff = pole.dyson_ao / np.sqrt(pole.strength)
    filename = output / (
        f"{calculation.method}_{pole.sector}_MO{pole.target}_normalized.cube"
    )
    cubegen.orbital(
        mol, str(filename), shape_coeff,
        nx=100, ny=100, nz=100, margin=5.0,
    )
```

This normalization loop is for the single-method `calculation` and `poles` from section 2. In a multi-method loop, use its `method` variable in the filename. Compare orbitals with the **same normalization, isovalue, geometry orientation, and grid settings**. Report PS separately when plotting normalized shapes.

In the AO basis, the norm is `c.T @ S @ c`, where `S = mf.get_ovlp()`; it is not generally `c.T @ c`. A finite cube grid only approximates the full-space normalization.

## 6. BD-T1 and electron attachment

For BD-T1, replace the constructor with:

```python
calculation = EPT(mf, "BD-T1", frozen=1, sector="ip")
```

The export code stays the same. `pole.dyson_ao` already incorporates the optimized Brueckner orbitals. Do not transform it again with the RHF MO coefficients. Target labels refer to the resulting Brueckner orbitals, which can change order or character relative to RHF. For the N2 fixed-core benchmark, use `frozen=2` and targets `[6, 4, 3]`; see the [BD-T1 example](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/examples/n2_bdt1.py).

For attachment, use `sector="ea"` and suitable virtual-orbital targets. Export the resulting `dyson_ao` in exactly the same way.

## 7. Grid quality and troubleshooting

| Situation | What to check |
|---|---|
| No surface appears | Confirm the cube loaded, then lower the isovalue. A weak or diffuse orbital may have amplitudes below the chosen level. |
| Jagged surface | Increase `nx`, `ny`, and `nz`, for example from 100 to 140. |
| Surface cut off at the box boundary | Increase `margin`; diffuse attachment orbitals may need substantially more padding. Increase grid counts too, to preserve spatial resolution. |
| Export is slow or files are large | Try 60–80 points along each axis for preview. Doubling all three grid counts increases the number of values eightfold. |
| Different methods have swapped colors | Overall orbital sign is arbitrary. Compare nodes and shape, not the color assignment alone. |
| Degenerate orbitals appear rotated | Individual orbitals within a degenerate subspace can mix or rotate. Compare the subspace and molecular orientation rather than treating that rotation as an error. |

`nx`, `ny`, and `nz` are grid-point counts. `margin` is in **bohr**, even when molecular coordinates were entered in ångströms. Increase the box and refine the grid until the features of interest are stable.

This guide plots the signed orbital amplitude. Plotting its squared magnitude gives a nonnegative field and hides the phase information; it is not the molecule's total electron density.

[PySCF cube-generation documentation](https://pyscf.org/_modules/pyscf/tools/cubegen.html)

[[Home]] | [[Examples]] | [[API Reference|API-Reference]]
