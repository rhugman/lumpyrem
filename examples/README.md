# Examples

Three notebooks, each building on the one before.

| notebook | covers |
|---|---|
| [01_single_site.ipynb](01_single_site.ipynb) | Forcing and its gap-filling rules, building and running a `Model`, reading results and choosing when they are reported, delays, the lower store and water-table conversion, convergence, and a first `ModelGrid` as a parameter sweep. |
| [02_spatial_parameters.ipynb](02_spatial_parameters.ipynb) | One model per cell of a synthetic valley: soil and land-use tables and continuous maps turned into per-cell parameters, a shared climate, water-budget maps, single cells, one- and two-store cells in one grid, and engine speed. |
| [03_spatial_forcing.ipynb](03_spatial_forcing.ipynb) | Forcing that differs by cell — a rainfall gradient, crop calendars by land use, irrigation districts — built as an xarray Dataset on `(time, y, x)`, written to and read from netCDF as a `GridForcing`, with results returned as gridded maps and stress-period rates. |

To run them:

```bash
pip install -e ".[fast,examples]"
cd examples
jupyter lab
```

`fast` brings Numba. Without it the notebooks still run, on the pure-Python
kernel, but notebooks 2 and 3 take minutes rather than seconds.

`herebedragons.py` holds what is not about the model: the synthetic
landscape, the example climate reader and the map plots. Notebook 3 writes
netCDF files to `output/`, which git ignores.
