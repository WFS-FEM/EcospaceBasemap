# dbSEABED raw grids (tracked copy)

Four percent-composition grids for the northern Gulf of Mexico from the dbSEABED project
("Data for Modellers"), at their native 1.2 arc-minute (0.02 degree) resolution, 891 x 383
cells, lower-left corner -98.19, 23.36. The header says `NODATA_value -9999` but missing
cells are coded `-99`; `fn.rasterize_dbseabed()` treats them as NA. Only the `.asc` grids are
tracked; the zips' auxiliary `.jpg`/`.avl`/`.url` files stay ignored.

| Folder | File | Variable | MD5 (LF form, as stored and checked out) |
|---|---|---|---|
| `Gmf_RCK/` | `gmf_RCK_val.asc` | percent rock outcrop (an areal fraction, separate from the grain-size triangle) | `e1332f6d7983bcae79a5e2b134021d21` |
| `Gmf_GVL/` | `gmf_GVL_val.asc` | percent gravel | `1242a3351fc554895c5a5c413d8cf5d8` |
| `Gmf_SND/` | `gmf_SND_val.asc` | percent sand | `8147b22640263554185940f712f74cf7` |
| `Gmf_MUD/` | `gmf_MUD_val.asc` | percent mud | `16e40884cbdc192cbdb03c3012d50adb` |

Source: https://csdms.colorado.edu/wiki/DBSEABED#Data_for_Modellers (zips
`https://csdms.colorado.edu/csdms_wiki/images/Gmf_<CLS>.zip`). Jenkins, C. (dbSEABED,
INSTAAR, University of Colorado). Cite dbSEABED when using these layers.

Provenance of this copy: downloaded 22 Sep 2026 by `fn.pull_dbseabed()` in
`R/dbSEABED_functions.R` and committed unchanged. The grids are tracked because the CSDMS
server was unreachable during the October 2026 GFISHER reproducibility review
(WFS-FEM/GFISHER#2), which would otherwise block a fresh clone of either repository.
`fn.pull_dbseabed()` remains the refresh path: remove the folder and `fn.pull_all()` fetches
it again. The `WFS-FEM/GFISHER` repository ships a byte-identical copy (its
`data/dbseabed/SOURCE.md` records the same MD5s) because its substrate-affinity stage reads
the raw grids directly; the two copies are reconciled under WFS-FEM/GFISHER#3 and
WFS-FEM/EcospaceBasemap#2.
