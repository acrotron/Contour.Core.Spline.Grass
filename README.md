# Contour.Core.Spline.Grass

A C# implementation of the **Completely Regularized Spline with Tension (CRST)** interpolation algorithm from [GRASS GIS `v.surf.rst`](https://grass.osgeo.org/grass-stable/manuals/v.surf.rst.html), plus marching-squares contour generation.

Converts scattered point data into a regular grid and produces contour lines/polygons suitable for noise mapping and similar applications.

## Pipeline

```
Scattered Points                    (WGS84 lon/lat with noise values)
    | LCC projection
Projected Points                    (meters, Lambert Conformal Conic)
    | SplineTensionGrass (GRASS v.surf.rst)
Regular Grid                        (CoordinateM[nCols, nRows] at cell centers)
    | RasterGrid.FromNodes (bilinear center interpolation)
Sub-Triangle Mesh                   (4 triangles per cell with adjacency)
    | MarchingSquaresContourLines (linear interpolation + tracing)
Contour Lines                       (Dictionary<interval, List<LineString>>)
    | MarchingSquaresContourPolygons
Contour Polygons                    (Dictionary<interval, MultiPolygon>)
```

## Algorithm

The interpolation finds a function *f(x, y)* that minimizes:

    Sum_i |f(xi, yi) - zi|^2 / sigma_i^2 + phi^2 * integral of curvature terms

The solution takes the form:

    f(x, y) = a0 + Sum_j lambda_j * R(rj)

where *R(r)* is the CRST radial basis function using the exponential integral, and coefficients are found by solving a linear system via LU decomposition with partial pivoting.

Processing follows GRASS's `v.surf.rst` flow:

1. **Quadtree construction** -- points inserted incrementally; leaves hold at most `segmax` points; density filter rejects points closer than `dmin`
2. **Coordinate translation** -- shift so minimum point is near zero for floating-point precision
3. **Per-segment interpolation** -- adaptive bisection gathers nearby points, normalizes coordinates, assembles and solves the linear system, evaluates grid cells at cell centers

## Default Parameters

| Parameter | Value | GRASS Name | Description |
|---|---|---|---|
| `tension` | 40.0 | `fi` / TENSION | Controls surface stiffness |
| `smoothing` | 0.1 | `rsm` / SMOOTH | Regularization on diagonal |
| `segmax` | 40 | KMAX / MAXSEGM | Max points per leaf before subdivision |
| `npmin` | 300 | KMIN / MINPOINTS | Base minimum points per segment |
| `kmax2` | 600 | KMAX2 (2x npmin) | Max points gathered per segment |
| `dmin` | cellSize/2 | dmin | Min inter-point distance |
| `cellSize` | min(W,H)/250 | ew_res/ns_res | Grid cell size from point extent |

## Contour Generation (Marching Squares)

After the spline produces a regular grid, contour lines are generated using a marching squares approach matching [Esri's contouring algorithm](https://pro.arcgis.com/en/pro-app/3.4/tool-reference/spatial-analyst/how-contouring-works.htm):

1. **Cell subdivision** -- each raster cell (2x2 corners) is split into 4 sub-triangles via a bilinear center point
2. **Contour tracing** -- linear interpolation along triangle edges finds contour crossings; adjacent intersections are connected into polylines
3. **Contour polygons** -- lines at adjacent intervals are combined into closed polygons representing bands between contour levels

## Performance Optimizations

Over a naive GRASS port:

1. **Parallel segment processing** via `Parallel.ForEach` on leaves
2. **Jagged array LU solver** (`double[][]`) for better cache locality
3. **Flat coordinate arrays** to avoid object dereference overhead in hot loops
4. **Leaf-bound evaluation** -- each leaf evaluates only its own grid cells

## Accuracy and Performance

Tested against GRASS GIS 8.4.2 with identical default parameters (tension=40, smooth=0.1, segmax=40, npmin=300) across six datasets ranging from 621 to 149,373 points:

- **RMSE**: 0.01 -- 0.13 dB across all test cases
- **Max difference**: under 1.1 dB in all cases
- **100% of grid cells within 1 dB** of the GRASS GIS reference in every test case
- **3--4x faster** than GRASS GIS for datasets from 1.7K to 149K points

Residual differences are caused by floating-point arithmetic across language runtimes (C compiled by GCC/MSVC vs C# .NET JIT), not algorithmic divergence. LU decomposition amplifies initial ~10^-15 rounding differences to ~10^-5 in spline coefficients. Density filtering and quadtree midpoint rounding at 1-ULP boundaries can also send individual points to different segments.

## GRASS C to C# Function Mapping

| GRASS C Source | C# Implementation | Description |
|---|---|---|
| `MT_insert()` | `QuadTree.Insert()` | Incremental point insertion with subdivision |
| `quad_add_data()` | `QuadTree.AddWithDensityCheck()` | Per-leaf density filter |
| `quad_compare()` | `QuadTree.GetQuadrant()` | Grid-based midpoint with odd-count rounding |
| `quad_divide_data()` | `QuadTree.Subdivide()` | Grid-based subdivision halving rows/cols |
| `quad_get_points()` | `QuadTree.FindPointsInRegion()` | Range search with strict inequality |
| `smallest_segment()` | `QuadTree.SmallestLeafWidth()` | X-width of smallest non-empty leaf |
| `translate_quad()` | `QuadTree.Translate()` | Shift all node bounds by offset |
| `IL_interp_segments_2d()` | `InterpolateToGrid()` per-leaf loop | Bisection search + solve + evaluate |
| `IL_crst()` | `CrstBasis()` | CRST radial basis function |
| `IL_crs_matrix_2d()` | `SolveCrstSystemFast()` | Matrix assembly + LU decomposition |
| `G_lubksb()` | `SolveLuJagged()` | In-place LU back-substitution |
| `IL_grid_calc_2d()` | Per-leaf grid evaluation | Cell-center evaluation with basis functions |

## Comparing Against GRASS GIS

A comparison script is included at `run_grass_comparison.py`. It takes a directory of CSV files (lon,lat,value per line, no header) and runs GRASS `v.surf.rst` with identical parameters, exporting grids for cell-by-cell comparison. Requires a [GRASS GIS](https://grass.osgeo.org/) installation:

```bash
grass --tmp-project EPSG:4326 --exec python run_grass_comparison.py example/
```

An example input file is included in `example/sample_points.csv` (102 points around New York). See the script header for full input format details.

## Scientific References

1. **Mitasova, H. and Mitas, L., 1993.** Interpolation by Regularized Spline with Tension: I. Theory and Implementation. *Mathematical Geology*, 25(6), pp. 641-655. [DOI: 10.1007/BF00893171](https://doi.org/10.1007/BF00893171)
2. **Mitasova, H. and Hofierka, J., 1993.** Interpolation by Regularized Spline with Tension: II. Application to Terrain Modeling and Surface Geometry Analysis. *Mathematical Geology*, 25(6), pp. 657-669. [DOI: 10.1007/BF00893172](https://doi.org/10.1007/BF00893172)
3. **Mitas, L. and Mitasova, H., 1988.** General Variational Approach to the Interpolation Problem. *Computers and Mathematics with Applications*, 16(12), pp. 983-992. [DOI: 10.1016/0898-1221(88)90255-6](https://doi.org/10.1016/0898-1221(88)90255-6)
4. **Talmi, A. and Gilat, G., 1977.** Method for Smooth Approximation of Data. *Journal of Computational Physics*, 23(2), pp. 93-123. [DOI: 10.1016/0021-9991(77)90115-2](https://doi.org/10.1016/0021-9991(77)90115-2)
5. **Wahba, G., 1990.** *Spline Models for Observational Data*. SIAM, Philadelphia. [DOI: 10.1137/1.9781611970128](https://doi.org/10.1137/1.9781611970128)
6. **Mitasova, H., Mitas, L., Brown, W.M., et al., 1995.** Modeling Spatially and Temporally Distributed Phenomena: New Methods and Tools for GRASS GIS. *International Journal of GIS*, 9(4), pp. 433-446. [DOI: 10.1080/02693799508902046](https://doi.org/10.1080/02693799508902046)
7. **Neteler, M. and Mitasova, H., 2008.** *Open Source GIS: A GRASS GIS Approach*, 3rd Edition. Springer. [DOI: 10.1007/978-0-387-68574-8](https://doi.org/10.1007/978-0-387-68574-8)

GRASS GIS 8.4 source: https://github.com/OSGeo/grass | Key directories: `lib/rst/interp_float/`, `lib/rst/qtree/`, `lib/rst/data/`, `vector/v.surf.rst/`

## License

See [LICENSE](LICENSE).
