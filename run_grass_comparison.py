"""
Run GRASS GIS v.surf.rst on CSV point data and export grids for comparison
with the C# SplineTensionGrass implementation.

Prerequisites:
  - GRASS GIS 8.x installed (https://grass.osgeo.org/)
  - Must be run inside a GRASS session with a suitable CRS

Usage:
  grass --tmp-project EPSG:4326 --exec python run_grass_comparison.py <input_dir> [output_dir]

Input:
  <input_dir>   Directory containing one or more CSV files. Each CSV must have
                three columns (no header row):

                    longitude,latitude,value

                Example:
                    -71.289010,42.469889,65.3
                    -71.288500,42.470100,62.1
                    ...

                Lines starting with non-numeric characters and lines containing
                "END" are skipped. Coordinates should be in the same CRS as
                the GRASS session (e.g. WGS84 lon/lat for EPSG:4326).

  [output_dir]  Optional. Directory for output grids. Defaults to
                "grass_comparison" next to this script.

Output:
  One CSV per input file named "<stem>_grass.csv" containing a metadata header
  followed by the raster values exported from GRASS. These can be compared
  cell-by-cell against grids produced by SplineTensionGrass.

Parameters:
  The script uses the same defaults as SplineTensionGrass:
    tension=40, smooth=0.1, segmax=40, npmin=300, cellSize=min(W,H)/250
"""

import os
import sys
import glob
import time
import tempfile

try:
    import grass.script as gs
except ImportError:
    print("ERROR: Must run within a GRASS GIS session.")
    print("Usage: grass --tmp-project EPSG:4326 --exec python run_grass_comparison.py <input_dir> [output_dir]")
    sys.exit(1)


def read_csv_points(csv_path):
    """Read a CSV of lon,lat,value (no header)."""
    points = []
    with open(csv_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.upper().startswith('END'):
                break
            parts = line.split(',')
            if len(parts) >= 3:
                try:
                    lon = float(parts[0].strip())
                    lat = float(parts[1].strip())
                    val = float(parts[2].strip())
                    points.append((lon, lat, val))
                except ValueError:
                    continue
    return points


def run_comparison(csv_path, output_dir):
    """Run v.surf.rst on a single CSV and export the grid."""
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    output_csv = os.path.join(output_dir, f"{stem}_grass.csv")

    print(f"\n{'=' * 60}")
    print(f"File: {csv_path}")
    print(f"{'=' * 60}")

    # Read points
    points = read_csv_points(csv_path)
    if not points:
        print("  SKIP: no valid points found")
        return False

    print(f"  Points: {len(points)}")

    # Write points to temp file for GRASS import
    point_file = os.path.join(tempfile.gettempdir(), f"grass_cmp_{stem}.csv")
    with open(point_file, 'w') as f:
        for p in points:
            f.write(f"{p[0]}|{p[1]}|{p[2]}\n")

    map_name = f"pts_{stem}"

    # Import as 3D points
    gs.run_command('v.in.ascii', input=point_file, output=map_name,
                   separator='pipe', flags='z',
                   x=1, y=2, z=3, skip=0, overwrite=True, quiet=True)

    # Compute grid parameters matching SplineTensionGrass defaults
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    width = max_x - min_x
    height = max_y - min_y
    cell_size = min(width, height) / 250.0

    # Expand by half cell (matching SplineTensionGrass)
    min_x -= 0.5 * cell_size
    min_y -= 0.5 * cell_size
    max_x += 0.5 * cell_size
    max_y += 0.5 * cell_size

    n_cols = round((max_x - min_x) / cell_size)
    n_rows = round((max_y - min_y) / cell_size)

    gs.run_command('g.region', n=max_y, s=min_y, e=max_x, w=min_x,
                   rows=n_rows, cols=n_cols, quiet=True)

    print(f"  Region: {n_cols}x{n_rows}, cell={cell_size:.10f}")
    print(f"  Extent: ({min_x:.10f}, {min_y:.10f}) - ({max_x:.10f}, {max_y:.10f})")

    # Run v.surf.rst with the same defaults as SplineTensionGrass
    raster_name = f"spline_{stem}"
    print(f"  Running v.surf.rst (tension=40, smooth=0.1, segmax=40, npmin=300)...")
    t0 = time.perf_counter()
    gs.run_command('v.surf.rst', input=map_name,
                   elevation=raster_name,
                   tension=40., smooth=0.1,
                   segmax=40, npmin=300,
                   overwrite=True)
    elapsed = time.perf_counter() - t0
    print(f"  v.surf.rst elapsed: {elapsed:.3f}s")

    # Get region info
    info_region = gs.parse_command('g.region', flags='g', quiet=True)
    r_n = float(info_region['n'])
    r_s = float(info_region['s'])
    r_e = float(info_region['e'])
    r_w = float(info_region['w'])
    r_rows = int(info_region['rows'])
    r_cols = int(info_region['cols'])
    r_ewres = float(info_region['ewres'])
    r_nsres = float(info_region['nsres'])

    # Export grid as headerless ASCII raster
    raster_data = gs.read_command('r.out.ascii', input=raster_name,
                                   flags='h', quiet=True)

    # Write CSV with header metadata
    os.makedirs(output_dir, exist_ok=True)
    with open(output_csv, 'w') as f:
        f.write(f"ncols,{r_cols}\n")
        f.write(f"nrows,{r_rows}\n")
        f.write(f"west,{r_w:.15f}\n")
        f.write(f"south,{r_s:.15f}\n")
        f.write(f"east,{r_e:.15f}\n")
        f.write(f"north,{r_n:.15f}\n")
        f.write(f"ewres,{r_ewres:.15f}\n")
        f.write(f"nsres,{r_nsres:.15f}\n")
        f.write("---\n")
        f.write(raster_data)

    print(f"  Output: {output_csv}")

    # Print stats
    info = gs.parse_command('r.univar', map=raster_name, flags='g', quiet=True)
    if info:
        print(f"  Stats: min={info.get('min','?')}, max={info.get('max','?')}, "
              f"mean={info.get('mean','?')}")

    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: grass --tmp-project EPSG:4326 --exec python run_grass_comparison.py <input_dir> [output_dir]")
        print()
        print("  <input_dir>   Directory containing CSV files (lon,lat,value per line)")
        print("  [output_dir]  Output directory (default: grass_comparison/)")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "grass_comparison")

    if not os.path.isdir(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)

    csv_files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
    if not csv_files:
        print(f"ERROR: No CSV files found in {input_dir}")
        sys.exit(1)

    print("GRASS GIS v.surf.rst Grid Export for Comparison")
    print(f"GRASS version: {gs.version()['version']}")
    print(f"Input:  {input_dir} ({len(csv_files)} CSV files)")
    print(f"Output: {output_dir}")

    for csv_path in csv_files:
        try:
            run_comparison(csv_path, output_dir)
        except Exception as e:
            print(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()

    print("\nDone.")


if __name__ == '__main__':
    main()
