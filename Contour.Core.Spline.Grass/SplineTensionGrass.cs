using Contour.Core.Interfaces;
using NetTopologySuite.Geometries;

namespace Contour.Core.Spline.Grass;

/// <summary>
/// Implements spline interpolation matching GRASS GIS v.surf.rst (Mitasova &amp; Mitas 1993).
/// Uses the Completely Regularized Spline with Tension (CRST) radial basis function,
/// rectangular segmentation with adaptive window sizing, dnorm normalization,
/// and optional smoothing.
/// Reference: https://github.com/OSGeo/grass/tree/main/lib/rst/interp_float
/// </summary>
public sealed class SplineTensionGrass : ISpline
{
    /// <summary>
    /// Tension parameter (fi). GRASS default is 40.
    /// Higher values produce surfaces that conform more tightly to data points.
    /// </summary>
    public double Tension { get; }

    /// <summary>
    /// Smoothing parameter (rsm). GRASS default is 0.1.
    /// Added to the matrix diagonal as -rsm, regularizing the system.
    /// Set to 0 for exact interpolation (no smoothing).
    /// </summary>
    public double Smoothing { get; }

    /// <summary>
    /// Maximum number of input points per segment before subdividing.
    /// GRASS default: 40 (MAXSEGM).
    /// </summary>
    public int SegMax { get; }

    /// <summary>
    /// Minimum number of points required per segment (KMIN/npmin).
    /// GRASS default: 300 (MINPOINTS).
    /// </summary>
    public int NpMin { get; }

    /// <summary>
    /// Maximum number of points to use per segment (KMAX2).
    /// GRASS default: 2 * npmin = 600.
    /// </summary>
    public int KMax2 { get; }

    /// <param name="tension">Tension parameter fi. Default 40 matches GRASS.</param>
    /// <param name="smoothing">Smoothing parameter rsm. Default 0.1 matches GRASS.</param>
    /// <param name="segMax">Max points per segment before subdividing. Default 40.</param>
    /// <param name="npMin">Minimum points per segment. Default 300.</param>
    /// <param name="kMax2">Maximum points per segment for solving. Default 2*npMin.</param>
    public SplineTensionGrass(double tension = 40.0, double smoothing = 0.1,
        int segMax = 40, int npMin = 300, int kMax2 = 0)
    {
        if (tension <= 0)
            throw new ArgumentOutOfRangeException(nameof(tension), "Tension must be positive.");
        if (segMax < 1)
            throw new ArgumentOutOfRangeException(nameof(segMax), "segMax must be at least 1.");
        if (npMin < 4)
            throw new ArgumentOutOfRangeException(nameof(npMin), "npMin must be at least 4.");

        Tension = tension;
        Smoothing = smoothing;
        SegMax = segMax;
        NpMin = npMin;
        KMax2 = kMax2 > 0 ? kMax2 : 2 * npMin;
    }

    /// <summary>
    /// Interpolates scattered WGS84 points onto a regular grid using LCC projection.
    /// </summary>
    public CoordinateM[,] InterpolateToGrid(List<CoordinateM> points, double airportLatitude, double airportLongitude)
    {
        if (points.Count < 4)
            throw new ArgumentException("Need at least 4 points for spline interpolation.", nameof(points));

        var lcc = new LambertConformalConic(airportLatitude, airportLongitude);

        var projected = new List<CoordinateM>(points.Count);
        foreach (var p in points)
        {
            var (easting, northing) = lcc.Forward(p.X, p.Y);
            projected.Add(new CoordinateM(easting, northing, p.M));
        }

        double cellSize = EstimateCellSize(projected);

        double minX = double.MaxValue, minY = double.MaxValue;
        double maxX = double.MinValue, maxY = double.MinValue;
        foreach (var p in projected)
        {
            if (p.X < minX) minX = p.X;
            if (p.Y < minY) minY = p.Y;
            if (p.X > maxX) maxX = p.X;
            if (p.Y > maxY) maxY = p.Y;
        }

        minX -= 0.5 * cellSize;
        minY -= 0.5 * cellSize;
        maxX += 0.5 * cellSize;
        maxY += 0.5 * cellSize;

        int nCols = (int)Math.Round((maxX - minX) / cellSize);
        int nRows = (int)Math.Round((maxY - minY) / cellSize);
        double yOriginTop = minY + nRows * cellSize;

        var projectedGrid = InterpolateToGrid(projected, cellSize, nCols, nRows, minX, yOriginTop);

        var grid = new CoordinateM[nCols, nRows];
        for (int col = 0; col < nCols; col++)
        {
            for (int row = 0; row < nRows; row++)
            {
                var cell = projectedGrid[col, row];
                var (lon, lat) = lcc.Inverse(cell.X, cell.Y);
                grid[col, row] = new CoordinateM(lon, lat, cell.M);
            }
        }

        return grid;
    }

    /// <summary>
    /// Interpolates scattered points onto a regular grid with auto-determined extent.
    /// </summary>
    public CoordinateM[,] InterpolateToGrid(List<CoordinateM> points, double cellSize)
    {
        if (points.Count < 4)
            throw new ArgumentException("Need at least 4 points for spline interpolation.", nameof(points));

        double minX = double.MaxValue, minY = double.MaxValue;
        double maxX = double.MinValue, maxY = double.MinValue;
        foreach (var p in points)
        {
            if (p.X < minX) minX = p.X;
            if (p.Y < minY) minY = p.Y;
            if (p.X > maxX) maxX = p.X;
            if (p.Y > maxY) maxY = p.Y;
        }

        minX -= 0.5 * cellSize;
        minY -= 0.5 * cellSize;
        maxX += 0.5 * cellSize;
        maxY += 0.5 * cellSize;

        int nCols = (int)Math.Round((maxX - minX) / cellSize);
        int nRows = (int)Math.Round((maxY - minY) / cellSize);
        double yOriginTop = minY + nRows * cellSize;

        return InterpolateToGrid(points, cellSize, nCols, nRows, minX, yOriginTop);
    }

    /// <summary>
    /// Interpolates scattered points onto a regular grid with explicit extent.
    /// Matches GRASS GIS v.surf.rst algorithm: grid-region quadtree with per-leaf
    /// density filtering, adaptive MINPTS per segment, bisection window search,
    /// dnorm normalization, CRST radial basis, LU decomposition, cell-center evaluation.
    /// </summary>
    public CoordinateM[,] InterpolateToGrid(List<CoordinateM> points, double cellSize,
        int nCols, int nRows, double xOrigin, double yOriginTop)
    {
        if (points.Count < 4)
            throw new ArgumentException("Need at least 4 points for spline interpolation.", nameof(points));

        var grid = new CoordinateM[nCols, nRows];

        // Initialize all cells with proper coordinates and NoData value.
        for (int row = 0; row < nRows; row++)
        {
            double y = yOriginTop - (row + 0.5) * cellSize;
            for (int col = 0; col < nCols; col++)
            {
                double x = xOrigin + (col + 0.5) * cellSize;
                grid[col, row] = new CoordinateM(x, y, -9999);
            }
        }

        // Pre-extract point coordinates into flat arrays and find min x,y
        int nPts = points.Count;
        var allX = new double[nPts];
        var allY = new double[nPts];
        var allZ = new double[nPts];
        double rawMinX = double.MaxValue, rawMinY = double.MaxValue;
        for (int i = 0; i < nPts; i++)
        {
            allX[i] = points[i].X;
            allY[i] = points[i].Y;
            allZ[i] = points[i].M;
            if (allX[i] < rawMinX) rawMinX = allX[i];
            if (allY[i] < rawMinY) rawMinY = allY[i];
        }

        // Translate all coordinates to near zero BEFORE tree construction.
        // This improves floating-point precision in midpoint calculations during
        // quadtree subdivision, ensuring boundary points are assigned to the same
        // quadrant as in GRASS GIS (which works with coordinates of similar magnitude).
        for (int i = 0; i < nPts; i++)
        {
            allX[i] -= rawMinX;
            allY[i] -= rawMinY;
        }

        // Grid region bounds (translated)
        double gridXMin = xOrigin - rawMinX;
        double gridYMax = yOriginTop - rawMinY;
        double gridXMax = xOrigin + nCols * cellSize - rawMinX;
        double gridYMin = yOriginTop - nRows * cellSize - rawMinY;

        // GRASS dmin: (min(ew_res, ns_res) / 2)^2 — squared for distance comparison
        double dminSq = (cellSize / 2.0) * (cellSize / 2.0);

        // --- Phase 1: Build quadtree on grid region with density filtering ---
        // Tree is built in translated coordinates for maximum floating-point precision
        // in midpoint calculations during subdivision.
        var tree = QuadTree.BuildFromGrid(allX, allY, nPts,
            gridXMin, gridYMin, gridXMax, gridYMax,
            nRows, nCols, SegMax, dminSq, out _);

        var leaves = tree.GetLeaves();

        // Collect accepted points and compute extent + zmin (in translated coords)
        double zMin = double.MaxValue;
        double ptMaxX = double.MinValue, ptMaxY = double.MinValue;
        double ptTransMinX = double.MaxValue, ptTransMinY = double.MaxValue;
        int acceptedCount = 0;
        foreach (var leaf in leaves)
        {
            if (leaf.PointIndices == null) continue;
            foreach (int i in leaf.PointIndices)
            {
                if (allZ[i] < zMin) zMin = allZ[i];
                if (allX[i] < ptTransMinX) ptTransMinX = allX[i];
                if (allY[i] < ptTransMinY) ptTransMinY = allY[i];
                if (allX[i] > ptMaxX) ptMaxX = allX[i];
                if (allY[i] > ptMaxY) ptMaxY = allY[i];
                acceptedCount++;
            }
        }

        if (acceptedCount < 4)
            throw new ArgumentException("Too few points survived density filtering.", nameof(points));

        // GRASS dnorm: sqrt(deltx * delty * KMIN / NPOINT) — using accepted point extent
        double deltX = ptMaxX - ptTransMinX;
        double deltY = ptMaxY - ptTransMinY;
        double dnorm = Math.Sqrt(deltX * deltY * NpMin / acceptedCount);
        if (dnorm < 1e-15) dnorm = 1.0;

        double smseg = tree.SmallestLeafWidth();

        // Global cell size
        double ewRes = (gridXMax - gridXMin) / nCols;
        double nsRes = (gridYMax - gridYMin) / nRows;

        // --- Phase 2: Process each leaf segment (parallelized) ---
        double fi = Tension;
        Parallel.ForEach(leaves, leaf =>
        {
            double leafXMin = leaf.XMin;
            double leafYMin = leaf.YMin;
            double leafXMax = leaf.XMax;
            double leafYMax = leaf.YMax;

            // GRASS: distx = n_cols * ew_res * 0.1 = leafWidth * 0.1
            double distx = (leafXMax - leafXMin) * 0.1;
            double disty = (leafYMax - leafYMin) * 0.1;
            double distxp = 0, distyp = 0;
            int MAXENC = 0;

            // GRASS adaptive MINPTS: based on segment size relative to smallest segment
            double pr = Math.Pow(2.0, (leafXMax - leafXMin) / smseg - 1.0);
            int minPts = (int)(NpMin * (pr / (1.0 + NpMin * pr / KMax2)));

            var regionIndices = new List<int>();
            tree.FindPointsInRegion(allX, allY,
                leafXMin - distx, leafYMin - disty,
                leafXMax + distx, leafYMax + disty,
                regionIndices, KMax2 + 1);

            // GRASS bisection loop: while (npt < MINPTS || npt > KMAX2), up to 70 iterations
            for (int iter = 0; iter < 70; iter++)
            {
                if (regionIndices.Count >= minPts && regionIndices.Count <= KMax2)
                    break;

                if (regionIndices.Count > KMax2)
                {
                    MAXENC = 1;
                    double temp1 = distxp;
                    distxp = distx;
                    distx = distxp - Math.Abs(distx - temp1) * 0.5;

                    double temp2 = distyp;
                    distyp = disty;
                    disty = distyp - Math.Abs(disty - temp2) * 0.5;
                }
                else
                {
                    double temp1 = distyp;
                    distyp = disty;
                    double temp2 = distxp;
                    distxp = distx;

                    if (MAXENC != 0)
                    {
                        disty = Math.Abs(disty - temp1) * 0.5 + distyp;
                        distx = Math.Abs(distx - temp2) * 0.5 + distxp;
                    }
                    else
                    {
                        distx += distx;
                        disty += disty;
                    }
                }

                regionIndices.Clear();
                tree.FindPointsInRegion(allX, allY,
                    leafXMin - distx, leafYMin - disty,
                    leafXMax + distx, leafYMax + disty,
                    regionIndices, KMax2 + 1);
            }

            if (regionIndices.Count == 0)
                return; // truly no data anywhere — skip

            if (regionIndices.Count > KMax2)
                regionIndices = regionIndices.GetRange(0, KMax2);

            int n = regionIndices.Count;

            // GRASS normalization: shift to leaf origin, divide by dnorm
            // (origin cancels in distance calculations, but matches GRASS convention)
            var pxn = new double[n];
            var pyn = new double[n];
            var pz = new double[n];
            for (int k = 0; k < n; k++)
            {
                int idx = regionIndices[k];
                pxn[k] = (allX[idx] - leafXMin) / dnorm;
                pyn[k] = (allY[idx] - leafYMin) / dnorm;
                pz[k] = allZ[idx] - zMin;
            }

            // Solve using in-place jagged-array LU decomposition
            var coeffs = SolveCrstSystemFast(pxn, pyn, pz, n, fi, Smoothing);

            // Evaluate grid cells for this leaf using its NRows/NCols
            // (matches GRASS: each segment evaluates exactly its own grid cells)
            int colStart = (int)Math.Round((leafXMin - gridXMin) / ewRes);
            int colEnd = colStart + leaf.NCols;
            int rowStart = (int)Math.Round((gridYMax - leafYMax) / nsRes);
            int rowEnd = rowStart + leaf.NRows;

            // Clamp to grid bounds (safety)
            colStart = Math.Max(0, colStart);
            colEnd = Math.Min(nCols, colEnd);
            rowStart = Math.Max(0, rowStart);
            rowEnd = Math.Min(nRows, rowEnd);

            if (coeffs != null)
            {
                double c0 = coeffs[0] + zMin;
                for (int row = rowStart; row < rowEnd; row++)
                {
                    // Translated grid cell y for normalization
                    double yt = gridYMax + nsRes * 0.5 - (row + 1) * nsRes;
                    double yn = (yt - leafYMin) / dnorm;

                    for (int col = colStart; col < colEnd; col++)
                    {
                        // Translated grid cell x for normalization
                        double xt = gridXMin + (col + 0.5) * ewRes;
                        double xn = (xt - leafXMin) / dnorm;

                        double z = c0;
                        for (int j = 0; j < n; j++)
                        {
                            double ddx = xn - pxn[j];
                            double ddy = yn - pyn[j];
                            z += coeffs[j + 1] * CrstBasis(ddx * ddx + ddy * ddy, fi);
                        }

                        // Output in original (untranslated) coordinates
                        double x = xOrigin + (col + 0.5) * cellSize;
                        double y = yOriginTop - (row + 0.5) * cellSize;
                        grid[col, row] = new CoordinateM(x, y, z);
                    }
                }
            }
            else
            {
                // IDW fallback (uses translated point coordinates via allX/allY)
                var rp = new List<CoordinateM>(n);
                foreach (int idx in regionIndices)
                    rp.Add(new CoordinateM(allX[idx], allY[idx], allZ[idx] + zMin));

                for (int row = rowStart; row < rowEnd; row++)
                {
                    double yt = gridYMax + nsRes * 0.5 - (row + 1) * nsRes;
                    for (int col = colStart; col < colEnd; col++)
                    {
                        double xt = gridXMin + (col + 0.5) * ewRes;
                        double x = xOrigin + (col + 0.5) * cellSize;
                        double y = yOriginTop - (row + 0.5) * cellSize;
                        grid[col, row] = new CoordinateM(x, y, FallbackIdw(rp, xt, yt));
                    }
                }
            }
        });

        return grid;
    }

    /// <summary>
    /// GRASS Completely Regularized Spline with Tension radial basis function (IL_crst).
    /// Computes R(r²) = -Ei(-fi²*r²/4) where Ei is the exponential integral.
    /// Input r is DISTANCE SQUARED.
    /// </summary>
    internal static double CrstBasis(double rSquared, double fi)
    {
        double x = fi * fi * rSquared / 4.0;

        if (x < 1.0)
        {
            // Power series: sum_{k=0}^{9} u[k] * x^(k+1) in Horner form
            return x * (1.0
                + x * (-0.25
                + x * (0.055555555555556
                + x * (-0.010416666666667
                + x * (0.00166666666666667
                + x * (-2.31481481481482e-04
                + x * (2.83446712018141e-05
                + x * (-3.10019841269841e-06
                + x * (3.06192435822065e-07
                + x * (-2.75573192239859e-08))))))))));
        }

        const double ce = 0.57721566; // Euler-Mascheroni

        double e1;
        if (x > 25.0)
        {
            e1 = 0.0;
        }
        else
        {
            // Rational approximation for E1(x) * x * exp(x)
            double ea = 0.2677737343 + x * (8.6347608925 + x * (18.0590169730 + x * (8.5733287401 + x)));
            double eb = 3.9584969228 + x * (21.0996530827 + x * (25.6329561486 + x * (9.5733223454 + x)));
            e1 = (ea / eb) / (x * Math.Exp(x));
        }

        return e1 + ce + Math.Log(x);
    }

    /// <summary>
    /// Builds and solves the GRASS spline system using jagged arrays and in-place LU.
    /// System is (n+1)x(n+1): row/col 0 = constant trend constraint,
    /// rows/cols 1..n = data points. Diagonal has -smoothing.
    /// </summary>
    private static double[]? SolveCrstSystemFast(double[] px, double[] py, double[] pz,
        int n, double fi, double smoothing)
    {
        int size = n + 1;

        // Build matrix directly as jagged array (better cache locality, no copy needed)
        var lu = new double[size][];
        for (int i = 0; i < size; i++)
            lu[i] = new double[size];

        var rhs = new double[size];

        // Row 0 / Col 0: trend constraint (sum of lambdas = 0)
        rhs[0] = 0.0;
        var row0 = lu[0];
        for (int i = 0; i < n; i++)
        {
            row0[i + 1] = 1.0;
            lu[i + 1][0] = 1.0;
            rhs[i + 1] = pz[i];
        }

        // Fill radial basis matrix (upper triangle, then mirror)
        for (int i = 0; i < n; i++)
        {
            var rowI = lu[i + 1];
            double pxi = px[i], pyi = py[i];
            for (int j = i + 1; j < n; j++)
            {
                double dx = pxi - px[j];
                double dy = pyi - py[j];
                double val = CrstBasis(dx * dx + dy * dy, fi);
                rowI[j + 1] = val;
                lu[j + 1][i + 1] = val;
            }
            rowI[i + 1] = -smoothing;
        }

        // In-place LU decomposition with partial pivoting
        return SolveLuJagged(lu, rhs, size);
    }

    /// <summary>
    /// In-place LU decomposition and solve on jagged arrays.
    /// Avoids the matrix copy and benefits from better cache locality.
    /// </summary>
    private static double[]? SolveLuJagged(double[][] lu, double[] b, int n)
    {
        const double Tiny = 1e-20;

        var indx = new int[n];
        var vv = new double[n];

        // Implicit scaling for each row
        for (int i = 0; i < n; i++)
        {
            var row = lu[i];
            double big = 0.0;
            for (int j = 0; j < n; j++)
            {
                double temp = Math.Abs(row[j]);
                if (temp > big) big = temp;
            }
            if (big == 0.0) return null;
            vv[i] = Tiny / big;
        }

        // Crout's algorithm with implicit pivoting
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < j; i++)
            {
                var row = lu[i];
                double sum = row[j];
                for (int k = 0; k < i; k++)
                    sum -= row[k] * lu[k][j];
                row[j] = sum;
            }

            double big = 0.0;
            int imax = j;
            for (int i = j; i < n; i++)
            {
                var row = lu[i];
                double sum = row[j];
                for (int k = 0; k < j; k++)
                    sum -= row[k] * lu[k][j];
                row[j] = sum;

                double temp = vv[i] * Math.Abs(sum);
                if (temp >= big) { big = temp; imax = i; }
            }

            if (imax != j)
            {
                (lu[imax], lu[j]) = (lu[j], lu[imax]);
                vv[imax] = vv[j];
            }
            indx[j] = imax;

            if (lu[j][j] == 0.0) lu[j][j] = Tiny;

            if (j < n - 1)
            {
                double dum = 1.0 / lu[j][j];
                for (int i = j + 1; i < n; i++)
                    lu[i][j] *= dum;
            }
        }

        // Forward/back substitution
        var x = new double[n];
        Array.Copy(b, x, n);

        int ii = -1;
        for (int i = 0; i < n; i++)
        {
            int ip = indx[i];
            double sum = x[ip];
            x[ip] = x[i];
            if (ii >= 0)
            {
                var row = lu[i];
                for (int j = ii; j < i; j++)
                    sum -= row[j] * x[j];
            }
            else if (sum != 0.0)
            {
                ii = i;
            }
            x[i] = sum;
        }

        for (int i = n - 1; i >= 0; i--)
        {
            var row = lu[i];
            double sum = x[i];
            for (int j = i + 1; j < n; j++)
                sum -= row[j] * x[j];
            x[i] = sum / row[i];
        }

        return x;
    }

    /// <summary>
    /// Estimates cell size from point extent: min(width, height) / 250.
    /// Matches Esri's formula used by AEDT.
    /// </summary>
    public static double EstimateCellSize(IList<CoordinateM> points)
    {
        double minX = double.MaxValue, minY = double.MaxValue;
        double maxX = double.MinValue, maxY = double.MinValue;

        foreach (var p in points)
        {
            if (p.X < minX) minX = p.X;
            if (p.Y < minY) minY = p.Y;
            if (p.X > maxX) maxX = p.X;
            if (p.Y > maxY) maxY = p.Y;
        }

        double width = maxX - minX;
        double height = maxY - minY;
        return Math.Min(width, height) / 250.0;
    }

    /// <summary>
    /// Inverse Distance Weighted interpolation fallback for singular spline systems.
    /// </summary>
    internal static double FallbackIdw(List<CoordinateM> points, double x, double y)
    {
        if (points.Count == 0)
            return 0.0;

        double numerator = 0, denominator = 0;
        foreach (var p in points)
        {
            double dx = x - p.X;
            double dy = y - p.Y;
            double d = Math.Sqrt(dx * dx + dy * dy);
            if (d < 1e-15) return p.M;
            double w = 1.0 / (d * d);
            numerator += w * p.M;
            denominator += w;
        }
        return numerator / denominator;
    }

    /// <summary>
    /// Diagnostic: returns the indices of points that survive density filtering
    /// for the given grid parameters. Used to compare filtering behavior against GRASS GIS.
    /// </summary>
    internal HashSet<int> GetAcceptedPointIndices(List<CoordinateM> points, double cellSize,
        int nCols, int nRows, double xOrigin, double yOriginTop)
    {
        int nPts = points.Count;
        var allX = new double[nPts];
        var allY = new double[nPts];
        double rawMinX = double.MaxValue, rawMinY = double.MaxValue;
        for (int i = 0; i < nPts; i++)
        {
            allX[i] = points[i].X;
            allY[i] = points[i].Y;
            if (allX[i] < rawMinX) rawMinX = allX[i];
            if (allY[i] < rawMinY) rawMinY = allY[i];
        }
        for (int i = 0; i < nPts; i++)
        {
            allX[i] -= rawMinX;
            allY[i] -= rawMinY;
        }

        double gridXMin = xOrigin - rawMinX;
        double gridYMax = yOriginTop - rawMinY;
        double gridXMax = xOrigin + nCols * cellSize - rawMinX;
        double gridYMin = yOriginTop - nRows * cellSize - rawMinY;
        double dminSq = (cellSize / 2.0) * (cellSize / 2.0);

        var tree = QuadTree.BuildFromGrid(allX, allY, nPts,
            gridXMin, gridYMin, gridXMax, gridYMax,
            nRows, nCols, SegMax, dminSq, out _);

        var accepted = new HashSet<int>();
        foreach (var leaf in tree.GetLeaves())
        {
            if (leaf.PointIndices == null) continue;
            foreach (int i in leaf.PointIndices)
                accepted.Add(i);
        }
        return accepted;
    }
}
