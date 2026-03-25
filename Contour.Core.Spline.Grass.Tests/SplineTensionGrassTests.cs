using AwesomeAssertions;
using NetTopologySuite.Geometries;

namespace Contour.Core.Spline.Grass.Tests;

[TestClass]
public class SplineTensionGrassTests
{
    // ──────────────────────────────────────────
    // Grid output sanity checks
    // ──────────────────────────────────────────

    [TestMethod]
    public void InterpolateToGrid_ProducesCorrectDimensions()
    {
        var points = CreateTestPoints();
        var grass = new SplineTensionGrass(npMin: 4, kMax2: 20);
        var grid = grass.InterpolateToGrid(points, cellSize: 0.5);

        grid.GetLength(0).Should().BeGreaterThan(0);
        grid.GetLength(1).Should().BeGreaterThan(0);
    }

    [TestMethod]
    public void InterpolateToGrid_ExplicitExtent_ProducesCorrectDimensions()
    {
        var points = CreateTestPoints();
        int nCols = 10, nRows = 8;

        var grass = new SplineTensionGrass(npMin: 4, kMax2: 20);
        var grid = grass.InterpolateToGrid(points, cellSize: 0.5, nCols, nRows,
            xOrigin: -0.5, yOriginTop: 3.0);

        grid.GetLength(0).Should().Be(nCols);
        grid.GetLength(1).Should().Be(nRows);
    }

    // ──────────────────────────────────────────
    // Interpolation at data points
    // ──────────────────────────────────────────

    [TestMethod]
    public void InterpolateToGrid_SingleSegment_CloseToDataPoints()
    {
        var points = new List<CoordinateM>
        {
            new(0, 0, 10.0),
            new(1, 0, 20.0),
            new(0, 1, 30.0),
            new(1, 1, 40.0),
            new(0.5, 0.5, 25.0),
        };

        // All points in one segment with smoothing=0.1 (not exact interpolation)
        var grass = new SplineTensionGrass(segMax: 1000, npMin: 4, kMax2: 500, smoothing: 0.1);
        var grid = grass.InterpolateToGrid(points, cellSize: 0.1);

        // Find grid cells closest to each data point
        foreach (var p in points)
        {
            double bestDist = double.MaxValue;
            double bestM = 0;
            for (int col = 0; col < grid.GetLength(0); col++)
            {
                for (int row = 0; row < grid.GetLength(1); row++)
                {
                    double dx = grid[col, row].X - p.X;
                    double dy = grid[col, row].Y - p.Y;
                    double dist = dx * dx + dy * dy;
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestM = grid[col, row].M;
                    }
                }
            }

            // With smoothing, values near data points should be close but not exact
            if (Math.Sqrt(bestDist) < 0.08)
            {
                bestM.Should().BeApproximately(p.M, 5.0,
                    $"Value near ({p.X},{p.Y}) should be close to {p.M}");
            }
        }
    }

    // ──────────────────────────────────────────
    // Smoothness check
    // ──────────────────────────────────────────

    [TestMethod]
    public void InterpolateToGrid_NoLargeJumpsAtSegmentBoundaries()
    {
        var points = new List<CoordinateM>();
        for (double x = 0; x <= 10; x += 1.0)
            for (double y = 0; y <= 10; y += 1.0)
                points.Add(new CoordinateM(x, y, 50.0 + 5.0 * Math.Sin(x) + 5.0 * Math.Cos(y)));

        // segMax=20 forces multiple segments with 121 points
        var grass = new SplineTensionGrass(segMax: 20, npMin: 30, kMax2: 200);
        var grid = grass.InterpolateToGrid(points, cellSize: 0.25);

        int cols = grid.GetLength(0);
        int rows = grid.GetLength(1);

        // Only check cells that were evaluated (non-zero coordinates)
        double maxJump = 0;
        for (int row = 0; row < rows; row++)
        {
            for (int col = 1; col < cols; col++)
            {
                if (grid[col, row].X == 0 && grid[col, row].Y == 0) continue;
                if (grid[col - 1, row].X == 0 && grid[col - 1, row].Y == 0) continue;
                double jump = Math.Abs(grid[col, row].M - grid[col - 1, row].M);
                if (jump > maxJump) maxJump = jump;
            }
        }

        maxJump.Should().BeLessThan(5.0,
            $"Adjacent cells should not have large jumps (got {maxJump:F3})");
    }

    // ──────────────────────────────────────────
    // All evaluated cells have finite values
    // ──────────────────────────────────────────

    [TestMethod]
    public void InterpolateToGrid_EvaluatedCellsHaveFiniteValues()
    {
        var points = CreateTestPoints();
        var grass = new SplineTensionGrass(npMin: 4, kMax2: 20);
        var grid = grass.InterpolateToGrid(points, cellSize: 0.5);

        for (int col = 0; col < grid.GetLength(0); col++)
        {
            for (int row = 0; row < grid.GetLength(1); row++)
            {
                var cell = grid[col, row];
                // Only check cells that were evaluated (covered by a leaf's window)
                if (cell.X == 0 && cell.Y == 0) continue;
                double.IsFinite(cell.M).Should().BeTrue(
                    $"Cell [{col},{row}] at ({cell.X},{cell.Y}) should have finite M value, got {cell.M}");
            }
        }
    }

    // ──────────────────────────────────────────
    // Values within reasonable range
    // ──────────────────────────────────────────

    [TestMethod]
    public void InterpolateToGrid_GridValuesWithinReasonableRange()
    {
        var points = CreateTestPoints();
        double minM = points.Min(p => p.M);
        double maxM = points.Max(p => p.M);
        double range = maxM - minM;

        var grass = new SplineTensionGrass(npMin: 4, kMax2: 20);
        var grid = grass.InterpolateToGrid(points, cellSize: 0.5);

        for (int col = 0; col < grid.GetLength(0); col++)
        {
            for (int row = 0; row < grid.GetLength(1); row++)
            {
                var cell = grid[col, row];
                if (cell.X == 0 && cell.Y == 0) continue;
                cell.M.Should().BeGreaterThan(minM - range);
                cell.M.Should().BeLessThan(maxM + range);
            }
        }
    }

    // ──────────────────────────────────────────
    // Constructor validation
    // ──────────────────────────────────────────

    [TestMethod]
    public void Constructor_NegativeTension_Throws()
    {
        Action act = () => new SplineTensionGrass(tension: -1);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [TestMethod]
    public void Constructor_ZeroSegMax_Throws()
    {
        Action act = () => new SplineTensionGrass(segMax: 0);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [TestMethod]
    public void Constructor_TooFewNpMin_Throws()
    {
        Action act = () => new SplineTensionGrass(npMin: 3);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [TestMethod]
    public void InterpolateToGrid_TooFewPoints_Throws()
    {
        var grass = new SplineTensionGrass();
        var points = new List<CoordinateM>
        {
            new(0, 0, 1),
            new(1, 0, 2),
            new(0, 1, 3),
        };

        Action act = () => grass.InterpolateToGrid(points, cellSize: 0.5);
        act.Should().Throw<ArgumentException>();
    }

    // ──────────────────────────────────────────
    // CrstBasis verification (matches GRASS IL_crst)
    // ──────────────────────────────────────────

    [TestMethod]
    public void CrstBasis_AtZero_ReturnsZero()
    {
        SplineTensionGrass.CrstBasis(0.0, 40.0).Should().Be(0.0);
    }

    [TestMethod]
    public void CrstBasis_ContinuousAtX1()
    {
        double fi = 2.0;
        double below = SplineTensionGrass.CrstBasis(0.999 / (fi * fi / 4.0), fi);
        double above = SplineTensionGrass.CrstBasis(1.001 / (fi * fi / 4.0), fi);
        double diff = Math.Abs(above - below);
        diff.Should().BeLessThan(0.002,
            "Basis function should be approximately continuous at x=1 boundary (GRASS has ~1.3e-3 gap)");
    }

    [TestMethod]
    public void CrstBasis_KnownValues()
    {
        // IL_crst(r²=1, fi=2): x = fi²*r²/4 = 4*1/4 = 1.0
        // Result ≈ E1(1) + gamma + ln(1) = 0.21938 + 0.57722 + 0 ≈ 0.7966
        double val = SplineTensionGrass.CrstBasis(1.0, 2.0);
        val.Should().BeApproximately(0.7966, 0.001);
    }

    [TestMethod]
    public void CrstBasis_LargeX_ApproachesGammaPlusLogX()
    {
        // For x > 25, E1 ≈ 0, so result ≈ gamma + ln(x)
        double fi = 20.0;
        double rSq = 1.0; // x = fi²*r²/4 = 100
        double val = SplineTensionGrass.CrstBasis(rSq, fi);
        double expected = 0.57721566 + Math.Log(fi * fi * rSq / 4.0);
        val.Should().BeApproximately(expected, 0.01);
    }

    // ──────────────────────────────────────────
    // GRASS defaults
    // ──────────────────────────────────────────

    [TestMethod]
    public void DefaultParameters_MatchGrassDefaults()
    {
        var grass = new SplineTensionGrass();
        grass.Tension.Should().Be(40.0);
        grass.Smoothing.Should().Be(0.1);
        grass.SegMax.Should().Be(40);
        grass.NpMin.Should().Be(300);
        grass.KMax2.Should().Be(600);
    }

    // ──────────────────────────────────────────
    // Helpers
    // ──────────────────────────────────────────

    private static List<CoordinateM> CreateTestPoints()
    {
        return
        [
            new(0.0, 0.0, 50.0),
            new(1.0, 0.0, 58.4),
            new(2.0, 0.0, 59.1),
            new(0.0, 1.0, 55.4),
            new(1.0, 1.0, 63.8),
            new(2.0, 1.0, 64.5),
            new(0.0, 2.0, 45.8),
            new(1.0, 2.0, 54.3),
            new(2.0, 2.0, 55.0),
            new(0.5, 0.5, 56.0),
            new(1.5, 1.5, 58.0),
        ];
    }
}
