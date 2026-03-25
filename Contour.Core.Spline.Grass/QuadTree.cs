using NetTopologySuite.Geometries;

namespace Contour.Core.Spline.Grass;

/// <summary>
/// Quadtree for point data matching GRASS GIS v.surf.rst's MT_tree/quaddata.
/// Supports incremental insertion with per-leaf spatial density filtering,
/// grid-based subdivision (halving rows/cols), and strict-inequality range search.
/// Reference: lib/rst/qtree/qtree.c, lib/rst/data/dataquad.c
/// </summary>
internal sealed class QuadTree
{
    public double XMin { get; private set; }
    public double YMin { get; private set; }
    public double XMax { get; private set; }
    public double YMax { get; private set; }

    /// <summary>Grid rows covered by this node (halved on subdivision, matching GRASS).</summary>
    public int NRows { get; }

    /// <summary>Grid columns covered by this node (halved on subdivision, matching GRASS).</summary>
    public int NCols { get; }

    public QuadTree?[]? Children { get; private set; }
    public List<int>? PointIndices { get; private set; }
    public bool IsLeaf => Children == null;

    // Quadrant indices: NW=0, NE=1, SW=2, SE=3
    private const int NW = 0;
    private const int NE = 1;
    private const int SW = 2;
    private const int SE = 3;

    private QuadTree(double xMin, double yMin, double xMax, double yMax,
        int nRows = 0, int nCols = 0)
    {
        XMin = xMin;
        YMin = yMin;
        XMax = xMax;
        YMax = yMax;
        NRows = nRows;
        NCols = nCols;
        PointIndices = new List<int>();
    }

    /// <summary>
    /// Builds a quadtree covering the grid region, inserting points one at a time
    /// with per-leaf density checking. Matches GRASS GIS MT_tree construction:
    /// tree covers the computational region, leaves hold at most kmax points,
    /// and points closer than sqrt(dminSq) to an existing point are rejected.
    /// </summary>
    public static QuadTree BuildFromGrid(double[] xs, double[] ys, int nPoints,
        double gridXMin, double gridYMin, double gridXMax, double gridYMax,
        int gridRows, int gridCols, int kmax, double dminSq, out int nFiltered)
    {
        var root = new QuadTree(gridXMin, gridYMin, gridXMax, gridYMax,
            gridRows, gridCols);
        nFiltered = 0;
        for (int i = 0; i < nPoints; i++)
        {
            double px = xs[i], py = ys[i];
            if (px < gridXMin || px > gridXMax || py < gridYMin || py > gridYMax)
            {
                nFiltered++;
                continue;
            }
            if (root.Insert(xs, ys, i, kmax, dminSq) == 0)
                nFiltered++;
        }
        return root;
    }

    /// <summary>
    /// Legacy build method for unit tests (batch partitioning, no density check).
    /// </summary>
    public static QuadTree Build(List<CoordinateM> points, int segMax)
    {
        double xMin = double.MaxValue, yMin = double.MaxValue;
        double xMax = double.MinValue, yMax = double.MinValue;
        foreach (var p in points)
        {
            if (p.X < xMin) xMin = p.X;
            if (p.Y < yMin) yMin = p.Y;
            if (p.X > xMax) xMax = p.X;
            if (p.Y > yMax) yMax = p.Y;
        }

        var indices = Enumerable.Range(0, points.Count).ToList();
        var root = new QuadTree(xMin, yMin, xMax, yMax);
        root.SubdivideLegacy(points, indices, segMax);
        return root;
    }

    private void SubdivideLegacy(List<CoordinateM> points, List<int> indices, int segMax)
    {
        if (indices.Count <= segMax)
        {
            PointIndices = indices;
            return;
        }

        double midX = (XMin + XMax) / 2.0;
        double midY = (YMin + YMax) / 2.0;

        if (Math.Abs(XMax - XMin) < 1e-15 && Math.Abs(YMax - YMin) < 1e-15)
        {
            PointIndices = indices;
            return;
        }

        Children = new QuadTree[4];
        Children[NW] = new QuadTree(XMin, midY, midX, YMax);
        Children[NE] = new QuadTree(midX, midY, XMax, YMax);
        Children[SW] = new QuadTree(XMin, YMin, midX, midY);
        Children[SE] = new QuadTree(midX, YMin, XMax, midY);

        var childIndices = new List<int>[4] { new(), new(), new(), new() };

        foreach (int idx in indices)
        {
            var p = points[idx];
            int quad;
            if (p.Y >= midY)
                quad = p.X < midX ? NW : NE;
            else
                quad = p.X < midX ? SW : SE;
            childIndices[quad].Add(idx);
        }

        for (int i = 0; i < 4; i++)
        {
            if (childIndices[i].Count > 0)
                Children[i]!.SubdivideLegacy(points, childIndices[i], segMax);
            else
                Children[i]!.PointIndices = new List<int>();
        }
    }

    /// <summary>
    /// Inserts a point into the tree. Returns 0 if rejected (too close to existing
    /// point or out of bounds), 1 if added. Matches GRASS MT_insert.
    /// </summary>
    private int Insert(double[] xs, double[] ys, int index, int kmax, double dminSq)
    {
        if (!IsLeaf)
        {
            // Internal node: navigate to correct child (GRASS quad_compare)
            int quad = GetQuadrant(xs[index], ys[index]);
            if (quad < 0) return 0;
            return Children![quad]!.Insert(xs, ys, index, kmax, dminSq);
        }

        // Leaf full (n_points >= kmax): subdivide then re-insert
        // Matches GRASS quad_division_check: returns 1 when n_points >= kmax
        if (PointIndices != null && PointIndices.Count >= kmax)
        {
            if (NCols <= 1 || NRows <= 1)
            {
                // Can't subdivide further — just add with density check
                return AddWithDensityCheck(xs, ys, index, dminSq);
            }

            Subdivide(xs, ys, kmax, dminSq);
            return Insert(xs, ys, index, kmax, dminSq);
        }

        // Leaf with room: add with density check
        return AddWithDensityCheck(xs, ys, index, dminSq);
    }

    /// <summary>
    /// Adds a point if it passes the density check. Matches GRASS quad_add_data:
    /// checks r = dx*dx + dy*dy against dmin (squared), using &lt;= comparison.
    /// </summary>
    private int AddWithDensityCheck(double[] xs, double[] ys, int index, double dminSq)
    {
        double px = xs[index], py = ys[index];
        if (PointIndices != null)
        {
            foreach (int existing in PointIndices)
            {
                double dx = xs[existing] - px;
                double dy = ys[existing] - py;
                if (dx * dx + dy * dy <= dminSq)
                    return 0; // too close — matches GRASS r <= dmin
            }
        }
        PointIndices ??= new List<int>();
        PointIndices.Add(index);
        return 1;
    }

    /// <summary>
    /// Returns the quadrant for a point. Matches GRASS quad_compare exactly:
    /// uses grid-based midpoint with odd-count rounding (cols1 = n_cols/2+1 for odd).
    /// </summary>
    private int GetQuadrant(double px, double py)
    {
        double ewRes = NCols > 0 ? (XMax - XMin) / NCols : (XMax - XMin);
        double nsRes = NRows > 0 ? (YMax - YMin) / NRows : (YMax - YMin);

        // GRASS: odd counts round up for the midpoint index
        int cols = (NCols % 2 == 0) ? NCols / 2 : NCols / 2 + 1;
        int rows = (NRows % 2 == 0) ? NRows / 2 : NRows / 2 + 1;

        double xMid = XMin + ewRes * cols;
        double yMid = YMin + nsRes * rows;

        // GRASS quad_compare: >= for all comparisons
        bool inX = px >= XMin;
        bool rightHalf = px >= xMid;
        bool inY = py >= YMin;
        bool upperHalf = py >= yMid;

        if (inX && inY)
        {
            if (rightHalf && upperHalf) return NE;
            if (rightHalf) return SE;
            if (upperHalf) return NW;
            return SW;
        }
        return -1; // out of bounds
    }

    /// <summary>
    /// Subdivides a leaf into 4 children. Matches GRASS quad_divide_data exactly:
    /// odd cols give cols1=n/2+1 to SW/NW, cols2=n/2 to NE/SE.
    /// Points are redistributed with density check (GRASS calls quad_add_data).
    /// </summary>
    private void Subdivide(double[] xs, double[] ys, int kmax, double dminSq)
    {
        double ewRes = NCols > 0 ? (XMax - XMin) / NCols : (XMax - XMin);
        double nsRes = NRows > 0 ? (YMax - YMin) / NRows : (YMax - YMin);

        int cols1, cols2, rows1, rows2;
        if (NCols % 2 == 0) { cols1 = NCols / 2; cols2 = cols1; }
        else { cols2 = NCols / 2; cols1 = cols2 + 1; }
        if (NRows % 2 == 0) { rows1 = NRows / 2; rows2 = rows1; }
        else { rows2 = NRows / 2; rows1 = rows2 + 1; }

        double xm = XMin + cols1 * ewRes;
        double ym = YMin + rows1 * nsRes;

        // GRASS child assignment:
        // NE: (xm, ym, xr, yr, rows2, cols2)
        // SW: (xl, yl, xm, ym, rows1, cols1)
        // SE: (xm, yl, xr, ym, rows1, cols2)
        // NW: (xl, ym, xm, yr, rows2, cols1)
        Children = new QuadTree[4];
        Children[NW] = new QuadTree(XMin, ym, xm, YMax, rows2, cols1);
        Children[NE] = new QuadTree(xm, ym, XMax, YMax, rows2, cols2);
        Children[SW] = new QuadTree(XMin, YMin, xm, ym, rows1, cols1);
        Children[SE] = new QuadTree(xm, YMin, XMax, ym, rows1, cols2);

        // Redistribute existing points (GRASS calls quad_add_data during redistribution)
        var oldPoints = PointIndices!;
        PointIndices = null; // no longer a leaf

        foreach (int idx in oldPoints)
        {
            int quad = GetQuadrant(xs[idx], ys[idx]);
            if (quad >= 0)
            {
                // GRASS calls quad_add_data with density check during redistribution.
                // Since all existing points already passed pairwise checks, this is
                // effectively a no-op, but we do it to match GRASS exactly.
                Children[quad]!.AddWithDensityCheck(xs, ys, idx, dminSq);
            }
        }
    }

    /// <summary>
    /// Collects all leaf nodes in traversal order (matching GRASS recursive traversal).
    /// </summary>
    public List<QuadTree> GetLeaves()
    {
        var leaves = new List<QuadTree>();
        CollectLeaves(leaves);
        return leaves;
    }

    private void CollectLeaves(List<QuadTree> leaves)
    {
        if (IsLeaf)
        {
            // Include all leaves (even empty ones) so every grid cell gets evaluated.
            // GRASS GIS evaluates all segments via bisection window search.
            if (NRows > 0 && NCols > 0)
                leaves.Add(this);
            return;
        }

        foreach (var child in Children!)
        {
            child?.CollectLeaves(leaves);
        }
    }

    /// <summary>
    /// Finds the X-width of the smallest non-empty leaf. Matches GRASS smallest_segment.
    /// </summary>
    public double SmallestLeafWidth()
    {
        if (IsLeaf)
            return (NRows > 0 && NCols > 0) ? XMax - XMin : double.MaxValue;

        double smallest = double.MaxValue;
        foreach (var child in Children!)
        {
            if (child != null)
            {
                double s = child.SmallestLeafWidth();
                if (s < smallest) smallest = s;
            }
        }
        return smallest;
    }

    /// <summary>
    /// Shifts all node bounds by (-dx, -dy). Matches GRASS translate_quad:
    /// after all points are inserted, coordinates are shifted so the minimum
    /// point is near zero, improving floating-point precision in segment processing.
    /// </summary>
    public void Translate(double dx, double dy)
    {
        XMin -= dx;
        YMin -= dy;
        XMax -= dx;
        YMax -= dy;

        if (Children != null)
        {
            foreach (var child in Children)
                child?.Translate(dx, dy);
        }
    }

    /// <summary>
    /// Finds all point indices within the given bounding box using strict inequality.
    /// Matches GRASS quad_get_points: point-&gt;x &gt; xmin &amp;&amp; point-&gt;x &lt; xmax.
    /// </summary>
    public void FindPointsInRegion(double[] xs, double[] ys, double qxMin, double qyMin,
        double qxMax, double qyMax, List<int> result, int maxPoints)
    {
        // GRASS quad_intersect: AABB overlap test
        if (qxMax <= XMin || qxMin >= XMax || qyMax <= YMin || qyMin >= YMax)
            return;
        if (result.Count >= maxPoints)
            return;

        if (IsLeaf)
        {
            if (PointIndices == null) return;
            foreach (int idx in PointIndices)
            {
                if (result.Count >= maxPoints) break;
                double px = xs[idx], py = ys[idx];
                // GRASS quad_get_points: strict inequality
                if (px > qxMin && px < qxMax && py > qyMin && py < qyMax)
                    result.Add(idx);
            }
            return;
        }

        foreach (var child in Children!)
        {
            child?.FindPointsInRegion(xs, ys, qxMin, qyMin, qxMax, qyMax, result, maxPoints);
        }
    }

    /// <summary>
    /// Finds all point indices within the given bounding box by traversing the tree.
    /// Uses non-strict inequality (legacy, for backward compatibility with unit tests).
    /// </summary>
    public void FindPointsInRegion(List<CoordinateM> points, double qxMin, double qyMin,
        double qxMax, double qyMax, List<int> result, int maxPoints)
    {
        if (qxMax < XMin || qxMin > XMax || qyMax < YMin || qyMin > YMax)
            return;
        if (result.Count >= maxPoints)
            return;

        if (IsLeaf)
        {
            if (PointIndices == null) return;
            foreach (int idx in PointIndices)
            {
                if (result.Count >= maxPoints) break;
                var p = points[idx];
                if (p.X >= qxMin && p.X <= qxMax && p.Y >= qyMin && p.Y <= qyMax)
                    result.Add(idx);
            }
            return;
        }

        foreach (var child in Children!)
        {
            child?.FindPointsInRegion(points, qxMin, qyMin, qxMax, qyMax, result, maxPoints);
        }
    }
}
