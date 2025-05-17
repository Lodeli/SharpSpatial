using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;
using SharpSpatial.Helpers;
using SharpSpatial.Model;

namespace SharpSpatial
{
    /// <summary>
    /// Represents a geography object, which is a geometry with a geodesic precision
    /// </summary>
    public class SharpGeography : BaseGeometry
    {
        /// <summary>
        /// Returns the Latitude coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double Lat => _geometry.Coordinate.Y;

        /// <summary>
        /// Returns the Longitude coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double Lon => _geometry.Coordinate.X;

        /// <summary>
        /// Returns or sets the SRID of this geography
        /// </summary>
        public int SRID
        {
            get => _geometry.SRID;
            set => _geometry.SRID = value;
        }

        /// <summary>
        /// Returns or sets the tolerance used for densifying this geography
        /// </summary>
        public double ToleranceUsedForDensify
        {
            get => _toleranceUsedForDensify;
            set => _toleranceUsedForDensify = value;
        }

        #region C'tors
        /// <summary>
        /// Creates an empty <see cref="SharpGeography"/> with the given SRID
        /// </summary>
        /// <param name="srid">The SRID of the geometry</param>
        public SharpGeography(int srid) : base(null, srid, 0)
        { }

        /// <summary>
        /// Creates a new <see cref="SharpGeography"/> with the given <see cref="Geometry"/> shape and the given SRID and eventually a tolerance already used for densification
        /// </summary>
        /// <param name="geometry">The <see cref="Geometry"/> shape for this Geometry</param>
        /// <param name="srid">The SRID of the geometry</param>
        /// <param name="toleranceUsedForDensify">The tolerance used for densifying the geometry (zero, if not already densified)</param>
        public SharpGeography(Geometry geometry, int srid, double toleranceUsedForDensify = 0) : base(geometry, srid, toleranceUsedForDensify)
        { }

        /// <summary>
        /// Creates a new <see cref="SharpGeography"/> from the given Well-Known Text (WKT) string with the given SRID and eventually a tolerance already used for densification
        /// </summary>
        /// <param name="wkt">WKT to parse</param>
        /// <param name="srid">The SRID of the geometry</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Value of tolerance to be used for geodesic densify, 0 to not densify result</param>
        public SharpGeography(string wkt, int srid, bool? makeValid, double toleranceForDensify = 0)
            : base(GeoHelper.CreateGeomFromWKT(wkt, makeValid, toleranceForDensify), srid, toleranceForDensify)
        { }
        #endregion

        #region Public methods
        /// <summary>
        /// Returns the <see cref="SharpGeography"/> point at the given index if it exists, otherwise null
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public SharpGeography? GetPointN(int index)
        {
            Coordinate? pointCoords = this.GetCoordinateN(index);
            if (pointCoords == null)
                return null;

            return GeoHelper.CreatePoint(pointCoords, this.SRID);
        }

        /// <summary>
        /// Returns the <see cref="SharpGeography"/> at the given index if this is a collection of geographies
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public SharpGeography GetGeographyN(int index) => GeoHelper.CreateGeography(_geometry.GetGeometryN(index), this.SRID);

        /// <summary>
        /// Returns the envelope (bounding box) of this geography
        /// </summary>
        /// <returns></returns>
        public SharpGeography GetEnvelope()
        {
            Geometry envelope = GeoHelper.GetEnvelope(_geometry);
            return GeoHelper.CreateGeography(envelope, this.SRID);
        }

        /// <summary>
        /// Returns the nearest points, and their distance, between this geography and another one with <see cref="GeodesicPrecision.Precise"/>" geodesic precision
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public NearestSolution? GetNearestPoints(SharpGeography other) => GeoHelper.GetNearestPoints(this, other, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns the nearest points, and their distance, between this geography and another one with the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="other"></param>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public NearestSolution? GetNearestPoints(SharpGeography other, GeodesicPrecision geodesicPrecision)
            => GeoHelper.GetNearestPoints(this, other, geodesicPrecision);

        /// <summary>
        /// Returns the shortest line between this geography and another one with <see cref="GeodesicPrecision.Precise"/>" geodesic precision
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography? GetShortestLineTo(SharpGeography other) => GeoHelper.GetShortestLine(this, other, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns the shortest line between this geography and another one with the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="other"></param>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public SharpGeography? GetShortestLineTo(SharpGeography other, GeodesicPrecision geodesicPrecision)
            => GeoHelper.GetShortestLine(this, other, geodesicPrecision);

        /// <summary>
        /// Returns geodesic distance, initial and final bearing between this geography and another one with <see cref="GeodesicPrecision.Precise"/>" geodesic precision
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public DistanceSolution? Distance(SharpGeography other) => GeoHelper.GetDistance(this, other, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns geodesic distance, initial and final bearing between this geography and another one with the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="other"></param>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public DistanceSolution? Distance(SharpGeography other, GeodesicPrecision geodesicPrecision)
            => GeoHelper.GetDistance(this, other, geodesicPrecision);

        /// <summary>
        /// Returns the geodesic buffer of this geography
        /// </summary>
        /// <param name="distance">The width of the buffer in meters</param>
        /// <returns></returns>
        public SharpGeography Buffer(double distance) => GeoHelper.GetBuffer(this, distance);

        /// <summary>
        /// Returns true if this geography intersects the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Intersects(SharpGeography other) => GeoHelper.Intersect(this, other);

        /// <summary>
        /// Returns the intersection between this geography and the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography Intersection(SharpGeography other) => GeoHelper.GetIntersection(this, other);

        /// <summary>
        /// Returns the union between this geography and the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography? Union(SharpGeography other) => GeoHelper.GetUnion(this, other);

        /// <summary>
        /// Returns the difference between this geography and the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography Difference(SharpGeography other) => GeoHelper.GetDifference(this, other);

        /// <summary>
        /// Returns the simplified geography of this one, using the given tolerance in meters
        /// </summary>
        /// <param name="tolerance">Tolerance to use, in meters</param>
        /// <returns></returns>
        public SharpGeography Simplify(double tolerance) => GeoHelper.GetSimplified(this, tolerance);

        /// <summary>
        /// Returns the densified <see cref="SharpGeography"/>, based on Mercator projected version of this one, using the given distance tolerance.
        /// Useful to force the resulting geography to follow the Mercator projection.
        /// </summary>
        /// <param name="distanceTolerance">The distance tolerance to use, in meters</param>
        /// <returns></returns>
        public SharpGeography DensifyProj(double distanceTolerance) => GeoHelper.GetDensifiedProj(this, distanceTolerance);

        /// <summary>
        /// Returns the densified geography of this one, using <see cref="GeoHelper.DEFAULT_BEARING_TOLERANCE"/> bearing tolerance
        /// </summary>
        /// <returns></returns>
        public SharpGeography Densify() => GeoHelper.GetDensified(this);

        /// <summary>
        /// Returns the length (perimeter for polygon, length for linear shapes) of this geography using the <see cref="GeodesicPrecision.Precise"/>" geodesic precision
        /// </summary>
        /// <returns></returns>
        public double GetLength() => GeoHelper.GetLength(this, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns the length (perimeter for polygon, length for linear shapes) of this geography using the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public double GetLength(GeodesicPrecision geodesicPrecision) => GeoHelper.GetLength(this, geodesicPrecision);

        /// <summary>
        /// Returns the geodesic approximate area of this geography (with Haversine formula)
        /// </summary>
        /// <returns></returns>
        public double GetArea() => GeoHelper.GetArea(this);

        /// <summary>
        /// Returns true if this geography crosses the International Date Line (IDL or Antimeridian)
        /// </summary>
        /// <returns></returns>
        public bool IsCrossingIDL() => GeoHelper.IsCrossingIDL(this);

        /// <summary>
        /// Returns a fixed geography for the International Date Line (IDL or Antimeridian) by shifting longitudes &lt; 0 by +360
        /// </summary>
        /// <returns></returns>
        public SharpGeography GetInternationalDateLineFix() => GeoHelper.GetInternationalDateLineFix(this);

        /// <summary>
        /// Computes a new geometry which has all component coordinate sequences in reverse order (opposite orientation) to this one
        /// </summary>
        /// <returns></returns>
        public SharpGeography Revert() => GeoHelper.Revert(this);

        /// <summary>
        /// Returns the fixed geography from this one
        /// </summary>
        /// <returns></returns>
        public SharpGeography? MakeValid() => GeoHelper.MakeValid(this);

        /// <summary>
        /// Returns a Mercator projected <see cref="SharpGeometry"/> from this geography
        /// </summary>
        /// <returns></returns>
        public SharpGeometry? Project() => GeoHelper.Project(this);
        #endregion
    }
}
