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
    public class SharpGeography : BaseGeometry
    {
        public int SRID
        {
            get => _geometry.SRID;
            set => _geometry.SRID = value;
        }

        public double ToleranceUsedForDensify
        {
            get => _toleranceUsedForDensify;
            set => _toleranceUsedForDensify = value;
        }

        #region C'tors
        /// <summary>
        /// Creates an empty <see cref="SharpGeography"/> with the given SRID
        /// </summary>
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
        #endregion

        #region Public methods
        /// <summary>
        /// Returns the point at the given index if it exists, otherwise null
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
        /// Returns the Latitude coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double Lat => _geometry.Coordinate.Y;

        /// <summary>
        /// Returns the Longitude coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double Lon => _geometry.Coordinate.X;

        /// <summary>
        /// Returns the envelope of this geography
        /// </summary>
        /// <returns></returns>
        public SharpGeography GetEnvelope()
        {
            Geometry envelope = GeoHelper.GetEnvelope(_geometry);
            return GeoHelper.CreateGeography(envelope, this.SRID);
        }

        /// <summary>
        /// Returns the geography at the given index
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public SharpGeography GetGeographyN(int index) => GeoHelper.CreateGeography(_geometry.GetGeometryN(index), this.SRID);

        /// <summary>
        /// Returns the geodesic buffer of this geography
        /// </summary>
        /// <param name="distance">The width of the buffer</param>
        /// <returns></returns>
        public SharpGeography Buffer(double distance) => GeoHelper.GetBuffer(this, distance);

        /// <summary>
        /// Returns the geodesic distance between this geography and another one with <see cref="GeodesicPrecision.Precise"/>" geodesic precision/>
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public DistanceSolution? Distance(SharpGeography other) => GeoHelper.GetDistance(this, other, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns the geodesic distance between this geography and another one with the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="other"></param>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public DistanceSolution? Distance(SharpGeography other, GeodesicPrecision geodesicPrecision)
            => GeoHelper.GetDistance(this, other, geodesicPrecision);

        /// <summary>
        /// Returns the nearest points between this geography and another one with <see cref="GeodesicPrecision.Precise"/>" geodesic precision/>
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public NearestSolution? GetNearestPoints(SharpGeography other) => GeoHelper.GetNearestPoints(this, other, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns the nearest points between this geography and another one with the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="other"></param>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public NearestSolution? GetNearestPoints(SharpGeography other, GeodesicPrecision geodesicPrecision)
            => GeoHelper.GetNearestPoints(this, other, geodesicPrecision);

        /// <summary>
        /// Returns the shortest line between this geography and another one with <see cref="GeodesicPrecision.Precise"/>" geodesic precision/>
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
        /// Returns the intersection between this geography and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography Intersection(SharpGeography other) => GeoHelper.GetIntersection(this, other);

        /// <summary>
        /// Returns true if this geography intersects with another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Intersects(SharpGeography other) => GeoHelper.Intersect(this, other);

        /// <summary>
        /// Returns the mercator projection of this geography
        /// </summary>
        /// <returns></returns>
        public SharpGeometry? Project() => GeoHelper.Project(this);

        /// <summary>
        /// Returns the reduced geography with the given tolerance
        /// </summary>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public SharpGeography Reduce(double tolerance) => GeoHelper.GetReduced(this, tolerance);

        /// <summary>
        /// Returns the densified geography, based on Mercator projected version of this geography, with the given distance tolerance.
        /// Useful to force the resulting geography to follow the Mercator projection.
        /// </summary>
        /// <param name="distanceTolerance"></param>
        /// <returns></returns>
        public SharpGeography DensifyProj(double distanceTolerance) => GeoHelper.GetDensifiedProj(this, distanceTolerance);

        /// <summary>
        /// Returns the union between this geography and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography? Union(SharpGeography other) => GeoHelper.GetUnion(this, other);

        /// <summary>
        /// Returns the difference between this geography and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeography Difference(SharpGeography other) => GeoHelper.GetDifference(this, other);

        /// <summary>
        /// Returns the length of this geography with <see cref="GeodesicPrecision.Precise"/>" geodesic precision/>
        /// </summary>
        /// <returns></returns>
        public double GetLength() => GeoHelper.GetLength(this, GeodesicPrecision.Precise);

        /// <summary>
        /// Returns the length of this geography with the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="geodesicPrecision"></param>
        /// <returns></returns>
        public double GetLength(GeodesicPrecision geodesicPrecision) => GeoHelper.GetLength(this, geodesicPrecision);

        /// <summary>
        /// Returns the geodesic approximate area of this geography (with the Haversine formula)
        /// </summary>
        /// <returns></returns>
        public double GetArea() => GeoHelper.GetArea(this);

        /// <summary>
        /// Returns the fixed geography: only if it crosses the International Date Line, it will be shifted to positive longitudes
        /// </summary>
        /// <returns></returns>
        public SharpGeography GetInternationalDateLineFix() => GeoHelper.GetInternationalDateLineFix(this);

        /// <summary>
        /// Returns true if this geography crosses the International Date Line
        /// </summary>
        /// <returns></returns>
        public bool IsCrossingIDL() => GeoHelper.IsCrossingIDL(this);

        public SharpGeography? MakeValid() => GeoHelper.MakeValid(this);

        public SharpGeography Revert() => GeoHelper.Revert(this);
        #endregion
    }
}
