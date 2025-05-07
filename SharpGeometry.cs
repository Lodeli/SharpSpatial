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
    /// Represents a geometry, tipically in the Mercator projection
    /// </summary>
    public class SharpGeometry : BaseGeometry
    {
        /// <summary>
        /// Returns the X coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double X => _geometry.Coordinate.X;

        /// <summary>
        /// Returns the Y coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double Y => _geometry.Coordinate.Y;

        #region C'tors
        /// <summary>
        /// Creates an empty <see cref="SharpGeometry"/>
        /// </summary>
        public SharpGeometry() : base(null, 0, 0)
        { }

        /// <summary>
        /// Creates a new <see cref="SharpGeometry"/> with the given <see cref="Geometry"/> shape
        /// </summary>
        /// <param name="geometry">The <see cref="Geometry"/> shape for this <see cref="SharpGeometry"/></param>
        public SharpGeometry(Geometry geometry) : base(geometry, 0, 0)
        { }

        /// <summary>
        /// Creates a new <see cref="SharpGeometry"/> from the given Well-Known Text (WKT) string
        /// </summary>
        /// <param name="wkt">WKT to parse</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        public SharpGeometry(string wkt, bool? makeValid) : this(GeoHelper.CreateGeom(wkt, makeValid, 0))
        { }
        #endregion

        #region Public methods
        /// <summary>
        /// Returns the <see cref="SharpGeometry"/> point at the given index if it exists, otherwise null
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public SharpGeometry? GetPointN(int index)
        {
            Coordinate? pointCoords = this.GetCoordinateN(index);
            if (pointCoords == null)
                return null;

            return GeoHelper.CreatePoint(pointCoords);
        }

        /// <summary>
        /// Returns the <see cref="SharpGeometry"/> at the given index if this is a collection of geometries
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public SharpGeometry GetGeometryN(int index) => GeoHelper.CreateGeometry(_geometry.GetGeometryN(index));

        /// <summary>
        /// Returns the envelope (bounding box) of this geometry
        /// </summary>
        /// <returns></returns>
        public SharpGeometry GetEnvelope()
        {
            Geometry envelope = GeoHelper.GetEnvelope(_geometry);
            return GeoHelper.CreateGeometry(envelope);
        }

        /// <summary>
        /// Returns the nearest points, and their distance, between this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public NearestSolution? GetNearestPoints(SharpGeometry other) => GeoHelper.GetNearestPoints(this, other);

        /// <summary>
        /// Returns the shortest line between this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry? GetShortestLineTo(SharpGeometry other) => GeoHelper.GetShortestLine(this, other);

        /// <summary>
        /// Returns distance, initial and final bearing between this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public DistanceSolution? Distance(SharpGeometry other) => GeoHelper.GetDistance(this, other);

        /// <summary>
        /// Returns the planar buffer of this geometry
        /// </summary>
        /// <param name="distance">The width of the buffer</param>
        /// <returns></returns>
        public SharpGeometry Buffer(double distance) => GeoHelper.GetBuffer(this, distance);

        /// <summary>
        /// Returns true if this geometry intersects the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Intersects(SharpGeometry other) => GeoHelper.Intersect(this, other);

        /// <summary>
        /// Returns the intersection between this geometry and the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry Intersection(SharpGeometry other) => GeoHelper.GetIntersection(this, other);

        /// <summary>
        /// Returns the union between this geometry and the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry? Union(SharpGeometry other) => GeoHelper.GetUnion(this, other);

        /// <summary>
        /// Returns the difference between this geometry and the given one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry Difference(SharpGeometry other) => GeoHelper.GetDifference(this, other);

        /// <summary>
        /// Returns the simplified geometry of this one, using the given tolerance
        /// </summary>
        /// <param name="tolerance">Tolerance to use</param>
        /// <returns></returns>
        public SharpGeometry Simplify(double tolerance) => GeoHelper.GetSimplified(this, tolerance);

        /// <summary>
        /// Returns the densified geometry of this one, using the given distance tolerance
        /// </summary>
        /// <param name="distanceTolerance"></param>
        /// <returns></returns>
        public SharpGeometry Densify(double distanceTolerance) => GeoHelper.GetDensified(this, distanceTolerance);

        /// <summary>
        /// Returns the length (perimeter for polygon, length for linear shapes) of this geometry.
        /// </summary>
        /// <returns></returns>
        public double GetLength() => GeoHelper.GetLength(this);

        /// <summary>
        /// Returns the area of this geometry
        /// </summary>
        /// <returns></returns>
        public double GetArea() => GeoHelper.GetArea(this);

        /// <summary>
        /// Returns true if this geometry crosses the International Date Line (IDL or Antimeridian)
        /// </summary>
        /// <returns></returns>
        public bool IsCrossingIDL() => GeoHelper.IsCrossingIDL(this);

        /// <summary>
        /// Returns a fixed geometry for the International Date Line (IDL or Antimeridian) by shifting longitudes &lt; 0 by +360
        /// </summary>
        /// <returns></returns>
        public SharpGeometry GetInternationalDateLineFix() => GeoHelper.GetInternationalDateLineFix(this);

        /// <summary>
        /// Computes a new geometry which has all component coordinate sequences in reverse order (opposite orientation) to this one
        /// </summary>
        /// <returns></returns>
        public SharpGeometry Revert() => GeoHelper.Revert(this);

        /// <summary>
        /// Returns the fixed geometry from this one
        /// </summary>
        /// <returns></returns>
        public SharpGeometry? MakeValid() => GeoHelper.MakeValid(this);

        /// <summary>
        /// Returns an unprojected <see cref="SharpGeography"/>, with the given SRID, from this (Mercator) geometry
        /// </summary>
        /// <param name="srid"></param>
        /// <returns></returns>
        public SharpGeography? Unproject(int srid) => GeoHelper.Unproject(this, srid);
        #endregion
    }
}
