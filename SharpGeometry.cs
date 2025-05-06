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
    public class SharpGeometry : BaseGeometry
    {
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
        #endregion

        /// <summary>
        /// Returns the point at the given index if it exists, otherwise null
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
        /// Returns the X coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double X => _geometry.Coordinate.X;

        /// <summary>
        /// Returns the Y coordinate of a vertex of the geometry (usually the first one)
        /// </summary>
        public double Y => _geometry.Coordinate.Y;

        /// <summary>
        /// Returns the envelope of this geometry
        /// </summary>
        /// <returns></returns>
        public SharpGeometry GetEnvelope()
        {
            Geometry envelope = GeoHelper.GetEnvelope(_geometry);
            return GeoHelper.CreateGeometry(envelope);
        }

        /// <summary>
        /// Returns the geometry at the given index
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public SharpGeometry GetGeometryN(int index) => GeoHelper.CreateGeometry(_geometry.GetGeometryN(index));

        /// <summary>
        /// Returns the planar buffer of this geometry
        /// </summary>
        /// <param name="distance">The width of the buffer</param>
        /// <returns></returns>
        public SharpGeometry Buffer(double distance) => GeoHelper.GetBuffer(this, distance);

        /// <summary>
        /// Returns the planar distance between this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public DistanceSolution? Distance(SharpGeometry other) => GeoHelper.GetDistance(this, other);

        /// <summary>
        /// Returns the nearest points, between this geometry and another one, and their distance
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
        /// Returns the intersection of this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry Intersection(SharpGeometry other) => GeoHelper.GetIntersection(this, other);

        /// <summary>
        /// Returns true if this geometry intersects another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Intersects(SharpGeometry other) => GeoHelper.Intersect(this, other);

        /// <summary>
        /// Returns the unprojected geography of this Mercator geometry
        /// </summary>
        /// <param name="srid"></param>
        /// <returns></returns>
        public SharpGeography? Unproject(int srid) => GeoHelper.Unproject(this, srid);

        /// <summary>
        /// Returns the reduced geometry of this geometry using the given tolerance
        /// </summary>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public SharpGeometry Reduce(double tolerance) => GeoHelper.GetReduced(this, tolerance);

        /// <summary>
        /// Returns the densified geometry of this geometry using the given distance tolerance
        /// </summary>
        /// <param name="distanceTolerance"></param>
        /// <returns></returns>
        public SharpGeometry Densify(double distanceTolerance) => GeoHelper.GetDensified(this, distanceTolerance);

        /// <summary>
        /// Returns the union between this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry? Union(SharpGeometry other) => GeoHelper.GetUnion(this, other);

        /// <summary>
        /// Returns the difference between this geometry and another one
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public SharpGeometry Difference(SharpGeometry other) => GeoHelper.GetDifference(this, other);

        /// <summary>
        /// Returns the length of this geometry
        /// </summary>
        /// <returns></returns>
        public double GetLength() => GeoHelper.GetLength(this);

        /// <summary>
        /// Returns the area of this geometry
        /// </summary>
        /// <returns></returns>
        public double GetArea() => GeoHelper.GetArea(this);

        /// <summary>
        /// Returns the fixed geometry: only if it crosses the International Date Line, it will be shifted to positive longitudes
        /// </summary>
        /// <returns></returns>
        public SharpGeometry GetInternationalDateLineFix() => GeoHelper.GetInternationalDateLineFix(this);

        /// <summary>
        /// Returns true if this geometry crosses the International Date Line
        /// </summary>
        /// <returns></returns>
        public bool IsCrossingIDL() => GeoHelper.IsCrossingIDL(this);

        public SharpGeometry? MakeValid() => GeoHelper.MakeValid(this);

        public SharpGeometry Revert() => GeoHelper.Revert(this);
    }
}
