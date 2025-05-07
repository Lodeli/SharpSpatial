using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;
using NetTopologySuite.Precision;
using SharpSpatial.Helpers;

namespace SharpSpatial.Model
{
    /// <summary>
    /// Base class for all geometries/geographies
    /// </summary>
    public abstract class BaseGeometry
    {
        /// <summary>
        /// Private field for the internal <see cref="Geometry"/> shape for this <see cref="BaseGeometry"/>
        /// </summary>
        protected readonly Geometry _geometry;

        /// <summary>
        /// The tolerance used for densifying the geometry, used to avoid densifying when not needed
        /// </summary>
        protected double _toleranceUsedForDensify;

        /// <summary>
        /// The internal <see cref="Geometry"/> shape for this entity
        /// </summary>
        public Geometry Geo => _geometry;

        /// <summary>
        /// Returns a vertex of this geometry (usually, but not necessarily, the first one), or <c>null</c> if the geometry is empty
        /// </summary>
        public Coordinate? Coordinate => _geometry?.Coordinate;

        /// <summary>
        /// Returns the total number of points in this geometry
        /// </summary>
        public int NumPoints => _geometry?.Coordinates.Length ?? 0;

        /// <summary>
        /// Tests whether the set of points covered in this <c>Geometry</c> is empty.
        /// <para/>
        /// Note this test is for topological emptiness, not structural emptiness.<br/>
        /// A collection containing only empty elements is reported as empty.<br/>
        /// To check structural emptiness use <see cref="NumGeometries"/>
        /// </summary>
        public bool IsEmpty => _geometry.IsEmpty;

        /// <summary>
        /// Returns the number of Geometryes in a GeometryCollection, or 1 if the geometry is not a collection
        /// </summary>
        public int NumGeometries => _geometry.NumGeometries;

        #region C'tors
        /// <summary>
        /// Creates a new <see cref="BaseGeometry"/> with the given <see cref="Geometry"/> shape and the given <see cref="GeoHelper"/> to execute spatial calculations
        /// </summary>
        /// <param name="geometry">The <see cref="Geometry"/>Shape for this Geometry</param>
        /// <param name="srid">SRID of the geometry</param>
        /// <param name="toleranceUsedForDensify"> The tolerance used for densifying the given geometry</param>
        public BaseGeometry(Geometry? geometry, int srid, double toleranceUsedForDensify)
        {
            _geometry = geometry?.Copy() ?? GeoHelper.CreateEmptyGeo();
            _geometry.SRID = srid;
            _toleranceUsedForDensify = toleranceUsedForDensify;
        }
        #endregion

        #region Public methods
        /// <summary>
        /// Returns the point at the given index if it exists, otherwise null
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public Coordinate? GetCoordinateN(int index) => (index < _geometry.Coordinates.Length) ? _geometry.Coordinates[index] : null;

        /// <summary>
        /// Returns the first point of the geometry if it exists, otherwise null
        /// </summary>
        /// <returns></returns>
        public Coordinate? GetStartCoordinate() => this.GetCoordinateN(0);

        /// <summary>
        /// Returns the last point of the geometry if it exists, otherwise null
        /// </summary>
        /// <returns></returns>
        public Coordinate? GetEndCoordinate() => (this.NumPoints > 0) ? this.GetCoordinateN(this.NumPoints - 1) : null;

        /// <summary>
        /// Returns the Well-Known Text (WKT) representation of the geometry
        /// </summary>
        /// <param name="maxDigits">The maximum number of decimal places for returned values</param>
        /// <returns></returns>
        public string ToWKT(int? maxDigits = null)
        {
            if (maxDigits < 1) throw new ArgumentOutOfRangeException(nameof(maxDigits), "maxDigits must be greater than 0");

            var wktWriter = new NetTopologySuite.IO.WKTWriter();
            if (maxDigits.HasValue)
                wktWriter.PrecisionModel = new(scale: Math.Pow(10, maxDigits.Value - 1)); ;
            return wktWriter.Write(_geometry);
        }

        /// <summary>
        /// Returns the Well-Known Binary (WKB) representation of the geometry
        /// </summary>
        /// <returns></returns>
        public byte[] ToWKB()
        {
            var wkbWriter = new NetTopologySuite.IO.WKBWriter();
            return wkbWriter.Write(_geometry);
        }

        /// <summary>
        /// Returns the text representation of this entity
        /// </summary>
        /// <returns></returns>
        public override string ToString() => this.ToWKT();

        /// <summary>
        /// Returns true if this shape is well oriented (CCW)
        /// </summary>
        /// <returns></returns>
        public bool IsWellOriented() => GeoHelper.IsWellOriented(this.Geo);
        #endregion
    }
}
