using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;
using SharpSpatial.Helpers;

namespace SharpSpatial.Model
{
    public abstract class BaseGeometry
    {
        protected readonly Geometry _geometry;
        protected double _toleranceUsedForDensify;

        public Geometry Geo => _geometry;
        public Coordinate? Coordinate => _geometry?.Coordinate;
        public int NumPoints => _geometry?.Coordinates.Length ?? 0;
        public bool IsEmpty => _geometry.IsEmpty;
        public int NumGeometries => _geometry.NumGeometries;

        #region C'tors
        /// <summary>
        /// Creates a new <see cref="BaseGeometry"/> with the given <see cref="Geometry"/> shape and the given <see cref="GeoHelper"/> to execute spatial calculations
        /// </summary>
        /// <param name="geometry">The <see cref="Geometry"/>Shape for this Geometry</param>
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

        public Coordinate? GetStartCoordinate() => this.GetCoordinateN(0);

        public Coordinate? GetEndCoordinate() => (this.NumPoints > 0) ? this.GetCoordinateN(this.NumPoints - 1) : null;

        public string ToWKT()
        {
            var wktWriter = new NetTopologySuite.IO.WKTWriter();
            return wktWriter.Write(_geometry);
        }

        public byte[] ToWKB()
        {
            var wkbWriter = new NetTopologySuite.IO.WKBWriter();
            return wkbWriter.Write(_geometry);
        }

        public override string ToString() => this.ToWKT();

        public bool IsWellOriented() => GeoHelper.IsWellOriented(this.Geo);
        #endregion
    }
}
