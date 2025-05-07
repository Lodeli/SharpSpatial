using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;

namespace SharpSpatial.Model
{
    /// <summary>
    /// Represents the solution of a destination calculation
    /// </summary>
    public class DestinationSolution
    {
        /// <summary>
        /// The precision used for the calculation
        /// </summary>
        public GeodesicPrecision? Precision { get; set; }

        /// <summary>
        /// The coordinate of the destination point
        /// </summary>
        public Coordinate Coordinate { get; set; }

        /// <summary>
        /// The bearing resulting at the destination point
        /// </summary>
        public double Bearing { get; set; }

        /// <summary>
        /// Creates a new <see cref="DestinationSolution"/> with the given precision, coordinate, and bearing
        /// </summary>
        /// <param name="precision"></param>
        /// <param name="coordinate"></param>
        /// <param name="bearing"></param>
        public DestinationSolution(GeodesicPrecision? precision, Coordinate coordinate, double bearing)
        {
            this.Precision = precision;
            this.Coordinate = new(coordinate);
            this.Bearing = bearing;
        }
    }
}
