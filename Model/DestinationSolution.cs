using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;

namespace SharpSpatial.Model
{
    public class DestinationSolution
    {
        public GeodesicPrecision? Precision { get; set; }
        public Coordinate Coordinate { get; set; }
        public double Bearing { get; set; }

        public DestinationSolution(GeodesicPrecision? precision, Coordinate coordinate, double bearing)
        {
            this.Precision = precision;
            this.Coordinate = new(coordinate);
            this.Bearing = bearing;
        }
    }
}
