using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;

namespace SharpSpatial.Model
{
    public class NearestSolution
    {
        public Coordinate Point1 { get; set; }
        public Coordinate Point2 { get; set; }
        public double Distance { get; set; }

        public NearestSolution(Coordinate point1, Coordinate point2, double distance)
        {
            this.Point1 = point1;
            this.Point2 = point2;
            this.Distance = distance;
        }
    }
}
