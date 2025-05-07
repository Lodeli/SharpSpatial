using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;

namespace SharpSpatial.Model
{
    /// <summary>
    /// Represents a solution for the nearest points calculation
    /// </summary>
    public class NearestSolution
    {
        /// <summary>
        /// The first point, the one on the first geometry
        /// </summary>
        public Coordinate Point1 { get; set; }

        /// <summary>
        /// The second point, the one on the second geometry
        /// </summary>
        public Coordinate Point2 { get; set; }

        /// <summary>
        /// The distance between the two points
        /// </summary>
        public double Distance { get; set; }

        /// <summary>
        /// Creates a new <see cref="NearestSolution"/> with the given points and distance
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="distance"></param>
        public NearestSolution(Coordinate point1, Coordinate point2, double distance)
        {
            this.Point1 = point1;
            this.Point2 = point2;
            this.Distance = distance;
        }
    }
}
