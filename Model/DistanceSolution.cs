using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SharpSpatial.Model
{
    /// <summary>
    /// Represents the solution of a distance calculation
    /// </summary>
    public class DistanceSolution
    {
        /// <summary>
        /// Represents a solution with zero distance
        /// </summary>
        public static readonly DistanceSolution Zero = new(null, 0, 0, 0);

        /// <summary>
        /// Represents a non-valid solution
        /// </summary>
        public static readonly DistanceSolution NaN = new(null, double.NaN, double.NaN, double.NaN);

        /// <summary>
        /// The precision used for the calculation
        /// </summary>
        public GeodesicPrecision? Precision { get; set; }

        /// <summary>
        /// The distance between the two points
        /// </summary>
        public double Distance { get; private set; }

        /// <summary>
        /// The initial bearing at the starting point
        /// </summary>
        public double InitialBearing { get; private set; }

        /// <summary>
        /// The final bearing at the destination point
        /// </summary>
        public double DestinationBearing { get; private set; }

        /// <summary>
        /// Creates a new <see cref="DistanceSolution"/> with the given precision, distance, initial bearing, and final bearing
        /// </summary>
        /// <param name="precision"></param>
        /// <param name="distance"></param>
        /// <param name="initialBearing"></param>
        /// <param name="finalBearing"></param>
        public DistanceSolution(GeodesicPrecision? precision, double distance, double initialBearing, double finalBearing)
        {
            this.Precision = precision;
            initialBearing = initialBearing % 360;
            finalBearing = finalBearing % 360;
            this.Distance = distance;
            this.InitialBearing = (double.IsNaN(initialBearing) || initialBearing >= 0)
                ? initialBearing
                : initialBearing + 360;
            this.DestinationBearing = (double.IsNaN(finalBearing) || finalBearing >= 0)
                ? finalBearing
                : finalBearing + 360;
        }

        /// <summary>
        /// Checks if this <see cref="DistanceSolution"/> is a valid solution
        /// </summary>
        /// <returns></returns>
        public bool IsNaN()
        {
            return DistanceSolution.IsNaN(this);
        }

        /// <summary>
        /// Checks if the given <see cref="DistanceSolution"/> is a valid solution
        /// </summary>
        /// <param name="distanceSolution"></param>
        /// <returns></returns>
        public static bool IsNaN(DistanceSolution distanceSolution)
        {
            return double.IsNaN(distanceSolution.Distance)
                || double.IsNaN(distanceSolution.InitialBearing)
                || double.IsNaN(distanceSolution.DestinationBearing);
        }
    }
}
