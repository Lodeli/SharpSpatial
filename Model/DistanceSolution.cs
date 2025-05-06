using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SharpSpatial.Model
{
    public class DistanceSolution
    {
        public static readonly DistanceSolution Zero = new(null, 0, 0, 0);
        public static readonly DistanceSolution NaN = new(null, double.NaN, double.NaN, double.NaN);

        public GeodesicPrecision? Precision { get; set; }
        public double Distance { get; private set; }
        public double InitialBearing { get; private set; }
        public double DestinationBearing { get; private set; }

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

        public bool IsNaN()
        {
            return DistanceSolution.IsNaN(this);
        }

        public static bool IsNaN(DistanceSolution distanceSolution)
        {
            return double.IsNaN(distanceSolution.Distance)
                || double.IsNaN(distanceSolution.InitialBearing)
                || double.IsNaN(distanceSolution.DestinationBearing);
        }
    }
}
