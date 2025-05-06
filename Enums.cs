using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SharpSpatial
{
    /// <summary>
    /// Enumeration of the different types of geodesic distance calculations.
    /// </summary>
    public enum GeodesicPrecision
    {
        /// <summary>
        /// Geodesic distances are calculated using the Haversine formula: less accurate but faster.
        /// </summary>
        Fast,
        /// <summary>
        /// Geodesic distances are calculated using the Vincenty formula: more accurate but slower.
        /// </summary>
        Precise,
    }
}
