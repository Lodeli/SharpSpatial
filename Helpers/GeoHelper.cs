using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using GeoAPI.CoordinateSystems;
using GeoAPI.CoordinateSystems.Transformations;
using NetTopologySuite;
using NetTopologySuite.Features;
using NetTopologySuite.Geometries;
using NetTopologySuite.Geometries.Utilities;
using NetTopologySuite.IO;
using NetTopologySuite.Operation.Polygonize;
using NetTopologySuite.Simplify;
using ProjNet.CoordinateSystems;
using ProjNet.CoordinateSystems.Transformations;
using SharpSpatial.Model;

namespace SharpSpatial.Helpers
{
    /// <summary>
    /// Class with static methods to create geometries and perform spatial calculations
    /// </summary>
    public static class GeoHelper
    {
        #region Consts
        /// <summary>
        /// Default tolerance for the bearing variation when densifying.
        /// </summary>
        public const double DEFAULT_BEARING_TOLERANCE = 0.05;

        /// <summary>
        /// Radius of Earth in meters for Haversine formula
        /// </summary>
        const double HAVERSINE_EARTH_R_M = 6371008.8; // Mean radius of Earth in meters
        /// <summary>
        /// Major semi-axes of Earth WGS84 ellipsoid for Vincenty formula
        /// </summary>
        private const double VINCENTY_ELLIPSOID_WGS84_A = 6378137;
        /// <summary>
        /// Minor semi-axes of Earth WGS84 ellipsoid for Vincenty formula
        /// </summary>
        private const double VINCENTY_ELLIPSOID_WGS84_B = 6356752.314245;
        /// <summary>
        /// Flattening of Earth WGS84 ellipsoid for Vincenty formula
        /// </summary>
        private const double VINCENTY_ELLIPSOID_WGS84_F = 1 / 298.257223563;
        #endregion

        #region Fields
        private static readonly ICoordinateTransformationFactory s_ctFact;
        private static readonly WKTReader s_wktReader;
        private static readonly WKBReader s_wkbReader;
        private static readonly GeoJsonReader s_geoJsonReader;
        private static readonly GeometryFactory s_geometryFactory;
        private static readonly ICoordinateTransformation s_ctProj;
        private static readonly ICoordinateTransformation s_ctUnproj;
        #endregion

        /// <summary>
        /// Class with static methods to create geometries and perform spatial calculations
        /// </summary>
        static GeoHelper()
        {
            s_ctFact = new CoordinateTransformationFactory();
            s_geometryFactory = new GeometryFactory(new PrecisionModel(PrecisionModels.Floating));
            NtsGeometryServices ntsGeometryServices = new(new PrecisionModel(PrecisionModels.Floating), 4326);
            s_wktReader = new(ntsGeometryServices);
            s_wkbReader = new(ntsGeometryServices);
            s_geoJsonReader = new();

            IGeographicCoordinateSystem sourceCS = GeographicCoordinateSystem.WGS84;
            IProjectedCoordinateSystem targetCS = ProjectedCoordinateSystem.WebMercator;
            s_ctProj = s_ctFact.CreateFromCoordinateSystems(sourceCS, targetCS);
            s_ctUnproj = s_ctFact.CreateFromCoordinateSystems(targetCS, sourceCS);
        }

        #region Public methods
        #region Creation methods
        /// <summary>
        /// Creates an empty <see cref="Geometry"/> object.
        /// </summary>
        /// <returns></returns>
        public static Geometry CreateEmptyGeo() => s_geometryFactory.CreateLineString();

        /// <summary>
        /// Creates an empty <see cref="SharpGeometry"/> object.
        /// </summary>
        /// <returns></returns>
        public static SharpGeometry CreateEmpty() => new();

        /// <summary>
        /// Creates an empty <see cref="SharpGeometry"/> object with the given SRID.
        /// </summary>
        /// <param name="srid"></param>
        /// <returns></returns>
        public static SharpGeography CreateEmpty(int srid) => new(srid);

        /// <summary>
        /// Creates a <see cref="SharpGeometry"/> Point from a <see cref="Coordinate"/> object.
        /// </summary>
        /// <param name="coordinate"></param>
        /// <returns></returns>
        public static SharpGeometry CreatePoint(Coordinate coordinate) => new(s_geometryFactory.CreatePoint(coordinate));

        /// <summary>
        /// Creates a <see cref="SharpGeography"/> Point from a <see cref="Coordinate"/> object with the given SRID.
        /// </summary>
        /// <param name="coordinate"></param>
        /// <param name="srid"></param>
        /// <returns></returns>
        public static SharpGeography CreatePoint(Coordinate coordinate, int srid) => new(s_geometryFactory.CreatePoint(coordinate), srid);

        /// <summary>
        /// Creates a <see cref="SharpGeometry"/> LineString from an array of <see cref="Coordinate"/> objects.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static SharpGeometry CreateLineString(Coordinate[] points) => new(s_geometryFactory.CreateLineString(points));

        /// <summary>
        /// Creates a <see cref="SharpGeography"/> LineString from an array of <see cref="Coordinate"/> objects with the given SRID.
        /// </summary>
        /// <param name="srid"></param>
        /// <param name="points"></param>
        /// <returns></returns>
        public static SharpGeography CreateLineString(int srid, params Coordinate[] points) => new(s_geometryFactory.CreateLineString(points), srid);

        /// <summary>
        /// Creates a <see cref="SharpGeometry"/> Polygon from an array of <see cref="Coordinate"/> objects.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static SharpGeometry CreatePolygon(Coordinate[] points) => new(s_geometryFactory.CreatePolygon(points));

        /// <summary>
        /// Creates a <see cref="SharpGeography"/> Polygon from an array of <see cref="Coordinate"/> objects with the given SRID.
        /// </summary>
        /// <param name="srid"></param>
        /// <param name="points"></param>
        /// <returns></returns>
        public static SharpGeography CreatePolygon(int srid, params Coordinate[] points) => new(s_geometryFactory.CreatePolygon(points), srid);

        /// <summary>
        /// Creates a <see cref="SharpGeometry"/> from a <see cref="Geometry"/> object.
        /// </summary>
        /// <param name="geometry"></param>
        /// <returns></returns>
        public static SharpGeometry CreateGeometry(Geometry geometry) => new(geometry);

        /// <summary>
        /// Creates a <see cref="SharpGeography"/> from a <see cref="Geometry"/> object with the given SRID.
        /// </summary>
        /// <param name="geometry"></param>
        /// <param name="srid"></param>
        /// <returns></returns>
        public static SharpGeography CreateGeography(Geometry geometry, int srid) => new(geometry, srid);

        /// <summary>
        /// Creates a <see cref="Geometry"/> from a Well-Known Text (WKT) string
        /// </summary>
        /// <param name="wkt">The WKT to parse</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Value of tolerance to be used for geodesic densify, 0 to not densify result</param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static Geometry CreateGeomFromWKT(string wkt, bool? makeValid, double toleranceForDensify = 0)
        {
            Geometry result = s_wktReader.Read(wkt);
            return FixParsedGeometry(result, makeValid, toleranceForDensify);
        }

        /// <summary>
        /// Creates a <see cref="SharpGeometry"/> from a Well-Known Text (WKT) string
        /// </summary>
        /// <param name="wkt">The WKT to parse</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <returns></returns>
        public static SharpGeometry CreateGeometryFromWKT(string wkt, bool? makeValid)
        {
            Geometry geo = CreateGeomFromWKT(wkt, makeValid, 0);
            return new SharpGeometry(geo);
        }

        /// <summary>
        /// Creates a new <see cref="SharpGeography"/> from the given Well-Known Text (WKT) string with the given SRID and eventually a tolerance to use for densification
        /// </summary>
        /// <param name="wkt">The WKT to parse</param>
        /// <param name="srid">The SRID to assign to the resulting geography</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static SharpGeography CreateGeographyFromWKT(string wkt, int srid, bool? makeValid, double toleranceForDensify = DEFAULT_BEARING_TOLERANCE)
        {
            Geometry result = CreateGeomFromWKT(wkt, makeValid, toleranceForDensify);
            return new SharpGeography(result, srid, toleranceForDensify);
        }

        /// <summary>
        /// Creates a <see cref="Geometry"/> from a Well-Known Binary (WKB) byte array
        /// </summary>
        /// <param name="wkb">The WKB to parse</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Value of tolerance to be used for geodesic densify, 0 to not densify result</param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        /// <returns></returns>
        public static Geometry CreateGeomFromWKB(byte[] wkb, bool? makeValid, double toleranceForDensify = 0)
        {
            Geometry result = s_wkbReader.Read(wkb);
            return FixParsedGeometry(result, makeValid, toleranceForDensify);
        }

        /// <summary>
        /// Creates a <see cref="SharpGeometry"/> from a Well-Known Binary (WKB) byte array
        /// </summary>
        /// <param name="wkb">The WKB to parse</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <returns></returns>
        public static SharpGeometry CreateGeometryFromWKB(byte[] wkb, bool? makeValid)
        {
            Geometry geo = CreateGeomFromWKB(wkb, makeValid, 0);
            return new SharpGeometry(geo);
        }

        /// <summary>
        /// Creates a new <see cref="SharpGeography"/> from the given Well-Known Binary (WKB) byte array with the given SRID and eventually a tolerance to use for densification
        /// </summary>
        /// <param name="wkb">The WKB to parse</param>
        /// <param name="srid">The SRID to assign to the resulting geography</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static SharpGeography CreateGeographyFromWKB(byte[] wkb, int srid, bool? makeValid, double toleranceForDensify = DEFAULT_BEARING_TOLERANCE)
        {
            Geometry result = CreateGeomFromWKB(wkb, makeValid, toleranceForDensify);
            return new SharpGeography(result, srid, toleranceForDensify);
        }

        /// <summary>
        /// Creates a list of <see cref="Geometry"/> from the given GeoJSON string
        /// </summary>
        /// <param name="geoJson"></param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static List<Geometry> CreateGeomFromGeoJson(string geoJson, bool? makeValid, double toleranceForDensify = 0)
        {
            List<Geometry> result = [];
            FeatureCollection features = s_geoJsonReader.Read<FeatureCollection>(geoJson);
            foreach (Feature feature in features)
            {
                Geometry geom = feature.Geometry;
                if (geom != null)
                {
                    result.Add(FixParsedGeometry(geom, makeValid, toleranceForDensify));
                }
            }
            return result;
        }

        /// <summary>
        /// Creates a list of <see cref="SharpGeometry"/> from the given GeoJSON string
        /// </summary>
        /// <param name="geoJson"></param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <returns></returns>
        public static List<SharpGeometry> CreateGeometryFromGeoJson(string geoJson, bool? makeValid)
        {
            List<SharpGeometry> result = [];
            List<Geometry> geometries = CreateGeomFromGeoJson(geoJson, makeValid, 0);
            foreach (Geometry geometry in geometries)
            {
                result.Add(new SharpGeometry(geometry));
            }

            return result;
        }

        /// <summary>
        /// Creates a list of <see cref="SharpGeography"/> from the given GeoJSON string with the given SRID and eventually a tolerance to  use for densification
        /// </summary>
        /// <param name="geoJson"></param>
        /// <param name="srid">The SRID to assign to the resulting geography</param>
        /// <param name="makeValid">False: do not perform any check, True: force a MakeValid on the geometry, Null: raise an exception if invalid</param>
        /// <param name="toleranceForDensify">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static List<SharpGeography> CreateGeographyFromGeoJson(string geoJson, int srid, bool? makeValid, double toleranceForDensify = DEFAULT_BEARING_TOLERANCE)
        {
            List<SharpGeography> result = [];
            List<Geometry> geometries = CreateGeomFromGeoJson(geoJson, makeValid, toleranceForDensify);
            foreach (Geometry geometry in geometries)
            {
                result.Add(new SharpGeography(geometry, srid));
            }
            return result;
        }
        #endregion

        /// <summary>
        /// Returns the envelope (bounding box) of the given <see cref="Geometry"/>
        /// </summary>
        /// <param name="geometry"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static Geometry GetEnvelope(Geometry geometry)
        {
            if (geometry == null || geometry.IsEmpty)
                throw new ArgumentException("Geometry must not be null or empty.", nameof(geometry));

            // Compute the envelope (bounding box)
            Geometry envelope = geometry.Envelope;

            // Ensure correct orientation if it's a polygon
            if (envelope.OgcGeometryType == OgcGeometryType.Polygon)
            {
                return IsWellOriented(envelope) ? envelope : envelope.Reverse();
            }

            // For non-polygon geometries, return as-is
            return envelope;
        }

        /// <summary>
        /// Returns the nearest points between the given two <see cref="Geometry"/> and their distance
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <param name="geodesicPrecision">Eventual geodesic precision to apply</param>
        /// <returns></returns>
        public static NearestSolution? GetNearestPoints(Geometry geom1, Geometry geom2, GeodesicPrecision? geodesicPrecision)
        {
            NearestSolution? result = null;

            Coordinate[] shell1 = GetExternalShell(geom1);
            Coordinate[] shell2 = GetExternalShell(geom2);

            if (shell1.Length < 1 || shell2.Length < 1) return result;

            double minDistance = double.MaxValue;
            Coordinate? nearestPoint1 = null;
            Coordinate? nearestPoint2 = null;

            // Check all points of geom1 against all segments of geom2
            for (int i = 0; i < shell1.Length; i++)
            {
                Coordinate point = shell1[i];
                for (int j = 0; j < shell2.Length - 1; j++)
                {
                    Coordinate segStart = shell2[j];
                    Coordinate segEnd = shell2[j + 1];

                    Coordinate projection = GetNearestPointOnSegment(point, segStart, segEnd);
                    double dist = geodesicPrecision == null
                        ? GetPlanarDistance(point.X, point.Y, projection.X, projection.Y).Distance
                        : GetHaversineDistance(point.Y, point.X, projection.Y, projection.X).Distance; // Useless using Vincenty here

                    if (dist < minDistance)
                    {
                        minDistance = dist;
                        nearestPoint1 = point;
                        nearestPoint2 = projection;
                    }
                }
            }

            // Check all points of geom2 against all segments of geom1
            for (int i = 0; i < shell2.Length; i++)
            {
                Coordinate point = shell2[i];
                for (int j = 0; j < shell1.Length - 1; j++)
                {
                    Coordinate segStart = shell1[j];
                    Coordinate segEnd = shell1[j + 1];

                    Coordinate projection = GetNearestPointOnSegment(point, segStart, segEnd);
                    double dist = geodesicPrecision == null
                        ? GetPlanarDistance(point.X, point.Y, projection.X, projection.Y).Distance
                        : GetHaversineDistance(point.Y, point.X, projection.Y, projection.X).Distance; // Useless using Vincenty here

                    if (dist < minDistance)
                    {
                        minDistance = dist;
                        nearestPoint1 = projection;
                        nearestPoint2 = point;
                    }
                }
            }

            return nearestPoint1 != null && nearestPoint2 != null
                ? new NearestSolution(nearestPoint1, nearestPoint2, minDistance)
                : null;
        }

        /// <summary>
        /// Returns the nearest points, and their distance, between the given two <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <returns></returns>
        public static NearestSolution? GetNearestPoints(SharpGeometry geom1, SharpGeometry geom2)
        {
            return GetNearestPoints(geom1.Geo, geom1.Geo, null);
        }

        /// <summary>
        /// Returns the nearest points, and their distance, between the given two <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geog1"></param>
        /// <param name="geog2"></param>
        /// <param name="geodesicPrecision">Geodesic precision to apply</param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static NearestSolution? GetNearestPoints(SharpGeography geog1, SharpGeography geog2,
            GeodesicPrecision geodesicPrecision, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            ThrowOnInvalidGeogArgs(geog1, geog2);

            geog1 = Densify(geog1, densifyWithBearingTolerance);
            geog2 = Densify(geog2, densifyWithBearingTolerance);

            return GetNearestPoints(geog1.Geo, geog2.Geo, geodesicPrecision);
        }

        /// <summary>
        /// Returns the shortest line between two <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <returns></returns>
        public static SharpGeometry? GetShortestLine(SharpGeometry geom1, SharpGeometry geom2)
        {
            NearestSolution? nearestPoints = GetNearestPoints(geom1, geom2);
            if (nearestPoints == null) return null;

            return CreateLineString([nearestPoints.Point1, nearestPoints.Point2]);
        }

        /// <summary>
        /// Returns the shortest line between two <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geog1"></param>
        /// <param name="geog2"></param>
        /// <param name="geodesicPrecision">Geodesic precision to apply</param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static SharpGeography? GetShortestLine(SharpGeography geog1, SharpGeography geog2,
            GeodesicPrecision geodesicPrecision, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            ThrowOnInvalidGeogArgs(geog1, geog2);

            geog1 = Densify(geog1, densifyWithBearingTolerance);
            geog2 = Densify(geog2, densifyWithBearingTolerance);

            NearestSolution? nearestSolution = GetNearestPoints(geog1.Geo, geog2.Geo, geodesicPrecision);
            if (nearestSolution == null) return null;

            return CreateLineString(geog1.SRID, [nearestSolution.Point1, nearestSolution.Point2]);
        }

        /// <summary>
        /// Returns distance, initial and final bearing between two points
        /// </summary>
        /// <param name="x_lon_1">X (or Longitude) coordinate of point 1</param>
        /// <param name="y_lat_1">Y (or Latitude) coordinate of point 1</param>
        /// <param name="x_lon_2">X (or Longitude) coordinate of point 2</param>
        /// <param name="y_lat_2">Y (or Latitude) coordinate of point 2</param>
        /// <param name="geodesicPrecision">Eventual geodesic precision to apply</param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public static DistanceSolution? GetDistance(double x_lon_1, double y_lat_1, double x_lon_2, double y_lat_2, GeodesicPrecision? geodesicPrecision)
        {
            if (geodesicPrecision == null)
                return GetPlanarDistance(x_lon_1, y_lat_1, x_lon_2, y_lat_2);

            switch (geodesicPrecision)
            {
                case GeodesicPrecision.Fast:
                    return GetHaversineDistance(y_lat_1, x_lon_1, y_lat_2, x_lon_2);
                case GeodesicPrecision.Precise:
                    return GetVincentyDistance(y_lat_1, x_lon_1, y_lat_2, x_lon_2);
                default:
                    throw new NotImplementedException(geodesicPrecision.ToString());
            }
        }

        /// <summary>
        /// Returns distance, initial and final bearing between two <see cref="Geometry"/>
        /// </summary>
        /// <param name="geo1"></param>
        /// <param name="geo2"></param>
        /// <param name="geodesicPrecision">Eventual geodesic precision to apply</param>
        /// <returns></returns>
        private static DistanceSolution? GetDistance(Geometry geo1, Geometry geo2, GeodesicPrecision? geodesicPrecision)
        {
            DistanceSolution? result = null;

            if (!geo1.Intersects(geo2))
            {
                NearestSolution? nearestSolution = GetNearestPoints(geo1, geo2, geodesicPrecision);
                if (nearestSolution != null)
                    result = GetDistance(nearestSolution.Point1, nearestSolution.Point2, geodesicPrecision);
            }
            else
            {
                result = new(geodesicPrecision, 0, 0, 0);
            }

            return result;
        }

        /// <summary>
        /// Returns distance, initial and final bearing between two <see cref="Coordinate"/>
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="geodesicPrecision">Eventual geodesic precision to apply</param>
        /// <returns></returns>
        public static DistanceSolution? GetDistance(Coordinate point1, Coordinate point2, GeodesicPrecision? geodesicPrecision)
        {
            return GetDistance(point1.X, point1.Y, point2.X, point2.Y, geodesicPrecision);
        }

        /// <summary>
        /// Returns distance, initial and final bearing between two <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <returns></returns>
        public static DistanceSolution? GetDistance(SharpGeometry geom1, SharpGeometry geom2)
        {
            return GetDistance(geom1.Geo, geom2.Geo, null);
        }

        /// <summary>
        /// Returns distance, initial and final bearing between two <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geog1"></param>
        /// <param name="geog2"></param>
        /// <param name="geodesicPrecision">Geodesic precision to apply</param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static DistanceSolution? GetDistance(SharpGeography geog1, SharpGeography geog2, GeodesicPrecision geodesicPrecision,
            double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            ThrowOnInvalidGeogArgs(geog1, geog2);

            geog1 = Densify(geog1, densifyWithBearingTolerance);
            geog2 = Densify(geog2, densifyWithBearingTolerance);

            return GetDistance(geog1.Geo, geog2.Geo, geodesicPrecision);
        }

        /// <summary>
        /// Returns the destination point from a given point, initial bearing and distance.
        /// </summary>
        /// <param name="startPoint">Coordinate of the starting point</param>
        /// <param name="initialBearing">Initial bearing (direction)</param>
        /// <param name="distance">Distance from starting point (in meters when geodesicPrecision is given)</param>
        /// <param name="geodesicPrecision">Eventual geodesic precision to apply</param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public static DestinationSolution? GetDestination(Coordinate startPoint, double initialBearing, double distance,
            GeodesicPrecision? geodesicPrecision)
        {
            if (geodesicPrecision == null)
                return GetPlanarDestination(startPoint, initialBearing, distance);

            switch (geodesicPrecision)
            {
                case GeodesicPrecision.Fast:
                    return GetHaversineDestination(startPoint, initialBearing, distance);
                case GeodesicPrecision.Precise:
                    return GetVincentyDestination(startPoint, initialBearing, distance);
                default:
                    throw new NotImplementedException(geodesicPrecision.ToString());
            }
        }

        /// <summary>
        /// Returns the planar buffer of the given <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom"></param>
        /// <param name="distance">The width of the buffer</param>
        /// <returns></returns>
        public static SharpGeometry GetBuffer(SharpGeometry geom, double distance)
        {
            return new SharpGeometry(geom.Geo.Buffer(distance));
        }

        /// <summary>
        /// Returns the geodesic buffer of the given <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geog"></param>
        /// <param name="distance">The width of the buffer in meters</param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static SharpGeography GetBuffer(SharpGeography geog, double distance, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            return new SharpGeography(CalculateBuffer(geog.Geo, distance, densifyWithBearingTolerance), geog.SRID, densifyWithBearingTolerance);
        }

        /// <summary>
        /// Returns true if the two given <see cref="SharpGeometry"/> intersect
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <returns></returns>
        public static bool Intersect(SharpGeometry geom1, SharpGeometry geom2)
        {
            return geom1.Geo.Intersects(geom2.Geo);
        }

        /// <summary>
        /// Returns true if the two given <see cref="SharpGeography"/> intersect
        /// </summary>
        /// <param name="geog1"></param>
        /// <param name="geog2"></param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static bool Intersect(SharpGeography geog1, SharpGeography geog2, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            ThrowOnInvalidGeogArgs(geog1, geog2);

            SharpGeography geog1Densified = Densify(geog1, densifyWithBearingTolerance);
            SharpGeography geog2Densified = Densify(geog2, densifyWithBearingTolerance);

            return geog1Densified.Geo.Intersects(geog2Densified.Geo);
        }

        /// <summary>
        /// Returns the intersection between the two given <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <returns></returns>
        public static SharpGeometry GetIntersection(SharpGeometry geom1, SharpGeometry geom2)
        {
            return new SharpGeometry(geom1.Geo.Intersection(geom2.Geo));
        }

        /// <summary>
        /// Returns the intersection between the two given <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geog1"></param>
        /// <param name="geog2"></param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static SharpGeography GetIntersection(SharpGeography geog1, SharpGeography geog2, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            ThrowOnInvalidGeogArgs(geog1, geog2);

            SharpGeography geog1Densified = Densify(geog1, densifyWithBearingTolerance);
            SharpGeography geog2Densified = Densify(geog2, densifyWithBearingTolerance);

            Geometry intersection = geog1Densified.Geo.Intersection(geog2Densified.Geo);

            return new SharpGeography(intersection, geog1.SRID, densifyWithBearingTolerance);
        }

        /// <summary>
        /// Returns the union of all the given <see cref="Geometry"/>
        /// </summary>
        /// <param name="geometries"></param>
        /// <returns></returns>
        public static Geometry? GetUnion(params Geometry[] geometries)
        {
            if (geometries == null || geometries.Length == 0)
                return null;

            // Get target SRID from the first non-null geometry
            int? targetSRID = geometries.FirstOrDefault(g => g != null)?.SRID;

            List<Geometry> validGeometries = [];

            foreach (Geometry geom in geometries)
            {
                if (geom == null || geom.IsEmpty)
                    continue;

                for (int i = 0; i < geom.NumGeometries; i++)
                {
                    Geometry current = geom.GetGeometryN(i);

                    // Optionally reproject or reject mismatched SRIDs
                    if (targetSRID.HasValue && current.SRID != targetSRID.Value)
                    {
                        // Example: force alignment (custom SRID transform could be applied here)
                        current = current.Copy();
                        current.SRID = targetSRID.Value;
                    }

                    // Attempt to fix invalid geometries if needed
                    if (!current.IsValid)
                    {
                        Geometry fixedGeom = MakeValid(current);
                        if (fixedGeom != null && !fixedGeom.IsEmpty)
                        {
                            current = fixedGeom;
                        }
                        else
                        {
                            continue; // Skip if still invalid
                        }
                    }

                    validGeometries.Add(current);
                }
            }

            if (validGeometries.Count == 0)
                return null;

            Geometry unionGeometry;

            if (validGeometries.Count == 1)
            {
                unionGeometry = validGeometries[0].Copy();
            }
            else
            {
                // More efficient way to combine multiple geometries
                unionGeometry = GeometryCombiner.Combine(validGeometries).Union();
            }

            unionGeometry.Normalize();
            return unionGeometry;
        }

        /// <summary>
        /// Returns the union of all the given <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geometries"></param>
        /// <returns></returns>
        public static SharpGeometry? GetUnion(params SharpGeometry[] geometries)
        {
            Geometry? unionGeometry = GetUnion(geometries.Select(g => g.Geo).ToArray());
            if (unionGeometry == null) return null;

            return new SharpGeometry(unionGeometry);
        }

        /// <summary>
        /// Returns the union of all the given <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geographies"></param>
        /// <returns></returns>
        public static SharpGeography? GetUnion(params SharpGeography[] geographies)
        {
            Geometry? unionGeometry = GetUnion(geographies.Select(g => g.Geo).ToArray());
            if (unionGeometry == null) return null;

            return new SharpGeography(unionGeometry, unionGeometry.SRID);
        }

        /// <summary>
        /// Returns the difference between the two given <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom1"></param>
        /// <param name="geom2"></param>
        /// <returns></returns>
        public static SharpGeometry GetDifference(SharpGeometry geom1, SharpGeometry geom2)
        {
            return new SharpGeometry(geom1.Geo.Difference(geom2.Geo));
        }

        /// <summary>
        /// Returns the difference between the two given <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geog1"></param>
        /// <param name="geog2"></param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static SharpGeography GetDifference(SharpGeography geog1, SharpGeography geog2, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            ThrowOnInvalidGeogArgs(geog1, geog2);

            SharpGeography geog1Densified = Densify(geog1, densifyWithBearingTolerance);
            SharpGeography geog2Densified = Densify(geog2, densifyWithBearingTolerance);
            Geometry difference = geog1Densified.Geo.Difference(geog2Densified.Geo);

            return new SharpGeography(difference, geog1.SRID);
        }

        /// <summary>
        /// Returns the simplified <see cref="SharpGeometry"/> of the given one, using the given tolerance
        /// </summary>
        /// <param name="geom"></param>
        /// <param name="tolerance">Tolerance to use</param>
        /// <returns></returns>
        public static SharpGeometry GetSimplified(SharpGeometry geom, double tolerance)
        {
            Geometry simplified = DouglasPeuckerSimplifier.Simplify(geom.Geo, tolerance);
            return new SharpGeometry(simplified);
        }

        /// <summary>
        /// Returns the simplified <see cref="SharpGeography"/> of the given one, using the given tolerance
        /// </summary>
        /// <param name="geog"></param>
        /// <param name="tolerance">Tolerance to use, in meters</param>
        /// <returns></returns>
        public static SharpGeography GetSimplified(SharpGeography geog, double tolerance)
        {
            double metersTolerance = tolerance / 111111.111111;
            Geometry simplified = Simplify(geog.Geo, metersTolerance);

            return new SharpGeography(simplified, geog.SRID);
        }

        /// <summary>
        /// Returns the densified <see cref="SharpGeometry"/> of the given one, using the given distance tolerance
        /// </summary>
        /// <param name="geom"></param>
        /// <param name="distanceTolerance">The distance tolerance to use</param>
        /// <returns></returns>
        public static SharpGeometry GetDensified(SharpGeometry geom, double distanceTolerance)
        {
            return new SharpGeometry(NetTopologySuite.Densify.Densifier.Densify(geom.Geo, distanceTolerance));
        }

        /// <summary>
        /// Returns the densified <see cref="SharpGeography"/> of the given one, using the given bearing tolerance
        /// </summary>
        /// <param name="geog"></param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static SharpGeography GetDensified(SharpGeography geog, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            return Densify(geog, densifyWithBearingTolerance);
        }

        /// <summary>
        /// Returns the densified <see cref="SharpGeography"/>, based on Mercator projected version of the given one, using the given distance tolerance in meters.
        /// Useful to force the resulting geography to follow the Mercator projection.
        /// </summary>
        /// <param name="geog"></param>
        /// <param name="distanceTolerance">The distance tolerance to use, in meters</param>
        /// <returns></returns>
        public static SharpGeography GetDensifiedProj(SharpGeography geog, double distanceTolerance)
        {
            List<Coordinate> resultCoordinates = [geog.Geo.Coordinates[0]];

            for (int p = 1; p < geog.Geo.Coordinates.Length; p++)
            {
                Coordinate point1 = geog.Geo.Coordinates[p - 1];
                Coordinate point2 = geog.Geo.Coordinates[p];
                SharpGeometry? projSegment = Project(CreateLineString(geog.SRID, point1, point2));

                if (projSegment != null)
                {
                    double dist = GetHaversineDistance(point1.Y, point1.X, point2.Y, point2.X).Distance;

                    int stepCount = (int)Math.Ceiling(dist / distanceTolerance);
                    if (stepCount > 1)
                    {
                        SharpGeometry? segP1 = projSegment.GetPointN(0);
                        if (segP1?.Coordinate != null)
                        {
                            SharpGeometry? segP2 = projSegment.GetPointN(1);
                            if (segP2?.Coordinate != null)
                            {
                                Coordinate point1Proj = segP1.Coordinate;
                                Coordinate point2Proj = segP2.Coordinate;
                                double stepX = (point2Proj.X - point1Proj.X) / stepCount;
                                double stepY = (point2Proj.Y - point1Proj.Y) / stepCount;
                                for (int s = 1; s < stepCount; s++)
                                {
                                    Coordinate newPointProj = new(point1Proj.X + stepX * s, point1Proj.Y + stepY * s);
                                    Coordinate newPoint = Unproject(newPointProj.X, newPointProj.Y);
                                    resultCoordinates.Add(newPoint);
                                }
                            }
                        }
                    }
                }

                resultCoordinates.Add(point2);
            }

            return CreateLineString(geog.SRID, resultCoordinates.ToArray());
        }

        /// <summary>
        /// Returns the length (perimeter for polygon, length for linear shapes) of the given <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom"></param>
        /// <returns></returns>
        public static double GetLength(SharpGeometry geom)
        {
            return geom.Geo.Length;
        }

        /// <summary>
        /// Returns the area of the given <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geom"></param>
        /// <returns></returns>
        public static double GetArea(SharpGeometry geom)
        {
            return geom.Geo.Area;
        }

        /// <summary>
        /// Returns the length (perimeter for polygon, length for linear shapes) of the given <see cref="SharpGeography"/> using the given <see cref="GeodesicPrecision"/>
        /// </summary>
        /// <param name="geog"></param>
        /// <param name="geodesicPrecision">The geodesic precision model to use</param>
        /// <param name="densifyWithBearingTolerance">Tolerance for the bearing variation when densifying, higher values reduce precision but improve speed, 0 to not densify result. Default value is <see cref="DEFAULT_BEARING_TOLERANCE"/></param>
        /// <returns></returns>
        public static double GetLength(SharpGeography geog, GeodesicPrecision geodesicPrecision, double densifyWithBearingTolerance = DEFAULT_BEARING_TOLERANCE)
        {
            double result = 0;

            if (geog.Geo.Coordinates.Length > 1)
            {
                for (int i = 1; i < geog.Geo.Coordinates.Length; i++)
                {
                    SharpGeography? p1 = geog.GetPointN(i - 1);
                    SharpGeography? p2 = geog.GetPointN(i);
                    if (p1 == null || p2 == null)
                        break;
                    result += GetDistance(p1, p2, geodesicPrecision, densifyWithBearingTolerance)?.Distance ?? 0;
                }
            }

            return result;
        }

        /// <summary>
        /// Returns the geodesic approximate area of this geography (with Haversine formula)
        /// </summary>
        /// <param name="geog"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static double GetArea(SharpGeography geog)
        {
            if (geog.Geo == null || geog.Geo.IsEmpty)
                return 0;

            double totalArea = 0;

            switch (geog.Geo.OgcGeometryType)
            {
                case OgcGeometryType.Polygon:
                    totalArea += ComputePolygonArea((Polygon)geog.Geo);
                    break;

                case OgcGeometryType.MultiPolygon:
                    var multi = (MultiPolygon)geog.Geo;
                    for (int i = 0; i < multi.NumGeometries; i++)
                    {
                        var poly = (Polygon)multi.GetGeometryN(i);
                        totalArea += ComputePolygonArea(poly);
                    }
                    break;

                default:
                    throw new ArgumentException("Geodesic area is only supported for Polygon or MultiPolygon geometries.");
            }

            return Math.Abs(totalArea); // in m²
        }

        #region Area support methods
        private static double ComputePolygonArea(Polygon polygon)
        {
            double area = ComputeRingArea(polygon.ExteriorRing);

            for (int i = 0; i < polygon.NumInteriorRings; i++)
            {
                LineString hole = polygon.GetInteriorRingN(i);
                area -= ComputeRingArea(hole);
            }

            return area;
        }

        private static double ComputeRingArea(LineString ring)
        {
            Coordinate[] coords = ring.Coordinates;

            if (!coords.First().Equals2D(coords.Last()))
                coords = coords.Append(coords.First()).ToArray();

            double total = 0;
            Coordinate origin = coords[0];

            for (int i = 1; i < coords.Length - 1; i++)
            {
                total += GetTriangleArea(origin, coords[i], coords[i + 1]);
            }

            return total;
        }

        private static double GetTriangleArea(Coordinate a, Coordinate b, Coordinate c)
        {
            double lat1 = ToRadians(a.Y), lon1 = ToRadians(a.X);
            double lat2 = ToRadians(b.Y), lon2 = ToRadians(b.X);
            double lat3 = ToRadians(c.Y), lon3 = ToRadians(c.X);

            double A = GetCentralAngle(lat1, lon1, lat2, lon2);
            double B = GetCentralAngle(lat2, lon2, lat3, lon3);
            double C = GetCentralAngle(lat3, lon3, lat1, lon1);

            double s = (A + B + C) / 2;
            double tanE = Math.Tan(s / 2) * Math.Tan((s - A) / 2) * Math.Tan((s - B) / 2) * Math.Tan((s - C) / 2);
            double E = 4 * Math.Atan(Math.Sqrt(Math.Abs(tanE)));

            return E * HAVERSINE_EARTH_R_M * HAVERSINE_EARTH_R_M;
        }

        private static double GetCentralAngle(double lat1, double lon1, double lat2, double lon2)
        {
            return Math.Acos(
                Math.Sin(lat1) * Math.Sin(lat2) +
                Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon2 - lon1)
            );
        }
        #endregion

        /// <summary>
        /// Returns true if the given <see cref="Geometry"/> crosses the International Date Line (IDL or Antimeridian).
        /// </summary>
        /// <param name="geometry"></param>
        /// <returns></returns>
        public static bool IsCrossingIDL(Geometry geometry)
        {
            for (int i = 0; i < geometry.NumPoints - 1; i++)
            {
                Coordinate point = geometry.Coordinates[i];
                Coordinate nextPoint = geometry.Coordinates[i + 1];

                if (IsCrossingIDL(point.X, nextPoint.X))
                    return true;
            }

            return false;
        }

        /// <summary>
        /// Returns true if the given <see cref="SharpGeography"/> crosses the International Date Line (IDL or Antimeridian).
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static bool IsCrossingIDL(params SharpGeography[] points)
        {
            if (points?.Any() == true)
            {
                for (int i = 0; i < points.Length - 1; i++)
                {
                    if (IsCrossingIDL(points[i].Lon, points[i + 1].Lon))
                        return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Returns true if the given <see cref="SharpGeometry"/> crosses the International Date Line (IDL or Antimeridian)
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static bool IsCrossingIDL(params SharpGeometry[] points)
        {
            if (points?.Any() == true)
            {
                for (int i = 0; i < points.Length - 1; i++)
                {
                    if (IsCrossingIDL(points[i].X, points[i + 1].X))
                        return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Returns a fixed <see cref="Geometry"/> for the International Date Line (IDL or Antimeridian) by shifting longitudes &lt; 0 by +360
        /// </summary>
        /// <param name="geometry"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="InvalidOperationException"></exception>
        public static Geometry GetInternationalDateLineFix(Geometry geometry)
        {
            if (geometry == null || geometry.IsEmpty)
                throw new ArgumentException("Geometry must not be null or empty.", nameof(geometry));

            // Ensure geometry is valid
            Geometry valid = MakeValid(geometry)?.Copy()
                ?? throw new InvalidOperationException("Cannot make the given geometry valid.");

            if (!IsWellOriented(valid))
                valid = valid.Reverse();

            if (!IsCrossingIDL(valid))
                return valid;

            // Shift all longitudes < 0 by +360
            Geometry fixedGeom = valid.Copy();
            fixedGeom.Apply(new CoordinateSequenceFilterAdapter(seq =>
            {
                for (int i = 0; i < seq.Count; i++)
                {
                    double x = seq.GetX(i);
                    if (x < 0)
                        seq.SetX(i, x + 360);
                }
            }));

            return fixedGeom;
        }

        /// <summary>
        /// Returns a fixed <see cref="SharpGeometry"/> for the International Date Line (IDL or Antimeridian) by shifting longitudes &lt; 0 by +360
        /// </summary>
        /// <param name="geom"></param>
        /// <returns></returns>
        public static SharpGeometry GetInternationalDateLineFix(SharpGeometry geom)
        {
            Geometry fixedGeo = GetInternationalDateLineFix(geom.Geo);
            return new SharpGeometry(fixedGeo);
        }

        /// <summary>
        /// Returns a fixed <see cref="SharpGeography"/> for the International Date Line (IDL or Antimeridian) by shifting longitudes &lt; 0 by +360
        /// </summary>
        /// <param name="geog"></param>
        /// <returns></returns>
        public static SharpGeography GetInternationalDateLineFix(SharpGeography geog)
        {
            Geometry fixedGeo = GetInternationalDateLineFix(geog.Geo);
            return new SharpGeography(fixedGeo, geog.SRID, geog.ToleranceUsedForDensify);
        }

        /// <summary>
        /// Returns true if the given <see cref="Geometry"/> is well oriented (CCW).
        /// </summary>
        /// <param name="geometry"></param>
        /// <returns></returns>
        public static bool IsWellOriented(Geometry geometry)
        {
            if (geometry.NumPoints < 4) return true;

            if (geometry.NumGeometries > 1)
            {
                for (int g = 0; g < geometry.NumGeometries; g++)
                {
                    if (!IsWellOriented(geometry.GetGeometryN(g)))
                        return false;
                }
                return true;
            }
            else
            {
                return NetTopologySuite.Algorithm.Orientation.IsCCW(geometry.Coordinates);
            }
        }

        /// <summary>
        /// Returns true if the given <see cref="SharpGeometry"/> is well oriented (CCW)
        /// </summary>
        /// <param name="geom"></param>
        /// <returns></returns>
        public static bool IsWellOriented(SharpGeometry geom)
        {
            return IsWellOriented(geom.Geo);
        }

        /// <summary>
        /// Returns true if the given <see cref="SharpGeography"/> is well oriented (CCW).
        /// </summary>
        /// <param name="geog"></param>
        /// <returns></returns>
        public static bool IsWellOriented(SharpGeography geog)
        {
            return IsWellOriented(geog.Geo);
        }

        /// <summary>
        /// Computes a new <see cref="SharpGeography"/> which has all component coordinate sequences
        /// in reverse order (opposite orientation) to this one.
        /// </summary>
        /// <param name="geog"></param>
        /// <returns></returns>
        public static SharpGeography Revert(SharpGeography geog)
        {
            SharpGeography result = CreateGeography(geog.Geo.Reverse(), geog.SRID);
            result.ToleranceUsedForDensify = geog.ToleranceUsedForDensify;
            return result;
        }

        /// <summary>
        /// Computes a new <see cref="SharpGeometry"/> which has all component coordinate sequences
        /// in reverse order (opposite orientation) to this one
        /// </summary>
        /// <param name="geom"></param>
        /// <returns></returns>
        public static SharpGeometry Revert(SharpGeometry geom)
        {
            return CreateGeometry(geom.Geo.Reverse());
        }

        /// <summary>
        /// Returns the fixed <see cref="SharpGeography"/> from the given one
        /// </summary>
        /// <param name="geog"></param>
        /// <returns></returns>
        public static SharpGeography? MakeValid(SharpGeography geog)
        {
            Geometry? validGeo = MakeValid(geog.Geo);
            if (validGeo == null) return null;

            return new SharpGeography(validGeo, geog.SRID, geog.ToleranceUsedForDensify);
        }

        /// <summary>
        /// Returns the fixed <see cref="SharpGeometry"/> from the given one
        /// </summary>
        /// <param name="geom"></param>
        /// <returns></returns>
        public static SharpGeometry? MakeValid(SharpGeometry geom)
        {
            Geometry? validGeo = MakeValid(geom.Geo);
            if (validGeo == null) return null;

            return new SharpGeometry(validGeo);
        }

        /// <summary>
        /// Returns the Mercator projected point from a given latitude and longitude WGS84 coordinates.
        /// </summary>
        /// <param name="lat"></param>
        /// <param name="lon"></param>
        /// <returns></returns>
        public static Coordinate Project(double lat, double lon)
        {
            GeoAPI.Geometries.Coordinate pointCoordinate = new(lon, lat);
            GeoAPI.Geometries.Coordinate result = s_ctProj.MathTransform.Transform(pointCoordinate);
            return new Coordinate(result.X, result.Y);
        }

        /// <summary>
        /// Returns the unprojected point from a given x and y Mercator coordinates.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static Coordinate Unproject(double x, double y)
        {
            GeoAPI.Geometries.Coordinate pointCoordinate = new(x, y);
            GeoAPI.Geometries.Coordinate result = s_ctUnproj.MathTransform.Transform(pointCoordinate);
            return new Coordinate(result.X, result.Y);
        }

        /// <summary>
        /// Returns a Mercator projected <see cref="SharpGeometry"/> from a <see cref="SharpGeography"/>
        /// </summary>
        /// <param name="geography"></param>
        /// <returns></returns>
        public static SharpGeometry? Project(SharpGeography? geography)
        {
            Geometry? projected = Project(geography?.Geo);
            return projected != null ? new SharpGeometry(projected) : null;
        }

        /// <summary>
        /// Returns an unprojected <see cref="SharpGeography"/>, with the given SRID, from a Mercator <see cref="SharpGeometry"/>
        /// </summary>
        /// <param name="geometry"></param>
        /// <param name="srid"></param>
        /// <returns></returns>
        public static SharpGeography? Unproject(SharpGeometry? geometry, int srid)
        {
            Geometry? unprojected = Unproject(geometry?.Geo);
            return unprojected != null ? new SharpGeography(unprojected, srid) : null;
        }
        #endregion

        #region Private methods
        private static Geometry FixAntimeridianCrossing(Geometry geometry)
        {
            if (geometry is Polygon)
            {
                return FixAntimeridianForPolygon((Polygon)geometry);
            }
            else if (geometry is MultiPolygon multiPolygon)
            {
                var fixedPolygons = new List<Geometry>();

                for (int i = 0; i < multiPolygon.NumGeometries; i++)
                {
                    var poly = multiPolygon.GetGeometryN(i) as Polygon;
                    if (poly != null)
                    {
                        fixedPolygons.Add(FixAntimeridianForPolygon(poly));
                    }
                }

                return s_geometryFactory.CreateGeometryCollection(fixedPolygons.ToArray()).Union();
            }

            return geometry;
        }

        private static Geometry FixAntimeridianForPolygon(Polygon polygon)
        {
            Coordinate[] coords = polygon.Coordinates.Select(c => new Coordinate(c)).ToArray();

            bool? splitPositive = null;

            for (int i = 1; i < coords.Length; i++)
            {
                double delta = coords[i].X - coords[i - 1].X;
                if (Math.Abs(delta) > 180)
                {
                    splitPositive = coords[i].X < 0;
                    if (splitPositive == true)
                        coords[i].X += 360;
                    else
                        coords[i].X -= 360;
                }
            }

            if (!splitPositive.HasValue)
                return polygon;

            // Crea maschere
            Geometry maskW = s_geometryFactory.CreatePolygon([
                new Coordinate(0, 89), new Coordinate(0, -89),
                new Coordinate(180, -89), new Coordinate(180, 89),
                new Coordinate(0, 89)]
            );

            Geometry maskE = s_geometryFactory.CreatePolygon([
                new Coordinate(180, 89), new Coordinate(180, -89),
                new Coordinate(360, -89), new Coordinate(360, 89),
                new Coordinate(180, 89)]
            );

            // Copy and apply traslation
            var adjusted = (Polygon)polygon.Copy();

            adjusted.Apply(new CoordinateSequenceFilterAdapter(seq =>
            {
                for (int i = 0; i < seq.Count; i++)
                {
                    double x = seq.GetX(i);
                    if (splitPositive == true && x < 0)
                        seq.SetX(i, x + 360);
                    else if (splitPositive == false && x > 0)
                        seq.SetX(i, x - 360);
                }
            }));

            // Intersect with masks
            Geometry west = adjusted.Intersection(maskW);
            Geometry east = adjusted.Intersection(maskE);

            // Set long to original
            if (splitPositive == true)
            {
                east.Apply(new CoordinateSequenceFilterAdapter(seq =>
                {
                    for (int i = 0; i < seq.Count; i++)
                        seq.SetX(i, seq.GetX(i) - 360);
                }));
            }
            else
            {
                west.Apply(new CoordinateSequenceFilterAdapter(seq =>
                {
                    for (int i = 0; i < seq.Count; i++)
                        seq.SetX(i, seq.GetX(i) + 360);
                }));
            }

            return east.Union(west);
        }

        private static Geometry FixParsedGeometry(Geometry geometry, bool? makeValid, double toleranceForDensify = 0)
        {
            Geometry result = geometry;

            if (result is Polygon or MultiPolygon)
            {
                result = FixAntimeridianCrossing(result);
            }

            if (makeValid != false)
            {
                if (!result.IsValid)
                {
                    if (makeValid == null)
                        throw new ArgumentException("Geometry is invalid and 'makeValid' is not explicitly allowed.");

                    result = MakeValid(result)
                        ?? throw new InvalidOperationException("Geometry is invalid and cannot be repaired.");
                }

                if (!IsWellOriented(result))
                {
                    if (makeValid == null)
                        throw new ArgumentException("Geometry orientation is incorrect and 'makeValid' is not explicitly allowed.");

                    result = result.Reverse();
                }
            }

            if (toleranceForDensify > 0)
            {
                result = Densify(result, toleranceForDensify);
            }

            return result;
        }

        private static void ThrowOnInvalidGeogArgs(params SharpGeography[] geographies)
        {
            if (geographies.Length > 0)
            {
                int srid = geographies[0].SRID;
                if (geographies.Any(g => g.SRID != srid)) throw new ArgumentException("Geographies must have the same SRID");
            }
        }

        private static DistanceSolution GetPlanarDistance(double x1, double y1, double x2, double y2)
        {
            double deltaX = x2 - x1; // East component
            double deltaY = y2 - y1; // North component

            Coordinate point1 = new(x1, y1);
            Coordinate point2 = new(x2, y2);
            double distance = point1.Distance(point2);

            // Bearing in radians, measured clockwise from north
            double bearingRad = Math.Atan2(deltaX, deltaY); // Note: Atan2(x, y), not (y, x)
            double bearingDeg = (ToDegrees(bearingRad) + 360) % 360;

            return new DistanceSolution(null, distance, bearingDeg, bearingDeg);
        }

        private static DistanceSolution GetHaversineDistance(double lat1Deg, double lon1Deg, double lat2Deg, double lon2Deg)
        {
            // Convert degrees to radians
            double lat1Rad = ToRadians(lat1Deg);
            double lon1Rad = ToRadians(lon1Deg);
            double lat2Rad = ToRadians(lat2Deg);
            double lon2Rad = ToRadians(lon2Deg);

            // Differences
            double deltaLat = lat2Rad - lat1Rad;
            double deltaLon = lon2Rad - lon1Rad;

            // Haversine formula
            double a = Math.Pow(Math.Sin(deltaLat / 2), 2) +
                       Math.Cos(lat1Rad) * Math.Cos(lat2Rad) * Math.Pow(Math.Sin(deltaLon / 2), 2);
            double centralAngle = 2 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1 - a));
            double distanceMeters = HAVERSINE_EARTH_R_M * centralAngle;

            // Initial bearing calculation
            double y = Math.Sin(deltaLon) * Math.Cos(lat2Rad);
            double x = Math.Cos(lat1Rad) * Math.Sin(lat2Rad) -
                       Math.Sin(lat1Rad) * Math.Cos(lat2Rad) * Math.Cos(deltaLon);
            double initialBearingRad = Math.Atan2(y, x);
            double initialBearingDeg = (ToDegrees(initialBearingRad) + 360) % 360;

            // Final bearing calculation (from point2 to point1)
            y = Math.Sin(-deltaLon) * Math.Cos(lat1Rad);
            x = Math.Cos(lat2Rad) * Math.Sin(lat1Rad) -
                Math.Sin(lat2Rad) * Math.Cos(lat1Rad) * Math.Cos(-deltaLon);
            double finalBearingRad = Math.Atan2(y, x);
            double finalBearingDeg = (ToDegrees(finalBearingRad) + 360) % 360;

            return new DistanceSolution(GeodesicPrecision.Fast, distanceMeters, initialBearingDeg, finalBearingDeg);
        }

        private static DistanceSolution GetVincentyDistance(double latitude1, double longitude1, double latitude2, double longitude2)
        {
            // Convert from degrees to radians
            double lat1Rad = ToRadians(latitude1);
            double lon1Rad = ToRadians(longitude1);
            double lat2Rad = ToRadians(latitude2);
            double lon2Rad = ToRadians(longitude2);

            // Vincenty formula constants
            double flattening = VINCENTY_ELLIPSOID_WGS84_F;
            double a = VINCENTY_ELLIPSOID_WGS84_A;
            double b = VINCENTY_ELLIPSOID_WGS84_B;

            double tanU1 = (1 - flattening) * Math.Tan(lat1Rad);
            double tanU2 = (1 - flattening) * Math.Tan(lat2Rad);
            double U1 = Math.Atan(tanU1);
            double U2 = Math.Atan(tanU2);

            double sinU1 = Math.Sin(U1), cosU1 = Math.Cos(U1);
            double sinU2 = Math.Sin(U2), cosU2 = Math.Cos(U2);

            double lambda = lon2Rad - lon1Rad;
            double prevLambda;
            int maxIterations = 20;
            double sinLambda, cosLambda;
            double sigma, sinSigma, cosSigma;
            double sinAlpha, cosSqAlpha, cos2SigmaM;
            double C;

            do
            {
                sinLambda = Math.Sin(lambda);
                cosLambda = Math.Cos(lambda);

                double tmp1 = cosU2 * sinLambda;
                double tmp2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
                sinSigma = Math.Sqrt(tmp1 * tmp1 + tmp2 * tmp2);

                if (sinSigma == 0)
                    return DistanceSolution.Zero; // coincident points

                cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
                sigma = Math.Atan2(sinSigma, cosSigma);

                sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
                cosSqAlpha = 1 - sinAlpha * sinAlpha;

                cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / (cosSqAlpha != 0 ? cosSqAlpha : 1); // equatorial line: cosSqAlpha=0, avoid div by zero
                if (double.IsNaN(cos2SigmaM)) cos2SigmaM = 0;

                C = flattening / 16 * cosSqAlpha * (4 + flattening * (4 - 3 * cosSqAlpha));
                prevLambda = lambda;
                lambda = lon2Rad - lon1Rad + (1 - C) * flattening * sinAlpha *
                         (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));

            } while (Math.Abs(lambda - prevLambda) > 1e-12 && --maxIterations > 0);

            if (maxIterations == 0)
                return DistanceSolution.NaN; // formula failed to converge

            double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
            double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
            double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

            double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 *
                (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
                 B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) *
                 (-3 + 4 * cos2SigmaM * cos2SigmaM)));

            double distance = b * A * (sigma - deltaSigma);

            double initialBearingRad = Math.Atan2(cosU2 * sinLambda,
                                                  cosU1 * sinU2 - sinU1 * cosU2 * cosLambda);
            double finalBearingRad = Math.Atan2(cosU1 * sinLambda,
                                                -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);

            return new DistanceSolution(GeodesicPrecision.Precise, distance, ToDegrees(initialBearingRad), ToDegrees(finalBearingRad));
        }

        private static DestinationSolution GetPlanarDestination(Coordinate startPoint, double initialBearingDeg, double distance)
        {
            // Normalize bearing between 0° and 360°
            double normalizedBearing = (initialBearingDeg % 360 + 360) % 360;

            // Convert bearing to radians (measured clockwise from north, like compass bearing)
            double bearingRad = ToRadians(normalizedBearing);

            // In 2D Cartesian: X = East, Y = North
            double deltaX = distance * Math.Sin(bearingRad);
            double deltaY = distance * Math.Cos(bearingRad);

            double finalX = startPoint.X + deltaX;
            double finalY = startPoint.Y + deltaY;

            // On a flat plane, the final bearing equals the initial one
            return new DestinationSolution(null, new Coordinate(finalX, finalY), normalizedBearing);
        }

        private static DestinationSolution GetHaversineDestination(Coordinate startPoint, double initialBearingDegrees, double distanceMeters)
        {
            const double EarthRadius = 6371000; // Earth radius in meters

            // Convert inputs to radians
            double lat1 = ToRadians(startPoint.Y);
            double lon1 = ToRadians(startPoint.X);
            double bearingRad = ToRadians(initialBearingDegrees);

            double angularDistance = distanceMeters / EarthRadius;

            // Destination latitude
            double lat2 = Math.Asin(
                Math.Sin(lat1) * Math.Cos(angularDistance) +
                Math.Cos(lat1) * Math.Sin(angularDistance) * Math.Cos(bearingRad)
            );

            // Destination longitude
            double lon2 = lon1 + Math.Atan2(
                Math.Sin(bearingRad) * Math.Sin(angularDistance) * Math.Cos(lat1),
                Math.Cos(angularDistance) - Math.Sin(lat1) * Math.Sin(lat2)
            );

            // Normalize longitude to [-180, 180]
            lon2 = (lon2 + 3 * Math.PI) % (2 * Math.PI) - Math.PI;

            // Final bearing calculation (from destination point back toward origin)
            double finalBearingRad = Math.Atan2(
                Math.Sin(-bearingRad) * Math.Sin(angularDistance) * Math.Cos(lat2),
                Math.Cos(angularDistance) - Math.Sin(lat2) * Math.Sin(lat1)
            );

            // Normalize to [0, 360]
            double finalBearingDeg = (ToDegrees(finalBearingRad) + 360) % 360;

            return new DestinationSolution(GeodesicPrecision.Fast, new Coordinate(ToDegrees(lon2), ToDegrees(lat2)), finalBearingDeg);
        }

        private static DestinationSolution? GetVincentyDestination(Coordinate startPoint, double initialBearingDegrees, double distanceMeters)
        {
            Coordinate destination = new(startPoint);
            double finalBearingDegrees = initialBearingDegrees;

            while (true)
            {
                // Normalize bearing to [0, 360]
                initialBearingDegrees %= 360;
                if (initialBearingDegrees < 0)
                    initialBearingDegrees += 360;

                // Convert inputs to radians
                double lat1Rad = ToRadians(startPoint.Y);
                double lon1Rad = ToRadians(startPoint.X);
                double alpha1Rad = ToRadians(initialBearingDegrees);

                // Ellipsoid constants
                double f = VINCENTY_ELLIPSOID_WGS84_F;
                double a = VINCENTY_ELLIPSOID_WGS84_A;
                double b = VINCENTY_ELLIPSOID_WGS84_B;

                // Pre-compute trigonometric values
                double sinAlpha1 = Math.Sin(alpha1Rad);
                double cosAlpha1 = Math.Cos(alpha1Rad);

                double tanU1 = (1 - f) * Math.Tan(lat1Rad);
                double cosU1 = 1 / Math.Sqrt(1 + tanU1 * tanU1);
                double sinU1 = tanU1 * cosU1;

                // Initial values
                double sigma1 = Math.Atan2(tanU1, cosAlpha1);
                double sinAlpha = cosU1 * sinAlpha1;
                double cosSqAlpha = 1 - sinAlpha * sinAlpha;

                double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
                double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
                double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

                double sigma = distanceMeters / (b * A);
                double sigmaPrev, deltaSigma;
                double cos2SigmaM, sinSigma, cosSigma;

                int iterations = 0;
                do
                {
                    cos2SigmaM = Math.Cos(2 * sigma1 + sigma);
                    sinSigma = Math.Sin(sigma);
                    cosSigma = Math.Cos(sigma);

                    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (
                        cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
                        B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) *
                        (-3 + 4 * cos2SigmaM * cos2SigmaM)));

                    sigmaPrev = sigma;
                    sigma = distanceMeters / (b * A) + deltaSigma;

                } while (Math.Abs(sigma - sigmaPrev) > 1e-12 && ++iterations < 100);

                if (iterations >= 100)
                    return null; // convergence failed

                double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;
                double lat2Rad = Math.Atan2(
                    sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
                    (1 - f) * Math.Sqrt(sinAlpha * sinAlpha + tmp * tmp)
                );

                double lambda = Math.Atan2(
                    sinSigma * sinAlpha1,
                    cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1
                );

                double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
                double L = lambda - (1 - C) * f * sinAlpha * (
                    sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM))
                );

                double lon2Rad = lon1Rad + L;

                destination = new Coordinate(ToDegrees(lon2Rad), ToDegrees(lat2Rad));

                // Calculate final bearing (reverse azimuth)
                double finalAlphaRad = Math.Atan2(
                    cosU1 * sinAlpha1,
                    -sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1
                );
                finalBearingDegrees = (ToDegrees(finalAlphaRad) + 360) % 360;

                // Check for large shift due to polar crossing
                if (Math.Abs(startPoint.X - destination.X) > 45)
                {
                    DistanceSolution vincentyCheck = GetVincentyDistance(
                        startPoint.Y, startPoint.X,
                        destination.Y, destination.X
                    );

                    double bearingDiff = Math.Abs(vincentyCheck.DestinationBearing - initialBearingDegrees);
                    if (bearingDiff > 90)
                    {
                        distanceMeters -= 1000;
                        continue;
                    }
                }

                break;
            }

            return new DestinationSolution(GeodesicPrecision.Precise, destination, finalBearingDegrees);
        }

        private static double ToRadians(double value) => Math.PI * value / 180.0;

        private static double ToDegrees(double value) => value * 180 / Math.PI;

        private static Coordinate[] GetExternalShell(Geometry geometry, int deepLevel = 0)
        {
            if (deepLevel > 10) throw new InvalidOperationException("Max deepLevel of 10 reached, possible future Stack Overflow exception");

            if (geometry == null) return [];

            if (geometry is Polygon polygon)
                return polygon.Shell.Coordinates;

            if (geometry is GeometryCollection geomCollection)
            {
                List<Coordinate> coordinates = [];
                foreach (Geometry geom in geomCollection.Geometries)
                {
                    coordinates.AddRange(GetExternalShell(geom, deepLevel + 1));
                    return [.. coordinates];
                }
            }

            if (geometry is MultiPolygon multiPolygon)
            {
                List<Coordinate> coordinates = [];
                foreach (Geometry geom in multiPolygon.Geometries)
                {
                    coordinates.AddRange(GetExternalShell(geom, deepLevel + 1));
                    return [.. coordinates];
                }
            }

            return geometry.Coordinates;
        }

        private static Coordinate GetNearestPointOnSegment(Coordinate point, Coordinate segStart, Coordinate segEnd)
        {
            // Convert to 2D Cartesian (approx): lat ≈ Y, lon ≈ X
            double px = point.X, py = point.Y;
            double x1 = segStart.X, y1 = segStart.Y;
            double x2 = segEnd.X, y2 = segEnd.Y;

            double dx = x2 - x1;
            double dy = y2 - y1;
            if (dx == 0 && dy == 0)
                return new Coordinate(x1, y1); // segment is a point

            double t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);
            t = Math.Max(0, Math.Min(1, t));
            return new Coordinate(x1 + t * dx, y1 + t * dy);
        }

        private static Geometry? BufferSegment(Coordinate p1, Coordinate p2, double distance, double degreesStep,
            double densifyWithBearingTolerance)
        {
            DistanceSolution solution = GetVincentyDistance(p1.Y, p1.X, p2.Y, p2.X);

            if (solution.Distance <= 0)
                return null;

            // Curve around P1
            double startDegP1 = solution.InitialBearing - 90;
            double endDegP1 = startDegP1 - 180;
            List<Coordinate> curveP1 = GetBufferCurvePoints(p1, startDegP1, endDegP1, degreesStep, distance, out double p1MinX, out double p1MaxX);

            // Curve around P2
            double startDegP2 = solution.DestinationBearing + 90;
            double endDegP2 = startDegP2 - 180;
            List<Coordinate> curveP2 = GetBufferCurvePoints(p2, startDegP2, endDegP2, degreesStep, distance, out double p2MinX, out double p2MaxX);

            // Combine coordinates
            curveP1.AddRange(curveP2);
            curveP1.Add(curveP1[0]); // close ring

            double minX = Math.Min(p1MinX, p2MinX);
            double maxX = Math.Max(p1MaxX, p2MaxX);

            Geometry result = Math.Abs(maxX - minX) < 170 && curveP1.Count > 3
                ? s_geometryFactory.CreatePolygon(curveP1.ToArray())
                : s_geometryFactory.CreateLineString(new[] { p1, p2 });

            return Densify(result, densifyWithBearingTolerance);
        }

        private static List<Coordinate> GetBufferCurvePoints(Coordinate center, double fromDeg, double toDeg, double stepDeg,
            double distance, out double minX, out double maxX)
        {
            minX = double.MaxValue;
            maxX = double.MinValue;

            List<Coordinate> result = [];
            double maxXDistance = distance / 111_111.1111 * 3; // max approx 3 km longitudinal deviation

            for (double deg = fromDeg; deg >= toDeg; deg -= stepDeg)
            {
                DestinationSolution? point = GetVincentyDestination(center, deg, distance);
                if (point == null) continue;

                double xDiff = Math.Abs(point.Coordinate.X - center.X);
                if (xDiff < maxXDistance)
                {
                    result.Add(point.Coordinate);
                    minX = Math.Min(minX, point.Coordinate.X);
                    maxX = Math.Max(maxX, point.Coordinate.X);
                }
            }

            return result;
        }

        private static Geometry CalculateBuffer(Geometry geometry, double distance, double densifyWithBearingTolerance)
        {
            if (distance == 0 || geometry.Coordinates.Length == 0)
                return geometry;

            if (geometry.Coordinates.Length == 1)
            {
                return BufferSinglePoint(geometry.Coordinates[0], distance, 1);
            }

            List<Geometry> buffers = [geometry];
            for (int i = 0; i < geometry.Coordinates.Length - 1; i++)
            {
                try
                {
                    Geometry? segmentBuffer = BufferSegment(geometry.Coordinates[i], geometry.Coordinates[i + 1], distance, 1, densifyWithBearingTolerance);
                    if (segmentBuffer != null)
                        buffers.Add(segmentBuffer);
                }
                catch (Exception ex)
                {
                    // Log error?
                    Debug.WriteLine($"BufferSegment error: {ex.Message}");
                }
            }

            return NetTopologySuite.Operation.Union.CascadedPolygonUnion.Union(buffers);
        }

        private static Geometry BufferSinglePoint(Coordinate point, double distance, double stepDeg)
        {
            var coords = new List<Coordinate>();
            for (double deg = 0; deg <= 360; deg += stepDeg)
            {
                DestinationSolution? dest = GetVincentyDestination(point, deg, distance);
                if (dest != null)
                    coords.Add(dest.Coordinate);
            }

            coords.Add(coords[0]); // close ring
            return s_geometryFactory.CreatePolygon(coords.ToArray());
        }

        private static Geometry Densify(Geometry geometry, double bearingTolerance)
        {
            if (geometry == null || geometry.IsEmpty)
                return CreateEmptyGeo();

            // Manage multi-part geometries
            if (geometry.NumGeometries > 1)
            {
                Geometry[] densifiedParts = new Geometry[geometry.NumGeometries];
                for (int i = 0; i < geometry.NumGeometries; i++)
                {
                    densifiedParts[i] = Densify(geometry.GetGeometryN(i), bearingTolerance);
                }
                return s_geometryFactory.BuildGeometry(densifiedParts);
            }

            return geometry.OgcGeometryType switch
            {
                OgcGeometryType.LineString => DensifyLineString((LineString)geometry, bearingTolerance),
                OgcGeometryType.Polygon => DensifyPolygon((Polygon)geometry, bearingTolerance),
                _ => geometry
            };
        }

        private static LineString DensifyLineString(LineString line, double tolerance)
        {
            List<Coordinate> coords = DensifyCoordinates(line.Coordinates, tolerance);
            return s_geometryFactory.CreateLineString(coords.ToArray());
        }

        private static Polygon DensifyPolygon(Polygon polygon, double tolerance)
        {
            List<Coordinate> shellCoords = DensifyCoordinates(polygon.ExteriorRing.Coordinates, tolerance, closeRing: true);
            LinearRing shell = s_geometryFactory.CreateLinearRing(shellCoords.ToArray());

            var holes = new LinearRing[polygon.NumInteriorRings];
            for (int i = 0; i < polygon.NumInteriorRings; i++)
            {
                List<Coordinate> holeCoords = DensifyCoordinates(polygon.GetInteriorRingN(i).Coordinates, tolerance, closeRing: true);
                holes[i] = s_geometryFactory.CreateLinearRing(holeCoords.ToArray());
            }

            return s_geometryFactory.CreatePolygon(shell, holes);
        }

        private static List<Coordinate> DensifyCoordinates(Coordinate[] coords, double tolerance, bool closeRing = false)
        {
            var result = new List<Coordinate>();
            int len = coords.Length;

            for (int i = 0; i < len - 1; i++)
            {
                Coordinate start = coords[i];
                Coordinate end = coords[i + 1];
                result.Add(start);

                AddIntermediatePoints(start, end, tolerance, result);
            }

            result.Add(coords[len - 1]);

            if (closeRing && !result.First().Equals2D(result.Last()))
                result.Add(result.First());

            return result;
        }

        private static void AddIntermediatePoints(Coordinate start, Coordinate end, double tolerance, List<Coordinate> result)
        {
            DistanceSolution vincenty = GetVincentyDistance(start.Y, start.X, end.Y, end.X);
            double bearingDiff = NormalizeBearingDifference(vincenty.InitialBearing, vincenty.DestinationBearing);

            if (bearingDiff > tolerance && vincenty.Distance > 1) // avoid points too close
            {
                double halfDistance = vincenty.Distance / 2.0;
                Coordinate? mid = GetVincentyDestination(start, vincenty.InitialBearing, halfDistance)?.Coordinate;
                if (mid != null)
                {
                    AddIntermediatePoints(start, mid, tolerance, result);
                    result.Add(mid);
                    AddIntermediatePoints(mid, end, tolerance, result);
                }
            }
        }

        private static double NormalizeBearingDifference(double bearing1, double bearing2)
        {
            double diff = Math.Abs(bearing1 - bearing2) % 360;
            return diff > 180 ? 360 - diff : diff;
        }

        private static SharpGeography Densify(SharpGeography geography, double bearingTolerance)
        {
            if (bearingTolerance <= 0) return geography;
            if (geography.ToleranceUsedForDensify > 0 && geography.ToleranceUsedForDensify < bearingTolerance)
                return geography;

            Geometry densifiedGeo = Densify(geography.Geo, bearingTolerance);

            return new SharpGeography(densifiedGeo, geography.SRID, bearingTolerance);
        }

        private static Geometry? Project(Geometry? geometry)
        {
            if (geometry == null)
                return null;

            if (geometry.IsEmpty)
                return s_geometryFactory.CreateGeometryCollection();

            switch (geometry.OgcGeometryType)
            {
                case OgcGeometryType.Point:
                    return s_geometryFactory.CreatePoint(ProjectCoordinate(geometry.Coordinate));

                case OgcGeometryType.MultiPoint:
                case OgcGeometryType.MultiLineString:
                case OgcGeometryType.MultiPolygon:
                case OgcGeometryType.GeometryCollection:
                    return s_geometryFactory.BuildGeometry(
                        [.. Enumerable.Range(0, geometry.NumGeometries)
                                  .Select(i => Project(geometry.GetGeometryN(i)))
                                  .Where(g => g != null)]
                    );

                case OgcGeometryType.LineString:
                    return s_geometryFactory.CreateLineString(
                        geometry.Coordinates.Select(ProjectCoordinate).ToArray()
                    );

                case OgcGeometryType.Polygon:
                    Polygon polygon = (Polygon)geometry;
                    LinearRing shell = s_geometryFactory.CreateLinearRing(
                        polygon.Shell.Coordinates.Select(ProjectCoordinate).ToArray()
                    );

                    LinearRing[] holes = polygon.Holes.Select(h =>
                        s_geometryFactory.CreateLinearRing(
                            h.Coordinates.Select(ProjectCoordinate).ToArray()
                        )).ToArray();

                    return s_geometryFactory.CreatePolygon(shell, holes);

                default:
                    throw new NotSupportedException($"Projection not supported for geometry type: {geometry.OgcGeometryType}");
            }
        }

        private static Coordinate ProjectCoordinate(Coordinate c)
        {
            return Project(c.Y, c.X); // Assuming this.Project(lat, lon)
        }

        private static Geometry? Unproject(Geometry? geometry)
        {
            if (geometry == null || geometry.IsEmpty)
                return s_geometryFactory.CreateGeometryCollection();

            switch (geometry.OgcGeometryType)
            {
                case OgcGeometryType.Point:
                    return s_geometryFactory.CreatePoint(UnprojectCoordinate(geometry.Coordinate));

                case OgcGeometryType.MultiPoint:
                case OgcGeometryType.MultiLineString:
                case OgcGeometryType.MultiPolygon:
                case OgcGeometryType.GeometryCollection:
                    return s_geometryFactory.BuildGeometry(
                        [.. Enumerable.Range(0, geometry.NumGeometries).Select(i => Unproject(geometry.GetGeometryN(i))).Where(g => g != null)]
                    );

                case OgcGeometryType.LineString:
                    return s_geometryFactory.CreateLineString(
                        geometry.Coordinates.Select(UnprojectCoordinate).ToArray()
                    );

                case OgcGeometryType.Polygon:
                    Polygon polygon = (Polygon)geometry;

                    LinearRing shell = s_geometryFactory.CreateLinearRing(
                        polygon.Shell.Coordinates.Select(UnprojectCoordinate).ToArray()
                    );

                    LinearRing[] holes = polygon.Holes.Select(h =>
                        s_geometryFactory.CreateLinearRing(
                            h.Coordinates.Select(UnprojectCoordinate).ToArray()
                        )).ToArray();

                    return s_geometryFactory.CreatePolygon(shell, holes);

                default:
                    throw new NotSupportedException($"Unprojection not supported for geometry type: {geometry.OgcGeometryType}");
            }
        }

        private static Coordinate UnprojectCoordinate(Coordinate c)
        {
            return Unproject(c.X, c.Y); // Assuming Unproject(x, y) returns a Coordinate in lat/lon
        }

        private static bool IsCrossingIDL(double lon1, double lon2)
        {
            double delta = Math.Abs(NormalizeLongitude(lon1) - NormalizeLongitude(lon2));
            return delta > 180;
        }

        private static double NormalizeLongitude(double lon)
        {
            lon %= 360;
            if (lon > 180) lon -= 360;
            if (lon < -180) lon += 360;
            return lon;
        }

        private static Geometry MakeValid(Geometry geometry)
        {
            return geometry.IsValid ? geometry : GeometryFixer.Fix(geometry);
        }

        private static Geometry Simplify(Geometry geometry, double tolerance)
        {
            if (geometry.IsEmpty)
                return geometry;

            switch (geometry.OgcGeometryType)
            {
                case OgcGeometryType.MultiPolygon:
                case OgcGeometryType.GeometryCollection:
                    {
                        var parts = new List<Geometry>();
                        for (int i = 0; i < geometry.NumGeometries; i++)
                        {
                            Geometry simplified = Simplify(geometry.GetGeometryN(i), tolerance);
                            if (simplified != null && !simplified.IsEmpty)
                                parts.Add(simplified);
                        }

                        return s_geometryFactory.BuildGeometry(parts).Union(); // single union at end
                    }

                case OgcGeometryType.LineString:
                    {
                        if (geometry.NumPoints <= 3)
                            return geometry;

                        Geometry simplified = DouglasPeuckerSimplifier.Simplify(geometry, tolerance);
                        return simplified.NumPoints >= 2 ? simplified : geometry;
                    }

                case OgcGeometryType.Polygon:
                    {
                        Polygon polygon = (Polygon)geometry;

                        var simplifiedShell = DouglasPeuckerSimplifier.Simplify(polygon.ExteriorRing, tolerance) as LineString;
                        if (simplifiedShell == null || simplifiedShell.NumPoints < 4)
                            return geometry;

                        LinearRing shell = s_geometryFactory.CreateLinearRing(simplifiedShell.CoordinateSequence);
                        var holes = new List<LinearRing>();

                        foreach (LineString ring in polygon.InteriorRings)
                        {
                            var simplifiedHole = DouglasPeuckerSimplifier.Simplify(ring, tolerance) as LineString;

                            if (simplifiedHole != null && simplifiedHole.NumPoints >= 4)
                                holes.Add(s_geometryFactory.CreateLinearRing(simplifiedHole.CoordinateSequence));
                        }

                        Geometry simplifiedPolygon = s_geometryFactory.CreatePolygon(shell, holes.ToArray());
                        return simplifiedPolygon.IsValid ? simplifiedPolygon : geometry;
                    }

                default:
                    return geometry; // fallback for unsupported types
            }
        }
        #endregion
    }
}
