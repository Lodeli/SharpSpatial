# SharpSpatial

[![Sponsor](https://img.shields.io/badge/Sponsor-❤-ff69b4)](https://github.com/sponsors/Lodeli)

**SharpSpatial** is a lightweight geospatial library for .NET, built as a modern alternative to `SqlGeography` and `SqlGeometry`, with full support for .NET Core and custom geodesic calculations written entirely in C#.

It combines the power of **NetTopologySuite**, **GeoAPI**, and **ProjNET4GeoAPI** for planar operations, while offering accurate geodesic calculations (Haversine and Vincenty) and basic projection capabilities (Spherical Mercator).

---

## ✨ Why use SharpSpatial?

- ✅ **Pure .NET Core** – no dependency on SQL Server or legacy types.
- ✅ **Custom geodesic engine** – Vincenty & Haversine distance, buffering, and more.
- ✅ **Planar geometry support** – thanks to NetTopologySuite.
- ✅ **Mercator projection** – convert from/to spherical Mercator (EPSG:3857).
- ✅ **Developer-friendly** – clean API, easy integration.

---

## ⚠️ Limitations

While useful in many real-world applications, SharpSpatial is not a full GIS engine.  
Here’s what it doesn’t do (yet):

- ❌ No support for arbitrary coordinate reference systems (CRS)
- ❌ Limited projections (only WGS84 and Mercator)
- ❌ Simplified geometry model (common types only)

Still, for many server-side or data-processing tasks, it offers exactly what’s needed — with minimal overhead.

---

## 📦 Installation

Via NuGet:

```bash
dotnet add package SharpSpatial
```

---

## 🚀 Quick Start

```csharp
// TODO
```

---

## 🧭 Use Cases

- Accurate distance and area calculations over the ellipsoid
- Buffering and manipulating spatial shapes
- Migrating from `SqlGeography` in .NET Core apps
- Lightweight geospatial features for APIs and services

---

## ❤️ Support SharpSpatial

I develop and maintain **SharpSpatial** in my free time, focusing on practical tools I use myself in the real world.

If you find this project useful, please consider [sponsoring me on GitHub](https://github.com/sponsors/Lodeli).  
Even a small donation helps me keep the library maintained, improved, and open for everyone.

---

## 🛣 Roadmap

- [ ] GeoJSON input/output utilities
- [ ] Unit testing improvements

Have an idea or feedback? [Open an issue](https://github.com/Lodeli/SharpSpatial/issues) or start a discussion.

---

## 📄 License

Licensed under the MIT License – free to use, modify, and distribute.

### Third-party licenses

SharpSpatial uses the following third-party libraries:

- [NetTopologySuite](https://github.com/NetTopologySuite/NetTopologySuite), licensed under the BSD-2-Clause license.
- [ProjNET4GeoAPI](https://github.com/NetTopologySuite/ProjNet4GeoAPI), licensed under the BSD-2-Clause license.

These libraries retain their respective licenses.

---

> Accurate geospatial logic for .NET, without the overhead.
