// Licensed to the .NET Foundation under one or more agreements.
// The .NET Foundation licenses this file to you under the MIT license.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NetTopologySuite.Geometries;

namespace SharpSpatial.Helpers
{
    internal class CoordinateSequenceFilterAdapter : ICoordinateSequenceFilter
    {
        private readonly Action<CoordinateSequence> _action;

        public CoordinateSequenceFilterAdapter(Action<CoordinateSequence> action)
        {
            _action = action;
        }

        public void Filter(CoordinateSequence seq, int i) { } // Not needed for full-sequence ops

        public bool Done => false;
        public bool GeometryChanged => true;

        public void Filter(CoordinateSequence seq)
        {
            _action(seq);
        }
    }

}
