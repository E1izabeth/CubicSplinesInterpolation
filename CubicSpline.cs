using System;

namespace CubicSplinesInterpolation
{
    /// <summary>
    /// Cubic spline interpolation.
    /// Call Fit (or use the corrector constructor) to compute spline coefficients, then Eval to evaluate the spline at other X coordinates.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is implemented based on the wikipedia article:
    /// http://en.wikipedia.org/wiki/Spline_interpolation
    /// </para>
    /// </remarks>
    public class CubicSpline
    {

        // N-1 spline coefficients for N points
        private double[] _a;
        private double[] _b;

        // Save the original x and y for Eval
        private double[] _xOrig;
        private double[] _yOrig;

        /// <summary>
        /// Default ctor.
        /// </summary>
        public CubicSpline()
        {
        }

        /// <summary>
        /// Construct and call Fit.
        /// </summary>
        /// <param name="x">Input. X coordinates to fit.</param>
        /// <param name="y">Input. Y coordinates to fit.</param>
        public CubicSpline(double[] x, double[] y)
        {
            this.Fit(x, y);
        }


        /// <summary>
        /// Throws if Fit has not been called.
        /// </summary>
        private void CheckAlreadyFitted()
        {
            if (_a == null) throw new Exception("Fit must be called before you can evaluate.");
        }

        private int _lastIndex = 0;

        /// <summary>
        /// Find where in xOrig the specified x falls, by simultaneous traverse.
        /// This allows xs to be less than x[0] and/or greater than x[n-1]. So allows extrapolation.
        /// This keeps state, so requires that x be sorted and xs called in ascending order, and is not multi-thread safe.
        /// </summary>
        private int GetNextXIndex(double x)
        {
            if (x < _xOrig[_lastIndex])
            {
                throw new ArgumentException("The X values to evaluate must be sorted.");
            }

            while ((_lastIndex < _xOrig.Length - 2) && (x > _xOrig[_lastIndex + 1]))
            {
                _lastIndex++;
            }

            return _lastIndex;
        }

        /// <summary>
        /// Evaluate the specified x value using the specified spline.
        /// </summary>
        /// <param name="x">The x value.</param>
        /// <param name="j">Which spline to use.</param>
        /// <returns>The y value.</returns>
        private double EvalSpline(double x, int j)
        {
            double dx = _xOrig[j + 1] - _xOrig[j];
            double t = (x - _xOrig[j]) / dx;
            double y = (1 - t) * _yOrig[j] + t * _yOrig[j + 1] + t * (1 - t) * (_a[j] * (1 - t) + _b[j] * t); // equation 9
            return y;
        }


        /// <summary>
        /// Fit x,y and then eval at points xs and return the corresponding y's.
        /// This does the "natural spline" style for ends.
        /// This can extrapolate off the ends of the splines.
        /// You must provide points in X sort order.
        /// </summary>
        /// <param name="x">Input. X coordinates to fit.</param>
        /// <param name="y">Input. Y coordinates to fit.</param>
        /// <param name="xs">Input. X coordinates to evaluate the fitted curve at.</param>
        /// <returns>The computed y values for each xs.</returns>

        public double[] FitAndEval(double[] x, double[] y, double[] xs)
        {
            this.Fit(x, y);
            return this.Eval(xs);
        }

        /// <summary>
        /// Compute spline coefficients for the specified x,y points.
        /// This does the "natural spline" style for ends.
        /// This can extrapolate off the ends of the splines.
        /// You must provide points in X sort order.
        /// </summary>
        /// <param name="x">Input. X coordinates to fit.</param>
        /// <param name="y">Input. Y coordinates to fit.</param>
        public void Fit(double[] x, double[] y)
        {
            // Save x and y for eval
            this._xOrig = x;
            this._yOrig = y;

            int n = x.Length;
            double[] r = new double[n]; // the right hand side numbers: wikipedia page overloads b

            TriDiagonalMatrixF m = new TriDiagonalMatrixF(n);
            double dx1, dx2, dy1, dy2;

            // First row is different (equation 16 from the article)
            dx1 = x[1] - x[0];
            m.C[0] = 1.0f / dx1;
            m.B[0] = 2.0f * m.C[0];
            r[0] = 3 * (y[1] - y[0]) / (dx1 * dx1);

            // Body rows (equation 15 from the article)
            for (int i = 1; i < n - 1; i++)
            {
                dx1 = x[i] - x[i - 1];
                dx2 = x[i + 1] - x[i];

                m.A[i] = 1.0f / dx1;
                m.C[i] = 1.0f / dx2;
                m.B[i] = 2.0f * (m.A[i] + m.C[i]);

                dy1 = y[i] - y[i - 1];
                dy2 = y[i + 1] - y[i];
                r[i] = 3 * (dy1 / (dx1 * dx1) + dy2 / (dx2 * dx2));
            }

            // Last row also different (equation 17 from the article)

            dx1 = x[n - 1] - x[n - 2];
            dy1 = y[n - 1] - y[n - 2];
            m.A[n - 1] = 1.0f / dx1;
            m.B[n - 1] = 2.0f * m.A[n - 1];
            r[n - 1] = 3 * (dy1 / (dx1 * dx1));

            // k is the solution to the matrix
            double[] k = m.Solve(r);

            // a and b are each spline's coefficients
            this._a = new double[n - 1];
            this._b = new double[n - 1];

            for (int i = 1; i < n; i++)
            {
                dx1 = x[i] - x[i - 1];
                dy1 = y[i] - y[i - 1];
                _a[i - 1] = k[i - 1] * dx1 - dy1; // equation 10 from the article
                _b[i - 1] = -k[i] * dx1 + dy1; // equation 11 from the article
            }
        }


        /// <summary>
        /// Evaluate the spline at the specified x coordinates.
        /// This can extrapolate off the ends of the splines.
        /// You must provide X's in ascending order.
        /// The spline must already be computed before calling this, meaning you must have already called Fit() or FitAndEval().
        /// </summary>
        /// <param name="x">Input. X coordinates to evaluate the fitted curve at.</param>
        /// <returns>The computed y values for each x.</returns>
        public double[] Eval(double[] x)
        {
            this.CheckAlreadyFitted();

            int n = x.Length;
            double[] y = new double[n];
            _lastIndex = 0; // Reset simultaneous traversal in case there are multiple calls

            for (int i = 0; i < n; i++)
            {
                // Find which spline can be used to compute this x (by simultaneous traverse)
                int j = this.GetNextXIndex(x[i]);

                // Evaluate using j'th spline
                y[i] = this.EvalSpline(x[i], j);
            }

            return y;
        }
        
        /// <summary>
        /// Fit the input x,y points using the parametric approach, so that y does not have to be an explicit
        /// function of x, meaning there does not need to be a single value of y for each x.
        /// </summary>
        /// <param name="x">Input x coordinates.</param>
        /// <param name="y">Input y coordinates.</param>
        /// <param name="nOutputPoints">How many output points to create.</param>
        /// <param name="xs">Output (interpolated) x values.</param>
        /// <param name="ys">Output (interpolated) y values.</param>
        public void FitParametric(double[] x, double[] y, int nOutputPoints, out double[] xs, out double[] ys)
        {
            xs = new double[nOutputPoints];
            var step = Math.Abs(x[0] - x[x.Length - 1]) / nOutputPoints;
            for (int i = 0; i < nOutputPoints; i++)
            {
                xs[i] = x[0] + i * step;
            }
            ys = this.FitAndEval(x, y, xs);
        }
    }
}
