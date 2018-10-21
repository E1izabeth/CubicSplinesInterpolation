using MathLibrary;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace TestMySpline
{
	class Program
	{
		static void Main(string[] args)
		{
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            var vis = new CubicSplinesInterpolation();
            Application.Run(vis);
            //TestSplineExample();
		}
        
		/// <summary>
		/// This is the Wikipedia "Spline Interpolation" article example.
		/// </summary>
		private static void TestSplineExample()
		{
            CubicSpline spl = new CubicSpline();
            // Create the test points.
            double[] x = { -5, -4, -3, -2, -1, 0, 1, 2};
            Func<double, double> func = Math.Sin;


            double[] y = new double[x.Length];

            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Math.Sin(x[i]);
            }

            double[] xs, ys;
            spl.FitParametric(x, y, 100, out xs, out ys);

            double[] xO = new double[x.Length * 10];
            double[] yO = new double[x.Length * 10];
            var step = Math.Abs(x[0] - x[x.Length - 1]) / x.Length / 10;
            for (int i = 0; i < x.Length * 10; i++)
            {
                xO[i] = x[0] + step * i;
                yO[i] = Math.Sin(xO[i]);
            }
            

            // Plot
            string path = @"..\..\spline.png";
			//var chart = PlotSplineSolution("Cubic Spline Interpolation", x, y, xs, ys, xO, yO, path);
            
            //if (File.Exists(path))
            //{
            //    File.Delete(path);
            //}

            //using (FileStream fs = new FileStream(path, FileMode.CreateNew))
            //{
            //    chart.SaveImage(fs, ChartImageFormat.Png);
            //}
        }
	}
}
