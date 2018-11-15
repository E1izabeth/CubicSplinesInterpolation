using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace CubicSplinesInterpolation
{
	class Program
	{
        static void Main(string[] args)
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            var vis = new CubicSplinesInterpolation();
            Application.Run(vis);
        }
        
	}
}
