using System;
using System.Collections.Generic;
using static System.Math;
namespace DSP_Model.Windows
{
    public class HannWindow : WindowWeightAbstract
    {
        public HannWindow(List<double> w, int n, WindowTypes windowType, SymmetryType symmeticType)
            : base(w, n, windowType, symmeticType)
        {
        }

        public override List<double> Weight()
        {
            double x = 0.0;
            double y;

            switch (symmeticType)
            {
                case SymmetryType.Symmetric:
                    x = (2 * PI) / (double)(n - 1);
                    break;

                case SymmetryType.Periodic:
                    x = (2 * PI) / (double)n;
                    break;
            }

            y = 0.0;
            for (int i = 0; i < n; i++)
            {
                w[i] = 0.5 * (1 - Cos(y));
                y += x;
            }

            return w;
        }
    }
}
