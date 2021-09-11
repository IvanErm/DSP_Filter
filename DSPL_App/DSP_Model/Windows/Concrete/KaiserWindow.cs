using System;
using System.Collections.Generic;
using static System.Math;
using static DSP_Model.HelperMethods.MathMethods;

namespace DSP_Model.Windows
{
    public class KaiserWindow : WindowWeightAbstract
    {
        public KaiserWindow(List<double> w, int n, WindowTypes windowType, SymmetryType symmeticType, double param)
            : base(w, n, windowType, symmeticType)
        {
            this.param = param;
        }

        public override List<double> Weight()
        {
            double x, y;
            double L = 0;

            switch (symmeticType)
            {
                case SymmetryType.Symmetric:
                    L = (double)(n - 1) / 2.0;
                    break;

                case SymmetryType.Periodic:
                    L = (double)n / 2.0;
                    break;
            }

            double den = Bessel_i0(param, 1);

            for (int i = 0; i < n; i++)
            {
                x = 2.0 * ((double)i - L) / (double)n;
                y = param * Sqrt(1.0 - x * x);
                double num = Bessel_i0(y, 1);

                w[i] = num / den;
            }

            return w;
        }
    }
}
