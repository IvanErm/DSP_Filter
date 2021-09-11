using System.Collections.Generic;
using static System.Math;

namespace DSP_Model.Windows
{
    public class BlackmanWindow : WindowWeightAbstract
    {
        public BlackmanWindow(List<double> w, int n, WindowTypes windowType, SymmetryType symmeticType)
            : base(w, n, windowType, symmeticType)
        {
        }

        public override List<double> Weight()
        {
            double x = 0.0;

            switch (symmeticType)
            {
                case SymmetryType.Symmetric:
                    x = (2 * PI) / (double)(n - 1);
                    break;

                case SymmetryType.Periodic:
                    x = (2 * PI) / (double)n;
                    break;
            }

            double y = 0.0;
            for (int i = 0; i < n; i++)
            {
                w[i] = 0.42 - 0.5 * Cos(y) + 0.08 * Cos(2.0 * y);
                y += x;
            }

            return w;
        }
    }
}
