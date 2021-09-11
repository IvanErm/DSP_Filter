using System.Collections.Generic;
using static System.Math;

namespace DSP_Model.Windows
{
    public class HammingWindow : WindowWeightAbstract
    {
        public HammingWindow(List<double> w, int n, WindowTypes windowType, SymmetryType symmeticType)
            : base(w, n, windowType, symmeticType)
        {
        }

        public override List<double> Weight()
        {
            double x = 0.0;
            double y;
            int i;

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
            for (i = 0; i < n; i++)
            {
                w[i] = 0.54 - 0.46 * Cos(y);
                y += x;
            }

            return w;
        }
    }
}
