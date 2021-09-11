using System;
using System.Collections.Generic;
using System.Text;

namespace DSP_Model.Windows
{
    public class RectWindow : WindowWeightAbstract
    {
        public RectWindow(List<double> w, int n, WindowTypes windowType, SymmetryType symmeticType, double param)
            : base(w, n, windowType, symmeticType)
        {
        }


        public override List<double> Weight()
        {
            for (int i = 0; i < n; i++)
                w[i] = 1.0;

            return w;
        }
    }
}
