using System.Collections.Generic;

namespace DSP_Model.Windows
{
    public abstract class WindowWeightAbstract
    {
        protected int n;
        protected List<double> w;
        protected WindowTypes windowType;
        protected SymmetryType symmeticType;
        protected double param = 0;

        public WindowWeightAbstract(List<double> w, int n, WindowTypes windowType, SymmetryType symmeticType)
        {
            this.w = w;
            this.n = n;
            this.windowType = windowType;
            this.symmeticType = symmeticType;
        }

        public abstract List<double> Weight();

    }
}
