using System.Collections.Generic;
using static DSP_Model.HelperMethods.MathMethods;

namespace DSP_Model.Filters
{
    public abstract class FilterBase
    {
        protected List<double> w = new List<double>();
        public List<double> W
        {
            get { return w; }
            set { w = value; }
        }

        public abstract List<double> GetMAG(int n, double min, double max);

        public abstract List<double> GetPHI(int n, double min, double max);

        public List<double> GetW(int n, double min, double max)
        {

            List<double> w = new List<double>();
            for (int i = 0; i < n; i++)
            {
                w.Add(0);
            }
            return Linspace(min, max, n, SymmetryType.Symmetric, w);
        }

    }
}
