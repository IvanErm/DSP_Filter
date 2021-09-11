using System;
using System.Numerics;
using static DSP_Model.HelperMethods.MathMethods;
using static System.Math;

namespace DSP_Model.Filters.IIRFilters
{

    public class ButterworthFilter : IIRFilterBase
    {
        public ButterworthFilter(FilterTypes filterType, double wp, double ws, double rp, double rs)
           : base(filterType, wp, ws, rp, rs)
        {

        }


        public override double GetCutoffFrenq(double wp = 0, double ws = 0)
        {
            if (wp == 0 && ws == 0) 
            {
                wp = this.wp;
                ws = this.ws;
            }
            double Wp = Tan(PI * wp / 2);
            double Ws = Tan(PI * ws / 2);
            double Wa = Ws / Wp;

            double W0 = Wa / Pow(Pow(10, 0.1 * Abs(rs)) - 1, 1 / (2.0 * CalculateOrder()));
            double Wn = W0 * Wp;
            return (2 / PI) * Atan(Wn);
        }

        public override int CalculateOrder()
        {
            double Wp1 = Tan(PI * Wp / 2);
            double Ws1 = Tan(PI * Ws / 2);
            double Wa = Ws1 / Wp1;

            return (int)Math.Abs(Ceiling(Log10((Pow(10, 0.1 * Abs(rs)) - 1) / (Pow(10, 0.1 * Abs(rp)) - 1)) / (2 * Log10(Wa))));
        }

        public override void Calculate_AP()
        {
            Calculate_AZ();
            Filter_zp2ab(this);
            B[0] = A[0];
        }

        public override void Calculate_AZ()
        {
            double alpha;
            double theta;
            double eps = Sqrt(Pow(10.0, Rp * 0.1) - 1.0);
            int r = Order % 2;
            int L = (int)((Order - r) / 2);

            alpha = Pow(eps, -1.0 / (double)Order);

            if (r != 0)
            {
                P.Add(new Complex(-alpha, 0.0));
            }

            for (int k = 0; k < L; k++)
            {
                theta = PI * (double)(2 * k + 1) / (double)(2 * Order);
                P.Add(new Complex(-alpha * Sin(theta), alpha * Cos(theta)));
                P.Add(new Complex(-alpha * Sin(theta), -alpha * Cos(theta)));
            }
        }



    }

}
