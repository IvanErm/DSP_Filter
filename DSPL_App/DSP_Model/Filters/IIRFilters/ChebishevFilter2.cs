using DSP_Model.HelperMethods;
using System.Collections.Generic;
using System.Numerics;
using static DSP_Model.HelperMethods.MathMethods;
using static System.Math;

namespace DSP_Model.Filters.IIRFilters
{
    public class ChebishevFilter2 : IIRFilterBase
    {

        public ChebishevFilter2(FilterTypes filterType, double wp, double ws, double rp, double rs)
           : base(filterType, wp, ws, rp, rs)
        {

        }

        public override int CalculateOrder()
        {
            double WPA = Tan(PI * Wp / 2);
            double WSA = Tan(PI * Ws / 2);
            double WA = WSA / WPA;

            return (int)Ceiling(Acosh(Sqrt(Pow(10, (0.1 * Abs(rs))) - 1) / (Pow(10, (0.1 * Abs(rp))) - 1)) / Acosh(WA));
        }

        public override double GetCutoffFrenq(double wp = 0, double ws = 0)
        {
            if (ws == 0 && wp == 0)
            {
                ws = this.ws;
            }
            return ws;
        }

        public override void Calculate_AP()
        {
            z = new List<Complex>();
            p = new List<Complex>();

            Calculate_AZ();

            Filter_zp2ab(this);

            double norm = a[0] / b[0];

            for (int k = 0; k < order + 1; k++)
                b[k] *= norm;
        }

        public override void Calculate_AZ()
        {
            double esp = Sqrt(Pow(10.0, Rs * 0.1) - 1.0);
            int r = order % 2;
            int L = (int)((order - r) / 2);

            double beta = MathMethods.Asinh(esp) / (double)order;
            double chb = Cosh(beta);
            double shb = Sinh(beta);


            if (r != 0)
            {
                p.Add(new Complex(-1.0 / Sinh(beta), 0.0));
            }

            double alpha;
            double ssh2, cch2;
            double sa, ca;
            for (int k = 0; k < L; k++)
            {
                alpha = PI * (double)(2 * k + 1) / (double)(2 * order);
                sa = Sin(alpha);
                ca = Cos(alpha);
                ssh2 = sa * shb;
                ssh2 *= ssh2;

                cch2 = ca * chb;
                cch2 *= cch2;

                Complex newZ = new Complex(0.0, 1.0 / ca);
                z.Add(newZ);
                z.Add(new Complex(0.0, -newZ.Imaginary));

                Complex newP = new Complex(-sa * shb / (ssh2 + cch2), ca * chb / (ssh2 + cch2));
                p.Add(newP);
                p.Add(new Complex(-sa * shb / (ssh2 + cch2), -newP.Imaginary));
            }
        }
    }
}
