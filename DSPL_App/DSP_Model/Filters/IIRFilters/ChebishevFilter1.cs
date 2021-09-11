using System;
using System.Collections.Generic;
using System.Numerics;
using static System.Math;
using static DSP_Model.HelperMethods.MathMethods;
using  DSP_Model.HelperMethods;


namespace DSP_Model.Filters.IIRFilters
{
    public class ChebishevFilter1 : IIRFilterBase
    {
        public ChebishevFilter1(FilterTypes filterType, double wp, double ws, double rp, double rs)
           : base(filterType, wp, ws, rp, rs)
        {

        }


        public override int CalculateOrder()
        {
            double Wp = Tan(PI * wp / 2);
            double Ws = Tan(PI * ws / 2);
            double Wa = Ws / Wp;

            return (int)Ceiling(Acosh(Sqrt(Pow(10, (0.1 * Abs(Rs))) - 1) / (Pow(10, (0.1 * Abs(Rp))) - 1)) / Acosh(Wa));
        }

        public override double GetCutoffFrenq(double wp = 0, double ws = 0)
        {
            if (wp == 0 && ws == 0)
            {
                wp = this.wp;
            }

            return wp;
        }


        public override void Calculate_AP()
        {
            p = new List<Complex>();
            z = new List<Complex>();
            Complex h0 = new Complex(1.0, 0.0);

            Calculate_AZ();

            Filter_zp2ab(this);

            if (!(order % 2 != 0))
                h0 = new Complex(1.0 / Pow(10.0, Rp * 0.05), 0.0);

            for (int k = 0; k < p.Count; k++)
            {
                h0 = new Complex(CMRE(h0, p[k]), CMIM(h0, p[k]));
            }

            b[0] = Math.Abs(h0.Real);

        }

        public override void Calculate_AZ()
        {
            double eps = Sqrt(Pow(10.0, Rp * 0.1) - 1.0);
            int r = order % 2;
            int L = (int)((order - r) / 2);


            double beta = MathMethods.Asinh(1.0 / eps) / (double)order;
            double chbeta = Cosh(beta);
            double shbeta = Sinh(beta);

            if (r != 0)
            {
                p.Add(new Complex(-shbeta, 0.0));
            }

            double theta;
            for (int k = 0; k < L; k++)
            {
                theta = PI * (double)(2 * k + 1) / (double)(2 * order);
                Complex newP = new Complex(-shbeta * Sin(theta), chbeta * Cos(theta));
                p.Add(newP);
                p.Add(new Complex(-shbeta * Sin(theta), -newP.Imaginary));
            }
        }
    }
}
