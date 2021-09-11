using System.Collections.Generic;
using System.Numerics;
using static DSP_Model.HelperMethods.MathMethods;
using static System.Math;

namespace DSP_Model.Filters.IIRFilters
{
    public class EllipFilter : IIRFilterBase
    {

        public EllipFilter(FilterTypes filterType, double wp, double ws, double rp, double rs)
           : base(filterType, wp, ws, rp, rs)
        {

        }


        public override int CalculateOrder()
        {
            double WP = Tan(PI * wp / 2);
            double WS = Tan(PI * ws / 2);

            double WA = WS / WP;
            double epsilon = Sqrt(Pow(10, 0.1 * rp) - 1);
            double k1 = epsilon / Sqrt(Pow(10, 0.1 * rs) - 1);
            double k = 1 / WA;

            return (int)Ceiling((Ellipke(k * k) * Ellipke(1 - k1 * k1)) / (Ellipke(1 - k * k) * Ellipke(k1 * k1)));
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
            z = new List<Complex>();
            p = new List<Complex>();

            Calculate_AZ();
            Filter_zp2ab(this);

            double g0 = 1.0;
            if (!(order % 2 != 0))
            {
                g0 = 1.0 / Pow(10.0, Rp * 0.05);
            }

            double norm = g0 * a[0] / b[0];

            for (int k = 0; k < order + 1; k++)
                b[k] *= norm;
        }

        public override void Calculate_AZ()
        {
            double ke, u, t;
            Complex tc, v0, jv0;

            double es = Sqrt(Pow(10.0, Rs * 0.1) - 1.0);
            double ep = Sqrt(Pow(10.0, Rp * 0.1) - 1.0);
            ke = ep / es;

            int r = order % 2;
            int L = (int)((order - r) / 2);

            double k = Ellip_modulareq(order, Rp, Rs);

            // v0
            tc = new Complex(0.0, 1.0 / ep);

            v0 = Ellip_asn_cmplx(tc, 1, ke);
            t = v0.Real;

            v0 = new Complex(v0.Imaginary / (double)order, -t / (double)order);
            jv0 = new Complex(-v0.Imaginary, v0.Real);

            if (r != 0)
            {
                tc = Ellip_sn_cmplx(jv0, 1, k, tc);
                p.Add(new Complex(-tc.Imaginary, tc.Real));

            }

            for (int n = 0; n < L; n++)
            {
                u = (double)(2 * n + 1) / (double)order;

                t = Ellip_cd(u, 1, k);

                z.Add(new Complex(0.0, 1.0 / (k * t)));
                z.Add(new Complex(0.0, -1.0 / (k * t)));

                tc = new Complex(u - jv0.Real, -jv0.Imaginary);

                Complex nextP2 = Ellip_cd_cmplx(tc, 1, k, new Complex());

                Complex nextP1 = new Complex(-nextP2.Imaginary, nextP2.Real);

                nextP2 = new Complex(nextP1.Real, -nextP1.Imaginary);
                p.Add(nextP1);
                p.Add(nextP2);
            }
        }
    }
}
