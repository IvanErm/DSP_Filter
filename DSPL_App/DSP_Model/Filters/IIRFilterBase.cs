using System;
using System.Collections.Generic;
using System.Numerics;
using DSP_Model.Filters.IIRFilters;
using static System.Math;
using static DSP_Model.HelperMethods.MathMethods;


namespace DSP_Model.Filters
{
    abstract public class IIRFilterBase : FilterBase
    {
        //Задаваемые пользователем параметры
        protected FilterTypes filterType;
        protected double rp;
        protected double rs;
        protected double w0;
        protected double w1;
        protected double wp;
        protected double ws;


        //Получаемые в результате расчета параметры фильтра 
        protected int order;
        protected List<Complex> z = new List<Complex>();
        protected List<Complex> p = new List<Complex>();
        protected List<double> a = new List<double>();
        protected List<double> b = new List<double>();
        protected List<double> bz = new List<double>();
        protected List<double> az = new List<double>();
        protected List<double> at = new List<double>();
        protected List<double> bt = new List<double>();

        public FilterTypes FilterType
        {
            get { return filterType; }
            set { filterType = value; }
        }

        public double W0
        {
            get { return w0; }
            set { w0 = value;}
        }  
        
        public double W1
        {
            get { return w1; }
            set { w1 = value;}
        }

        public double Wp
        {
            get
            {
                return wp;
            }
        }

        public double Ws
        {
            get
            {
                return ws;
            }
        }

        public List<double> At
        {
            get
            {
                if (at.Count == 0)
                {
                    for (int i = 0; i < order + 1; i++)
                    {
                        at.Add(0);
                    }
                }
                return at;
            }

            set { at = value; }
        }

        public List<double> Bt
        {
            get
            {
                if (bt.Count == 0)
                {
                    for (int i = 0; i < order + 1; i++)
                    {
                        bt.Add(0);
                    }
                }
                return bt;
            }

            set { bt = value; }
        }

        public double Rp
        {
            get { return rp; }
            set { rp = value; }
        }

        public double Rs
        {
            get { return rs; }
            set { rs = value; }
        }

        public int Order
        {
            get { return order; }
            set { order = value; }
        }

        public List<Complex> Z
        {
            get { return z; }
            set { z = value; }
        }
        public List<Complex> P
        {
            get { return p; }
            set { p = value; }
        }
        public List<double> A
        {
            get { return a; }
        }
        public List<double> B
        {
            get { return b; }
        }
        public List<double> Az
        {
            get
            {
                if (az.Count == 0)
                {
                    for (int i = 0; i < order + 1; i++)
                        az.Add(0);
                }
                return az;
            }
            set { az = value; }
        }
        public List<double> Bz
        {
            get
            {
                if (bz.Count == 0)
                {
                    for (int i = 0; i < order + 1; i++)
                        bz.Add(0);
                }
                return bz;
            }
            set { bz = value; }
        }


        public IIRFilterBase(FilterTypes filterType,
                         double wp, double ws, double rp, double rs)
        {
            this.filterType = filterType;
            this.rp = rp;
            this.rs = rs;
            this.wp = wp;
            this.ws = ws;
        }

        public void DesignFilter()
        {
            Order = CalculateOrder();
            this.W0 = GetCutoffFrenq();
            //this.W1 = Filter_ws1(order, Rp, Rs);

            double wa0, wa1, ws;
            int ord_ap = order;

            if (filterType == FilterTypes.Bandpass || filterType == FilterTypes.Bandstop)
            {
                Order *= 2;
                if (order % 2 != 0)
                {
                    order++;
                }
                else
                {
                    ord_ap = order / 2;
                }
            }

            Calculate_AP();

            wa0 = Tan(W0 * PI * 0.5);
            wa1 = Tan(W1 * PI * 0.5);

            switch (filterType)
            {
                case FilterTypes.Lowpass:
                    Low2low(A, B, ord_ap, 1.0, wa0, At, Bt);
                    break;

                case FilterTypes.Highpass:
                    ws = Filter_ws1(ord_ap, Rp, Rs, GetTypeToStr());
                    Low2low(A, B, ord_ap, 1.0, 1.0 / ws, A, B);

                    Low2high(A, B, ord_ap, 1.0, wa0, At, Bt);
                    break;

                case FilterTypes.Bandpass:
                    Low2bp(A, B, ord_ap, 1.0, wa0, wa1, At, Bt);
                    break;

                case FilterTypes.Bandstop:
                    ws = Filter_ws1(ord_ap, Rp, Rs, GetTypeToStr());
                    Low2low(A, B, ord_ap, 1.0, 1.0 / ws, A, B);
                    Low2bs(A, B, ord_ap, 1.0, wa0, wa1, At, Bt);
                    break;

                default:
                    throw new Exception("Недопустимый тип фильтра!");
            }

            Bilinear();

            for (int i = 1; i <= ord_ap; i++)
            {
                Az[i] = Az[i] / Az[0];
                Bz[i] = Bz[i] / Az[0];
            }

            Bz[0] = Bz[0] / Az[0];
            Az[0] = 1.0;

        }

        public string GetTypeToStr()
        {
            if (this is ButterworthFilter)
                return "Butter";
            if (this is EllipFilter)
                return "Ellip";
            else
                return "Cheb";
        }

        public abstract int CalculateOrder();

        public abstract double GetCutoffFrenq(double wp = 0, double ws = 0);
        public abstract void Calculate_AP();
        public abstract void Calculate_AZ();


        //Частотное преобразование ФНЧ-ФНЧ
        public void Low2low(List<double> a, List<double> b, int ord, double w0, double w1, List<double> alpha, List<double> beta)
        {
            double[] num = { 0.0, 1.0 };
            double[] den = { 0.0, 0.0 };


            den[0] = w1 / w0;


            Ratcompos(a, b, num, den, 1, ord, alpha, beta);
        }

        //Частотное преобразование ФНЧ-ФВЧ
        public void Low2high(List<double> a, List<double> b, int ord, double w0, double w1, List<double> alpha, List<double> beta)
        {
            double[] num = { 0.0, 0.0 };
            double[] den = { 0.0, 1.0 };

            num[0] = w1 / w0;

            Ratcompos(a, b, num, den, 1, ord, alpha, beta);
        }

        //ФНЧ - ПФ
        public void Low2bp(List<double> a, List<double> b, int ord, double w0, double wpl, double wph, List<double> alpha, List<double> beta)
        {
            double[] num = { 0.0, 0.0, 1.0 };
            double[] den = { 0.0, 0.0, 0.0 };


            num[0] = (wph * wpl) / (w0 * w0);
            den[1] = (wph - wpl) / w0;

            Ratcompos(a, b, num, den, 2, ord, alpha, beta);
        }

        //ФНЧ - РФ
        public void Low2bs(List<double> a, List<double> b, int ord, double w0, double wsl, double wsh, List<double> alpha, List<double> beta)
        {
            double[] den = { 0.0, 0.0, 1.0 };
            double[] num = { 0.0, 0.0, 0.0 };


            den[0] = (wsh * wsl) / (w0 * w0);
            num[1] = (wsh - wsl) / w0;

            Ratcompos(a, b, num, den, 2, ord, alpha, beta);
        }

        //Функция рассчитывает композицию вида Y(s) = (H \circ F)(s) = H(F(s))
        private void Ratcompos(List<double> a1, List<double> b1, double[] c, double[] d, int p, int ord, List<double> a2, List<double> b2)
        {
            int pd, ld;

            int k2 = (ord * p) + 1;
            int nk2s = (ord + 1) * k2;

            double[] num = new double[nk2s];
            double[] den = new double[nk2s];
            double[] ndn = new double[nk2s];
            double[] ndd = new double[nk2s];

            Memset(num, 0, nk2s, 0);
            Memset(den, 0, nk2s, 0);
            Memset(ndn, 0, nk2s, 0);
            Memset(ndd, 0, nk2s, 0);


            num[0] = 1.0;
            den[0] = 1.0;
            int pn = 0;
            int ln = 1;
            for (int i = 1; i < ord + 1; i++)
            {

                Conv(num, pn, ln, c, 0, p + 1, num, pn + k2);

                Conv(den, pn, ln, d, 0, p + 1, den, pn + k2);

                pn += k2;
                ln += p;
            }

            pn = 0;
            pd = ord * k2;
            ln = 1;
            ld = k2;

            for (int i = 0; i < ord + 1; i++)
            {
                Conv(num, pn, ln, den, pd, ld, ndn, i * k2);

                ln += p;
                ld -= p;
                pn += k2;
                pd -= k2;
            }

            for (int i = 0; i < ord + 1; i++)
            {
                for (int k = 0; k < k2; k++)
                {
                    ndd[i * k2 + k] = ndn[i * k2 + k] * a1[i];
                    ndn[i * k2 + k] *= b1[i];
                }
            }



            Memset(a2, 0, k2, 0);
            Memset(b2, 0, k2, 0);

            for (int k = 0; k < k2; k++)
            {
                for (int i = 0; i < ord + 1; i++)
                {
                    b2[k] += ndn[i * k2 + k];
                    a2[k] += ndd[i * k2 + k];
                }
            }
        }

        //Линейная свертка двух вещественных векторов
        private void Conv(double[] a1, int a1Start, int na, double[] b1, int b1Start, int nb, double[] c1, int c1Start)
        {
            double[] t;

            int bufsize = na + nb - 1;
            if ((a1 != c1) && (b1 != c1))
            {
                int len = c1.Length;
                t = c1;

                t = new double[len];
            }
            else
            {
                t = new double[bufsize];
                Memset(t, 0, bufsize, 0);
            }

            for (int k = 0; k < na; k++)
                for (int n = 0; n < nb; n++)
                    t[k + n] += a1[k + a1Start] * b1[n + b1Start];

            if (t != c1)
            {
                Memcpy(c1, t, bufsize, c1Start);
            }
        }

        //Билинейное преобразование
        public void Bilinear()
        {
            double[] c = { 1.0, -1.0 };
            double[] d = { 1.0, 1.0 };


            Ratcompos(At, Bt, c, d, 1, order, Az, Bz);
        }


        public override List<double> GetMAG(int n, double min, double max)
        {
            W = GetW(n, min, max);

            return CalclMAG(Az, Bz, Order, W, n);
        }

        public override List<double> GetPHI(int n, double min, double max)
        {
            W = GetW(n, min, max);

            return CalcPHI(Az, Bz, Order, W, n);
        }


    }
}
