using DSP_Model.Windows;
using System;
using System.Collections.Generic;
using static System.Math;
using static DSP_Model.HelperMethods.MathMethods;

namespace DSP_Model.Filters
{
    public class FIRFilter : FilterBase
    {
        private int order;
        private double w0;
        private double w1;
        private WindowTypes winowType;
        private double winParameter;
        private WindowWeightAbstract winWeighting;
        private List<double> h = new List<double>();
        private FilterTypes filterType;
        private SymmetryType symmetricType;
        private double rp;
        private double rs;
        private double wp;
        private double ws;

        public double Rp
        {
            get { return rp; }
        }

        public double Rs
        {
            get { return rs; }
        }

        public double Wp
        {
            get { return wp; }
        }

        public double Ws
        {
            get { return ws; }
        }

        public int Order
        {
            get { return order; }
            set { order = value; }
        }

        public double W0
        {
            get { return w0; }
            set { w0 = value; }
        }

        public double W1
        {
            get { return w1; }
            set { w1 = value; }
        }

        public WindowTypes WinType
        {
            get { return winowType; }
            set { winowType = value; }
        }

        public double WinParameter
        {
            get { return winParameter; }
            set { winParameter = value; }
        }

        public List<double> H
        {
            get { return h; }
            set { h = value; }
        }

        public FilterTypes FilterType
        {
            get { return filterType; }
            set { filterType = value; }
        }

        public FIRFilter(int order, FilterTypes type, double rp,
            double rs, double wp, double ws)
        {
            this.Order = order;
            this.filterType = type;
            this.rp = rp;
            this.rs = rs;
            this.wp = wp;
            this.ws = ws;
        }

        public void DesignFilter(WindowTypes winType)
        {
            this.winowType = winType;
            this.W0 = Wp;
            this.W1 = Ws;

            h = Fir_linphase(W0, W1, filterType, winType);
        }

        //Расчет коэффициентов линейно-фазового КИХ-фильтра методом оконного взвешивания
        public List<double> Fir_linphase(double w0l, double w11, FilterTypes fType, WindowTypes winType, List<double> H = null)
        {
            if (H == null)
            {
                H = new List<double>();
                for (int i = 0; i < order + 1; i++)
                {
                    H.Add(0);
                }
            }
            switch (fType)
            {
                case FilterTypes.Lowpass:
                    Fir_linphase_lpf(w0l, winType, H);
                    break;

                case FilterTypes.Highpass:
                    Fir_linphase_lpf(1.0 - w0l, winType, H);

                    for (int n = 0; n < order + 1; n += 2)
                        H[n] = -H[n];

                    break;

                case FilterTypes.Bandpass:
                    if (w11 < w0l)
                        throw new Exception("w0 должно быть больше w1");


                    double wc = (w0l + w11) * 0.5;
                    double b = w11 - w0l;
                    Fir_linphase_lpf(b * 0.5, winType, H);


                    double del = 0.5 * (double)order;
                    for (int n = 0; n < order + 1; n++)
                        H[n] *= 2.0 * Cos(PI * ((double)n - del) * wc);

                    break;

                /* BandStop FIR coefficients calculation */
                case FilterTypes.Bandstop:
                    {

                        if (order % 2 != 0)
                            order *= 2;
                            //throw new Exception("Порядок фильтра должен быть четным");


                        List<double> h0 = new List<double>();
                        for (int i = 0; i < order + 1; i++)
                            h0.Add(0);

                        Fir_linphase(w0l, 0.0, FilterTypes.Lowpass, winType, h0);

                        Fir_linphase(w11, 0.0, FilterTypes.Highpass, winType, H);

                        for (int n = 0; n < order + 1; n++)
                            H[n] += h0[n];

                        break;
                    }
            }
            return H;
        }

        // Расчет ФНЧ
        public void Fir_linphase_lpf(double wp, WindowTypes winType, List<double> h)
        {
            w = CreateInitializeList(order + 1);
            w = Linspace(-(double)order * 0.5, (double)order * 0.5, order + 1, SymmetryType.Symmetric, w);

            Sinc(w, order + 1, PI * wp, h);

            winWeighting = CreateWinWeigthing(w, order + 1);
            w = winWeighting.Weight();

            for (int n = 0; n < order + 1; n++)
                h[n] *= w[n] * wp;

        }

        // Расчет функции оконного взвешивания
        private WindowWeightAbstract CreateWinWeigthing(List<double> w, int n)
        {
            return new HammingWindow(w, n, winowType, symmetricType);
        }

        public List<double> CreateInitializeList(int n)
        {
            List<double> list = new List<double>();
            for (int i = 0; i < n; i++)
            {
                list.Add(0);
            }
            return list;
        }


        public override List<double> GetMAG(int n, double min, double max)
        {
            W = GetW(n, min, max);

            return CalclMAG(null, H, Order, W, n);
        }

        public override List<double> GetPHI(int n, double min, double max)
        {
            W = GetW(n, min, max);

            return CalcPHI(null, H, Order, W, n);
        }
    }
}





