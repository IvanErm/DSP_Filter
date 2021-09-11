using DSP_Model.Filters;
using System;
using System.Collections.Generic;
using System.Numerics;
using static System.Math;

namespace DSP_Model.HelperMethods
{
    public static class MathMethods
    {
        private static int ELLIP_ITER = 16;
        //Расчет ареакосинуса
        public static double Acosnh(double x)
        {
            return Log10(x + Sqrt(Pow(x, 2) - 1));
        }

        //Расчет ареасинуса
        public static double Asinh(double x)
        {
            return Log(x + Sqrt(Pow(x, 2) + 1));
        }

        // Многочлен Чебышева первого рода порядка ord
        public static double[] Cheby_poly1(double x, int n, int ord, double[] y)
        {
            double[] t = new double[2];

            if (ord == 0)
            {
                for (int k = 0; k < n; k++)
                {
                    y[0] = 1.0;
                }
            }

            else if (ord == 1)
            {
                y[0] = x;
            }
            else
            {
                int m;
                for (int k = 0; k < n; k++)
                {
                    m = 2;
                    t[1] = x;
                    t[0] = 1.0;
                    while (m <= ord)
                    {
                        y[k] = 2.0 * x * t[1] - t[0];
                        t[0] = t[1];
                        t[1] = y[k];
                        m++;
                    }
                }
            }

            return y;
        }

        // Модифицированная функция Бесселя первого рода I_0(x)
        public static double Bessel_i0(double x, int n)
        {
            double[] P16 = { 1.0000000000000000000000801e+00,
                       2.4999999999999999999629693e-01,
                       2.7777777777777777805664954e-02,
                       1.7361111111111110294015271e-03,
                       6.9444444444444568581891535e-05,
                       1.9290123456788994104574754e-06,
                       3.9367598891475388547279760e-08,
                       6.1511873265092916275099070e-10,
                       7.5940584360755226536109511e-12,
                       7.5940582595094190098755663e-14,
                       6.2760839879536225394314453e-16,
                       4.3583591008893599099577755e-18,
                       2.5791926805873898803749321e-20,
                       1.3141332422663039834197910e-22,
                       5.9203280572170548134753422e-25,
                       2.0732014503197852176921968e-27,
                       1.1497640034400735733456400e-29};

            double[] P22 = { 3.9894228040143265335649948e-01,
                       4.9867785050353992900698488e-02,
                       2.8050628884163787533196746e-02,
                       2.9219501690198775910219311e-02,
                       4.4718622769244715693031735e-02,
                       9.4085204199017869159183831e-02,
                      -1.0699095472110916094973951e-01,
                       2.2725199603010833194037016e+01,
                      -1.0026890180180668595066918e+03,
                       3.1275740782277570164423916e+04,
                      -5.9355022509673600842060002e+05,
                       2.6092888649549172879282592e+06,
                       2.3518420447411254516178388e+08,
                      -8.9270060370015930749184222e+09,
                       1.8592340458074104721496236e+11,
                      -2.6632742974569782078420204e+12,
                       2.7752144774934763122129261e+13,
                      -2.1323049786724612220362154e+14,
                       1.1989242681178569338129044e+15,
                      -4.8049082153027457378879746e+15,
                       1.3012646806421079076251950e+16,
                      -2.1363029690365351606041265e+16,
                       1.6069467093441596329340754e+16};


            double x2;
            double y = 0.0;
            for (int k = 0; k < n; k++)
            {

                if (x < 7.75)
                {
                    x2 = x * x * 0.25;
                    y = Polyval(P16, 16, x2, 1);
                    y = x2 * y + 1.0;
                }
                else
                {
                    x2 = 1.0 / x;
                    y = Polyval(P22, 22, x2, 1);
                    y *= Exp(x) / Sqrt(x);
                }
            }

            return y;
        }

        //Расчет вещественного полинома
        private static double Polyval(double[] a, int ord, double x, int n)
        {
            double y = 0;

            for (int k = 0; k < n; k++)
            {
                y = a[ord];
                for (int m = ord - 1; m > -1; m--)
                    y = y * x + a[m];
            }

            return y;
        }

        // Функция рассчитывает значения функции для вещественного вектора x
        public static void Sinc(List<double> x, int n, double a, List<double> y)
        {
            for (int k = 0; k < n; k++)
                y[k] = (x[k] == 0.0) ? 1.0 : Sin(a * x[k]) / (a * x[k]);
        }

        //Копирует элементы из одного массива комплексных чисел в другой
        public static void Memcpy(List<Complex> a1, List<Complex> b1, int n)
        {
            for (int i = 0; i < n; i++)
            {
                a1[i] = new Complex(b1[i].Real, b1[i].Imaginary);
            }
        }

        //Копирует элементы из одного массива чисел с плавающей запятой в другой
        public static void Memcpy(double[] a1, double[] b1, int n, int start)
        {
            for (int i = start, j = 0; j < n && i < a1.Length; i++, j++)
            {
                a1[i] = b1[j];
            }
        }

        //Копирует элементы из одного массива чисел с плавающей запятой в другой
        public static void Memcpy(List<double> a1, List<double> b1, int n)
        {
            for (int i = 0; i < n; i++)
            {
                a1[i] = b1[i];
            }
        }

        // CMIM(a,b) = (RE(a)) * (IM(b)) + (IM(a)) * (RE(b))
        public static double CMIM(Complex a1, Complex b1)
        {
            return (a1.Real * b1.Imaginary) + (a1.Imaginary * b1.Real);
        }

        // CMRE(a,b) = (RE(a)) * (RE(b)) - (IM(a)) * (IM(b))
        public static double CMRE(Complex a1, Complex b1)
        {
            return (a1.Real * b1.Real) - (a1.Imaginary * b1.Imaginary);
        }

        //Заполняет массив комплексных указанными значения
        public static void Memset(List<Complex> o, double c, int size, int start)
        {
            for (int i = start; i < start + size; i++)
            {
                o[i] = new Complex(c, c);
            }
        }

        //Заполняет массив чисел  плавающей запятой указанными значения
        public static void Memset(double[] o, double c, int size, int start)
        {
            for (int i = start; i < start + size; i++)
            {
                o[i] = c;
            }
        }

        //Заполняет массив комплексных указанными значения
        public static void Memset(List<double> o, double c, int size, int start)
        {
            for (int i = start; i < start + size; i++)
            {
                if (o.Count >= size)
                    o[i] = c;
                else
                    o.Add(c);
            }
        }

        //Преобразование массива комплексных данных в два массива вещественных данных, 
        //содержащих реальную и мнимую части исходного массива
        public static void Cmplx2re(List<Complex> x, int n, List<double> re, List<double> im)
        {
            if (re != null)
            {
                for (int k = 0; k < n; k++)
                    re.Add(x[k].Real);
            }
            if (im != null)
            {
                for (int k = 0; k < n; k++)
                    im.Add(x[k].Imaginary);
            }
        }

        public static void Poly_z2a_cmplx(List<Complex> w, int ord, List<Complex> q)
        {
            List<Complex> x = new List<Complex>();
            x.Add(new Complex());
            x.Add(new Complex(1.0, 0.0));


            Memset(q, 0, ord + 1, 0);
            q[0] = new Complex(1.0, 0.0);

            int ind = 1;
            for (int k = 0; k < w.Count; k++)
            {
                x[0] = new Complex(-w[k].Real, -w[k].Imaginary);
                q = Conv_cmplx(q, ind, x, 2, q);
                ind++;
            }
        }

        //Функция рассчитывает линейную свертку двух векторов
        public static List<Complex> Conv_cmplx(List<Complex> a1, int na, List<Complex> b1, int nb, List<Complex> c)
        {
            List<Complex> t;

            if ((a1 != c) && (b1 != c))
                t = c;
            else
                t = new List<Complex>();

            for (int i = 0; i < na + nb - 1; i++)
                t.Add(new Complex());

            Memset(t, 0, na + nb - 1, 0);


            for (int k = 0; k < na; k++)
            {
                for (int n = 0; n < nb; n++)
                {
                    t[k + n] = new Complex(
                        real: t[k + n].Real + CMRE(a1[k], b1[n]),
                        imaginary: t[k + n].Imaginary + CMIM(a1[k], b1[n]));
                }
            }

            if (t != c)
            {
                Memcpy(c, t, na + nb - 1);
            }
            return c;
        }

        //Функция, преобразующая нули и полюса в коэффициенты цифрового фильтра
        public static void Filter_zp2ab(IIRFilterBase filter)
        {
            List<Complex> acc = new List<Complex>();
            for (int i = 0; i < filter.Order + 1; i++)
                acc.Add(new Complex());

            Poly_z2a_cmplx(filter.Z, filter.Order, acc);

            Cmplx2re(acc, filter.Order + 1, filter.B, null);

            Poly_z2a_cmplx(filter.P, filter.Order, acc);

            Cmplx2re(acc, filter.Order + 1, filter.A, null);

        }

        //Обратная эллиптическая функция Якоби
        public static Complex Ellip_acd_cmplx(Complex w, int n, double k)
        {
            double[] lnd = new double[ELLIP_ITER];
            double t;
            Complex tmp0 = new Complex();
            Complex tmp1 = new Complex();

            lnd = Ellip_landen(k, ELLIP_ITER);

            Complex u = new Complex();
            for (int m = 0; m < n; m++)
            {
                u = new Complex(w.Real, w.Imaginary);

                for (int i = 1; i < ELLIP_ITER; i++)
                {
                    tmp0 = new Complex(lnd[i - 1] * u.Real, lnd[i - 1] * u.Imaginary);
                    tmp1 = new Complex(1.0 - CMRE(tmp0, tmp0), -CMIM(tmp0, tmp0));


                    tmp0 = Sqrt_cmplx(tmp1, 1);
                    tmp0 = new Complex(tmp0.Real + 1.0, tmp0.Imaginary);

                    tmp1 = new Complex(tmp0.Real * (1.0 + lnd[i]), tmp0.Imaginary * (1.0 + lnd[i]));

                    t = 2.0 / Abssqr(tmp1);
                    tmp0 = new Complex(t * CMCONJRE(u, tmp1), t * CMCONJIM(u, tmp1));

                    u = new Complex(tmp0.Real, tmp0.Imaginary);

                }

                Acos_cmplx(tmp0, 1, u);

                t = 2.0 / PI;
                u = new Complex(u.Real * t, u.Imaginary * t);
            }

            return u;
        }

        // Арккосинус комплексного аргумента x
        public static Complex Acos_cmplx(Complex x, int n, Complex y)
        {
            double pi2 = 0.5 * PI;

            y = Asin_cmplx(x, n, y);

            for (int k = 0; k < n; k++)
            {
                y = new Complex(pi2 - y.Real, -y.Imaginary);

            }
            return y;
        }

        // Полный эллептический интегралл
        public static double Ellipke(double m, double tol = 0.1)
        {
            double a0 = 1;
            double b0 = Sqrt(1 - m);
            double s0 = m;
            double i1 = 0;
            double mm = Double.MaxValue;
            double a1 = 0;

            while (mm > tol)
            {
                a1 = (a0 + b0) / 2;
                double b1 = Sqrt(a0 * b0);
                double c1 = (a0 - b0) / 2;
                i1 = i1 + 1;
                double w1 = Pow(2, i1) * Pow(c1, 1);

                mm = Min(mm, w1);

                s0 = s0 + w1;
                a0 = a1;
                b0 = b1;
            }

            return PI / (2 * a1);
        }

        public static double Ellip_modulareq(int ord, double rp, double rs)
        {
            double t, sn = 0.0;

            double ep = Sqrt(Pow(10.0, rp * 0.1) - 1.0);
            double es = Sqrt(Pow(10.0, rs * 0.1) - 1.0);
            double ke = ep * 1.0 / es;
            ke = Sqrt(1.0 - ke * ke);

            int r = ord % 2;
            int L = (ord - r) / 2;

            double kp = 1.0;

            for (int i = 0; i < L; i++)
            {
                t = (double)(2 * i + 1) / (double)ord;
                sn = Ellip_sn(t, 1, ke)[0];
                sn *= sn;
                kp *= sn * sn;
            }
            kp *= Pow(ke, (double)ord);

            return Sqrt(1.0 - kp * kp);
        }

        //Эллиптическая функция Якоби
        public static double[] Ellip_sn(double u, int n, double k)
        {

            double[] lnd = new double[ELLIP_ITER];
            double[] y = new double[n];
            lnd = Ellip_landen(k, ELLIP_ITER);

            for (int m = 0; m < n; m++)
            {
                y[m] = Sin(u * PI * 0.5);
                for (int i = ELLIP_ITER - 1; i > 0; i--)
                {
                    y[m] = (1.0 + lnd[i]) / (1.0 / y[m] + lnd[i] * y[m]);
                }
            }
            return y;
        }

        //Расчет коэффициентов ряда полного эллиптического интеграла
        public static double[] Ellip_landen(double k, int n)
        {
            double[] y = new double[n];
            y[0] = k;

            for (int i = 1; i < n; i++)
            {
                y[i] = y[i - 1] / (1.0 + Sqrt(1.0 - y[i - 1] * y[i - 1]));
                y[i] *= y[i];
            }

            return y;
        }

        //Функция рассчитывает значения значения обратной эллиптической функции\
        public static Complex Ellip_asn_cmplx(Complex w, int n, double k)
        {
            double t;
            double[] lnd = new double[ELLIP_ITER];
            Complex tmp0 = new Complex();
            Complex tmp1;

            lnd = Ellip_landen(k, ELLIP_ITER);

            Complex u = new Complex();
            for (int m = 0; m < n; m++)
            {
                u = new Complex(w.Real, w.Imaginary);

                for (int i = 1; i < ELLIP_ITER; i++)
                {
                    tmp0 = new Complex(lnd[i - 1] * u.Real, lnd[i - 1] * u.Imaginary);
                    tmp1 = new Complex(1.0 - CMRE(tmp0, tmp0), -CMIM(tmp0, tmp0));

                    tmp0 = Sqrt_cmplx(tmp1, 1);

                    tmp0 = new Complex(tmp0.Real + 1.0, tmp0.Imaginary);
                    tmp1 = new Complex(tmp0.Real * (1.0 + lnd[i]), tmp0.Imaginary * (1.0 + lnd[i]));

                    t = 2.0 / Abssqr(tmp1);

                    tmp0 = new Complex(t * CMCONJRE(u, tmp1), t * CMCONJIM(u, tmp1));

                    u = new Complex(tmp0.Real, tmp0.Imaginary);
                }

                u = Asin_cmplx(tmp0, 1, Complex.Add(u, m));
                t = 2.0 / PI;

                u = new Complex(u.Real * t, u.Imaginary * t);
                ;
            }
            return u;
        }

        //Квадратный корень из комплексного вектора x (поэлементный)
        public static Complex Sqrt_cmplx(Complex x, int n)
        {
            int k;
            double r, zr, at;
            Complex t;
            Complex y = new Complex();

            for (k = 0; k < n; k++)
            {
                r = x.Magnitude;
                if (Math.Round(r, 4) == 0.0000)
                {
                    y = new Complex(0.0, 0.0);
                }
                else
                {
                    t = new Complex(x.Real + r, x.Imaginary);
                    at = t.Magnitude;
                    if (Math.Round(at, 4) == 0.0000)
                    {
                        y = new Complex(0.0, Sqrt(r));
                    }
                    else
                    {
                        zr = 1.0 / t.Magnitude;
                        r = Sqrt(r);
                        y = new Complex(t.Real * zr * r, t.Imaginary * zr * r);
                    }
                }
            }
            return y;
        }

        // Возвращает квадрат модуля комплексного числа x:
        // SQR(RE(x)) + SQR(IM(x))
        public static double Abssqr(Complex x)
        {
            return Pow(x.Real, 2) + Pow(x.Imaginary, 2);
        }

        // ((IM(a)) * (RE(b)) - (RE(a)) * (IM(b)))
        public static double CMCONJIM(Complex a, Complex b)
        {
            return (a.Imaginary * b.Real) - (a.Real * b.Imaginary);
        }

        // ((RE(a)) * (RE(b)) + (IM(a)) * (IM(b)))
        public static double CMCONJRE(Complex a, Complex b)
        {
            return (a.Real * b.Real) + (a.Imaginary * b.Imaginary);
        }

        // Арксинус комплексного аргумента x
        public static Complex Asin_cmplx(Complex x, int n, Complex y)
        {
            int k;
            Complex tmp;

            for (k = 0; k < n; k++)
            {
                tmp = new Complex(1.0 - CMRE(x, x), -CMIM(x, x));
                y = Sqrt_cmplx(tmp, 1);

                y = new Complex(y.Real - x.Imaginary, y.Imaginary + x.Real);

                tmp = Log_cmplx(Complex.Add(y, k), 1);
                y = new Complex(tmp.Imaginary, -tmp.Real);
            }
            return y;
        }

        // Натуральный логарифм комплексного аргумента x
        public static Complex Log_cmplx(Complex x, int n)
        {
            Complex y = new Complex();
            for (int k = 0; k < n; k++)
            {
                y = new Complex(0.5 * Log(Abssqr(x)), Atan2(x.Imaginary, x.Real));
            }
            return y;
        }

        //Эллиптическая функция Якоби
        public static Complex Ellip_sn_cmplx(Complex u, int n, double k, Complex y)
        {
            double[] lnd = new double[ELLIP_ITER];
            double t;
            Complex tmp;

            lnd = Ellip_landen(k, ELLIP_ITER);


            for (int m = 0; m < n; m++)
            {
                tmp = new Complex(u.Real * PI * 0.5, u.Imaginary * PI * 0.5);

                y = Sin_cmplx(tmp, 1, Complex.Add(y, m));

                for (int i = ELLIP_ITER - 1; i > 0; i--)
                {
                    t = 1.0 / Abssqr(y);
                    tmp = new Complex(y.Real * t + y.Real * lnd[i], -y.Imaginary * t + y.Imaginary * lnd[i]);
                    t = (1.0 + lnd[i]) / Abssqr(tmp);
                    y = new Complex(tmp.Real * t, -tmp.Imaginary * t);
                }
            }
            return y;
        }

        // Синус комплексного аргумента x
        public static Complex Sin_cmplx(Complex x, int n, Complex y)
        {
            double ep, em, sx, cx; ;

            for (int k = 0; k < n; k++)
            {
                ep = Exp(x.Imaginary);
                em = Exp(-x.Imaginary);
                sx = 0.5 * Sin(x.Real);
                cx = 0.5 * Cos(x.Real);
                y = new Complex(sx * (em + ep), cx * (ep - em));

            }
            return y;
        }

        // Эллиптическая функция Якоб
        public static double Ellip_cd(double u, int n, double k)
        {
            double[] lnd = new double[ELLIP_ITER];

            lnd = Ellip_landen(k, ELLIP_ITER);
            double y = 0;
            for (int m = 0; m < n; m++)
            {
                y = Cos(u * PI * 0.5);
                for (int i = ELLIP_ITER - 1; i > 0; i--)
                {
                    y = (1.0 + lnd[i]) / (1.0 / y + lnd[i] * y);
                }
            }
            return y;
        }

        //Эллиптическая функция Якоби
        public static Complex Ellip_cd_cmplx(Complex u, int n, double k, Complex y)
        {
            double t;
            Complex tmp;

            double[] lnd = Ellip_landen(k, ELLIP_ITER);

            for (int m = 0; m < n; m++)
            {
                tmp = new Complex(u.Real * PI * 0.5, u.Imaginary * PI * 0.5);

                y = Cos_cmplx(tmp, 1, Complex.Add(y, m));

                for (int i = ELLIP_ITER - 1; i > 0; i--)
                {
                    t = 1.0 / Abssqr(y);
                    tmp = new Complex(y.Real * t + y.Real * lnd[i], -y.Imaginary * t + y.Imaginary * lnd[i]);

                    t = (1.0 + lnd[i]) / Abssqr(tmp);

                    y = new Complex(tmp.Real * t, -tmp.Imaginary * t);
                }
            }
            return y;
        }

        //Косинус комплексного аргумента x
        public static Complex Cos_cmplx(Complex x, int n, Complex y)
        {
            double ep, em, sx, cx;

            for (int k = 0; k < n; k++)
            {
                ep = Exp(x.Imaginary);
                em = Exp(-x.Imaginary);
                sx = 0.5 * Sin(x.Real);
                cx = 0.5 * Cos(x.Real);
                y = new Complex(cx * (em + ep), sx * (em - ep));

            }
            return y;
        }

        //Расчет частоты Ws
        public static double Filter_ws1(int ord, double Rp, double Rs, string IIRType)
        {
            double es2 = Pow(10.0, Rs * 0.1) - 1.0;
            double ep2 = Pow(10.0, Rp * 0.1) - 1.0;
            double gs2 = 1.0 / (1.0 + es2);

            double x = (1.0 - gs2) / (gs2 * ep2);
            double ws;

            switch (IIRType)
            {
                case "Butter":
                    ws = Pow(x, 0.5 / (double)ord);
                    break;

                case "Cheb":
                    x = Sqrt(x) + Sqrt(x - 1.0);
                    x = Log(x) / (double)ord;
                    ws = 0.5 * (Exp(-x) + Exp(x));
                    break;
                case "Ellip":
                    Complex y, z;
                    double k = Sqrt(ep2 / es2);
                    double k1 = Ellip_modulareq(ord, Rp, Rs);
                    z = new Complex(Sqrt(x), 0.0);
                    y = Ellip_acd_cmplx(z, 1, k);
                    y = new Complex(y.Real / (double)ord, y.Imaginary / (double)ord);
                    z = Ellip_cd_cmplx(y, 1, k1, z);
                    ws = z.Real;
                    break;
                default:
                    ws = -1.0;
                    break;
            }


            return ws;
        }

        // Преобразование массива вещественных данных в массив комплексных
        public static void Re2cmplx(List<double> x, int n, List<Complex> y)
        {
            for (int k = 0; k < n; k++)
            {
                y[k] = new Complex(x[k], 0.0);
            }
        }

        public static void Unwrap(List<double> phi, int n, double lev, double mar)
        {
            double[] a = { 0.0, 0.0 };
            double d;
            int k;
            int flag = 1;

            double th = mar * lev;
            while (flag != 0)
            {
                flag = 0;
                a[0] = a[1] = 0.0;
                for (k = 0; k < n - 1; k++)
                {
                    d = phi[k + 1] - phi[k];
                    if (d > th)
                    {
                        a[0] -= lev;
                        flag = 1;
                    }
                    if (d < -th)
                    {
                        a[0] += lev;
                        flag = 1;
                    }
                    phi[k] += a[1];
                    a[1] = a[0];
                }
                phi[n - 1] += a[1];
            }

        }

        // Расчет комплексного коэффициента передачи
        public static void Freqz(List<double> b, List<double> a, int ord, List<double> w, int n, List<Complex> h)
        {
            Complex jw;
            List<Complex> bc = new List<Complex>();
            List<Complex> ac = new List<Complex>();
            List<Complex> num = new List<Complex>();
            num.Add(new Complex(0, 0));

            List<Complex> den = new List<Complex>();
            den.Add(new Complex(0, 0));

            double mag;
            int k;


            for (int i = 0; i < ord + 1; i++)
            {
                bc.Add(new Complex());
            }

            Re2cmplx(b, ord + 1, bc);

            if (a != null)
            {
                /* IIR filter if a != NULL */
                for (int i = 0; i < ord + 1; i++)
                {
                    ac.Add(new Complex(0, 0));
                }

                Re2cmplx(a, ord + 1, ac);

                for (k = 0; k < n; k++)
                {
                    jw = new Complex(Cos(w[k]), -Sin(w[k]));

                    Polyval_cmplx(bc, ord, jw, 1, num);
                    Polyval_cmplx(ac, ord, jw, 1, den);

                    mag = Abssqr(den[0]);
                    mag = 1.0 / mag;
                    h[k] = new Complex(CMCONJRE(num[0], den[0]) * mag, CMCONJIM(num[0], den[0]) * mag);
                }
            }
            else
            {
                /* FIR filter if a == NULL */
                for (k = 0; k < n; k++)
                {
                    jw = new Complex(Cos(w[k]), -Sin(w[k]));

                    Polyval_cmplx(bc, ord, jw, 1, h, k);
                }
            }

        }

        // Расчет комплексного полинома
        public static void Polyval_cmplx(List<Complex> a, int ord, Complex x, int n, List<Complex> y, int start = 0)
        {
            int k;
            Complex t = new Complex();

            for (k = start; k < n + start; k++)
            {
                y[k] = new Complex(a[ord].Real, a[ord].Imaginary);

                for (int m = ord - 1; m > -1; m--)
                {
                    t = new Complex(CMRE(y[k], x), CMIM(y[k], x));
                    y[k] = new Complex(t.Real + a[m].Real, t.Imaginary + a[m].Imaginary);
                }
            }

        }

        // Расчет амплитудно-частотной (АЧХ) характеристики
        public static List<double> CalclMAG(List<double> a, List<double> b, int ord, List<double> w, int n)
        {
            List<Complex> hc = new List<Complex>();
            List<double> mag = new List<double>();
            for (int i = 0; i < n; i++)
            {
                hc.Add(0);
                mag.Add(0);
            }

            Freqz(b, a, ord, w, n, hc);


            for (int k = 0; k < n; k++)
                mag[k] = Sqrt(Abssqr(hc[k]));
            

            return mag;
        }


        // Расчет фазочастотной характеристики (ФЧХ)
        public static List<double> CalcPHI(List<double> a, List<double> b, int ord, List<double> w, int n)
        {
            List<Complex> hc = new List<Complex>();
            List<double> phi = new List<double>();
            for (int i = 0; i < n; i++)
            {
                hc.Add(0);
                phi.Add(0);
            }

            Freqz(b, a, ord, w, n, hc);


            for (int k = 0; k < n; k++)
                phi[k] = Atan2(hc[k].Imaginary, hc[k].Real);

            return phi;
        }


        //Функция заполняет массив линейно-нарастающими, равноотстоящим значениями от x0 до x1
        public static List<double> Linspace(double x0, double x1, int n, SymmetryType symmetricType, List<double> x)
        {
            if (n < 2)
                throw new Exception("Количество точек массива n должно быть больше 1!");
            if (x.Count == 0)
                throw new Exception("Массивы не может быть пустым!");

            double dx;
            switch (symmetricType)
            {
                case SymmetryType.Symmetric:
                    dx = (x1 - x0) / (double)(n - 1);
                    x[0] = x0;
                    for (int k = 1; k < n; k++)
                        x[k] = x[k - 1] + dx;
                    break;
                case SymmetryType.Periodic:
                    dx = (x1 - x0) / (double)n;
                    x[0] = x0;
                    for (int k = 1; k < n; k++)
                        x[k] = x[k - 1] + dx;
                    break;
            }

            return x;
        }

    }
}
