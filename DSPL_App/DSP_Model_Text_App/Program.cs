using System;
using DSP_Model;
using DSP_Model.Filters;
using DSP_Model.Filters.IIRFilters;
using static System.Math;
using DSP_Model.Windows;
using System.Collections.Generic;

namespace DSP_Model_Text_App
{
    class Program
    {
        static void Main(string[] args)
        {
            double rp = 3;
            double rs = 60;
            double wp = 0.3;
            double ws = 0.7;
            double w0 = 0.3;
            double w1 = 0.7;
            ButterworthFilter filter = new ButterworthFilter(FilterTypes.Bandpass, wp, ws, rp, rs);
            // EllipFilter filter = new EllipFilter(FilterTypes.Bandpass, wp, ws, rp, rs);
            filter.DesignFilter(w0, w1);
            int N = 1000;
            List<double> mag = filter.GetMAG(N, 0, PI);
            List<double> phi = filter.GetPHI(N, 0, PI);


            ////for (int i = 0; i < mag.Count; i++)
            //{
            //    Console.WriteLine($"{mag[i]:0.00000}");
            //}
            //Console.WriteLine("\n\n\n");

            //for (int i = 0; i < mag.Count; i++)
            //{
            //    Console.WriteLine($"{phi[i]:0.00000}");
            //}
            //Console.WriteLine("\n\n\n");

            //for (int i = 0; i < filter.W.Count; i++)
            //{
              // Console.WriteLine($"{filter.W[i]:0.000}");
            //}

            //PrintFilterInfo(filter);

            Console.ReadKey();
        }
        private static void PrintFilterInfo(FIRFilter filter)
        {
            Console.WriteLine("Порядок фильтра: " + filter.Order);

            Console.WriteLine("H : ");
            foreach (var h in filter.H)
            {
                Console.WriteLine($"{h:0.0000}");
            }             
            
            Console.WriteLine("\nW : ");
            foreach (var w in filter.W)
            {
                Console.WriteLine($"{w:0.0000}");
            }            

        }

        private static void PrintFilterInfo(IIRFilterBase filter)
        {
            Console.WriteLine("Порядок фильтра: " + filter.Order);

            Console.WriteLine("Z : ");
            foreach (var z in filter.Z)
            {
                string signRe = (z.Real < 0) ? "-" : " ";
                string signIm = (z.Imaginary < 0) ? "-" : " ";
                Console.WriteLine($"\t{signRe}" + Math.Round(Math.Abs(z.Real), 4) + $"\t\t{signIm}" + Math.Round(Math.Abs(z.Imaginary), 4) + "i");

            }

            Console.WriteLine("\n\n\n");
            Console.WriteLine("P : ");
            Console.WriteLine("\t");
            foreach (var p in filter.P)
            {
                string signRe = (p.Real < 0) ? "-" : " ";
                string signIm = (p.Imaginary < 0) ? "-" : " ";
                Console.WriteLine($"\t{signRe}" + Math.Round(Math.Abs(p.Real), 4) + $"\t\t{signIm}" + Math.Round(Math.Abs(p.Imaginary), 4) + "i");
            }



            Console.WriteLine("\n\n\n");
            Console.WriteLine("B : ");
            Console.Write("\t");
            foreach (var b in filter.B)
            {
                Console.Write(Math.Round(b, 4) + "\t\t");
            }

            Console.WriteLine("\n\n\n");
            Console.WriteLine("A : ");
            Console.Write("\t");
            foreach (var a in filter.A)
            {
                Console.Write(Math.Round(a, 4) + "\t\t");
            }

            Console.WriteLine("\n\n\n");
            Console.WriteLine("Az : ");
            Console.Write("\t");
            foreach (var a in filter.Az)
            {
                Console.Write(Math.Round(a, 4) + "\t\t");
            }

            Console.WriteLine("\n\n\n");
            Console.WriteLine("Bz : ");
            Console.Write("\t");
            foreach (var b in filter.Bz)
            {
                Console.Write(Math.Round(b, 4) + "\t\t");
            }
        }


    }
}
