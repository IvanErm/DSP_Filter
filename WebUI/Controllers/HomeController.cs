using DSP_Model;
using DSP_Model.Filters;
using DSP_Model.Filters.IIRFilters;
using DSP_Model.Windows;
using Microsoft.AspNetCore.Mvc;
using System;
using System.Collections.Generic;
using WebUI.ViewModel;

namespace WebUI.Controllers
{
    public class HomeController : Controller
    {
        List<double> MFC = new List<double>(); // АЧХ
        List<double> PFC = new List<double>(); // ФЧХ

        [HttpGet]
        public IActionResult Index()
        {
            return View();
        }      
        
        [HttpPost]
        public IActionResult Index(string filterClass, double Rp, double Rs, string WpStr, string WsStr, FilterTypes filterType, string curvedType = "", int order = 1)
        {
            int N = 100;
            string[] WpStrArr = WpStr.Split(' ');
            string[] WsStrArr = WsStr.Split(' ');
            double Wp = Convert.ToDouble(WpStrArr[0]);
            double Ws = Convert.ToDouble(WsStrArr[0]);
            WindowTypes winType = WindowTypes.Hamming;
            FilterBaseViewModel filterViewModel;


            // Логика создания и использования фильтра
            if (filterClass.Equals("IIR"))
            {
                IIRFilterBase IIRfilter = CreateIIRFilter(Rp, Rs, Ws, Wp, (FilterTypes)filterType, curvedType);
                if(IIRfilter.FilterType == FilterTypes.Bandpass || IIRfilter.FilterType == FilterTypes.Bandstop)
                {
                    IIRfilter.W1 = IIRfilter.GetCutoffFrenq(Convert.ToDouble(WpStrArr[1]), Convert.ToDouble(WsStrArr[1]));
                }
                IIRfilter.DesignFilter();
                MFC = IIRfilter.GetMAG(N, 0, Math.PI);
                PFC = IIRfilter.GetPHI(N, 0, Math.PI);

                filterViewModel = new IIRFilterViewModel();
                filterViewModel.Order = IIRfilter.Order;
                filterViewModel.Az = IIRfilter.Az;
                filterViewModel.Bz = IIRfilter.Bz;
                filterViewModel.MFC = MFC;
                filterViewModel.PFC = PFC;
                filterViewModel.W = IIRfilter.W;
                (filterViewModel as IIRFilterViewModel).As = IIRfilter.A;
                (filterViewModel as IIRFilterViewModel).Bs = IIRfilter.B;
                (filterViewModel as IIRFilterViewModel).Z = IIRfilter.Z;
                (filterViewModel as IIRFilterViewModel).P = IIRfilter.P;
            }
            else
            {
                FIRFilter FIRfilter = new FIRFilter(order, (FilterTypes)filterType, Rp, Rs, Wp, Ws);
                FIRfilter.DesignFilter(winType);
                MFC = FIRfilter.GetMAG(N, 0, Math.PI);
                PFC = FIRfilter.GetPHI(N, 0, Math.PI);
                filterViewModel = new FilterBaseViewModel(); ;
                filterViewModel.Order = FIRfilter.Order;
                filterViewModel.Bz = FIRfilter.H;
                filterViewModel.MFC = MFC;
                filterViewModel.PFC = PFC;
                filterViewModel.W = FIRfilter.W;
            }


            return View("Result", filterViewModel);
        }

        private IIRFilterBase CreateIIRFilter(double rp, double rs, double ws, double wp, FilterTypes filterType, string curvedType)
        {
            switch (curvedType)
            {
                case "Batterrfort":
                    return new ButterworthFilter(filterType, wp, ws, rp, rs);

                case "Cheby1":
                    return new ChebishevFilter1(filterType, wp, ws, rp, rs);

                case "Cheby2":
                    return new ChebishevFilter2(filterType, wp, ws, rp, rs);

                default:
                    return new EllipFilter(filterType, wp, ws, rp, rs);
            }
        }

        public IActionResult Author()
        {
            return View();
        }      
        
        public IActionResult About()
        {
            return View();
        }

        public ActionResult _GetIIRResult(IIRFilterViewModel model)
        {
            return PartialView("_GetIIRResult", model);
        }

        public ActionResult _GetFIRResult(FilterBaseViewModel model)
        {
            return PartialView("_GetFIRResult", model);
        }
    }
}
