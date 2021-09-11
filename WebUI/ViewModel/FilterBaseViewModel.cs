using System.Collections.Generic;

namespace WebUI.ViewModel
{
    public class FilterBaseViewModel
    {
        private ChartViewModel mfc_chart = null;
        private ChartViewModel ffc_chart = null;

        public int Order { get; set; }
        public List<double> Az { get; set; }
        public List<double> Bz { get; set; }
        public List<double> W { get; set; }
        public List<double> MFC { get; set; }  // АЧХ
        public List<double> PFC { get; set; }  // ФЧХ

        public ChartViewModel MFC_Chart
        {
            get
            {
                if (mfc_chart == null)
                {
                    mfc_chart = new ChartViewModel(W, MFC);
                }
                return mfc_chart;
            }

        }
        public ChartViewModel PFC_Chart
        {
            get
            {
                if(ffc_chart == null)
                {
                    ffc_chart = new ChartViewModel(W, PFC);
                }

                return ffc_chart;
            }
        }
    } 
}
