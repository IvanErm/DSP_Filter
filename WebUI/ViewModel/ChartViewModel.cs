using System.Collections.Generic;

namespace WebUI.ViewModel
{
    public class ChartViewModel
    {
        public List<double> XAxial { get; set; }
        public List<double> YAxial { get; set; }

        public ChartViewModel(List<double> x, List<double> y)
        {
            XAxial = x;
            YAxial = y;
        }
    }
}
