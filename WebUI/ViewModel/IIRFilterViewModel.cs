using System.Collections.Generic;
using System.Numerics;

namespace WebUI.ViewModel
{
    public class IIRFilterViewModel : FilterBaseViewModel
    {
        public List<Complex> Z { get; set; }
        public List<Complex> P { get; set; }

        public List<double> As { get; set; }
        public List<double> Bs { get; set; }
    }
}
