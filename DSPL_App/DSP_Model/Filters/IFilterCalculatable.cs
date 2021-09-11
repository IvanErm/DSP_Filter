using System;
using System.Collections.Generic;
using System.Text;

namespace DSP_Model.Filters
{
    interface IFilterCalculatable
    {
        void Calculate_AP();
        void Calculate_AZ();
        int CalculateOrder();

        void DesignFilter();

    }
}
