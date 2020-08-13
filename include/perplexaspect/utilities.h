/*
 * Copyright (C) 2011 - 2020 by the authors of the ASPECT code.
 * Copyright (C) 2020 Connor Ward
 *
 * This file is part of PerpleX-ASPECT.
 * PerpleX-ASPECT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PerpleX-ASPECT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with PerpleX-ASPECT. If not, see <https://www.gnu.org/licenses/>.
 */


#ifndef _perplexaspect_utilities_h
#define _perplexaspect_utilities_h


#include <deal.II/base/parameter_handler.h>
#include <perplexcpp/base.h>


namespace aspect
{
  namespace PerplexUtils
  {
    using namespace dealii;


    /**
     * TODO
     */
    void
    declare_parameters(ParameterHandler &prm);


    /**
     * TODO
     */
    void
    parse_parameters(ParameterHandler &prm);


    /**
     * Given the porosity (melt volume fraction) return the melt composition
     * scaled to match it.
     */
    void
    put_melt_composition(const perplexcpp::MinimizeResult &result,
	                 const double porosity,
	                 std::vector<double> &melt_composition);


    /**
     * Restrict the pressure to lie within the bounds accepted by Perple_X.
     */
    double
    limit_pressure(const double pressure);


    /**
     * Restrict the temperature to lie within the bounds accepted by Perple_X.
     */
    double
    limit_temperature(const double temperature);
  }
}


#endif
