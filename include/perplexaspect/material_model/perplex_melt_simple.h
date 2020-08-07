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


#ifndef _perplexaspect_material_model_perplex_melt_simple_h
#define _perplexaspect_material_model_perplex_melt_simple_h


#include <aspect/material_model/melt_simple.h>
#include <perplexcpp/base.h>


namespace aspect
{
  namespace MaterialModel
  {
    /**
     * TODO
     */
    template <int dim>
    class PerplexMeltSimple : public MeltSimple<dim>
    {
      public:

        void evaluate(const MaterialModelInputs<dim> &in,
                      MaterialModelOutputs<dim> &out) const override;


        static
        void
        declare_parameters (ParameterHandler &prm);


        void
        parse_parameters (ParameterHandler &prm) override;

      private:

	void
	get_melt_composition(const std::vector<double> &initial_composition,
	                     const perplexcpp::MinimizeResult &result,
	                     std::vector<double> &melt_composition) const;


	void
	load_perplex_composition_from_fields(const std::vector<double> &aspect_composition,
					     std::vector<double> &perplex_composition) const;


	void
	update_reaction_rates(const std::vector<double> &initial_composition,
			      const std::vector<double> &final_composition,
			      const std::vector<double> &melt_composition,
			      std::vector<double> &reaction_rates) const;
    };
  }
}


#endif
