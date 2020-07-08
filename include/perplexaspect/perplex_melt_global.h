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


#include <aspect/material_model/melt_global.h>


namespace aspect
{
  namespace MaterialModel
  {

    /**
     * ???
     */
    template <int dim>
    class PerpleXMelt : public MeltGlobal<dim>
    {
      public:

	void
	initialize() override;

	/**
	 * ???
	 */
	void
	evaluate(const MaterialModelInputs<dim> &in, 
		 MaterialModelOutputs<dim> &out) const override;


	/**
	 * ???
	 */
	static
	void
	declare_parameters(ParameterHandler &prm);


	/**
	 * ???
	 */
	void
	parse_parameters(ParameterHandler &prm) override;


      protected:

	std::vector<double>
	make_composition(const std::vector<double>& compositional_fields) const;


	/**
	 * ???
	 */
	std::vector<double>
	get_composition(const std::vector<double>& comp_fields, 
	                const std::string& name) const;

	double melting_time_scale;
	bool include_melting_and_freezing;
    };
  }
}

