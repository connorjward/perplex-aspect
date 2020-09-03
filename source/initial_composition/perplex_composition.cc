/*
  Copyright (C) 2017 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <perplexaspect/initial_composition/perplex_composition.h>

#include <aspect/initial_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <perplexaspect/material_model/perplex_melt_simple.h>
#include <perplexaspect/utilities.h>
#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace InitialComposition
  {
    using namespace perplexcpp;


    template <int dim>
    double
    PerplexComposition<dim>::
    initial_composition(const Point<dim> &position,
                        const unsigned int compositional_index) const
    {
      AssertThrow(Plugins::plugin_type_matches
	          <const MaterialModel::PerplexMeltSimple<dim>>(this->get_material_model()),
	          ExcMessage("The material model is not `perplex melt simple'."));

      const auto& px = Wrapper::get_instance();

      for (unsigned int c = 0; c < px.n_composition_components; c++) {
	const std::string cname = px.composition_component_names[c];

	AssertThrow(this->introspection().compositional_name_exists("melt_"+cname),
		    ExcMessage("The initial composition plugin `perplex composition' did not find a "
			       "compositional field called `melt_"+cname+"' to initialize. Please add a "
			       "compositional field with this name."));

	AssertThrow(this->introspection().compositional_name_exists("residue_"+cname),
		    ExcMessage("The initial composition plugin `perplex composition' did not find a "
			       "compositional field called `residue_"+cname+"' to initialize. Please add a "
			       "compositional field with this name."));

	const unsigned int melt_idx = 
	  this->introspection().compositional_index_for_name("melt_"+cname);
	const unsigned int residue_idx = 
	  this->introspection().compositional_index_for_name("residue_"+cname);


	if (compositional_index == melt_idx || compositional_index == residue_idx) {
          const double pressure = 
	    PerplexUtils::limit_pressure(this->get_adiabatic_conditions().pressure(position));
          const double temperature = 
	    PerplexUtils::
	    limit_temperature(this->get_initial_temperature_manager().initial_temperature(position));

	  // Use the bulk composition from the file.
	  const MinimizeResult result = px.minimize(pressure, temperature);

	  const Phase melt = find_phase(result.phases, "liquid");

	  if (compositional_index == melt_idx)
	    return melt.composition_ratio[c] * melt.n_moles;
	  else
	    return std::max(result.composition[c] - melt.composition_ratio[c]*melt.n_moles, 0.0);
        }
      }
      return 0.0;
    }
  }
}


// Explicit instantiations.
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(
      PerplexComposition,
      "perplex composition",
      "Initial composition plugin for the `perplex melt simple' material model."
    )
  }
}

