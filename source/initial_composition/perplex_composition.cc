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
#include <perplexaspect/material_model/perplex_melt.h>
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
	          <const MaterialModel::MeltFractionModel<dim>>(this->get_material_model()),
	          ExcMessage("The material model is not derived from the 'MeltFractionModel' class, "
	                     "and therefore does not support computing equilibrium melt fractions. "
	                     "This is incompatible with the `perplex composition' "
	                     "initial composition plugin, which needs to compute these melt fractions."));

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
          double pressure = this->get_adiabatic_conditions().pressure(position);
          double temperature = this->get_initial_temperature_manager().initial_temperature(position);

	  if (pressure < px.min_pressure)
	    pressure = px.min_pressure;
	  if (pressure > px.max_pressure)
	    pressure = px.max_pressure;
	  if (temperature < px.min_temperature)
	    temperature = px.min_temperature;
	  if (temperature > px.max_temperature)
	    temperature = px.max_temperature;

	  // Use the bulk composition from the file.
	  const MinimizeResult result { px.minimize(pressure, temperature) };
	  const Phase melt = find_phase(result.phases, "liquid");

	  // Get the composition of the melt.
	  std::vector<double> melt_composition(px.n_composition_components);
	  {
	    const unsigned int porosity_idx = 
	      this->introspection().compositional_index_for_name("porosity");
	    const double porosity = 
	      this->get_initial_composition_manager().initial_composition(position, porosity_idx);
	    PerplexUtils::put_melt_composition(result, porosity, melt_composition);
	  }

	  if (compositional_index == melt_idx)
	    return melt_composition[c];
	  else
	    return std::max(result.composition[c] - melt_composition[c], 0.0);
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
      "???"
    )
  }
}

