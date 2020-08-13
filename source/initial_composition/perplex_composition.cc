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
          const double pressure = 
	    PerplexUtils::limit_pressure(this->get_adiabatic_conditions().pressure(position));
          const double temperature = 
	    PerplexUtils::
	    limit_temperature(this->get_initial_temperature_manager().initial_temperature(position));

	  // Use the bulk composition from the file.
	  const MinimizeResult result { px.minimize(pressure, temperature) };
	  
	  // Get the composition of the melt.
	  const double porosity = this->calc_porosity(position);
	  std::vector<double> melt_composition(px.n_composition_components);
	  PerplexUtils::put_melt_composition(result, porosity, melt_composition);

	  if (compositional_index == melt_idx)
	    return melt_composition[c];
	  else
	    return std::max(result.composition[c] - melt_composition[c], 0.0);
        }
      }
      return 0.0;
    }


    template <int dim>
    double
    PerplexComposition<dim>::
    calc_porosity(const Point<dim> position) const
    {
      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());

      in.position[0] = position;
      in.temperature[0] = this->get_initial_temperature_manager().initial_temperature(position);
      in.pressure[0] = this->get_adiabatic_conditions().pressure(position);
      in.pressure_gradient[0] = 0.0;
      in.velocity[0] = 0.0;

      for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
	in.composition[0][i] = 0.0;
      
      in.strain_rate[0] = SymmetricTensor<2,dim>();

      std::vector<double> melt_fraction(1);
      
      Plugins::get_plugin_as_type<const MaterialModel::MeltFractionModel<dim>>
	(this->get_material_model()).melt_fractions(in, melt_fraction);

      return melt_fraction[0];
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

