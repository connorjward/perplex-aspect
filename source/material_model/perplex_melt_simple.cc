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


#include <perplexaspect/material_model/perplex_melt_simple.h>

#include <perplexaspect/utilities.h>
#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace MaterialModel
  { 
    template <int dim>
    void
    PerplexMeltSimple<dim>::evaluate(const MaterialModelInputs<dim> &in,
	                             MaterialModelOutputs<dim> &out) const
    {
      MeltSimple<dim>::evaluate(in, out);

      const auto& px = perplexcpp::Wrapper::get_instance();

      ReactionRateOutputs<dim> *reaction_rate_out = 
	out.template get_additional_output<ReactionRateOutputs<dim>>();

      if (reaction_rate_out != nullptr &&
	  this->get_timestep_number() > 0 &&
	  in.requests_property(MaterialProperties::reaction_terms))
      {
	for (unsigned int q = 0; q < in.n_evaluation_points(); q++)
	{
	  const double pressure = PerplexUtils::limit_pressure(in.pressure[q]);

	  // 200K is added to make sure that Perple_X reports melt being present 
	  // when the ASPECT parametrisation reports that there is.
	  const double temperature = PerplexUtils::limit_temperature(in.temperature[q]) + 200;

	  std::vector<double> composition(px.n_composition_components);
	  this->load_perplex_composition_from_fields(in.composition[q], composition);

	  const auto result = px.minimize(pressure, temperature, composition);

	  std::vector<double> melt_composition(px.n_composition_components);
	  {
	    const unsigned int porosity_idx = 
	      this->introspection().compositional_index_for_name("porosity");
	    const double porosity =
	      std::max(in.composition[q][porosity_idx], 0.0);
	    PerplexUtils::put_melt_composition(result, porosity, melt_composition);
	  }

	  this->update_reaction_rates(in.composition[q], 
	                              result.composition, 
				      melt_composition, 
				      reaction_rate_out->reaction_rates[q]);
	}
      }
    }


    template <int dim>
    void
    PerplexMeltSimple<dim>::declare_parameters(ParameterHandler &prm)
    {
      MeltSimple<dim>::declare_parameters(prm);
      PerplexUtils::declare_parameters(prm);
    }



    template <int dim>
    void
    PerplexMeltSimple<dim>::parse_parameters(ParameterHandler &prm)
    {
      MeltSimple<dim>::parse_parameters(prm);
      PerplexUtils::parse_parameters(prm);
    }


    template <int dim>
    void
    PerplexMeltSimple<dim>::
    load_perplex_composition_from_fields(const std::vector<double> &aspect_composition,
					 std::vector<double> &perplex_composition) const
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      Assert(perplex_composition.size() == px.n_composition_components,
	     ExcInternalError("The composition is the wrong size."));

      for (unsigned int c = 0; c < px.n_composition_components; c++)
      {
	const std::string cname = px.composition_component_names[c];

	const unsigned int cmelt_idx = 
	  this->introspection().compositional_index_for_name("melt_" + cname);
	const unsigned int cres_idx = 
	  this->introspection().compositional_index_for_name("residue_" + cname);

	perplex_composition[c] = 
	  aspect_composition[cmelt_idx] + aspect_composition[cres_idx];
      }
    }


    template <int dim>
    void
    PerplexMeltSimple<dim>::
    update_reaction_rates(const std::vector<double> &initial_composition,
	                  const std::vector<double> &final_composition,
			  const std::vector<double> &melt_composition,
	                  std::vector<double> &reaction_rates) const
    {
      const auto &px = perplexcpp::Wrapper::get_instance();

      for (unsigned int c = 0; c < px.n_composition_components; c++) 
      {
	const std::string cname = px.composition_component_names[c];
	const unsigned int cmelt_idx = 
	  this->introspection().compositional_index_for_name("melt_" + cname);
	const unsigned int cres_idx = 
	  this->introspection().compositional_index_for_name("residue_" + cname);

	reaction_rates[cmelt_idx] = 
	  (melt_composition[c] - initial_composition[cmelt_idx]) 
	  / this->get_parameters().reaction_time_step;

	reaction_rates[cres_idx] = 
	  (final_composition[c] - melt_composition[c] - initial_composition[cres_idx])
	  / this->get_parameters().reaction_time_step;
      }
    }
  }
}


namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(
      PerplexMeltSimple,
      "perplex melt simple",
      "Description here."
    )
  }
}
