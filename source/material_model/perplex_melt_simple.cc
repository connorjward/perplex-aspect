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
	  double pressure = in.pressure[q];
	  double temperature = in.temperature[q];

	  // 200K is added to the temperature to make the results more consistent
	  // between Perple_X and ASPECT. This is a hack. In order to avoid this
	  // the porosity value should really be calculated from Perple_X too.
	  temperature += 200;

	  if (pressure < px.min_pressure)
	    pressure = px.min_pressure;
	  else if (pressure > px.max_pressure)
	    pressure = px.max_pressure;

	  if (temperature < px.min_temperature)
	    temperature = px.min_temperature;
	  if (temperature > px.max_temperature)
	    temperature = px.max_temperature;

	  std::vector<double> composition(px.n_composition_components);
	  this->load_perplex_composition_from_fields(in.composition[q], composition);

	  const auto result = px.minimize(pressure, temperature, composition);

	  std::vector<double> melt_composition(px.n_composition_components);
	  this->get_melt_composition(in.composition[q], result, melt_composition);

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
    get_melt_composition(const std::vector<double> &initial_composition,
			 const perplexcpp::MinimizeResult &result,
	                 std::vector<double> &melt_composition) const
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      const unsigned int porosity_idx = 
	this->introspection().compositional_index_for_name("porosity");
      const double porosity = std::max(initial_composition[porosity_idx], 0.0);

      const perplexcpp::Phase melt = perplexcpp::find_phase(result.phases, "liquid");

      // If no melt is present return a vector of zeros.
      if (porosity < 1e-8 || melt.n_moles < 1e-8)
      {
	for (unsigned int c = 0; c < px.n_composition_components; c++)
	  melt_composition[c] = 0.0;
	return;
      }

      // Initial calculations for readability later on.
      double phase_densities_sum = 0.0;
      for (auto p : result.phases)
	phase_densities_sum += p.density;

      double melt_mol_mass = 0.0;
      for (unsigned int c = 0; c < px.n_composition_components; ++c)
	melt_mol_mass += melt.composition_ratio[c] * px.composition_molar_masses[c];

      double melt_comp_ratio_sum = 0.0;
      for (double melt_comp_ratio : melt.composition_ratio)
	melt_comp_ratio_sum += melt_comp_ratio;

      double bulk_weight = 0.0;
      for (unsigned int c = 0; c < px.n_composition_components; ++c)
	bulk_weight += result.composition[c] * px.composition_molar_masses[c];

      const double melt_weight_frac = 
	porosity * melt.density / phase_densities_sum;

      // Populate the melt composition vector.
      for (unsigned int c = 0; c < px.n_composition_components; c++)
      {
	const double cmelt_weight_frac = 
	  melt.composition_ratio[c] * px.composition_molar_masses[c] / melt_mol_mass;

	// The weight fraction of the composition component in the melt is the weight fraction 
	// of the melt multiplied by the weight fraction of the component inside the melt.
	const double cbulk_weight_frac = melt_weight_frac * cmelt_weight_frac;

	// Turn the weight fraction into an actual weight.
	const double cbulk_weight = cbulk_weight_frac * bulk_weight;

	// Figure out the number of moles of the component.
	const double cbulk_mol = cbulk_weight / px.composition_molar_masses[c];

	melt_composition[c] = cbulk_mol;
      }
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
