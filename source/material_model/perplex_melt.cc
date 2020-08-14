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


#include <perplexaspect/material_model/perplex_melt.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <perplexaspect/utilities.h>
#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace MaterialModel
  { 
    template <int dim>
    double
    PerplexMelt<dim>::reference_darcy_coefficient() const
    {
      return this->permeability * std::pow(0.01,3.0) / this->fluid_viscosity;
    }


    template <int dim>
    bool
    PerplexMelt<dim>::is_compressible() const
    {
      return this->base_model->is_compressible();
    }


    template <int dim>
    double
    PerplexMelt<dim>::reference_viscosity() const
    {
      return this->base_model->reference_viscosity();
    }


    template <int dim>
    void
    PerplexMelt<dim>::
    melt_fractions(const MaterialModelInputs<dim> &in, 
	           std::vector<double> &melt_fractions) const
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      for (unsigned int q = 0; q < in.n_evaluation_points(); q++) {
	const double pressure = PerplexUtils::limit_pressure(in.pressure[q]);
	const double temperature = PerplexUtils::limit_temperature(in.temperature[q]);

	std::vector<double> composition(px.n_composition_components);
	this->load_perplex_composition_from_fields(in.composition[q], composition);

	const perplexcpp::MinimizeResult result = 
	  px.minimize(pressure, temperature, composition);

	const perplexcpp::Phase melt = 
	  perplexcpp::find_phase(result.phases, "liquid");

	melt_fractions[q] = melt.volume_frac;
      }
    }


    template <int dim>
    void
    PerplexMelt<dim>::evaluate(const MaterialModelInputs<dim> &in,
	                             MaterialModelOutputs<dim> &out) const
    {
      this->base_model->evaluate(in, out);

      const auto& px = perplexcpp::Wrapper::get_instance();

      ReactionRateOutputs<dim> *reaction_rate_out = 
	out.template get_additional_output<ReactionRateOutputs<dim>>();

      if (reaction_rate_out != nullptr &&
	  this->get_timestep_number() > 0 &&
	  in.requests_property(MaterialProperties::reaction_terms))
      {
	for (unsigned int q = 0; q < in.n_evaluation_points(); q++) {
	  const double pressure = PerplexUtils::limit_pressure(in.pressure[q]);
	  const double temperature = PerplexUtils::limit_temperature(in.temperature[q]);

	  std::vector<double> composition(px.n_composition_components);
	  this->load_perplex_composition_from_fields(in.composition[q], composition);

	  const perplexcpp::MinimizeResult result = 
	    px.minimize(pressure, temperature, composition);

	  this->put_reaction_rates(in.composition[q], result, reaction_rate_out->reaction_rates[q]);
	}
      }

      // Fill melt outputs if they exist.
      MeltOutputs<dim> *melt_out = 
	out.template get_additional_output<MeltOutputs<dim>>();

      if (melt_out != nullptr)
	this->fill_melt_outputs(in, melt_out);
    }


    template <int dim>
    void
    PerplexMelt<dim>::declare_parameters(ParameterHandler &prm)
    {
      PerplexUtils::declare_parameters(prm);
    }


    template <int dim>
    void
    PerplexMelt<dim>::parse_parameters(ParameterHandler &prm)
    {
      /* base_model.reset(create_material_model<dim>(prm.get("Base model"))); */
      base_model.reset(create_material_model<dim>("simple"));
      if (Plugins::plugin_type_matches<SimulatorAccess<dim>>(*base_model))
	Plugins::
	  get_plugin_as_type<SimulatorAccess<dim>>(*base_model)
	  .initialize_simulator(this->get_simulator());

      PerplexUtils::parse_parameters(prm);

      this->base_model->parse_parameters(prm);
      this->model_dependence = this->base_model->get_model_dependence();
    }


    template <int dim>
    void
    PerplexMelt<dim>::
    create_additional_named_outputs(MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting && out.template get_additional_output<ReactionRateOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }


    template <int dim>
    void
    PerplexMelt<dim>::
    load_perplex_composition_from_fields(const std::vector<double> &aspect_composition,
	                                 std::vector<double> &perplex_composition) const
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      Assert(perplex_composition.size() == px.n_composition_components,
	     ExcInternalError("The composition is the wrong size."));

      for (unsigned int c = 0; c < px.n_composition_components; c++) {
	const std::string cname = px.composition_component_names[c];

	const unsigned int cmelt_idx = 
	  this->introspection().compositional_index_for_name("melt_" + cname);
	const unsigned int cres_idx = 
	  this->introspection().compositional_index_for_name("residue_" + cname);

	perplex_composition[c] =
	  std::max(aspect_composition[cmelt_idx], 0.0)
	  + std::max(aspect_composition[cres_idx], 0.0);
      }
    }


    template <int dim>
    void
    PerplexMelt<dim>::
    put_reaction_rates(const std::vector<double> &initial_composition,
	               const perplexcpp::MinimizeResult &result,
	               std::vector<double> &reaction_rates) const
    {
      const auto &px = perplexcpp::Wrapper::get_instance();

      const perplexcpp::Phase melt =
	perplexcpp::find_phase(result.phases, "liquid");

      // Store porosity.
      const unsigned int porosity_idx = 
	this->introspection().compositional_index_for_name("porosity");
      reaction_rates[porosity_idx] = 
	(melt.volume_frac - initial_composition[porosity_idx]) 
	/ this->get_parameters().reaction_time_step;

      // Store composition.
      for (unsigned int c = 0; c < px.n_composition_components; c++) {
	const std::string cname = px.composition_component_names[c];
	const unsigned int cmelt_idx = 
	  this->introspection().compositional_index_for_name("melt_"+cname);
	const unsigned int cres_idx = 
	  this->introspection().compositional_index_for_name("residue_"+cname);

	reaction_rates[cmelt_idx] = 
	  (melt.composition_ratio[c]*melt.n_moles - initial_composition[cmelt_idx]) 
	  / this->get_parameters().reaction_time_step;

	reaction_rates[cres_idx] = 
	  (result.composition[c] - melt.composition_ratio[c]*melt.n_moles - initial_composition[cres_idx])
	  / this->get_parameters().reaction_time_step;
      }
    }


    template <int dim>
    void
    PerplexMelt<dim>::
    fill_melt_outputs(const MaterialModelInputs<dim> &in, 
	              MeltOutputs<dim> *melt_out) const
    {
      /* for (unsigned int q = 0; q < in.n_evaluation_points(); q++) { */
	/* melt_out->fluid_viscosities[q] = this->fluid_viscosity; */
	/* melt_out->permeabilities[q] = this->permeability; */
	/* melt_out->fluid_density_gradients[q] = Tensor<1,dim>(); */
	/* melt_out->fluid_densities[q] = this->fluid_density; */
	/* melt_out->compaction_viscosities[q] = this->compaction_viscosity; */
      /* } */


      double reference_permeability = 1e-8;
      double eta_f = 10;
      double reference_rho_f = 2500;
      double thermal_expansivity = 2e-5;
      double xi_0 = 1e22;
      double melt_bulk_modulus_derivative = 0.0;
      double reference_T = 293;
      double melt_compressibility = 0.0;
      double thermal_bulk_viscosity_exponent = 0.0;



      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
      {
	double porosity = std::max(in.composition[i][porosity_idx],0.0);

	melt_out->fluid_viscosities[i] = eta_f;
	melt_out->permeabilities[i] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);

	// first, calculate temperature dependence of density
	double temperature_dependence = 1.0;
	if (this->include_adiabatic_heating ())
	{
	  // temperature dependence is 1 - alpha * (T - T(adiabatic))
	  temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
	    * thermal_expansivity;
	}
	else
	  temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

	// the fluid compressibility includes two parts, a constant compressibility, and a pressure-dependent one
	// this is a simplified formulation, experimental data are often fit to the Birch-Murnaghan equation of state
	const double fluid_compressibility = melt_compressibility / (1.0 + in.pressure[i] * melt_bulk_modulus_derivative * melt_compressibility);

	melt_out->fluid_densities[i] = reference_rho_f * std::exp(fluid_compressibility * (in.pressure[i] - this->get_surface_pressure()))
	  * temperature_dependence;

	melt_out->fluid_density_gradients[i] = melt_out->fluid_densities[i] * melt_out->fluid_densities[i]
	  * fluid_compressibility
	  * this->get_gravity_model().gravity_vector(in.position[i]);

	const double phi_0 = 0.05;
	porosity = std::max(std::min(porosity,0.995),1e-4);
	melt_out->compaction_viscosities[i] = xi_0 * phi_0 / porosity;

	double visc_temperature_dependence = 1.0;
	if (this->include_adiabatic_heating ())
	{
	  const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
	  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
	}
	else
	{
	  const double delta_temp = in.temperature[i]-reference_T;
	  const double T_dependence = (thermal_bulk_viscosity_exponent == 0.0
	      ?
	      0.0
	      :
	      thermal_bulk_viscosity_exponent*delta_temp/reference_T);
	  visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence),1e4),1e-4);
	}
	melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
      }
    }
  }
}


// Explicit instantiations.
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(
      PerplexMelt,
      "perplex melt",
      "Description here."
    )
  }
}
