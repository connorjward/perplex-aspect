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
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    PerplexMelt<dim>::initialize()
    {
      AssertThrow(this->include_melt_transport(),
	          ExcMessage("Material model 'Perple_X melt' on works if melt "
		             "transport is enabled."));
      AssertThrow(this->introspection().compositional_name_exists("porosity"),
		  ExcMessage("Material model 'Perple_X melt' only works if "
		             "there is a compositional field called 'porosity'."));
    }



    template <int dim>
    bool
    PerplexMelt<dim>::is_compressible() const
    {
      return this->equation_of_state.is_compressible();
    }



    template <int dim>
    double
    PerplexMelt<dim>::reference_viscosity() const
    {
      return this->constant_rheology.compute_viscosity();
    }



    template <int dim>
    void
    PerplexMelt<dim>::evaluate(const MaterialModelInputs<dim> &in,
                               MaterialModelOutputs<dim> &out) const
    {
      EquationOfStateOutputs<dim> eos_outputs(1);

      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
      {
	equation_of_state.evaluate(in, i, eos_outputs);

	out.viscosities[i] = constant_rheology.compute_viscosity();
	out.densities[i] = eos_outputs.densities[0];
	out.thermal_expansion_coefficients[i] = eos_outputs.thermal_expansion_coefficients[0];
	out.specific_heat[i] = eos_outputs.specific_heat_capacities[0];
	out.thermal_conductivities[i] = k_value;
	out.compressibilities[i] = eos_outputs.compressibilities[0];

	// Added to stop it breaking (same as melt_global).
	out.entropy_derivative_pressure[i] = 0.0;
	out.entropy_derivative_temperature[i] = 0.0;

	for (unsigned int c=0; c<in.composition[i].size(); ++c)
	  out.reaction_terms[i][c] = 0.0;

	    /* // start ADDED section */
	    /* auto& perplex_wrapper = perplexcpp::Wrapper::get_instance(); */

	    /* std::vector<double> composition; */
	    /* { */
	      /* auto melt_composition = this->get_composition(in.composition[i], "melt"); */
	      /* auto residue_composition = this->get_composition(in.composition[i], "residue"); */

	      /* for (unsigned int c = 0; */ 
		   /* c < perplexcpp::Wrapper::get_instance().n_composition_components; */ 
		   /* ++c) */
	      /* { */
		/* composition.push_back(melt_composition[c] + residue_composition[c]); */
	      /* } */
	    /* } */

	    /* double pressure = in.pressure[i]; */
	    /* double temperature = in.temperature[i]; */

	    /* // 200K is added to the temperature to make the results more consistent */
	    /* // between Perple_X and ASPECT. This is a hack. In order to avoid this */
	    /* // the porosity value should really be calculated from Perple_X too. */
	    /* temperature += 200; */

	    /* if (pressure < perplex_wrapper.min_pressure) */
	      /* pressure = perplex_wrapper.min_pressure; */
	    /* else if (pressure > perplex_wrapper.max_pressure) */
	      /* pressure = perplex_wrapper.max_pressure; */

	    /* if (temperature < perplex_wrapper.min_temperature) */
	      /* temperature = perplex_wrapper.min_temperature; */
	    /* if (temperature > perplex_wrapper.max_temperature) */
	      /* temperature = perplex_wrapper.max_temperature; */

	    /* const auto result = perplex_wrapper.minimize(pressure, */ 
							 /* temperature, */
							 /* composition); */
	    
	    /* // Make sure that the melt amount is non-negative. */
	    /* const double melt_frac = old_porosity > 0 ? old_porosity : 0.0; */

	    /* // Get the melt composition. */
	    /* const std::vector<double> melt_composition = this->calc_melt_composition(melt_frac, result); */

	    /* // Determine the composition changes and alter the reaction rates accordingly. */
	    /* for (unsigned int c = 0; c < perplex_wrapper.n_composition_components; ++c) */ 
	    /* { */
	      /* const std::string comp_name = perplex_wrapper.composition_component_names[c]; */

	      /* const unsigned int melt_comp_idx = */ 
		/* this->introspection().compositional_index_for_name("melt_" + comp_name); */

	      /* const unsigned int residue_comp_idx = */ 
		/* this->introspection().compositional_index_for_name("residue_" + comp_name); */

	      /* const double melt_comp_change = */ 
		/* melt_composition[c] - in.composition[i][melt_comp_idx]; */

	      /* const double residue_comp_change = */ 
		/* result.composition[c] - melt_composition[c] */
		/* - in.composition[i][residue_comp_idx]; */

	      /* // Populate the reaction rates. */
	      /* if (reaction_rate_out != nullptr && this->get_timestep_number() > 0) */
	      /* { */
		/* reaction_rate_out->reaction_rates[i][melt_comp_idx] = */ 
		  /* melt_comp_change / this->melting_time_scale; */

		/* reaction_rate_out->reaction_rates[i][residue_comp_idx] = */ 
		  /* residue_comp_change / this->melting_time_scale; */
	      /* } */
	    /* } */
	    /* // end ADDED section */

	// Fill melt outputs if they exist.
	MeltOutputs<dim> *melt_out = 
	  out.template get_additional_output<MeltOutputs<dim>>();

	if (melt_out != nullptr)
	  this->fill_melt_outputs(in, melt_out);
      }
    }



    template <int dim>
    void
    PerplexMelt<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Perple_X melt");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters(prm);

          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in the density formula. Units: $\\si{K}$.");
          
	  prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");

          Rheology::ConstantViscosity::declare_parameters(prm, 5e24);

	  //meltstuff
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0.),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference melt viscosity", "10.",
                             Patterns::Double (0.),
                             "The value of the constant melt viscosity $\\eta_f$. "
			     "Units: $Pa \\, s$.");

          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");

          prm.declare_entry ("Reference melt viscosity", "10.",
                             Patterns::Double (0.),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa \\, s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27.",
                             Patterns::Double (0.),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Melt compressibility", "0.0",
                             Patterns::Double (0.),
                             "The value of the compressibility of the melt. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Reference melt density", "2500.",
                             Patterns::Double (0.),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. Units: $kg/m^3$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    PerplexMelt<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Perple_X melt");
        {
          equation_of_state.parse_parameters (prm);

          this->reference_T                = prm.get_double ("Reference temperature");
          this->k_value                    = prm.get_double ("Thermal conductivity");

          constant_rheology.parse_parameters(prm);

	  this->xi_0 = prm.get_double ("Reference bulk viscosity");
          this->eta_f                      = prm.get_double ("Reference melt viscosity");
          this->reference_permeability     = prm.get_double ("Reference permeability");
	  this->thermal_expansivity        = prm.get_double ("Thermal expansion coefficient");
	  this->melt_compressibility = prm.get_double ("Melt compressibility");
	  this->alpha_phi = prm.get_double ("Exponential melt weakening factor");
	  this->reference_rho_f = prm.get_double ("Reference melt density");
	  this->thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
	  this->thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables.
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }



    template <int dim>
    void 
    PerplexMelt<dim>::
    melt_fractions(const MaterialModelInputs<dim> &in,
		   std::vector<double> &melt_fractions) const
    {
      // TODO: Implement this properly.
      melt_fractions = std::vector<double>(in.n_evaluation_points(), 0.0);
    }



    template <int dim>
    double
    PerplexMelt<dim>::reference_darcy_coefficient() const
    {
      // The same as melt_simple.cc and melt_global.cc.
      return this->reference_permeability * std::pow(0.01,3.0) / this->eta_f;
    }



    // Note: This function is semantically identical to the relevant part of
    // evaluate() in melt_global.cc.
    template <int dim>
    void
    PerplexMelt<dim>::
    fill_melt_outputs(const MaterialModelInputs<dim> &in, 
	              MeltOutputs<dim> *melt_out) const
    {
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
      {
	double porosity = std::max(in.composition[i][porosity_idx],0.0);

	melt_out->fluid_viscosities[i] = this->eta_f;
	melt_out->permeabilities[i] = 
	  this->reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);
	melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

	// temperature dependence of density is 1 - alpha * (T - T(adiabatic))
	double temperature_dependence = 1.0;
	if (this->include_adiabatic_heating ())
	  temperature_dependence -= 
	    (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i])) 
	    * this->thermal_expansivity;
	else
	  temperature_dependence -= 
	    (in.temperature[i] - reference_T) * thermal_expansivity;

	melt_out->fluid_densities[i] = 
	  reference_rho_f * temperature_dependence
	  * std::exp(melt_compressibility * (in.pressure[i] - this->get_surface_pressure()));

	melt_out->compaction_viscosities[i] = this->xi_0 * exp(-this->alpha_phi * porosity);

	double visc_temperature_dependence = 1.0;
	if (this->include_adiabatic_heating ())
	{
	  const double delta_temp = 
	    in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
	  visc_temperature_dependence = 
	    std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp 
		                       / this->get_adiabatic_conditions().temperature(in.position[i])),
		  1e4),
		1e-4);
	}
	else if (thermal_viscosity_exponent != 0.0)
	{
	  const double delta_temp = in.temperature[i]-reference_T;
	  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
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
    ASPECT_REGISTER_MATERIAL_MODEL(PerplexMelt,
                                   "perplex melt",
                                   "A material model that has constant values "
                                   "except for density, which depends linearly on temperature: "
                                   "\\begin{align}"
                                   "  \\rho(p,T) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0."
                                   "\\end{align}"
                                   "\n\n"
                                   "\\note{This material model fills the role the ``simple'' material "
                                   "model was originally intended to fill, before the latter acquired "
                                   "all sorts of complicated temperature and compositional dependencies.}")
  }
}