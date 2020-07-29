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


#ifndef _perplexaspect_material_model_perplex_melt_h
#define _perplexaspect_material_model_perplex_melt_h


#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/constant_viscosity.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/melt.h>
#include <aspect/simulator_access.h>
#include<perplexcpp/wrapper.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class PerplexMelt : public MeltInterface<dim>, 
                        public MeltFractionModel<dim>,
	                public ::aspect::SimulatorAccess<dim>
    {
      public:

	void 
	initialize() override;


	bool 
	is_compressible() const override;


        double 
	reference_viscosity() const override;


	void 
	evaluate(const MaterialModelInputs<dim> &in,
                 MaterialModelOutputs<dim> &out) const override;


	static
	void
	declare_parameters(ParameterHandler &prm);


	void
	parse_parameters(ParameterHandler &prm) override;


	void 
	melt_fractions(const MaterialModelInputs<dim> &in,
                       std::vector<double> &melt_fractions) const override;


        double 
	reference_darcy_coefficient() const override;


      private:

        double reference_T;
        double k_value;
        Rheology::ConstantViscosity constant_rheology;
        EquationOfState::LinearizedIncompressible<dim> equation_of_state;

	// melt stuff
	double eta_f;
	double reference_permeability;
	double thermal_expansivity;
	double melt_compressibility;
	double alpha_phi;
	double reference_rho_f;
	double thermal_bulk_viscosity_exponent;
	double thermal_viscosity_exponent;
	double xi_0;


	std::vector<double>
	get_bulk_composition(const MaterialModelInputs<dim> &in,
	                     const unsigned int q) const;

	std::vector<double>
    get_composition(const std::vector<double>& comp_fields, 
	            const std::string& name) const;

	std::vector<double>
	calc_melt_composition(const double porosity, 
	                      const perplexcpp::MinimizeResult &result) const;


	void
	fill_melt_outputs(const MaterialModelInputs<dim> &in, 
			  MeltOutputs<dim> *melt_out) const;


    };
  }
}


#endif
