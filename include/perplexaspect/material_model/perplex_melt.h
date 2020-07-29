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

#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

#include <perplexcpp/base.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * ???
     */
    template <int dim>
    class PerplexMelt : public MaterialModel::MeltInterface<dim>, 
	                public ::aspect::SimulatorAccess<dim>, 
			public MaterialModel::MeltFractionModel<dim>
    {
      public:

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;


        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;


        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;


        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;


        /**
         * @name Reference quantities
         * @{
         */
        double reference_viscosity () const override;


        double reference_darcy_coefficient () const override;


        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;


        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        double reference_rho_s;
        double reference_rho_f;
        double reference_T;
        double eta_0;
        double xi_0;
        double eta_f;
        double thermal_viscosity_exponent;
        double thermal_bulk_viscosity_exponent;
        double thermal_expansivity;
        double reference_specific_heat;
        double thermal_conductivity;
        double reference_permeability;
        double alpha_phi;
        double extraction_depth;
        double compressibility;
        double melt_compressibility;
        double melt_bulk_modulus_derivative;
        double depletion_density_change;
        double depletion_solidus_change;
        bool model_is_compressible;
        bool fractional_melting;
        double freezing_rate;
        double melting_time_scale;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        // entropy change upon melting
        double peridotite_melting_entropy_change;

        /**
         * Percentage of material that is molten for a given @p temperature and
         * @p pressure (assuming equilibrium conditions). Melting model after Katz,
         * 2003, for dry peridotite.
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;

        /**
         * Compute the change in entropy due to melting for a given @p temperature
         * and @p pressure, and under the assumption that a fraction
         * @p maximum_melt_fraction of the material has already been molten
         * previously. The entropy change is computed with respect to temperature
         * or pressure, depending on @p dependence.
         * This is needed to calculate the latent heat of melt.
         */
        virtual
        double
        entropy_change (const double temperature,
                        const double pressure,
                        const double maximum_melt_fraction,
                        const NonlinearDependence::Dependence dependence) const;

	/**
	 * ADDED
	 * Retrieve the composition from the input compositional fields given
	 * the phase name.
	 */
	std::vector<double>
        get_composition(const std::vector<double>& comp_fields, 
	                              const std::string& name) const;


	/**
	 * ADDED
	 * Calculate the melt composition (number of moles of each component) 
	 * from the porosity (melt volume fraction) and melt composition ratio.
	 */
	std::vector<double> 
	calc_melt_composition(const double porosity, const perplexcpp::MinimizeResult &result)
	const;
    };

  }
}


#endif
