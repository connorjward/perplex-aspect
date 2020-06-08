/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_visco_plastic_h
#define _aspect_material_model_visco_plastic_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/strain_dependent.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/material_model/rheology/constant_viscosity_prefactors.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
#include <aspect/material_model/rheology/elasticity.h>

#include<deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for the plastic parameters weakened (or hardened)
     * by strain to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class PlasticAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        PlasticAdditionalOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Cohesions at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> cohesions;

        /**
         * Internal angles of friction at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> friction_angles;

        /**
         * The area where the viscous stress exceeds the plastic yield stress,
         * and viscosity is rescaled back to the yield envelope.
         */
        std::vector<double> yielding;
    };

    /**
     * A material model combining viscous and plastic deformation, with
     * the option to also include viscoelastic deformation.
     *
     * Viscous deformation is defined by a viscous flow law describing
     * dislocation and diffusion creep:
     *   $ v = \frac{1}{2}   A^{-\frac{1}{n}} d^{\frac{m}{n}}
     *               \dot{\varepsilon}_{ii}^{\frac{1-n}{n}}
     *               \exp\left(\frac{E + PV}{nRT}\right) $
     * where
     *   where $A$ is the prefactor, $n$ is the stress exponent,
     *   $\dot{\varepsilon}_{ii}$ is the square root of the deviatoric
     *   strain rate tensor second invariant, $d$ is grain size,
     *   $m$ is the grain size exponent, $E$ is activation energy,
     *   $V$ is activation volume, $P$ is pressure, $R$ is the gas
     *   exponent and $T$ is temperature.
     *
     * One may select to use the diffusion ($v_{diff}$; $n=1$, $m!=0$),
     * dislocation ($v_{disl}$, $n>1$, $m=0$) or composite
     * $\frac{v_{diff}v_{disl}}{v_{diff}+v_{disl}}$ equation form.
     *
     * Viscous stress is limited by plastic deformation, which follows
     * a Drucker Prager yield criterion:
     *  $\sigma_y = C\cos(\phi) + P\sin(\phi)$  (2D)
     * or in 3D
     *  $\sigma_y = \frac{6C\cos(\phi) + 2P\sin(\phi)}{\sqrt{3}(3+\sin(\phi))}$
     * where
     *   $\sigma_y$ is the yield stress, $C$ is cohesion, $phi$ is the angle
     *   of internal friction and $P$ is pressure.
     * If the viscous stress ($2v{\varepsilon}_{ii})$) exceeds the yield
     * stress ($\sigma_{y}$), the viscosity is rescaled back to the yield
     * surface: $v_{y}=\sigma_{y}/(2{\varepsilon}_{ii})$
     *
     * When included, the viscoelastic rheology takes into account the elastic shear
     * strength (e.g., shear modulus), while the tensile and volumetric
     * strength (e.g., Young's and bulk modulus) are not considered. The model
     * is incompressible and allows specifying an arbitrary number of
     * compositional fields, where each field represents a different rock type
     * or component of the viscoelastic stress tensor. The symmetric stress tensor in
     * 2D and 3D, respectively, contains 3 or 6 components. The compositional fields
     * representing these components must be named and listed in a very specific
     * format, which is designed to minimize mislabeling stress tensor components
     * as distinct 'compositional rock types' (or vice versa). For 2D models,
     * the first three compositional fields must be labeled stress_xx, stress_yy
     * and stress_xy. In 3D, the first six compositional fields must be labeled
     * stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz.
     *
     * Combining this viscoelasticity implementation with non-linear viscous flow
     * and plasticity produces a constitutive relationship commonly referred to
     * as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community.
     * While extensively discussed and applied within the geodynamics literature,
     * notable references include:
     * Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497.
     * Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105.
     * Gerya (2010), Introduction to Numerical Geodynamic Modeling.
     * Kaus (2010), Tectonophysics, v. 484, p. 36-47.
     * Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444.
     * Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442.
     *
     * The overview below directly follows Moresi et al. (2003) eqns. 23-32.
     * However, an important distinction between this material model and
     * the studies above is the use of compositional fields, rather than
     * tracers, to track individual components of the viscoelastic stress
     * tensor. The material model will be updated when an option to track
     * and calculate viscoelastic stresses with tracers is implemented.
     * Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric
     * rate of deformation ($\hat{D}$) as the sum of elastic
     * ($\hat{D_{e}}$) and viscous ($\hat{D_{v}}$) components:
     * $\hat{D} = \hat{D_{e}} + \hat{D_{v}}$.
     * These terms further decompose into
     * $\hat{D_{v}} = \frac{\tau}{2\eta}$ and
     * $\hat{D_{e}} = \frac{\overset{\triangledown}{\tau}}{2\mu}$, where
     * $\tau$ is the viscous deviatoric stress, $\eta$ is the shear viscosity,
     * $\mu$ is the shear modulus and $\overset{\triangledown}{\tau}$ is the
     * Jaumann corotational stress rate. If plasticity is included the
     * deviatoric rate of deformation may be writted as:
     * $\hat{D} = \hat{D_{e}} + \hat{D_{v}} + \hat{D_{p}}$, where $\hat{D_{p}}$
     * is the plastic component. As defined in the second paragraph, $\hat{D_{p}}$
     * decomposes to $\frac{\tau_{y}}{2\eta_{y}}$, where $\tau_{y}$ is the yield
     * stress and $\eta_{y}$ is the viscosity rescaled to the yield surface.
     *
     * Above, the Jaumann corotational stress rate (eqn. 24) from the elastic
     * component contains the time derivative of the deviatoric stress ($\dot{\tau}$)
     * and terms that  account for material spin (e.g., rotation) due to advection:
     * $\overset{\triangledown}{\tau} = \dot{\tau} + {\tau}W -W\tau$.
     * Above, $W$ is the material spin tensor (eqn. 25):
     * $W_{ij} = \frac{1}{2} \left (\frac{\partial V_{i}}{\partial x_{j}} -
     * \frac{\partial V_{j}}{\partial x_{i}} \right )$.
     *
     * The Jaumann stress-rate can also be approximated using terms from the time
     * at the previous time step ($t$) and current time step ($t + \Delta t^{e}$):
     * $\smash[t]{\overset{\triangledown}{\tau}}^{t + \Delta t^{e}} \approx
     * \frac{\tau^{t + \Delta t^{e} - \tau^{t}}}{\Delta t^{e}} -
     * W^{t}\tau^{t} + \tau^{t}W^{t}$.
     * In this material model, the size of the time step above ($\\Delta t^{e}$)
     * can be specified as the numerical time step size or an independent fixed time
     * step. If the latter case is selected, the user has an option to apply a
     * stress averaging scheme to account for the differences between the numerical
     * and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time
     * step throughout the model run, an equal numerical and elastic time step can be
     * achieved by using CFL and maximum time step values that restrict the numerical
     * time step to the fixed elastic time step.
     *
     * The formulation above allows rewriting the total rate of deformation (eqn. 29) as
     * $\tau^{t + \Delta t^{e}} = \eta_{eff} \left (
     * 2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} +
     * \\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$.
     *
     * The effective viscosity (eqn. 28) is a function of the viscosity ($\eta$),
     * elastic time step size ($\Delta t^{e}$) and shear relaxation time
     * ($ \alpha = \frac{\eta}{\mu} $):
     * $\eta_{eff} = \eta \frac{\Delta t^{e}}{\Delta t^{e} + \alpha}$
     * The magnitude of the shear modulus thus controls how much the effective
     * viscosity is reduced relative to the initial viscosity.
     *
     * Elastic effects are introduced into the governing Stokes equations through
     * an elastic force term (eqn. 30) using stresses from the previous time step:
     * $F^{e,t} = -\frac{\eta_{eff}}{\mu \Delta t^{e}} \tau^{t}$.
     * This force term is added onto the right-hand side force vector in the
     * system of equations.
     *
     * Several model parameters (reference densities, thermal expansivities
     * thermal diffusivities, heat capacities and rheology parameters) can
     * be defined per-compositional field.
     * For each material parameter the user supplies a comma delimited list of
     * length N+1, where N is the number of compositional fields.  The
     * additional field corresponds to the value for background material.  They
     * should be ordered ``background, composition1, composition2...''
     *
     * If a list of values is given for the density, thermal expansivity,
     * thermal diffusivity and heat capacity, the volume weighted sum of the
     * values of each of the compositional fields is used in their place,
     * for example $\rho = \sum \left( \rho_i V_i \right)$
     *
     * The individual output viscosities for each compositional field are
     * also averaged. The user can choose from a range of options for this
     * viscosity averaging. If only one value is given for any of these parameters,
     * all compositions are assigned the same value.
     * The first value in the list is the value assigned to "background material"
     * (regions where the sum of the compositional fields is < 1.0).
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ViscoPlastic : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
        *
        * This material model is incompressible.
         */
        bool is_compressible () const override;

        double reference_viscosity () const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

        double get_min_strain_rate() const;

        /**
         * A function that returns whether the material is plastically yielding at
         * the given pressure, temperature, composition, and strain rate.
         */
        bool
        is_yielding ( const double &pressure,
                      const double &temperature,
                      const std::vector<double> &composition,
                      const SymmetricTensor<2,dim> &strain_rate) const;

      private:

        double min_strain_rate;
        double ref_strain_rate;
        double min_visc;
        double max_visc;
        double ref_visc;

        std::vector<double> thermal_diffusivities;

        /**
         * Whether to use user-defined thermal conductivites instead of thermal diffusivities.
         */
        bool define_conductivities;

        std::vector<double> thermal_conductivities;

        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;

        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         */
        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

        /**
         * Enumeration for selecting which type of viscous flow law to use.
         * Select between diffusion, dislocation or composite.
         */
        enum ViscosityScheme
        {
          diffusion,
          dislocation,
          composite
        } viscous_flow_law;

        /**
         * Enumeration for selecting which type of yield mechanism to use.
         * Select between Drucker Prager and stress limiter.
         */
        enum YieldScheme
        {
          stress_limiter,
          drucker_prager
        } yield_mechanism;

        /**
         * This function calculates viscosities assuming that all the compositional fields
         * experience the same strain rate (isostrain).
         */
        std::pair<std::vector<double>, std::vector<bool> >
        calculate_isostrain_viscosities ( const std::vector<double> &volume_fractions,
                                          const double &pressure,
                                          const double &temperature,
                                          const std::vector<double> &composition,
                                          const SymmetricTensor<2,dim> &strain_rate,
                                          const ViscosityScheme &viscous_type,
                                          const YieldScheme &yield_type) const;


        /**
         * A function that fills the plastic additional output in the
         * MaterialModelOutputs object that is handed over, if it exists.
         * Does nothing otherwise.
         */
        void fill_plastic_outputs (const unsigned int point_index,
                                   const std::vector<double> &volume_fractions,
                                   const bool plastic_yielding,
                                   const MaterialModel::MaterialModelInputs<dim> &in,
                                   MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * A function that fills the viscosity derivatives in the
         * MaterialModelOutputs object that is handed over, if they exist.
         * Does nothing otherwise.
         */
        void compute_viscosity_derivatives(const unsigned int point_index,
                                           const std::vector<double> &volume_fractions,
                                           const std::vector<double> &composition_viscosities,
                                           const MaterialModel::MaterialModelInputs<dim> &in,
                                           MaterialModel::MaterialModelOutputs<dim> &out) const;


        /**
         * A function that returns a ComponentMask that represents all compositional
         * fields that should be considered 'volumetric', that is representing a
         * physical proportion of the material, e.g. volume fraction of peridotite
         * (as opposed to non-volumetric quantities like the amount of finite-strain).
         */
        ComponentMask get_volumetric_composition_mask() const;

        std::vector<double> exponents_stress_limiter;

        /**
         * temperature gradient added to temperature used in the flow law.
         */
        double adiabatic_temperature_gradient_for_viscosity;

        Rheology::StrainDependent<dim> strain_rheology;

        /**
         * Objects for computing viscous creep viscosities.
         */
        Rheology::DiffusionCreep<dim> diffusion_creep;
        Rheology::DislocationCreep<dim> dislocation_creep;

        /**
         * Object for computing the viscosity multiplied by a constant prefactor.
         * This multiplication step is done just prior to calculating the effective
         * viscoelastic viscosity or plastic viscosity.
         */
        Rheology::ConstantViscosityPrefactors<dim> constant_viscosity_prefactors;

        /*
         * Objects for computing plastic stresses, viscosities, and additional outputs
         */
        Rheology::DruckerPrager<dim> drucker_prager_plasticity;

        /*
         * Input parameters for the drucker prager plasticity.
         */
        Rheology::DruckerPragerParameters drucker_prager_parameters;

        /**
         * Object that handles phase transitions.
         */
        MaterialUtilities::PhaseFunction<dim> phase_function;

        /**
         * Object for computing viscoelastic viscosities and stresses.
         */
        Rheology::Elasticity<dim> elastic_rheology;

        /**
         * Whether to include viscoelasticity in the constitutive formulation.
         */
        bool use_elasticity;
    };

  }
}

#endif
