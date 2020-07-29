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


#include <aspect/particle/property/interface.h>

#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * An enum enumerating the different possible properties that can be tracked.
       */
      enum class PhaseProperty
      {
	composition,
	n_moles,
	molar_fraction,
	volume_fraction,
	weight_fraction,
      };



      /**
       * A class that calculates phase and system properties using the thermodynamical 
       * code Perple_X.
       */
      template <int dim>
      class PerplexComposition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
	  void
	  initialize() override;


          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
           */
          void
          initialize_one_particle_property(const Point<dim> &position,
                                           std::vector<double> &properties) const override;


          /**
           * Update function. This function is called every time an update is
           * request by need_update() for every particle for every property.
           * It is obvious that
           * this function is called a lot, so its code should be efficient.
           * The interface provides a default implementation that does nothing,
           * therefore derived plugins that do not require an update do not
           * need to implement this function.
           *
           * @param [in] data_position An unsigned integer that denotes which
           * component of the particle property vector is associated with the
           * current property. For properties that own several components it
           * denotes the first component of this property, all other components
           * fill consecutive entries in the @p particle_properties vector.
           *
           * @param [in] position The current particle position.
           *
           * @param [in] solution The values of the solution variables at the
           * current particle position.
           *
           * @param [in] gradients The gradients of the solution variables at
           * the current particle position.
           *
           * @param [in,out] particle_properties The properties of the particle
           * that is updated within the call of this function.
           */
	  void
          update_one_particle_property(const unsigned int data_position,
                                       const Point<dim>&,
                                       const Vector<double> &solution,
                                       const std::vector<Tensor<1,dim>>&,
                                       const ArrayView<double> &properties) const override;


          /**
           * Returns an enum, which determines at what times particle properties
           * are updated. The default implementation returns update_never, which
           * signals that particle properties should never be updated.
           * See the documentation of UpdateTimeFlags for a list of possible
           * values and examples for their use. Every
           * plugin that implements this function should return the value
           * appropriate for its purpose, unless it does not need any update,
           * which is the default. This option saves considerable computation
           * time in cases, when no plugin needs to update particle properties
           * over time.
           */
	  UpdateTimeFlags 
	  need_update() const override;


          /**
           * Return which data has to be provided to update all properties.
           * Note that particle properties can only ask for update_default
           * (no data), update_values (solution values), and update_gradients
           * (solution gradients). All other update flags will have no effect.
           *
           * @return The necessary update flags for this particle property.
           */
	  UpdateFlags 
	  get_needed_update_flags () const override;


          /**
           * Set up the information about the names and number of components
           * this property requires. Derived classes need to implement this
           * function.
           *
           * The purpose of this function is to return a vector of pairs of a
           * property name and the number of components associated with this
           * name (e.g. 1 for scalar properties,
           * n for n-dimensional vectors).
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
	  std::vector<std::pair<std::string, unsigned int>>
	  get_property_information() const override;


          /**
           * Declare the parameters this class takes through input files.
           * Derived classes should overload this function if they actually do
           * take parameters; this class declares a fall-back function that
           * does nothing, so that property classes that do not take any
           * parameters do not have to do anything at all.
           *
           * This function is static (and needs to be static in derived
           * classes) so that it can be called without creating actual objects
           * (because declaring parameters happens before we read the input
           * file and thus at a time when we don't even know yet which
           * property objects we need).
           */
	  static
	  void
	  declare_parameters(ParameterHandler &prm);


          /**
           * Read the parameters this class declares from the parameter file.
           * The default implementation in this class does nothing, so that
           * derived classes that do not need any parameters do not need to
           * implement it.
           */
	  void
	  parse_parameters(ParameterHandler &prm) override;


	private:

	  /**
	   * A vector containing the names of the tracked Perple_X phases.
	   */
	  std::vector<std::string> phase_names;


	  /**
	   * A vector containing the selected phase properties to report.
	   */
	  std::vector<PhaseProperty> phase_properties;


	  /**
	   * Flag indicating whether or not to report the bulk composition.
	   */
	  bool track_bulk_composition;


	  /**
	   * Flag indicating whether or not to extract melt.
	   */
	  bool extract_melt;


	  /**
	   * The volume fraction of melt required in order to trigger melt extraction.
	   */
	  double melt_extraction_threshold;


	  /**
	   * The minimum amount of substance (in moles) necessary for a
	   * sensible call to Perple_X. The reason for this is because 
	   * Perple_X can produce strange results if called with negligible 
	   * amounts of substance.
	   */
	  double min_amount_of_substance;


	  /**
	   * Retrieve the composition in the correct order for Perple_X from
	   * the compositional fields.
	   */
	  std::vector<double>
	  get_composition(const Vector<double> &solution) const;


	  /**
	   * Retrieve the particle properties from the MinimizeResult in the
	   * correct order.
	   */
	  std::vector<double>
	  get_properties(const perplexcpp::MinimizeResult &result) const;


	  /**
	   * Return a zero array the same size as the required particle 
	   * properties.
	   */
	  std::vector<double>
	  get_zero_properties() const;
      };
    }
  }
}
