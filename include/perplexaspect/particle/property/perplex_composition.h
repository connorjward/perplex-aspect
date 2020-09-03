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
       * An enum listing the different possible properties that can be tracked.
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

          void
          initialize_one_particle_property(const Point<dim> &position,
                                           std::vector<double> &properties) const override;


	  void
          update_one_particle_property(const unsigned int data_position,
                                       const Point<dim>&,
                                       const Vector<double> &solution,
                                       const std::vector<Tensor<1,dim>>&,
                                       const ArrayView<double> &properties) const override;


	  UpdateTimeFlags 
	  need_update() const override;


	  UpdateFlags 
	  get_needed_update_flags () const override;


	  std::vector<std::pair<std::string, unsigned int>>
	  get_property_information() const override;


	  static
	  void
	  declare_parameters(ParameterHandler &prm);


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
	  void
	  put_properties(const perplexcpp::MinimizeResult &result, 
	                 ArrayView<double>::iterator &iter) const;


	  /**
	   * Return a zero array the same size as the required particle 
	   * properties.
	   */
	  void
	  put_zero_properties(ArrayView<double>::iterator &iter) const;


	  /**
	   * Save the properties of the extracted melt.
	   */
	  void
	  put_extracted_melt_properties(const perplexcpp::Phase &melt,
				        ArrayView<double>::iterator &iter) const;


	  /**
	   * Save the properties of the phases.
	   */
	  void
	  put_phase_properties(const perplexcpp::MinimizeResult &result,
			       const std::string &phase_name,
			       ArrayView<double>::iterator &it) const;
      };
    }
  }
}
