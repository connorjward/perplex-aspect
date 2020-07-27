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

/* #include <aspect/initial_composition/interface.h> */
/* #include <aspect/utilities.h> */
/* #include <deal.II/base/exceptions.h> */
/* #include <perplexcpp/wrapper.h> */


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
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
	  std::vector<std::string> tracked_phases;


	  /**
	   * A vector containing the selected phase properties to report.
	   */
	  std::vector<std::string> tracked_phase_properties;



	  /**
	   * Flag indicating whether or not to report the bulk composition.
	   */
	  bool track_bulk_composition;



	  /**
	   * ???
	   */
	  bool extract_melt;



	  /**
	   * ???
	   */
	  double melt_extraction_threshold;


	  /**
	   * ???
	   */
	  double min_amount_of_substance;



	  /**
	   * ???
	   */
	  std::vector<double>
	  get_composition(const Vector<double> &solution) const;
      };
    }
  }
}
