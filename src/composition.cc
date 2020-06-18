/*
 Copyright (C) 2015 - 2018 by the authors of the ASPECT code.
 Copyright (C) 2020 Connor Ward.

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

#include <aspect/particle/property/interface.h>

#include <perplex/state.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that initializes particle properties based on a
       * functional description provided in the input file.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class PhaseComposition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        private:
	  std::string perplex_dat_filename;

        public:
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
          initialize_one_particle_property(const Point<dim>&,
                                           std::vector<double> &particle_properties) const override
	  {
	    perplex::State state{perplex::State::get_instance()};
	    for (unsigned int i = 0; i < state.get_n_soln_phases() + 1; ++i)
	      for (unsigned int j = 0; j < state.get_n_composition_components(); ++j)
		particle_properties.push_back(0.0);
	  }

          /**
           * Update function. This function is called every time an update is
           * request by need_update() for every particle for every property.
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
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim>&,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim>>&,
                                        const ArrayView<double> &particle_properties) const override
	  {
	    perplex::State state{perplex::State::get_instance()};

	    const double pressure = solution[this->introspection().component_indices.pressure];
	    const double temperature = solution[this->introspection().component_indices.temperature];

	    state.minimize(pressure, temperature);

	    // store as property
	    unsigned int pos { data_position };

	    std::vector<double> bulk_composition{state.get_bulk_composition()};
	    for (unsigned int i = 0; i < state.get_n_composition_components(); ++i) {
	      particle_properties[pos+i] = bulk_composition[i];
	    }
	    pos += state.get_n_composition_components();

	    // phase compositions
	    std::vector<std::string> soln_names{state.get_abbr_soln_phase_names()};
	    for (unsigned int i = 0; i < state.get_n_soln_phases(); ++i) {
	      for (unsigned int j = 0; j < state.get_n_end_phases(); ++j) {
		std::string name{state.get_end_phase_name(j)};
		if (state.find_abbr_phase_name(name) == soln_names[i]) {
		  std::vector<double> composition{state.get_end_phase_composition(j)};
		  for (unsigned int k = 0; k < state.get_n_composition_components(); ++k)
		    particle_properties[pos+k] = composition[k];

		  pos += state.get_n_composition_components();
		  break;
		}
	      }
	    }
	  }

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           *         number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override
	  {
	    perplex::State state{perplex::State::get_instance()};
	    std::vector<std::pair<std::string,unsigned int>> property_information;

	    std::vector<std::string> comp_names(state.get_composition_component_names());
	    std::vector<std::string> phase_names(state.get_abbr_soln_phase_names());

	    // bulk composition
	    for (unsigned int i = 0; i < comp_names.size(); ++i)
	      property_information.push_back(std::make_pair("bulk::"+comp_names[i], 1));

	    // phase compositions
	    for (unsigned int i = 0; i < phase_names.size(); ++i)
	      for (unsigned int j = 0; j < comp_names.size(); ++j) 
		property_information.push_back(std::make_pair(phase_names[i]+"::"+comp_names[j], 1));
      
	    return property_information;
	  }

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters(ParameterHandler &prm)
	  {
	    prm.enter_subsection("Postprocess");
	    {
	      prm.enter_subsection("Particles");
	      {
		prm.enter_subsection("Phase composition");
		{
		  prm.declare_entry("PerpleX file name", "",
		                    Patterns::FileName(),
				    "The name of the PerpleX .dat file in use.");
		}
		prm.leave_subsection();
	      }
	      prm.leave_subsection();
	    }
	    prm.leave_subsection();
	  }

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters(ParameterHandler &prm) override
	  {
	    prm.enter_subsection("Postprocess");
	    {
	      prm.enter_subsection("Particles");
	      {
		prm.enter_subsection("Phase composition");
		  // TODO
		  // Add decent assertions here to catch errors.
		  // specify phases to measure, defaulting to 'all'
		  // have different ways of specifying the initial composition (default: from file)
		  perplex_dat_filename = prm.get("PerpleX file name");
		  perplex::State state{perplex::State::get_instance()};
		    
		  state.initialize(perplex_dat_filename);
		prm.leave_subsection();
	      }
	      prm.leave_subsection();
	    }
	    prm.leave_subsection();
	  }
      };
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(PhaseComposition,
                                        "phase composition",
                                        "Information about plugin.")
    }
  }
}
