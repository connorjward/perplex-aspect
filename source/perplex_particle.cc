/*
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

#include <aspect/utilities.h>
#include <deal.II/base/exceptions.h>
#include <perplexcpp/wrapper.h>

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
      class PerpleXParticle : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          void
          initialize_one_particle_property(const Point<dim>&,
                                           std::vector<double> &particle_properties) const override
	  {
	    auto wrapper = perplexcpp::Wrapper::get_instance();
	    for (unsigned int i = 0; i < wrapper.get_n_phases() + 1; ++i)
	      for (unsigned int j = 0; j < wrapper.get_n_composition_components(); ++j)
		particle_properties.push_back(0.0);
	  }

	  void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim>&,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim>>&,
                                        const ArrayView<double> &particle_properties) const override
	  {
	    auto wrapper = perplexcpp::Wrapper::get_instance();

	    const double pressure = solution[this->introspection().component_indices.pressure];
	    const double temperature = solution[this->introspection().component_indices.temperature];

	    wrapper.minimize(pressure, temperature);

	    // store as property
	    unsigned int pos { data_position };

	    std::vector<double> bulk_composition{wrapper.get_bulk_composition()};
	    for (unsigned int i = 0; i < wrapper.get_n_composition_components(); ++i) {
	      particle_properties[pos+i] = bulk_composition[i];
	    }
	    pos += wrapper.get_n_composition_components();

	    // phase compositions
	    auto compositions = wrapper.get_phase_compositions();
	    for (unsigned int i = 0; i < wrapper.get_n_phases(); ++i) {
	      for (unsigned int j = 0; j < wrapper.get_n_composition_components(); ++j)
		particle_properties[pos+j] = compositions[i][j];

	      pos += wrapper.get_n_composition_components();
	    }
	  }

	  UpdateTimeFlags need_update() const override
	  {
	    return update_output_step;
	  }

	  UpdateFlags get_needed_update_flags () const override
	  {
	    return update_values;
	  }

          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override
	  {
	    auto wrapper = perplexcpp::Wrapper::get_instance();
	    std::vector<std::pair<std::string,unsigned int>> property_information;

	    auto composition_names = wrapper.get_composition_component_names();
	    auto phase_names = wrapper.get_abbr_phase_names();

	    // bulk composition
	    for (unsigned int i = 0; i < composition_names.size(); ++i)
	      property_information.push_back(std::make_pair("bulk::"+composition_names[i], 1));

	    // phase compositions
	    for (unsigned int i = 0; i < phase_names.size(); ++i)
	      for (unsigned int j = 0; j < composition_names.size(); ++j) 
		property_information.push_back(std::make_pair(phase_names[i]+"::"+composition_names[j], 1));
      
	    return property_information;
	  }

          static
          void
          declare_parameters(ParameterHandler &prm)
	  {
	    prm.enter_subsection("Postprocess");
	    {
	      prm.enter_subsection("Particles");
	      {
		prm.enter_subsection("Perple_X particle");
		{
		  prm.declare_entry("Data directory", ".", Patterns::DirectoryName(),
				    "The location of the Perple_X data files.");
		  prm.declare_entry("Problem definition file", "", Patterns::FileName(),
				    "The name of the PerpleX .dat file in use.");
		}
		prm.leave_subsection();
	      }
	      prm.leave_subsection();
	    }
	    prm.leave_subsection();
	  }

          void
          parse_parameters(ParameterHandler &prm) override
	  {
	    prm.enter_subsection("Postprocess");
	    {
	      prm.enter_subsection("Particles");
	      {
		prm.enter_subsection("Perple_X particle");
		  auto data_dirname = prm.get("Data directory");
		  auto problem_filename = prm.get("Problem definition file");

		  AssertThrow(Utilities::fexists(data_dirname + "/" + problem_filename),
		              ExcMessage("The Perple_X problem file "
			                 "could not be found."));

		  // The Perple_X wrapper must be initialized in this method rather than 
		  // in the (optional) initialize() because initialize() is called after 
		  // get_property_information() which requires Perple_X to have been 
		  // already initialized.
		  perplexcpp::Wrapper::get_instance().initialize(problem_filename, 
		                                                 data_dirname);
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(PerpleXParticle,
                                        "perplex particle",
                                        "A plugin that calls MEEMUM from Perple_X in "
					"order to track properties such as phase "
					"compositions and amounts.")
    }
  }
}
