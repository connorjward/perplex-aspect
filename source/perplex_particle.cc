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
	  get_composition(const Vector<double> &solution) const
	  { 
	    // Retrieve the compositional fields from the solution.
	    std::vector<double> compositional_fields(this->introspection().n_compositional_fields);
	    solution.extract_subvector_to(this->introspection().component_indices.compositional_fields,
		                          compositional_fields);

	    // Load the composition in the correct order for Perple_X.
	    std::vector<double> composition;
	    for (std::string comp_name : perplexcpp::Wrapper::get_instance().composition_component_names) 
	    {
	      const unsigned int idx = 
		this->introspection().compositional_index_for_name(comp_name);

	      const double n_moles = compositional_fields[idx];

	      // Some values end up being just below zero. 
	      // Zeroing them allows Perple_X to work.
	      if (n_moles >= 0)
		composition.push_back(n_moles);
	      else
		composition.push_back(0.0);
	    }
	    return composition;
	  }


        public:

          void
          initialize_one_particle_property(const Point<dim>&,
                                           std::vector<double> &particle_properties) 
	  const override
	  {
	    auto& wrapper = perplexcpp::Wrapper::get_instance();

	    for (unsigned int i = 0; i < tracked_phases.size(); ++i) {
	      for (std::string property : tracked_phase_properties) {
		if (property == "composition") {
		  for (unsigned int j = 0; j < wrapper.n_composition_components; ++j)
		    particle_properties.push_back(0.0);
		} else {
		  particle_properties.push_back(0.0);
		}
	      }
	    }

	    if (track_bulk_composition)
	      for (unsigned int i = 0; i < wrapper.n_composition_components; ++i)
		particle_properties.push_back(0.0);
	  }



	  void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim>&,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim>>&,
                                        const ArrayView<double> &particle_properties) 
	  const override
	  {
	    auto& px = perplexcpp::Wrapper::get_instance();

	    double pressure = solution[this->introspection().component_indices.pressure];
	    double temperature = solution[this->introspection().component_indices.temperature];

	    if (pressure < px.min_pressure)
	      pressure = px.min_pressure;
	    else if (pressure > px.max_pressure)
	      pressure = px.max_pressure;
	    if (temperature < px.min_temperature)
	      temperature = px.min_temperature;
	    if (temperature > px.max_temperature)
	      temperature = px.max_temperature;

	    const std::vector<double> composition = get_composition(solution);

	    // Track the current data position.
	    unsigned int current_position = data_position;

	    if (std::accumulate(composition.begin(), composition.end(), 0) > this->min_amount_of_substance)
	    {
	      const perplexcpp::MinimizeResult result = 
		px.minimize(pressure, temperature, composition);

	      for (std::string name : tracked_phases) 
	      {
		perplexcpp::Phase phase = perplexcpp::find_phase(result.phases, name);

		for (std::string property : tracked_phase_properties) {
		  if (property == "composition") {
		    for (double comp : phase.composition_ratio) {
		      particle_properties[current_position] = comp;
		      current_position++;
		    }
		  } 
		  else if (property == "molar amount") {
		    particle_properties[current_position] = phase.n_moles;
		    current_position++;
		  } 
		  else if (property == "molar fraction") {
		    particle_properties[current_position] = phase.molar_frac;
		    current_position++;
		  } 
		  else if (property == "volume fraction") {
		    particle_properties[current_position] = phase.volume_frac;
	 	    current_position++;
		  } 
		  else if (property == "weight fraction") {
		    particle_properties[current_position] = phase.weight_frac;
		    current_position++;
		  } 
		  else {
		    Assert(false, ExcInternalError(property + " could not be found."));
		  }
		}
	      }

	      std::vector<double> residue_composition = result.composition;

	      if (this->extract_melt)
	      {
		const perplexcpp::Phase melt = perplexcpp::find_phase(result.phases, "liquid");

		if (melt.volume_frac > this->melt_extraction_threshold)
		{
		  for (unsigned int c = 0; c < px.n_composition_components; c++)
		  {
		    Assert(melt.composition_ratio[c] * melt.n_moles >= 0,
			   ExcInternalError("The extracted melt should be non-negative"));
		    residue_composition[c] -= melt.composition_ratio[c] * melt.n_moles;

		    if (residue_composition[c] < 0)
		      residue_composition[c] = 0.0;
		  }
		}
	      }

	      if (track_bulk_composition) {
		for (double comp : residue_composition) 
		{
		  Assert(comp >=0, ExcInternalError("The new composition should be non-negative"));
		  particle_properties[current_position] = comp;
		  current_position++;
		}
	      }
	    }
	    else
	    {
	      for (std::string name : tracked_phases) 
	      {
		for (std::string property : tracked_phase_properties) {
		  if (property == "composition") {
		    for (unsigned int c = 0; c < px.n_composition_components; c++) 
		    {
		      particle_properties[current_position] = 0.0;
		      current_position++;
		    }
		  } 
		  else if (property == "molar amount") {
		    particle_properties[current_position] = 0.0;
		    current_position++;
		  } 
		  else if (property == "molar fraction") {
		    particle_properties[current_position] = 0.0;
		    current_position++;
		  } 
		  else if (property == "volume fraction") {
		    particle_properties[current_position] = 0.0;
	 	    current_position++;
		  } 
		  else if (property == "weight fraction") {
		    particle_properties[current_position] = 0.0;
		    current_position++;
		  } 
		  else {
		    Assert(false, ExcInternalError(property + " could not be found."));
		  }
		}
	      }

	      if (track_bulk_composition) {
		for (unsigned int c = 0; c < px.n_composition_components; c++) {
		  particle_properties[current_position] = 0.0;
		  current_position++;
		}
	      }
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
	    std::vector<std::pair<std::string,unsigned int>> property_information;

	    const auto& wrapper = perplexcpp::Wrapper::get_instance();

	    // Iterate through the tracked phases and the selected phase properties for each.
	    for (std::string phase_name : tracked_phases) {
	      for (std::string property : tracked_phase_properties) {
		if (property == "composition") {
		  for (std::string composition_name : wrapper.composition_component_names)
		    property_information.push_back(std::make_pair(phase_name + " " 
			                                          + composition_name, 1));
		} 
		else if (property == "molar amount") {
		  property_information.push_back(std::make_pair("molar amount " +
			phase_name, 1));
		} 
		else if (property == "molar fraction") {
		  property_information.push_back(std::make_pair("molar fraction " +
			phase_name, 1));
		} 
		else if (property == "volume fraction") {
		  property_information.push_back(std::make_pair("volume fraction " +
			phase_name, 1));
		} 
		else if (property == "weight fraction") {
		  property_information.push_back(std::make_pair("weight fraction " +
			phase_name, 1));
		} 
		else {
		  Assert(false, ExcInternalError("Could not find phase property '" + property + "'."));
		}
	      }
	    }

	    if (track_bulk_composition)
	      for (std::string name : wrapper.composition_component_names)
		property_information.push_back(std::make_pair("bulk " + name, 1));

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
		  prm.declare_entry(
		    "Data directory", 
		    ".", 
		    Patterns::DirectoryName(),
		    "The location of the Perple_X data files."
		  );

		  prm.declare_entry(
		    "Problem definition file", 
		    "", 
		    Patterns::FileName(),
		    "The name of the PerpleX .dat file in use.",
		    true
		  );

		  prm.declare_entry(
		    "List of phases",
		    "all", 
		    Patterns::List(Patterns::Anything()),
		    "The phases that will be tracked during the simulation. "
		    "Possible options are: 'all' or any solution phase listed in "
		    "the Perple_X problem definition file (e.g. melt(HGP)). These can "
		    "be in any of the three formats supported by Perple_X (e.g. "
		    "melt can be written as 'melt(HGP)', 'Melt' or 'liquid')."
		  );

		  prm.declare_entry(
		    "List of phase properties",
		    "molar amount",
		    Patterns::MultipleSelection("composition|molar amount|"
		                                "molar fraction|volume fraction|"
		                                "weight fraction"),
		    "The phase properties to track during the simulation."
		  );

		  prm.declare_entry(
		    "Track bulk composition",
		    "false",
		    Patterns::Bool(),
		    "Whether to track the evolution of the bulk composition. In batch "
		    "melting this will be constant. Required for fractional melting."
		  );

		  prm.declare_entry(
		    "Extract melt",
		    "false",
		    Patterns::Bool(),
		    "Whether to extract the melt composition as the simulation "
		    "progresses or not."
		  );

		  prm.declare_entry(
		    "Melt extraction threshold",
		    "0.05",
		    Patterns::Double(0.0, 1.0),
		    "The volume fraction of the bulk material that must be melt "
		    "in order to trigger melt extraction."
		  );

		  prm.declare_entry(
		    "Minimum amount of substance",
		    "1.0",
		    Patterns::Double(0.0),
		    "The minimum number of moles of material necessary for a valid "
		    "Perple_X call."
		  );
		}
		prm.leave_subsection();
	      }
	      prm.leave_subsection();
	    }
	    prm.leave_subsection();
	  }



	  void
	  parse_parameters(ParameterHandler &prm) 
	  override
	  {
	    prm.enter_subsection("Postprocess");
	    {
	      prm.enter_subsection("Particles");
	      {
		prm.enter_subsection("Perple_X particle");
		{
		  const std::string data_dirname = prm.get("Data directory");
		  const std::string problem_filename = prm.get("Problem definition file");

		  AssertThrow(Utilities::fexists(data_dirname + "/" + problem_filename),
			      ExcMessage("The Perple_X problem file could not be found."));

		  perplexcpp::Wrapper::initialize(problem_filename, data_dirname);
		  auto& perplex_wrapper = perplexcpp::Wrapper::get_instance();

		  this->tracked_phases = Utilities::split_string_list(prm.get("List of phases"));

		  // See if 'all' was selected (only valid if no other phases are included). 
		  // If so simply replace the list with one that contains all the phases. If not,
		  // find the phase indices that correspond to the submitted phase names.
		  if (std::find(this->tracked_phases.begin(), this->tracked_phases.end(), "all") 
		      != this->tracked_phases.end()) 
		  {
		    AssertThrow(this->tracked_phases.size() == 1, 
				ExcMessage("'all' specified for the parameter " 
					   "'List of tracked phases' but other phases are also "
					   "specified. Please check your parameter file."));

		    this->tracked_phases.clear();
		    for (perplexcpp::PhaseName name : perplex_wrapper.phase_names)
		      this->tracked_phases.push_back(name.full);
		  }

		  // Check that all of the specified phase names actually exist.
		  for (std::string tracked_name : tracked_phases) 
		  {
		    bool found = false;

		    for (perplexcpp::PhaseName name : perplex_wrapper.phase_names) {
		      if (tracked_name == name.standard ||
			  tracked_name == name.abbreviated ||
			  tracked_name == name.full)
		      {
			found = true;
			break;
		      }
		    }

		    AssertThrow(found, 
			        ExcMessage("The phase '" + tracked_name 
				           + "' could not be found."));
		  }

		  this->tracked_phase_properties = Utilities::split_string_list(prm.get("List of phase "
										  "properties"));
		  AssertThrow(Utilities::has_unique_entries(tracked_phase_properties),
			      ExcMessage("The list of strings for the parameter "
					 "'Postprocess/Particles/Perple_X particle/Phase properties "
					 "contains entries more than once. "
					 "This is not allowed. Please check your parameter file."));

		  this->track_bulk_composition = prm.get_bool("Track bulk composition");

		  this->extract_melt = prm.get_bool("Extract melt");

		  if (this->extract_melt)
		    AssertThrow(this->track_bulk_composition,
			        ExcMessage("In order to be able to extract melt from the "
				           "simulation 'Track bulk composition' must be enabled."));

		  this->melt_extraction_threshold = prm.get_double("Melt extraction threshold");
		  this->min_amount_of_substance = prm.get_double("Minimum amount of substance");
		}
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
