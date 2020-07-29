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


#include <perplexaspect/particle/property/perplex_composition.h>

#include <aspect/initial_composition/interface.h>
#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      PerplexComposition<dim>::initialize()
      {
	if (this->extract_melt)
	  AssertThrow(this->track_bulk_composition,
		      ExcMessage("In order to be able to extract melt from the "
				 "simulation 'Track bulk composition' must be enabled."));
      }



      template <int dim>
      void
      PerplexComposition<dim>::
      initialize_one_particle_property(const Point<dim> &position,
				       std::vector<double> &properties) const 
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	for (std::string pname : phase_names) 
	{
	  for (PhaseProperty pprop : phase_properties) 
	  {
	    switch (pprop)
	    {
	      case PhaseProperty::composition:
		for (unsigned int c = 0; c < px.n_composition_components; c++)
		  properties.push_back(0.0);
		break;
	      default:
		properties.push_back(0.0);
	    }
	  }
	}

	// Load the initial composition specified in the parameter file.
	if (track_bulk_composition)
	{
	  for (std::string cname : px.composition_component_names)
	  {
	    const unsigned int idx = this->introspection().compositional_index_for_name(cname);
	    properties.push_back(
	      this->get_initial_composition_manager().initial_composition(position, idx)
	    );
	  }
	}
      }



      template <int dim>
      void
      PerplexComposition<dim>::
      update_one_particle_property(const unsigned int data_position,
				   const Point<dim>&,
				   const Vector<double> &solution,
				   const std::vector<Tensor<1,dim>>&,
				   const ArrayView<double> &particle_properties) const
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	double pressure = solution[this->introspection().component_indices.pressure];
	double temperature = solution[this->introspection().component_indices.temperature];

	// Adjust the pressure and temperature to lie within the bounds
	// permitted by Perple_X.
	if (pressure < px.min_pressure)
	  pressure = px.min_pressure;
	else if (pressure > px.max_pressure)
	  pressure = px.max_pressure;
	if (temperature < px.min_temperature)
	  temperature = px.min_temperature;
	if (temperature > px.max_temperature)
	  temperature = px.max_temperature;

	const std::vector<double> composition = get_composition(solution);

	const double n_moles = std::accumulate(composition.begin(), composition.end(), 0);

	std::vector<double> props;
	if (n_moles > this->min_amount_of_substance)
	{
	  const perplexcpp::MinimizeResult result = 
	    px.minimize(pressure, temperature, composition);

	  props = get_properties(result);
	}
	else
	{
	  props = get_zero_properties();
	}

	// Copy the props array into the solution.
	for (unsigned int i = 0; i < props.size(); i++)
	  particle_properties[data_position+i] = props[i];
      }



      template <int dim>
      UpdateTimeFlags 
      PerplexComposition<dim>::need_update() const
      {
	return update_output_step;
      }



      template <int dim>
      UpdateFlags 
      PerplexComposition<dim>::get_needed_update_flags () const
      {
	return update_values;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      PerplexComposition<dim>::get_property_information() const 
      {
	static const std::unordered_map<PhaseProperty,std::string> pprop_map = {
	  { PhaseProperty::composition, "composition" },
	  { PhaseProperty::n_moles, "n moles" },
	  { PhaseProperty::molar_fraction, "molar fraction" },
	  { PhaseProperty::volume_fraction, "volume fraction" },
	  { PhaseProperty::weight_fraction, "weight fraction" },
	};

	const auto& px = perplexcpp::Wrapper::get_instance();

	std::vector<std::pair<std::string,unsigned int>> prop_info;
	for (std::string pname : phase_names) 
	{
	  for (PhaseProperty pprop : phase_properties) 
	  {
	    switch (pprop)
	    {
	      case PhaseProperty::composition:
		for (std::string cname : px.composition_component_names)
		  // E.g. 'composition liquid SiO2'.
		  prop_info.push_back(std::make_pair(pprop_map.at(pprop) + " " + pname + " " + cname, 1));
		break;
	      default:
		// E.g. 'volume fraction olivine'.
		prop_info.push_back(std::make_pair(pprop_map.at(pprop) + " " + pname, 1));
	    }
	  }
	}

	if (track_bulk_composition)
	  for (std::string cname : px.composition_component_names)
	    // E.g. 'bulk SiO2'.
	    prop_info.push_back(std::make_pair("bulk " + cname, 1));

	return prop_info;
      }



      template <int dim>
      void
      PerplexComposition<dim>::declare_parameters(ParameterHandler &prm)
      {
	prm.enter_subsection("Perple_X configuration");
	{
	  prm.declare_entry(
	    "Data directory", 
	    ".", 
	    Patterns::DirectoryName(),
	    "The location of the Perple_X data files."
	  );

	  prm.declare_entry(
	    "Problem file", 
	    "", 
	    Patterns::FileName(),
	    "The name of the PerpleX .dat file in use (within the specified "
	    "directory).",
	    true
	  );

	  prm.declare_entry(
	    "Cache capacity",
	    "0",
	    Patterns::Integer(0),
	    "The number of results held in the cache."
	  );

	  prm.declare_entry(
	    "Cache tolerance",
	    "0.0",
	    Patterns::Double(0.0, 1.0),
	    "The relative tolerance accepted by the cache. A result will only "
	    "be returned if all of the input parameters (temperature, pressure "
	    "and composition) vary by less than this amount."
	  );
	}
	prm.leave_subsection();

	prm.enter_subsection("Postprocess");
	{
	  prm.enter_subsection("Particles");
	  {
	    prm.enter_subsection("Perple_X composition");
	    {
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
		"none",
		Patterns::MultipleSelection("none|composition|molar amount|"
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



      namespace
      {
	void
	initialize_perplex(const ParameterHandler &prm)
	{
	  const std::string dirname = prm.get("Data directory");
	  const std::string fname = prm.get("Problem file");
	  const unsigned int cache_capacity = prm.get_integer("Cache capacity");
	  const double cache_rtol = prm.get_double("Cache tolerance");

	  AssertThrow(Utilities::fexists(dirname + "/" + fname),
		      ExcMessage("The Perple_X problem file could not be found."));

	  perplexcpp::Wrapper::initialize(fname, dirname, cache_capacity, cache_rtol);
	}



	std::vector<std::string>
	parse_phase_names(const ParameterHandler &prm)
	{
	  const auto& px = perplexcpp::Wrapper::get_instance();

	  std::vector<std::string> phase_names = 
	    Utilities::split_string_list(prm.get("List of phases"));

	  // See if 'all' was selected.
	  if (std::find(phase_names.begin(), phase_names.end(), "all") != phase_names.end()) 
	  {
	    AssertThrow(phase_names.size() == 1, 
			ExcMessage("'all' specified for the parameter 'Postprocess/Particle/"
			           "Perple_X composition/List of tracked phases' but other phases "
				   "are also specified. This is not allowed. Please check your "
				   "parameter file."));

	    // Return a list containing all of the possible phase names.
	    phase_names.clear();
	    for (perplexcpp::PhaseName pname : px.phase_names)
	      phase_names.push_back(pname.full);
	    return phase_names;
	  }

	  // Check that all of the specified phase names actually exist.
	  for (std::string pname : phase_names) 
	  {
	    bool found = false;

	    for (perplexcpp::PhaseName pname_actual : px.phase_names) {
	      if (pname == pname_actual.standard ||
		  pname == pname_actual.abbreviated ||
		  pname == pname_actual.full)
	      {
		found = true;
		break;
	      }
	    }
	    AssertThrow(found, ExcMessage("The phase " + pname + " could not be found."));
	  }
	  return phase_names;
	}



	std::vector<PhaseProperty>
	parse_phase_properties(const ParameterHandler &prm)
	{
	  static const std::unordered_map<std::string,PhaseProperty> pprop_map = {
	    { "composition", PhaseProperty::composition },
	    { "molar amount", PhaseProperty::n_moles },
	    { "molar fraction", PhaseProperty::molar_fraction },
	    { "volume fraction", PhaseProperty::volume_fraction },
	    { "weight fraction", PhaseProperty::weight_fraction },
	  };

	  const std::vector<std::string> pprops_str = 
	    Utilities::split_string_list(prm.get("List of phase properties"));

	  AssertThrow(Utilities::has_unique_entries(pprops_str),
		      ExcMessage("The list of strings for the parameter "
				 "'Postprocess/Particles/Perple_X particle/Phase properties "
				 "contains entries more than once. "
				 "This is not allowed. Please check your parameter file."));

	  // Check for 'none' option.
	  if (std::find(pprops_str.begin(), pprops_str.end(), "none") != pprops_str.end())
	  {
	    AssertThrow(pprops_str.size() == 1,
			ExcMessage("The list of strings for the parameter 'Postprocess/"
				   "Particles/Perple_X composition/Phase properties' "
				   "contains the option 'none' alongside other options. "
				   "This is not allowed. Please check your parameter file."));

	    // Return an empty vector.
	    return std::vector<PhaseProperty>();
	  }

	  std::vector<PhaseProperty> pprops;
	  for (std::string pprop_str : pprops_str)
	  {
	    try {
	      pprops.push_back(pprop_map.at(pprop_str));
	    }
	    catch (const std::out_of_range &e) {
	      AssertThrow(false, ExcMessage(pprop_str + " is not a valid option. "
					    "Check your parameter file."));
	    }
	  }
	  return pprops;
	}
      }


      template <int dim>
      void
      PerplexComposition<dim>::parse_parameters(ParameterHandler &prm) 
      {
	prm.enter_subsection("Perple_X configuration");
	{
	  initialize_perplex(prm);
	}
	prm.leave_subsection();

	prm.enter_subsection("Postprocess");
	{
	  prm.enter_subsection("Particles");
	  {
	    prm.enter_subsection("Perple_X composition");
	    {
	      this->phase_names = parse_phase_names(prm);
	      this->phase_properties = parse_phase_properties(prm);
	      this->track_bulk_composition = prm.get_bool("Track bulk composition");
	      this->extract_melt = prm.get_bool("Extract melt");
	      this->melt_extraction_threshold = prm.get_double("Melt extraction threshold");
	      this->min_amount_of_substance = prm.get_double("Minimum amount of substance");
	    }
	    prm.leave_subsection();
	  }
	  prm.leave_subsection();
	}
	prm.leave_subsection();
      }



      template <int dim>
      std::vector<double>
      PerplexComposition<dim>::get_composition(const Vector<double> &solution) const
      { 
	const auto& px = perplexcpp::Wrapper::get_instance();

	// Retrieve the compositional fields from the solution.
	std::vector<double> compositional_fields(this->introspection().n_compositional_fields);
	solution.extract_subvector_to(this->introspection().component_indices.compositional_fields,
				      compositional_fields);

	// Load the composition in the correct order for Perple_X. Set any
	// negative values to zero.
	std::vector<double> composition;
	for (std::string cname : px.composition_component_names) 
	{
	  const unsigned int idx = this->introspection().compositional_index_for_name(cname);
	  composition.push_back(compositional_fields[idx] >= 0 
				? compositional_fields[idx] : 0.0);
	}
	return composition;
      }



      template <int dim>
      std::vector<double>
      PerplexComposition<dim>::
      get_properties(const perplexcpp::MinimizeResult &result) const
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	std::vector<double> props;
	for (std::string pname : this->phase_names) 
	{
	  const perplexcpp::Phase phase = perplexcpp::find_phase(result.phases, pname);

	  for (PhaseProperty pprop : this->phase_properties) 
	  {
	    switch (pprop)
	    {
	      case PhaseProperty::composition:
		for (double comp : phase.composition_ratio) 
		  props.push_back(comp);
		break;
	      case PhaseProperty::n_moles:
		props.push_back(phase.n_moles);
		break;
	      case PhaseProperty::molar_fraction:
		props.push_back(phase.molar_frac);
		break;
	      case PhaseProperty::volume_fraction:
		props.push_back(phase.volume_frac);
		break;
	      case PhaseProperty::weight_fraction:
		props.push_back(phase.weight_frac);
		break;
	    }
	  }
	}

	std::vector<double> residue_composition = result.composition;
	if (this->extract_melt)
	{
	  const perplexcpp::Phase melt = perplexcpp::find_phase(result.phases, "liquid");

	  if (melt.volume_frac > this->melt_extraction_threshold)
	    for (unsigned int c = 0; c < px.n_composition_components; c++)
	      residue_composition[c] -= melt.composition_ratio[c] * melt.n_moles;	    
	}

	if (this->track_bulk_composition) 
	  for (double camount : residue_composition) 
	    props.push_back(camount);

	return props;
      }



      template <int dim>
      std::vector<double>
      PerplexComposition<dim>::get_zero_properties() const
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	std::vector<double> props;

	for (std::string pname : this->phase_names) 
	{
	  for (PhaseProperty pprop : this->phase_properties)
	  {
	    switch (pprop)
	    {
	      case PhaseProperty::composition:
		for (unsigned int c = 0; c < px.n_composition_components; c++) 
		  props.push_back(0.0);
		break;
	      default:
		props.push_back(0.0);
	    }
	  }
	}

	if (this->track_bulk_composition)
	  for (unsigned int c = 0; c < px.n_composition_components; c++)
	    props.push_back(0.0);

	return props;
      }
    }
  }
}



// Explicit instantiations.
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(PerplexComposition,
                                        "perplex composition",
                                        "A plugin that calls MEEMUM from Perple_X in "
					"order to track properties such as phase "
					"compositions and amounts.")
    }
  }
}
