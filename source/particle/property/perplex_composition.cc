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
#include <perplexaspect/utilities.h>
#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      PerplexComposition<dim>::
      initialize_one_particle_property(const Point<dim> &,
				       std::vector<double> &properties) const 
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	// Bulk composition.
	for (double camount : px.initial_bulk_composition)
	  properties.push_back(camount);

	// Phase properties.
	for (std::string pname : this->phase_names) {
	  for (PhaseProperty pprop : this->phase_properties) {
	    switch (pprop) {
	      case PhaseProperty::composition:
		for (unsigned int c = 0; c < px.n_composition_components; c++)
		  properties.push_back(0.0);
		break;
	      default:
		properties.push_back(0.0);
	    }
	  }
	}

	// Extracted melt properties.
	if (this->extract_melt) {
	  for (PhaseProperty pprop : this->phase_properties) {
	    switch (pprop) {
	      case PhaseProperty::composition:
		for (unsigned int c = 0; c < px.n_composition_components; c++)
		  properties.push_back(0.0);
		break;
	      default:
		properties.push_back(0.0);
	    }
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

	std::vector<double> composition;
	for (unsigned int c = 0; c < px.n_composition_components; c++)
	  composition.push_back(particle_properties[data_position+c]);

	ArrayView<double>::iterator iter = particle_properties.begin() + data_position;

	const double n_moles = std::accumulate(composition.begin(), composition.end(), 0);
	if (n_moles > this->min_amount_of_substance) {
	  const perplexcpp::MinimizeResult result = 
	    px.minimize(pressure, temperature, composition);

	  this->put_properties(result, iter);
	} else {
	  this->put_zero_properties(iter);
	}
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

	std::vector<std::pair<std::string,unsigned int>> info;

	// Bulk composition (e.g. 'bulk composition SiO2'). 
	for (std::string cname : px.composition_component_names)
	  info.push_back(std::make_pair("bulk composition " + cname, 1));

	// Phase properties.
	for (std::string pname : this->phase_names) {
	  for (PhaseProperty pprop : this->phase_properties) {
	    switch (pprop) {
	      case PhaseProperty::composition:
		// E.g. 'liquid composition SiO2'.
		for (std::string cname : px.composition_component_names)
		  info.push_back(std::make_pair(pname+" "+pprop_map.at(pprop)+" "+cname, 1));
		break;
	      default:
		// E.g. 'olivine volume fraction'.
		info.push_back(std::make_pair(pname+" "+pprop_map.at(pprop), 1));
	    }
	  }
	}

	// Extracted melt properties.
	if (this->extract_melt) {
	  for (PhaseProperty pprop : this->phase_properties) {
	    switch (pprop) {
	      case PhaseProperty::composition:
		// E.g. 'extracted melt composition SiO2'.
		for (std::string cname : px.composition_component_names)
		  info.push_back(std::make_pair("extracted melt "+pprop_map.at(pprop)+" "+cname, 1));
		break;
	      default:
		// E.g. 'extracted melt volume fraction'.
		info.push_back(std::make_pair("extracted melt "+pprop_map.at(pprop), 1));
	    }
	  }
	}

	return info;
      }


      template <int dim>
      void
      PerplexComposition<dim>::declare_parameters(ParameterHandler &prm)
      {
	PerplexUtils::declare_parameters(prm);

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
	PerplexUtils::parse_parameters(prm);

	prm.enter_subsection("Postprocess");
	{
	  prm.enter_subsection("Particles");
	  {
	    prm.enter_subsection("Perple_X composition");
	    {
	      this->phase_names = parse_phase_names(prm);
	      this->phase_properties = parse_phase_properties(prm);
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
      void
      PerplexComposition<dim>::
      put_properties(const perplexcpp::MinimizeResult &result, 
	             ArrayView<double>::iterator &iter) const
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	std::vector<double> residue_composition = result.composition;
	for (double camount : residue_composition) {
	  *iter = camount;
	  iter++;
	}

	for (std::string pname : this->phase_names) 
	  this->put_phase_properties(result, pname, iter);

	if (this->extract_melt) {
	  const perplexcpp::Phase melt { perplexcpp::find_phase(result.phases, "liquid") };

	  this->put_extracted_melt_properties(melt, iter);

	  if (melt.volume_frac > this->melt_extraction_threshold)
	    for (unsigned int c = 0; c < px.n_composition_components; c++) {
	      residue_composition[c] -= melt.composition_ratio[c] * melt.n_moles;	    
	      if (residue_composition[c] < 0)
		residue_composition[c] = 0;
	    }
	}
      }


      template <int dim>
      void
      PerplexComposition<dim>::
      put_zero_properties(ArrayView<double>::iterator &iter) const
      {
	const auto& px = perplexcpp::Wrapper::get_instance();

	// Bulk composition.
	for (unsigned int c = 0; c < px.n_composition_components; c++) {
	  *iter = 0.0;
	  iter++;
	}

	// Phase properties.
	for (std::string pname : this->phase_names) {
	  for (PhaseProperty pprop : this->phase_properties) {
	    switch (pprop) {
	      case PhaseProperty::composition:
		for (unsigned int c = 0; c < px.n_composition_components; c++) {
		  *iter = 0.0;
		  iter++;
		}
		break;
	      default:
		*iter = 0.0;
		iter++;
	    }
	  }
	}

	// Extracted melt properties.
	if (this->extract_melt) {
	  for (PhaseProperty pprop : this->phase_properties) {
	    switch (pprop) {
	      case PhaseProperty::composition:
		for (unsigned int c = 0; c < px.n_composition_components; c++) {
		  *iter = 0.0;
		  iter++;
		}
		break;
	      default:
		*iter = 0.0;
		iter++;
	    }
	  }
	}
      }


      template <int dim>
      void 
      PerplexComposition<dim>::
      put_extracted_melt_properties(const perplexcpp::Phase &melt,
				    ArrayView<double>::iterator &it) const
      {
	for (PhaseProperty pprop : this->phase_properties) 
	{
	  switch (pprop)
	  {
	    case PhaseProperty::composition:
	      for (double camount : melt.composition_ratio) {
		*it += camount;
		it++;
	      }
	      break;
	    case PhaseProperty::n_moles:
	      *it += melt.n_moles;
	      it++;
	      break;
	    case PhaseProperty::molar_fraction:
	      *it += melt.molar_frac;
	      it++;
	      break;
	    case PhaseProperty::volume_fraction:
	      *it += melt.volume_frac;
	      it++;
	      break;
	    case PhaseProperty::weight_fraction:
	      *it += melt.weight_frac;
	      it++;
	      break;
	  }
	}
      }


      template <int dim>
      void 
      PerplexComposition<dim>::
      put_phase_properties(const perplexcpp::MinimizeResult &result,
			   const std::string &phase_name,
			   ArrayView<double>::iterator &it) const
      {
	const perplexcpp::Phase phase { perplexcpp::find_phase(result.phases, phase_name) };

	for (PhaseProperty pprop : this->phase_properties) {
	  switch (pprop) {
	    case PhaseProperty::composition:
	      for (double camount : phase.composition_ratio) {
		*it = camount;
		it++;
	      }
	      break;
	    case PhaseProperty::n_moles:
	      *it = phase.n_moles;
	      it++;
	      break;
	    case PhaseProperty::molar_fraction:
	      *it = phase.molar_frac;
	      it++;
	      break;
	    case PhaseProperty::volume_fraction:
	      *it = phase.volume_frac;
	      it++;
	      break;
	    case PhaseProperty::weight_fraction:
	      *it = phase.weight_frac;
	      it++;
	      break;
	  }
	}
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
