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


#include <perplexaspect/utilities.h>

#include <aspect/utilities.h>
#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace PerplexUtils
  {
    void
    declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Perple_X settings");
      {
	prm.declare_entry(
	  "Data directory", 
	  ".", 
	  Patterns::DirectoryName(),
	  "The location of the Perple_X data files."
	);

	prm.declare_entry(
	  "Problem filename", 
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
    }


    void 
    parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Perple_X settings");
      {
	const std::string dirname         = prm.get("Data directory");
	const std::string fname           = prm.get("Problem filename");
	const unsigned int cache_capacity = prm.get_integer("Cache capacity");
	const double cache_rtol           = prm.get_double("Cache tolerance");

	AssertThrow(Utilities::fexists(dirname + "/" + fname),
		    ExcMessage("The Perple_X problem file could not be found."));

	perplexcpp::Wrapper::initialize(fname, dirname, cache_capacity, cache_rtol);
      }
      prm.leave_subsection();
    }


    void
    put_melt_composition(const perplexcpp::MinimizeResult &result,
	                 const double porosity,
	                 std::vector<double> &melt_composition)
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      Assert(porosity >= 0, ExcInternalError("The porosity must be non-negative"));
      Assert(melt_composition.size() == px.n_composition_components, 
	     ExcInternalError("The melt composition vector is the wrong size."));

      const perplexcpp::Phase melt = perplexcpp::find_phase(result.phases, "liquid");

      // If no melt is present return a vector of zeros.
      if (porosity < 1e-8 || melt.n_moles < 1e-8) {
	for (unsigned int c = 0; c < px.n_composition_components; c++)
	  melt_composition[c] = 0.0;
	return;
      }

      // Initial calculations for readability later on.
      double phase_densities_sum = 0.0;
      for (auto p : result.phases)
	phase_densities_sum += p.density;

      double melt_mol_mass = 0.0;
      for (unsigned int c = 0; c < px.n_composition_components; ++c)
	melt_mol_mass += melt.composition_ratio[c] * px.composition_molar_masses[c];

      double melt_comp_ratio_sum = 0.0;
      for (double melt_comp_ratio : melt.composition_ratio)
	melt_comp_ratio_sum += melt_comp_ratio;

      double bulk_weight = 0.0;
      for (unsigned int c = 0; c < px.n_composition_components; ++c)
	bulk_weight += result.composition[c] * px.composition_molar_masses[c];

      const double melt_weight_frac = 
	porosity * melt.density / phase_densities_sum;

      // Populate the melt composition vector.
      for (unsigned int c = 0; c < px.n_composition_components; c++)
      {
	const double cmelt_weight_frac = 
	  melt.composition_ratio[c] * px.composition_molar_masses[c] / melt_mol_mass;

	// The weight fraction of the composition component in the melt is the weight fraction 
	// of the melt multiplied by the weight fraction of the component inside the melt.
	const double cbulk_weight_frac = melt_weight_frac * cmelt_weight_frac;

	// Turn the weight fraction into an actual weight.
	const double cbulk_weight = cbulk_weight_frac * bulk_weight;

	// Figure out the number of moles of the component.
	const double cbulk_mol = cbulk_weight / px.composition_molar_masses[c];

	melt_composition[c] = cbulk_mol;
      }
    }


    double
    limit_pressure(const double pressure)
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      if (pressure < px.min_pressure)
	return px.min_pressure;
      else if (pressure > px.max_pressure)
	return px.max_pressure;
      else
	return pressure;
    }


    double
    limit_temperature(const double temperature)
    {
      const auto& px = perplexcpp::Wrapper::get_instance();

      if (temperature < px.min_temperature)
	return px.min_temperature;
      else if (temperature > px.max_temperature)
	return px.max_temperature;
      else
	return temperature;
    }
  }
}
