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
  }
}
