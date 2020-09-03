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


#include <perplexaspect/postprocess/perplex_cache_statistics.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PerplexCacheStatistics<dim>::execute(TableHandler &statistics)
    {
      perplexcpp::ResultCache &cache = perplexcpp::Wrapper::get_instance().get_cache();

      AssertThrow(cache.capacity > 0, 
	          ExcMessage("You are trying to record the cache statistics "
		             "but the cache is empty."));

      // Compute the sum over all processors.
      const unsigned int total_hits = 
	Utilities::MPI::sum(cache.get_n_hits(), this->get_mpi_communicator());
      const unsigned int total_misses = 
	Utilities::MPI::sum(cache.get_n_misses(), this->get_mpi_communicator());
      const double hit_rate = (double) total_hits / (total_hits + total_misses);

      // Reset cache counters.
      cache.reset_counters();

      statistics.add_value("Cache hits", total_hits);
      statistics.add_value("Cache misses", total_misses);
      statistics.add_value("Cache hit rate", hit_rate);

      statistics.set_precision("Cache hit rate", 3);

      std::ostringstream output;
      output.precision(4);

      output << total_hits << ", " << total_misses << ", " << hit_rate;

      return std::pair<std::string,std::string>("Cache hits / misses / hit rate: ",
                                                 output.str());
    }
  }
}


// Explicit instantiations.
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PerplexCacheStatistics,
                                  "perplex cache statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the Perple_X wrapper cache usage.")
  }
}
