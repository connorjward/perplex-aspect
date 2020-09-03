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


#ifndef _perplexaspect_postprocess_perplex_cache_statistics_h
#define _perplexaspect_postprocess_perplex_cache_statistics_h


#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that compute some statistics about the Perple_X cache usage.
     */
    template <int dim>
    class PerplexCacheStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        std::pair<std::string,std::string>
        execute(TableHandler &statistics) override;
    };
  }
}


#endif
