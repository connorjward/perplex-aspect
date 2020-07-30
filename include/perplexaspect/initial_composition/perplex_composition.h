/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#ifndef _perplexaspect_initial_composition_perplex_composition_h
#define _perplexaspect_initial_composition_perplex_composition_h


#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * ???
     */
    template <int dim>
    class PerplexComposition : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double 
	initial_composition(const Point<dim> &position,
                            const unsigned int compositional_index) const override;
    };
  }
}


#endif

