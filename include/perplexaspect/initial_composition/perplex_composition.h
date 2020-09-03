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


#ifndef _perplexaspect_initial_composition_perplex_composition_h
#define _perplexaspect_initial_composition_perplex_composition_h


#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * Class that provides initial composition information for the `perplex melt simple'
     * material model.
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


      private:

	/**
	 * Calculate the porosity (volume fraction of melt) using the material model.
	 */
	double
	calc_porosity(const Point<dim> position) const;
    };
  }
}


#endif

