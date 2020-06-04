/*
   Copyright (C) 2011 - 2018 by the authors of the ASPECT code.
   Copyright (C) 2020 Connor Ward

   This file is part of PerplexASPECT.

   PerplexASPECT is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PerplexASPECT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with PerplexASPECT.  If not, see <https://www.gnu.org/licenses/>.
   */

#include "material_model.h"

#include <deal.II/base/multithread_info.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class PhaseAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
	PhaseAdditionalOutputs(const std::vector<std::string> &phase_names, const unsigned int n_points);

	std::vector<double> get_nth_output(const unsigned int idx) const override;

	std::vector<std::vector<double>> n_moles;
    };

    template <int dim>
    PhaseAdditionalOutputs<dim>::
    PhaseAdditionalOutputs(const std::vector<std::string> &phase_names, 
			   const unsigned int n_points)
    : 
    NamedAdditionalMaterialOutputs<dim>(phase_names) {
      for (unsigned int i = 0; i < phase_names.size(); i++)
	n_moles.push_back(std::vector<double>(n_points, 0));
    }

    template <int dim>
    std::vector<double>
    PhaseAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      /* AssertIndexRange(idx, n_moles.size()); */
      return n_moles[idx];
    }










    template <int dim>
    void
    MyPerpleXLookup<dim>::initialize()
    {
      std::cout << "entering initialize" << std::endl;
      AssertThrow(dealii::MultithreadInfo::is_running_single_threaded(),
                  ExcMessage("The PerpleXLookup MaterialModel only works in single threaded mode (do not use -j)!"));

      wrapper.init(perplex_file_name.c_str()); // this line initializes meemum

      // create phase index map
      const std::vector<std::string> names = wrapper.solution_phase_names();
      for (unsigned int i = 0; i < names.size(); i++)
	phase_idx_map[names[i]] = i;

      for (auto it = phase_idx_map.begin(); it != phase_idx_map.end(); ++it) {
	std::cout << it->first << ", " << it->second << '\n';
      }
      std::cout << "leaving initialize" << std::endl;
    }

    template <int dim>
    bool
    MyPerpleXLookup<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    double
    MyPerpleXLookup<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    void
    MyPerpleXLookup<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      std::cout << "entering evaluate" << std::endl;
      /* Instead of evaluating at every quadrature point per cell,
       * we here average the P, T and X values, and evaluate once.
       * This is much quicker than evaluating at all quadrature
       * points, and if the grid is fine, it should be a reasonable
       * approximation
       */

 //    std::vector<double> wtphases(p_size_phases);
 //    std::vector<double> cphases(p_size_phases * p_size_components);
 //    std::vector<char> namephases(p_size_phases * p_pname_len);
 //    std::vector<double> sysprop(p_size_sysprops);

      int phaseq_dbg = 0;

      unsigned int n_quad = in.n_evaluation_points(); // number of quadrature points in cell
      unsigned int n_comp = in.composition[0].size(); // number of components in rock

      /* const double average_temperature = std::min(max_temperature, */
      /*                                             std::max(min_temperature, */
      /*                                                      (accumulate( in.temperature.begin(), in.temperature.end(), 0.0) / */
      /*                                                       n_quad))); */
      /* const double average_pressure = std::min(max_pressure, */
      /*                                          std::max(min_pressure, */
      /*                                                   (accumulate( in.pressure.begin(), in.pressure.end(), 0.0) / */
      /*                                                    n_quad))); */

      std::vector<double> comp;
      comp.resize(n_comp);

      // average composition over all quadrature points
      for (unsigned int c=0; c<n_comp; ++c)
        {
          for (unsigned int i=0; i<n_quad; ++i)
            {
              comp[c] += in.composition[i][c];
              out.reaction_terms[i][c] = 0.0;
            }
          comp[c] /= (double)n_quad;
        }

      /* MinimizeResult res = wrapper.minimize(in.pressure[0], in.temperature[0], in.composition[0]); */
      MinimizeResult res = wrapper.minimize(20000, 1500, in.composition[0]);

      for (unsigned int i=0; i<n_quad; ++i)
      {
	out.viscosities[i] = eta;
	out.thermal_conductivities[i] = k_value;
	out.densities[i] = 1000;
	out.specific_heat[i] = 100;
	out.thermal_expansion_coefficients[i] = 100;
	out.compressibilities[i] = 100;
      }

      // store phase information
      PhaseAdditionalOutputs<dim> *phase_out = out.template get_additional_output<PhaseAdditionalOutputs<dim>>();

      if (phase_out != nullptr) {
	for (Phase phase : res.phases) {
	  unsigned int phase_idx;
	  try {
	    /* const unsigned int phase_idx = phase_idx_map.at(phase.name); */
	    phase_idx = phase_idx_map.at(phase.name);
	  } catch(std::out_of_range ex) {
	    std::cout << "out of range exception triggered" << std::endl;
	    std::cout << "looking for: " << phase.name << std::endl;
	    for (auto it = phase_idx_map.begin(); it != phase_idx_map.end(); ++it) {
	      std::cout << it->first << ", " << it->second << '\n';
	    }
	    exit(1);
	  }

	  for (unsigned int i = 0; i < n_quad; i++)
	    phase_out->get_nth_output(phase_idx)[i] = phase.n_moles;
	}
      }
      std::cout << "leaving evaluate" << std::endl;
    }


    template <int dim>
    void
    MyPerpleXLookup<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("PerpleX lookup model v2");
        {

          prm.declare_entry ("PerpleX input file name", "rock.dat",
                             Patterns::Anything (),
                             "The name of the PerpleX input file (should end with .dat).");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the viscosity $\\eta$. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Minimum material temperature", "0.",
                             Patterns::Double (0),
                             "The value of the minimum temperature used to query PerpleX. "
                             "Units: $\\si{K}$.");
          prm.declare_entry ("Maximum material temperature", "6000.",
                             Patterns::Double (0),
                             "The value of the maximum temperature used to query PerpleX. "
                             "Units: $\\si{K}$.");
          prm.declare_entry ("Minimum material pressure", "1.e5",
                             Patterns::Double (0),
                             "The value of the minimum pressure used to query PerpleX. "
                             "Units: $Pa$.");
          prm.declare_entry ("Maximum material pressure", "1.e12",
                             Patterns::Double (0),
                             "The value of the maximum pressure used to query PerpleX. "
                             "Units: $Pa$.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MyPerpleXLookup<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("PerpleX lookup model v2");
        {
          perplex_file_name   = prm.get ("PerpleX input file name");
          eta                 = prm.get_double ("Viscosity");
          k_value             = prm.get_double ("Thermal conductivity");
          min_temperature     = prm.get_double ("Minimum material temperature");
          max_temperature     = prm.get_double ("Maximum material temperature");
          min_pressure        = prm.get_double ("Minimum material pressure");
          max_pressure        = prm.get_double ("Maximum material pressure");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      this->model_dependence.density = NonlinearDependence::temperature
                                       | NonlinearDependence::pressure
                                       | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::temperature
                                               | NonlinearDependence::pressure
                                               | NonlinearDependence::compositional_fields;
      this->model_dependence.specific_heat = NonlinearDependence::temperature
                                             | NonlinearDependence::pressure
                                             | NonlinearDependence::compositional_fields;
    }

    template <int dim>
    void
    MyPerpleXLookup<dim>::
    create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      std::cout << "entering create_additional_named_outputs" << std::endl;
      if (out.template get_additional_output<PhaseAdditionalOutputs<dim>>() == nullptr) {
	const std::vector<std::string> names = wrapper.solution_phase_names();
	std::cout << names.size() << std::endl;
	/* std::cout << names[0] << std::endl; */
	const unsigned int n_points = out.n_evaluation_points();

	out.additional_outputs.push_back(std_cxx14::
	                                 make_unique<MaterialModel::
					 PhaseAdditionalOutputs<dim>>(names, n_points));
      }
      std::cout << "leaving create_additional_named_outputs" << std::endl;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MyPerpleXLookup,
                                   "my perplex lookup",
                                   "A material model that has constant values "
                                   "for viscosity and thermal conductivity, and "
                                   "calculates other properties on-the-fly using "
                                   "PerpleX meemum. Compositional fields correspond "
                                   "to the individual components in the order given "
                                   "in the PerpleX file.")
  }
}
