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


#include <perplexaspect/perplex_melt_global.h>

#include <perplexcpp/wrapper.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;


    namespace
    {
      double
      sum(const std::vector<double>& vec)
      {
	double sum = 0.0;
	for (double x : vec)
	  sum += x;
	return sum;
      }


      void
      normalize(std::vector<double>& vec)
      {
	const double sum_ = std::abs(sum(vec));
	for (double& x : vec)
	  x /= sum_;
      }
    }



    template <int dim>
    void
    PerpleXMelt<dim>::initialize()
    {
      AssertThrow(this->include_melt_transport(), ExcMessage("Need this param."));
      AssertThrow(this->get_parameters().use_operator_splitting,
		  ExcMessage("The material model ``Perple_X melt'' can only be used "
			     "with operator splitting."));
      AssertThrow(this->include_melting_and_freezing, ExcMessage("Need this param."));
    }


    template <int dim>
    void
    PerpleXMelt<dim>::evaluate(const MaterialModelInputs<dim> &in, 
	                  MaterialModelOutputs<dim> &out) const
    {
      // Retain the complicated logic of the MeltGlobal class without duplicating it.
      MeltGlobal<dim>::evaluate(in, out);
      
      ReactionRateOutputs<dim> *reaction_rate_out = 
	out.template get_additional_output<ReactionRateOutputs<dim>>();

      // TODO Explain why.
      if (reaction_rate_out == nullptr)
	return;

      // DEBUG
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
                  const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
      const double peridotite_change = reaction_rate_out->reaction_rates[0][peridotite_idx];
      const double porosity_change = reaction_rate_out->reaction_rates[0][porosity_idx];
      if (peridotite_change > 0 || porosity_change > 0) 
      {
	std::cout << "periodite change: " << peridotite_change << std::endl;
	std::cout << "porosity change: " << porosity_change << std::endl;
      }
      // END DEBUG


      /* if (this->get_timestep_number() == 0) */
	/* return; */

      auto& perplex_wrapper = perplexcpp::Wrapper::get_instance();

      if(in.requests_property(MaterialProperties::reaction_terms))
      {
	for (unsigned int i = 0; i < in.n_evaluation_points(); ++i) {
	  const double melt_frac = 
	    in.composition[i][this->introspection().compositional_index_for_name("porosity")];

	  Assert(melt_frac >= 0 && melt_frac <= 1, ExcInternalError("The porosity should be between 0 and 1"));

	  auto melt_composition = get_composition(in.composition[i], "melt");
	  auto residue_composition = get_composition(in.composition[i], "residue");

	  /* std::cout << "Melt composition" << std::endl; */
	  /* for (auto c : melt_composition) */
	  /*   std::cout << c << " "; */
	  /* std::cout << std::endl; */
	    
	  /* std::cout << "Residue composition" << std::endl; */
	  /* for (auto c : residue_composition) */
	  /*   std::cout << c << " "; */
	  /* std::cout << std::endl; */

	  /* std::cout << "melt composition" << std::endl; */
	  /* std::cout << "sum: " << sum(melt_composition) << std::endl; */
	  /* for (auto c : melt_composition) */
	  /*   std::cout << c << std::endl; */
	  /* std::cout << std::endl; */

	  /* std::cout << "residue composition" << std::endl; */
	  /* std::cout << "sum: " << sum(residue_composition) << std::endl; */
	  /* for (auto c : residue_composition) */
	  /*   std::cout << c << std::endl; */
	  /* std::cout << std::endl; */

	  if (std::abs(sum(melt_composition)) > 1e-8)
	    normalize(melt_composition);

	  if(std::abs(sum(residue_composition)) > 1e-8)
	    normalize(residue_composition);

	  // Determine the composition to pass into Perple_X.
	  std::vector<double> composition;
	  for (unsigned int c = 0; 
	       c < perplexcpp::Wrapper::get_instance().n_composition_components; 
	       ++c)
	  {
	    composition.emplace_back(melt_frac * melt_composition[c] +
				     (1-melt_frac) * residue_composition[c]);
	  }

	  /* std::cout << "Composition" << std::endl; */
	  /* for (auto c : composition) */
	  /*   std::cout << c << " "; */
	  /* std::cout << std::endl; */

	  /* if (sum(composition) < 0) { */
	  /*   for (auto c : composition) */
	  /*     std::cout << c << " " << std::endl; */
	  /*   exit(0); */ 
	  /* } */
      
	  /* std::cout << "min pressure: " << perplex_wrapper.min_pressure << std::endl; */
	  /* std::cout << "max pressure: " << perplex_wrapper.max_pressure << std::endl; */
	  /* std::cout << "min temperature: " << perplex_wrapper.min_temperature << std::endl; */
	  /* std::cout << "max temperature: " << perplex_wrapper.max_temperature << std::endl; */
	  /* std::cout << "sum comp: " << sum(composition) << std::endl; */


	  /* std::cout << "pressure: " << in.pressure[i] << std::endl; */
	  /* std::cout << "temperature: " << in.temperature[i] << std::endl; */

	  /* for (auto c : composition) */
	  /*   std::cout << c << " "; */
	  /* std::cout << std::endl; */
	    
	  if (in.pressure[i] < perplex_wrapper.min_pressure || 
	      in.pressure[i] > perplex_wrapper.max_pressure ||
	      in.temperature[i] < perplex_wrapper.min_temperature ||
	      in.temperature[i] > perplex_wrapper.max_temperature ||
	      sum(composition) < 1e-8)
	  {
	    continue;
	  }

	  const auto result = perplex_wrapper.minimize(in.pressure[i], 
	                                               in.temperature[i], 
						       composition);

	  Assert(result.n_moles > 1 - 1e8 && result.n_moles < 1 + 1e8,
	         ExcInternalError("The number of moles should be normalized to 1"));

	  perplexcpp::Phase melt = perplexcpp::find_phase(result.phases, "liquid");

	  for (unsigned int c = 0; c < perplex_wrapper.n_composition_components; ++c) 
	  {
	    const std::string comp_name = perplex_wrapper.composition_component_names[c];

	    const unsigned int melt_comp_idx = 
	      this->introspection().compositional_index_for_name("melt_" + comp_name);

	    const unsigned int residue_comp_idx = 
	      this->introspection().compositional_index_for_name("residue_" + comp_name);

	    /* std::cout << "melt.composition[" << c << "]: " << melt.composition[c] << std::endl; */
	    /* std::cout << "result.composition[" << c << "]: " << result.composition[c] << std::endl; */
	    /* std::cout << "in.comp[i][melt]: " << in.composition[i][melt_comp_idx] << std::endl; */
	    /* std::cout << "in.comp[i][res]: " << in.composition[i][residue_comp_idx] << std::endl; */

	    const double melt_comp_change = 
	      melt.composition[c] - in.composition[i][melt_comp_idx];
	    const double residue_comp_change =
	      result.composition[c] - melt.composition[c] - in.composition[i][residue_comp_idx];

	    /* std::cout << "melt_comp_change: " << melt_comp_change << std::endl; */
	    /* std::cout << "residue_comp_change: " << residue_comp_change << std::endl; */

	    /* std::cout << "melting_time_scale: " << this->melting_time_scale << std::endl; */


	    reaction_rate_out->reaction_rates[i][melt_comp_idx] = 
	      melt_comp_change / this->melting_time_scale;

	    reaction_rate_out->reaction_rates[i][residue_comp_idx] = 
	      residue_comp_change / this->melting_time_scale;

	    // DEBUG
	    if (melt_comp_change > 0 || residue_comp_change > 0) 
	    {
	      std::cout << "melt change: " << melt_comp_change / this->melting_time_scale << std::endl;
	      std::cout << "residue change: " << residue_comp_change / this->melting_time_scale << std::endl;
	    }
	    //END DEBUG
	  }
	}
      }
    } 


    // parse parameters
    // check that the compositional fields are named correctly (maybe do in an init 
    // step instead to make sure that ordering is good.
    template <int dim>
    void
    PerpleXMelt<dim>::declare_parameters(ParameterHandler &prm)
    {
      MeltGlobal<dim>::declare_parameters(prm);

      prm.enter_subsection("Material model");
      {
	/* prm.enter_subsection("Melt global"); */
	/* { */
	/*   prm.declare_entry("Melting time scale for operator splitting", "1e3", */
	/*                     Patterns::Double (0.), */ 
	/* 		    "See other one."); */
	/* } */
	/* prm.leave_subsection(); */

	prm.enter_subsection("Perple_X melt");
	{
	  prm.declare_entry("Data directory", 
			    ".", 
			    Patterns::DirectoryName(),
			    "The location of the Perple_X data files.");

	  prm.declare_entry("Problem definition file", 
			    "", 
			    Patterns::FileName(),
			    "The name of the PerpleX .dat file in use.",
			    true);
	}
	prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    PerpleXMelt<dim>::parse_parameters(ParameterHandler &prm)
    {
      MeltGlobal<dim>::parse_parameters(prm);

      prm.enter_subsection("Material model");
      {
	prm.enter_subsection("Melt global");
	{
	  // This needs to be included because the attribute is private and cannot otherwise be
	  // accessed by this class.
	  include_melting_and_freezing = prm.get_bool("Include melting and freezing");
	  melting_time_scale = prm.get_double("Melting time scale for operator splitting");
	}
	prm.leave_subsection();
	prm.enter_subsection("Perple_X melt");
	{
	      const auto data_dirname = prm.get("Data directory");
	      const auto problem_filename = prm.get("Problem definition file");

	      AssertThrow(Utilities::fexists(data_dirname + "/" + problem_filename),
			  ExcMessage("The Perple_X problem file could not be found."));

	      perplexcpp::Wrapper::initialize(problem_filename, data_dirname);
	}
	prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    std::vector<double>
    PerpleXMelt<dim>::get_composition(const std::vector<double>& comp_fields, 
	                              const std::string& name) const
    { 
      std::vector<double> composition;
      for (std::string comp_name : 
	   perplexcpp::Wrapper::get_instance().composition_component_names) 
      {
	const unsigned int idx = 
	  this->introspection().compositional_index_for_name(name + "_" + comp_name);

	composition.emplace_back(comp_fields[idx]);
      }
      return composition;
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PerpleXMelt,
                                   "perplex melt",
                                   "Description here.")
  }
}
