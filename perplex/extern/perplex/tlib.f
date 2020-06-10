c----------------------------------------------------------------------

c TLIB - a library of subprograms called by the PERPLEX programs.

c Copyright (C) 2018 James A D Connolly

c This file is part of Perple_X.

c Perple_X is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2, or (at your option)
c any later version.

c Perple_X is distributed in the hope that it will be useful,
c but without any warranty; without even the implied warranty of
c merchantability or fitness for a particular purpose.  See the
c GNU General Public License for more details.

c You should have received a copy of the GNU General Public License
c along with Perple_X (file license.txt). If not see
c <http://www.gnu.org/licenses/>.

c----------------------------------------------------------------------

      subroutine vrsion (n)
c----------------------------------------------------------------------
c a version stamp for each executable
c----------------------------------------------------------------------
      implicit none

      integer n

      write (n,'(/,a,//,a)') 
     *      'Perple_X version 6.9.0, source updated May 27, 2020.',

     *      'Copyright (C) 1986-2020 James A D Connolly '//
     *      '<www.perplex.ethz.ch/copyright.html>.'

      end

      logical function chksol (new)
c----------------------------------------------------------------------
c check that the version flag in the solution model file is consistent.
c This creates headaches, but is used to prevent old versions of Perple_X 
c from crashing while reading a new solution model format
c----------------------------------------------------------------------
      implicit none

      character*3 new

      if (new.eq.'682'.or.new.eq.'683'.or.new.eq.'688'.or.
     *    new.eq.'685'.or.new.eq.'687') then

         call error (3,0d0,0,new)

      else if (new.eq.'008'.or.new.eq.'011'.or.new.eq.'670'.or.
     *         new.eq.'672'.or.new.eq.'673'.or.new.eq.'674'.or.
     *         new.eq.'675'.or.new.eq.'676'.or.new.eq.'678'.or.
     *         new.eq.'679'.or.new.eq.'689'.or.new.eq.'690') then 

         chksol = .true.

      else 

         chksol = .false.

      end if 

      end 

      subroutine redop1 (output,opname)
c----------------------------------------------------------------------
c redop1 - redop1 looks for the perplex_option.dat file, if it finds
c the option file it scans the file for keywords and sets options
c accordingly. also sets machine precision dependent (wmach) and 
c fractional constants. also checks for goofball parameter choices.

c iam - indicates calling program 1 - vertex
c                                 2 - meemum
c                                 3 - werami
c                                 all other values no output

c special internal values for lopt, iopt, nopt

c               lop_28-30 
c               iop_28-30
c               nop_28-30
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, jer, i, loopx, loopy, ibeg, iend, ik1, ik21

      logical output

      character*3 key*22, val, nval1*12, nval2*12,
     *            nval3*12,opname*100,strg*40,strg1*40

      double precision dnan, res0, r2

      integer jtest,jpot
      common/ debug /jtest,jpot

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer iam
      common/ cst4 /iam

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      double precision wmach(9)
      common/ cstmch /wmach
c----------------------------------------------------------------------
c                                 periodic fractions
      r13 = 1d0/3d0
      r23 = 2d0/3d0
      r43 = 4d0/3d0
      r59 = 5d0/9d0
c                                 seismic speed conversion, this shouldn't be
c                                 hardwired.
      units = dsqrt(1d1)/1d1
c                                 -------------------------------------
c                                 loop to find machine precision (mainly
c                                 for nag)
      r1 = 1d-12
      r2 = 0d0

      do
         if (1d0+r1.eq.1d0) exit
         r2 = r1 
         r1 = r1/2d0
      end do 

      if (r2.eq.0d0) call error (72,r1,i,
     *                         'starting precision for r1 < zero')
c                                 wmach(1-2,5-6,9) do not have the 
c                                 standard BLAS values, additionally
c                                 BLAS routines may assume a dimension
c                                 of 15 for wmach.
c                                 as currently set wmach(1,2) are
c                                 function precision and and optimality
c                                 tolerance.
      wmach(1) = r2**0.9d0
      wmach(2) = r2**0.8d0
c                                 relative precision (eps)
      wmach(3) = r2
      wmach(4) = dsqrt(r2)
      wmach(5) = 1d0 + r2
c                                 largest number
      wmach(7) = huge(0d0)
      wmach(8) = dsqrt(wmach(7))
      wmach(9) = max(1d0/wmach(4),1d2)
c                                 solution composition zero and one
      r2 = 1d2 * r2

      zero = r2
      r1 = 1d0 + r2
      one = 1d0 - r2
c                                 -------------------------------------
c                                 default option values:
c                                 reserved for temporary use:
      do i = 28, 30
         nopt(i) = 1d0
         lopt(i) = .false.
         iopt(i) = 0
      end do 
c                                 -------------------------------------
c                                 minimum replicate label distance
      nopt(4) = 0.025
c                                 speciation_factor
      nopt(5) = 1d2
c                                 vrh_weighting keyword
      nopt(6) = 0.5d0
c                                 bad_number keyword
      nopt(7) = dnan()
c                                 zero_mode (<0 off)
      nopt(9) = 1d-6
c                                 set zero threshold for fractionation calculations
      if (icopt.eq.7.or.icopt.eq.9.or.icopt.eq.12) nopt(9) = zero
c                                 tolerance below which a component is considered to 
c                                 be zero during fractionation
      nopt(11) = 1d-6
c                                 quench temperature (K)
      nopt(12) = 0d0
c                                 initial resolution for adaptive 
c                                 refinement
      nopt(13) = 1d0/16d0
c                                 solvus_tolerance
      nopt(8) = 1.5*nopt(13)
c                                 solvus_tolerance_II
      nopt(25) = nopt(13)
c                                 compositional resolution for conformal
c                                 subdivision
      nopt(14) = 2d-3
c                                 perturbation to eliminate pseudocompound
c                                 degeneracies
      nopt(15) = 5d-3
c                                 poisson ratio to be use for calculating
c                                 missing shear moduli, see iopt(15)
      nopt(16) = 0.35d0
c                                 increase in resolution for adaptive minimization 
      nopt(17) = 3d0
c                                 increase in resolution for composition and mixed variable diagram calculations
      nopt(18) = 1d1
c                                 increase in resolution for Schreinemakers diagram calculations   
      nopt(19) = 3d0 
c                                 T_melt cutoff 
      nopt(20) = 873d0
c                                 iteration keyword 1
      nopt(21) = 2d0
c                                 global reach factor
      nopt(23) = 0d0
c                                 finite_difference_p threshold for finite difference estimates
      nopt(26) = 1d4
c                                 finite_difference_p fraction for first order difference estimates
      nopt(27) = 1d-3
c                                 fd_expansion_factor is the factor by which finite difference
c                                 increments increase for higher order derivatives
      nopt(31) = 2d0
c                                 fractionation_lower_threshold
      nopt(32) = 0d0
c                                 fractionation_upper_threshold
      nopt(33) = 0d0
c                                 aq_vapor_epsilon
      nopt(34) = 1d0
c                                 -------------------------------------
c                                 composition_phase
      iopt(2) = 0 
      valu(2) = 'mol'
c                                 porportions
      iopt(3) = 0 
      valu(3) = 'vol'
c                                 interpolation keyword
      iopt(4) = 2
      valu(4) = 'on '
c                                 autorefine, 2 - automatic, 1 - manual, 0 - no
      iopt(6) = 2
      valu(6) = 'aut'
c                                 subdivision model, 0 - solution model
c                                 1 - cartesian, 2 - stretch
      iopt(13) = 0 
      valu(13)  = 'off'
c                                 poisson ratio switch fpr calculating
c                                 missing shear moduli
      iopt(16) = 1
      valu(15) = 'on '
c                                 compare local and max disorder state for o-d models
      iopt(17) = 1
      valu(17) = 'on '
c                                 assume linear boundaries within a cell during gridded minimization
      iopt(18) = 1
      valu(18) = 'on '
c                                 seismic data output for WERAMI/MEEMUM, 0 - none, 1 - some, 2 - all
      iopt(14) = 1
      valu(14) = 'som'
c                                 reach_increment_switch 0 - off, 1 - on (for auto-refine), 2 - all
      iopt(20) = 1 
      valu(20) = 'on'
c                                 conditional for MEEMUM
      if (iam.eq.2) then 
         valu(20) = 'all'
         iopt(20) = 2
      end if 
c                                 speciation_max_it - for speciation calculations
      iopt(21) = 100
c                                 aq_bad_results 
c                                       0 - err - treat as bad result (optimization error)
c                                       1 - 101 - cease iteration at solute component saturation
c                                       2 - 102
c                                       3 - 103
c                                      99 - ign - ignore 102/103 conditions and cease as in 101.
      iopt(22) = 0
      valu(5) = 'err'
c                                 for infiltration/fractionation calculations set default to 101
      if (icopt.gt.6) then 
         iopt(22) = 1
         valu(5) = '101'
      end if 
c                                 solution_names 0 - model, 1 - abbreviation, 2 - full
      iopt(24) = 0
      valu(22) = 'mod'
c                                 hyb_h2o - eos to be used for pure h2o, 0-2, 4-5, 6-7
      iopt(25) = 4
c                                 hyb_co2 - eos to be used for pure co2, 0-4, 6-7
      iopt(26) = 4
c                                 hyb_ch4 - eos to be used for pure ch4, 0-1, 6-7
      iopt(27) = 0

c     iopt(28-30)                 reserved as debug options iop_28 - iop_30

c                                 refinement_points_II renamed refinement_points
      iopt(31) = 5
c                                 maximum number of aqueous species
      iopt(32) = 20
c                                 aq_lagged_iterations
      iopt(33) = 0
c                                 interim_results, 1 - auto, 0 - off, 2 - man
      iopt(34) = 1
      valu(34) = 'aut'
c                                 keep_max
      iopt(52) = 20000
c                                 -------------------------------------
c                                 closed or open compositional space
      lopt(1) = .true.
c                                 lopt(2) set by WERAMI
c     lopt(2)
c                                 hard_limits for solution model refinement
      lopt(3) = .false.
      valu(16) = 'off'
c                                 Anderson-Gruneisen correction
      lopt(4) = .false.
c                                 auto_exclude 
      lopt(5) = .true.
c                                 melt_is_fluid
      lopt(6) = .false.
c                                 set locally
c     lopt(7) 
c                                 approx_alpha
      lopt(8) = .true.
c                                 automatic solvus tolerance
      lopt(9) = .true.
c                                 pseudocompound_file
      lopt(10) = .false.
c                                 auto_refine_file
      lopt(11) = .true.
c                                 option_list_files
      lopt(12) = .false.
c                                 automatic solvus_tolerance_II
      lopt(13) = .true.
c                                 logarithimic P
      lopt(14) = .false.
c                                 spreadsheet format
      lopt(15) = .true.
c                                 averaging scheme
      lopt(16) = .true.
      valu(19) = 'VRH'
c                                 use explicit bulk modulus when available
      lopt(17) = .true.
c                                 refine_bad_nodes -> not used
      lopt(18) = .true. 
c                                 pause_on_error
      lopt(19) = .true.
c                                 poisson_test (reject results with poisson < 0)
      lopt(20) = .false.
c                                 species_output
      lopt(21) = .true. 
c                                 composition_constant
      lopt(22) = .false.
c                                 composition_system (true = wt)
      lopt(23) = .true. 
      valu(21) = 'wt '
c                                 output endmember gibbs energies (werami/meemum)
      lopt(24) = .false. 
c                                 aq_output output back-calculated solute speciation
      lopt(25) = .true.
c                                 aq_solvent_composition (true = molar)
      lopt(26) = .true.
      valu(26) = 'y'
c                                 aq_solute_composition (true = molal)
      lopt(27) = .true.
      valu(27) = 'm'
c                                 lop_28, lop_29, lop_30 temporary options
c     lopt(28-30)
c                                 warning_ver637, werami interpolation
      lopt(31) = .true.
c                                 aq_lagged_speciation
      lopt(32) = .false.
c                                 output_iteration_g
      lopt(33) = .false.
c                                 output_iteration_details
      lopt(34) = .false.
c                                 fractionation threshold, set by nopt(32)/nopt(33)
      lopt(35) = .false.
c                                 aq_oxide_components
      lopt(36) = .false.
c                                 allow null phases
      lopt(37) = .false.
c                                 non_linear_switch
      lopt(38) = .false.
c                                 refine_endmembers
      lopt(39) = .false.
c                                 automatic specification of refinement_points_II
      lopt(40) = .true.
c                                 absolute (amounts)
      lopt(41) = .false.
c                                 cumulative (amounts)
      lopt(42) = .false.
c                                 reject_negative_sites
      lopt(43) = .true.
c                                 aq_ion_H+ 
      lopt(44) = .true.
c                                 fancy_cumulative_modes
      lopt(45) = .false.
c                                 aq_solvent_solvus
      lopt(46) = .true.
c                                 sample_on_grid 
      lopt(48) = .true.
c                                 refinement_switch
      lopt(49) = .false.
c                                 seismic_data_file
      lopt(50) = .true.
c                                 structural formula options
      lopt(51) = .true.
c                                 keep_auto
      lopt(52) = .true.
c                                 final resolution, auto-refine stage
      rid(2,2) = 1d-3
c                                 final resolution, exploratory stage
      rid(2,1) = 1d-2
c                                 if meemum set auto-refine vale
      if (iam.eq.2) rid(2,1) = rid(2,2)
c                                 initialize mus flag lagged speciation
      mus = .false.
c                                 -------------------------------------
c                                 for gridded minimization:
c                                 # nodes in i direction
      grid(1,1) = 40 
      grid(1,2) = 40
c                                 # nodes in j direction
      grid(2,1) = 40 
      grid(2,2) = 40
c                                 # of levels
      grid(3,1) = 1
      grid(3,2) = 4
c                                 1d fractionation path
      grid(4,1) = 40 
      grid(4,2) = 150
c                                 -------------------------------------
c                                 for schreinemakers etc:
c                                 max variance 
      grid(5,1) = 1
      grid(5,2) = 99
c                                 default increment (relative)
      rid(1,1) = 0.1d0
      rid(1,2) = 0.025d0
c                                 reaction format
      ifull = 0 
      valu(7) = 'min'
c                                 reaction list 
      jtest = 0 
      valu(9) = 'off'
c                                 console msgs
      imsg = 0
      valu(8) = 'on '
c                                 efficiency level
      isec = 3
c                                 short print
      io3p = 1
      valu(10) = 'on '
c                                 print dependent potentials
      jpot = 0
      valu(11) = 'on '
c                                 -------------------------------------
c                                 look for file
      open (n8, file = opname, iostat = jer, status = 'old')
c                                 if no option file (jer.ne.0) use defaults
      ier = jer
c                                 read cards to end of 
c                                 option file
      do while (ier.eq.0) 

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 ier ne 0 = eof
         if (ier.ne.0) exit 
c                                 cycle on default specification
         if (strg.eq.'default') cycle
c                                 if here we should have a keyword and
c                                 value
         if (key.eq.'composition'.or.key.eq.'composition_phase') then 
c                                 phase composition key
            if (val.eq.'wt') then
               iopt(2) = 1
               valu(2) = val 
            end if 

         else if (key.eq.'aq_species') then

            read (strg,*) iopt(32)

         else if (key.eq.'auto_exclude') then 

            if (val.ne.'T') lopt(5) = .false.

         else if (key.eq.'aq_lagged_iterations') then

            read (strg,*) iopt(33)

         else if (key.eq.'aq_output') then

            if (val.ne.'T') lopt(25) = .false.

         else if (key.eq.'aq_lagged_speciation') then 

            if (val.eq.'T') lopt(32) = .true.

         else if (key.eq.'aq_oxide_components') then 

            if (val.eq.'T') lopt(36) = .true.

         else if (key.eq.'aq_ion_H+') then 

             if (val.eq.'F') lopt(44) = .false.

         else if (key.eq.'fancy_cumulative_modes') then 

             if (val.eq.'T') lopt(45) = .true.

         else if (key.eq.'null_phase') then 

            if (val.eq.'T') lopt(37) = .true.

         else if (key.eq.'non_linear_switch') then 

            if (val.eq.'T') lopt(38) = .true.

         else if (key.eq.'output_iteration_G') then 

            if (val.eq.'T') lopt(33) = .true.

         else if (key.eq.'output_iteration_detai') then 

            if (val.eq.'T') lopt(34) = .true.

         else if (key.eq.'aq_bad_results') then 

            if (val.eq.'err') then 
c                                 abort on any hint of trouble
               iopt(22) = 0
            else if (val.eq.'101') then 
c                                 continue on solute undersaturation (unwise)
               iopt(22) = 1
            else if (val.eq.'102') then 
c                                 continue if pure solvent coexists with immiscible impure solvent
               iopt(22) = 2
            else if (val.eq.'103') then
c                                 abort if pure solvent is stable
               iopt(22) = 3
            else if (val.eq.'ign') then 
               iopt(22) = 99
            end if

            valu(5) = val

         else if (key.eq.'refine_endmembers') then 

            if (val.eq.'T') lopt(39) = .true.

         else if (key.eq.'absolute') then 

            if (val.eq.'T') lopt(41) = .true.

         else if (key.eq.'cumulative') then 

            if (val.eq.'T') lopt(42) = .true.

         else if (key.eq.'reject_negative_sites') then 

            if (val.eq.'F') lopt(43) = .false.

         else if (key.eq.'aq_solvent_composition') then

            if (val.ne.'y') then
               lopt(26) = .false.
               valu(26) = 'm'
            end if 

         else if (key.eq.'aq_solute_composition') then

            if (val.ne.'m') then
               lopt(27) = .false.
               valu(27) = 'y'
            end if 

         else if (key.eq.'hybrid_EoS_H2O') then

            read (strg,*) iopt(25)

            if (iopt(25).lt.0.or.iopt(25).gt.7.or.iopt(25).eq.3) then 
               write (*,1180) strg,key
               call errpau
            end if 

         else if (key.eq.'hybrid_EoS_CO2') then

            read (strg,*) iopt(26)

            if (iopt(26).lt.0.or.(iopt(26).gt.4.and.iopt(26).ne.7)) then 
               write (*,1180) strg,key
               call errpau
            end if

         else if (key.eq.'hybrid_EoS_CH4') then

            read (strg,*) iopt(27)

            if (iopt(27).lt.0.or.(iopt(27).gt.1.and.iopt(27).ne.7)) then 
               write (*,1180) strg,key
               call errpau
            end if 

         else if (key.eq.'proportions') then 
c                                 phase proportion key
c                                 volume is default
            if (val.eq.'wt') then
               iopt(3) = 1
            else if (val.eq.'mol') then 
               iopt(3) = 2 
            end if
 
            valu(3) = val

         else if (key.eq.'interpolation') then 
c                                 interpolation key
            if (val.ne.'on ') then
               iopt(4) = 0
               valu(4) = val
            else 
               read (nval1,*) iopt(4)
            end if 

         else if (key.eq.'bounds') then 
c                                 
            if (val.eq.'HS'.or.val.eq.'hs')  then
              lopt(16) = .false.
              valu(19) = 'HS '
            end if 

         else if (key.eq.'vrh/hs_weighting'.or.
     *            key.eq.'vrh_weighting') then 
c                                 vrh/hs weighting key
            read (strg,*) nopt(6)

         else if (key.eq.'bad_number') then
c                                 bad number key 
            if (val.eq.'NaN'.or.val.eq.'nan') then
               nopt(7) = dnan()
            else 
               read (strg,*) nopt(7)
            end if 

         else if (key.eq.'solution_names') then 

              if (val.eq.'abb') then 
                 iopt(24) = 1
              else if (val.eq.'ful') then 
                 iopt(24) = 2
              end if 

              valu(22) = val

         else if (key.eq.'solvus_tolerance') then 

            if (val.ne.'aut') then 
               lopt(9) = .false.
               read (strg,*) nopt(8)
            end if 

         else if (key.eq.'solvus_tolerance_II') then 

            if (val.ne.'aut') then 
               lopt(13) = .false.
               read (strg,*) nopt(25)
            end if 

         else if (key.eq.'speciation_factor') then 

            read (strg,*) nopt(5)
            if (nopt(5).lt.10) nopt(5) = 1d1

         else if (key.eq.'zero_bulk') then
c                                 zero_bulk key
            read (strg,*) nopt(11)

         else if (key.eq.'aq_vapor_epsilon') then
c                                 "vapor" threshold
            read (strg,*) nopt(34)

         else if (key.eq.'aq_max_molality') then

c             obsolete
         else if (key.eq.'Tisza_test') then

c             to be implemented

         else if (key.eq.'warning_ver637') then
c                                  allow for solvent immiscisibiliy
            if (val.eq.'F') lopt(31) = .false.

         else if (key.eq.'aq_solvent_solvus') then
c                                  allow for solvent immiscisibiliy
            if (val.eq.'F') lopt(46) = .false.

         else if (key.eq.'interim_results') then
c                                  output interim results (VERTEX/PSSECT/WERAMI)
            if (val.eq.'off') then 
               iopt(34) = 0
            else if (val.eq.'man') then 
               iopt(34) = 2
            end if

            valu(34) = val 

         else if (key.eq.'sample_on_grid') then
c                                  sample on computational grid (WERAMI)
            if (val.eq.'F') lopt(48) = .false.

         else if (key.eq.'zero_mode') then
c                                 zero_mode key
            read (strg,*) nopt(9)

         else if (key.eq.'iteration'.or.
     *            key.eq.'resolution_factor') then
c                                 how fast resolution improves with iteration
            read (val,*) nopt(21)

         else if (key.eq.'initial_resolution') then
c                                 initial_resolution key
            read (strg1,'(40a)') chars(1:40)
            ibeg = 1
            call readfr (nopt(13),ibeg,iend,40,ier)
            if (ier.ne.0) call error (72,nopt(1),iopt(1),key//
     *                                 'has an invalid value.')
c                                 initial resolution
            if (nopt(13).ge.1d0.or.nopt(13).lt.0d0) then
               nopt(13) = 1d0/16d0
               write (*,1050)
            end if
c                                  686+ read second value as arf value
            ibeg = iend + 1

            call readfr (nopt(17),ibeg,iend,40,ier)

            if (ier.ne.0) then
c                                  special backward compatibility msg
               write (*,1010) opname
               call errpau

            end if 

            if (nopt(17).ge.1d0.or.nopt(17).lt.0d0) then
               nopt(17) = 1d0/48d0
               write (*,1050)
            end if

            nopt(17) = nopt(13)/nopt(17)
            nopt(18) = nopt(17)
            nopt(19) = nopt(17)

         else if (key.eq.'final_resolution') then
c                                 final_resolution keys 
            read (strg,*) rid(2,1)
            read (nval1,*,iostat=ier) rid(2,2)
            if (ier.ne.0) call error (72,r1,i,'the final_resolution '//
     *                               'option requires two values')

         else if (key.eq.'fd_expansion_factor') then 

            read (strg,*) nopt(31)

         else if (key.eq.'finite_difference_p') then 
c                                 p threshold
            read (strg,*) nopt(26)
c                                 p fraction
            read (nval1,*) nopt(27)

         else if (key.eq.'global_reach_increment') then
          
            read (strg,*) nopt(23) 

         else if (key.eq.'seismic_output') then 
c                                 seismic data output WERAMI/MEEMUM/FRENDLY
            valu(14) = val

            if (val.eq.'non') then 
               iopt(14) = 0
            else if (val.eq.'all') then
               iopt(14) = 2
            else
               valu(14) = 'som'
            end if

         else if (key.eq.'refinement_points_II'.or.
     *            key.eq.'refinement_points') then
c                                 refinement points
            if (val.ne.'aut') then 
               read (strg,*) iopt(31)
               lopt(40) = .false.
            end if

         else if (key.eq.'refinement_switch') then
c                                 during iteration allow metastable 
c                                 refinement points for stable solutions.
            if (val.ne.'F') lopt(49) = .true.

         else if (key.eq.'seismic_data_file') then
 
            if (val.ne.'T') lopt(50) = .false.

         else if (key.eq.'structural_formulae') then

            if (val.ne.'T') lopt(51) = .false. 

         else if (key.eq.'keep_auto') then

            if (val.ne.'T') lopt(52) = .false. 

         else if (key.eq.'keep_max') then 

            read (strg,*) iopt(52)

         else if (key.eq.'max_aq_species_out') then 
c                                 max number of aq species output for
c                                 back-calculated and lagged speciation
            read (strg,*) iopt(32)

         else if (key.eq.'reach_increment_switch') then 
c                                 reach_increment_switch
             valu(20) = val

             if (val.eq.'off') then 
                iopt(20) = 0
             else if (val.eq.'on ') then
                iopt(20) = 1
             else if (val.eq.'all') then 
                iopt(20) = 2
             end if 

         else if (key.eq.'stretch_factor') then
c                                 stretch_factor key = b - 1       
            read (strg,*) nopt(14)

         else if (key.eq.'subdivision_override') then 
c                                 subdivision override key
            valu(13) = val

            if (val.eq.'lin') then
               iopt(13) = 1
            else if (val.eq.'str') then
               iopt(13) = 2
            end if 

         else if (key.eq.'auto_refine') then
c                                 autorefine
            valu(6) = val

            if (val.eq.'off') then
               iopt(6) = 0
            else if (val.eq.'man') then
               iopt(6) = 1
            end if 

         else if (key.eq.'auto_refine_factor_I') then
c   
            read (strg,*) nopt(17)

         else if (key.eq.'auto_refine_factor_II') then
c   
            read (strg,*) nopt(18)

         else if (key.eq.'auto_refine_factor_III') then
c   
            read (strg,*) nopt(19)

         else if (key.eq.'speciation_max_it') then
c   
            read (strg,*) iopt(21)

         else if (key.eq.'x_nodes') then
c                                 number of x nodes at level 1 before autorefine
            read (strg,*) grid(1,1)
c                                 number of x nodes for autorefine
            read (nval1,*,iostat=ier) grid(1,2)

            if (ier.ne.0) call error (72,r1,i,'the x_nodes option '//
     *                               'requires two values')

         else if (key.eq.'y_nodes') then
c                                 number of y nodes at level 1
            read (strg,*) grid(2,1)
c                                 number of y nodes for autorefine
            read (nval1,*,iostat=ier) grid(2,2)

            if (ier.ne.0) call error (72,r1,i,'the y_nodes option '//
     *                               'requires two values')

         else if (key.eq.'grid_levels') then 
c                                 number of grid levels before autorefine
            read (strg,*) grid(3,1)
c                                 number of grid levels for autorefine
            read (nval1,*,iostat=ier) grid(3,2)

            if (ier.ne.0) call error (72,r1,i,'the grid_levels option'//
     *                               ' requires two values')

         else if (key.eq.'1d_path') then 
c                                 number of grid points for 1d path before autorefine
            read (strg,*) grid(4,1)
c                                 number of grid points for 1d path for autorefine
            read (nval1,*,iostat=ier) grid(4,2)

            if (ier.ne.0) call error (72,r1,i,'the 1d_path option '//
     *                               'requires two values')

         else if (key.eq.'variance') then 
c                                 max variance of traced equilibria before autorefine
            read (strg,*) grid(5,1)
c                                 max variance of traced equilibria for autorefine
            read (nval1,*,iostat=ier) grid(5,2)

            if (ier.ne.0) call error (72,r1,i,'the variance option '//
     *                               'requires two values')

         else if (key.eq.'increment') then 
c                                 default exploratory relative increment    
            read (strg,*) rid(1,1)  
c                                 default autorefine relative increment
            read (nval1,*,iostat=ier) rid(1,2)

            if (ier.ne.0) call error (72,r1,i,'the increment option '//
     *                               'requires two values')

         else if (key.eq.'reaction_format') then 

            valu(7) = val

            if (val.eq.'ful') then 
               ifull = 1
            else if (val.eq.'sto') then 
               ifull = 2
            else if (val.eq.'S+V') then 
               ifull = 3
            else if (val.eq.'eve') then
               ifull = 4
            end if 
  
         else if (key.eq.'console_messages') then 
            
            if (val.eq.'off') then 
               valu(8) = val
               imsg = 1
            end if 

         else if (key.eq.'reaction_list') then

            if (val.eq.'on ') then 
               valu(9) = val
               jtest = 3 
            end if

         else if (key.eq.'efficiency') then 

            read (strg,*) isec

            if (isec.lt.1.or.isec.gt.5) isec = 3 

         else if (key.eq.'short_print') then 

            if (val.eq.'off') then 
               io3p = 0
               valu(10) = val
            end if 

         else if (key.eq.'dependent_potentials') then 

            if (val.eq.'off') then 
               jpot = 1
               valu(11) = val
            end if

         else if (key.eq.'hard_limits') then 

            if (val.eq.'on ') then 
               lopt(3) = .true.
               valu(16) = val
            end if

         else if (key.eq.'T_stop') then 
c                                 equilibrium cutoff T (K)    
            read (strg,*) nopt(12)

         else if (key.eq.'T_melt') then 
c                                 cutoff T (K) for melt endmember stability    
            read (strg,*) nopt(20)

         else if (key.eq.'fractionation_hi_limit') then 
c                                 upper fractionation threshold
            read (strg,*) nopt(33)

         else if (key.eq.'fractionation_lo_limit') then 
c                                 lower fractionation threshold
            read (strg,*) nopt(32)

         else if (key.eq.'order_check') then 
c                                 compare local and max disorder state for o-d models
            if (val.eq.'off') then 
               iopt(17) = 0
               valu(17) = 'off'
            end if 

         else if (key.eq.'linear_model') then   
c                                 assume linear boundaries within a cell during gridded minimization
            if (val.eq.'off') then 
               iopt(18) = 0
               valu(18) = 'off'
            end if 

         else if (key.eq.'closed_c_space') then
 
            if (val.ne.'T') lopt(1) = .false. 

         else if (key.eq.'pause_on_error') then
 
            if (val.ne.'T') lopt(19) = .false. 

         else if (key.eq.'poisson_test') then
 
            if (val.ne.'F') lopt(20) = .true. 

         else if (key.eq.'species_output') then
 
            if (val.ne.'T') lopt(21) = .false. 

         else if (key.eq.'composition_constant') then
 
            if (val.ne.'F') lopt(22) = .true. 

         else if (key.eq.'composition_system') then
 
            if (val.ne.'wt') then
               lopt(23) = .false. 
               valu(21) = val
            end if

         else if (key.eq.'species_Gibbs_energies') then
 
            if (val.ne.'F') lopt(24) = .true.

         else if (key.eq.'logarithmic_p') then
 
            if (val.eq.'T') lopt(14) = .true.

         else if (key.eq.'spreadsheet') then
 
            if (val.eq.'F') lopt(15) = .false.

         else if (key.eq.'Anderson-Gruneisen') then

            if (val.eq.'T') lopt(4) = .true.

         else if (key.eq.'approx_alpha') then

            if (val.eq.'F') lopt(8) = .false.

         else if (key.eq.'melt_is_fluid') then 

            if (val.eq.'T') lopt(6) = .true.

         else if (key.eq.'pc_perturbation') then
c                                 perturbation to eliminate pseudocompound degeneracies  
            read (strg,*) nopt(15)

         else if (key.eq.'pseudocompound_file') then

            if (val.eq.'T') lopt(10) = .true.

         else if (key.eq.'auto_refine_file') then

            if (val.eq.'T') lopt(11) = .true.

         else if (key.eq.'option_list_files') then

            if (val.eq.'T') lopt(12) = .true.

         else if (key.eq.'explicit_bulk_modulus') then 

            if (val.eq.'F') lopt(17) = .false.

         else if (key.eq.'poisson_ratio') then 
c                                 handle missing shear moduli
            if (val.eq.'on ') then
               read (nval1,*) nopt(16)
            else if (val.eq.'off') then 
               valu(15) = val
               iopt(16) = 0
            else if (val.eq.'all') then 
               read (nval1,*) nopt(16)
               valu(15) = val
               iopt(16) = 2
            end if   
          
         else if (key.eq.'lop_28') then
c                                 reserved values for debugging, etc
            if (val.eq.'T') lopt(28) = .true.
         else if (key.eq.'lop_29') then
            if (val.eq.'T') lopt(29) = .true.
         else if (key.eq.'lop_30') then
            if (val.eq.'T') lopt(30) = .true.
         else if (key.eq.'iop_28') then
            read (strg,*) iopt(28) 
         else if (key.eq.'iop_29') then
            read (strg,*) iopt(29) 
         else if (key.eq.'iop_30') then
            read (strg,*) iopt(30) 
         else if (key.eq.'nop_28') then
            read (strg,*) nopt(28) 
         else if (key.eq.'nop_29') then
            read (strg,*) nopt(29) 
         else if (key.eq.'nop_30') then
            read (strg,*) nopt(30) 
         else if (key.ne.'|') then

            call error (72,nopt(1),iopt(1),key//' is not a valid Perpl'
     *                 //'e_X option file keyword and must be deleted '
     *                 //'or corrected.')

         end if

      end do 

      close (n8)
c                                 -------------------------------------
c                                 computation dependent options
c                                 -------------------------------------
c                                 automatic specification of metastable
c                                 refinement points
      if (lopt(40)) iopt(31) = icp + 2
c                                 turn dependent potentials on if auto_keep
      if (lopt(52)) then
         jpot = 0
         valu(11) = 'on'
      end if 
c                                 always allow null phases if not CONVEX
      if (iam.eq.15) lopt(37) = .true. 
c                                 write and optional file choices
      if (iam.ne.14) then 
         if (jer.ne.0) then 
            write (*,1120) opname
         else 
            write (*,1130) opname
         end if
      end if 

      if (iam.eq.1.or.iam.eq.15) then 
c                                 vertex only files:
         if (icopt.eq.1.or.icopt.eq.3) then 

            if (jtest.eq.3) then 
               call mertxt (tfname,prject,'_reaction_list.txt',0)
               open (n6,file=tfname)
            else 
               tfname = 'not requested'
            end if 

            write (*,1170) tfname

         end if 
c                                 auto refine summary
         if (iopt(6).ne.0.and.io9.eq.0) then
 
            if (lopt(11)) then 
               call mertxt (tfname,prject,'_auto_refine.txt',0)
            else 
               tfname = 'not requested'
            end if 

            write (*,1150) tfname

         end if

      end if

      if (iam.lt.3) then

         if (lopt(50)) then 
            call mertxt (tfname,prject,'_seismic_data.txt',0)
         else 
            tfname = 'not requested'
         end if 

         write (*,'(a)') 'Writing seismic data options to: '//tfname

      end if 
c                                 pseudocompound glossary
      if (io9.eq.0.and.(iam.lt.3.or.iam.eq.15)) then
 
         if (lopt(10)) then 
            call mertxt (tfname,prject,'_pseudocompound_glossary.txt',0)
         else 
            tfname = 'not requested'
         end if 

         write (*,1140) tfname

      end if 
c                                 computational options this is redundant
      if (iam.lt.4.or.iam.eq.15) then 
         if (lopt(12)) then 
            if (iam.eq.1) then           
               call mertxt (tfname,prject,'_VERTEX_options.txt',0)
            else if (iam.eq.2) then 
               call mertxt (tfname,prject,'_MEEMUM_options.txt',0)
            else if (iam.eq.3) then 
               call mertxt (tfname,prject,'_WERAMI_options.txt',0)
            else if (iam.eq.15) then 
               call mertxt (tfname,prject,'_CONVEX_options.txt',0)
            end if
         else 
            tfname = 'not requested'
         end if 
 
         write (*,1160) tfname

      end if 
c                                 -------------------------------------
c                                 dependent parameters and error traps:
c                                 fractionation theshold flag
      if (nopt(33).gt.nopt(32)) lopt(35) = .true. 

      if (nopt(21).lt.2d0) then 
         write (*,1040)
         nopt(21) = 2d0
      end if

      if (iopt(31).gt.k5+2.or.iopt(31).lt.1) then 
         write (*,1090) icp + 2
         iopt(31) = icp + 2
      end if 
c                                 stretching parameters
      if (nopt(14).le.0d0) then 
         write (*,1060)
         nopt(14) = 2d-3
      end if 
c                                 grid parameters
      do i = 1, 2

         if (grid(3,i).le.0.or.grid(3,i).gt.l8) grid(3,i) = 4

         loopy = (grid(2,i)-1) * 2**(grid(3,i)-1) + 1
         loopx = (grid(1,i)-1) * 2**(grid(3,i)-1) + 1

         if (loopy.gt.l7) then 
            call warn (92,nopt(1),loopy,'y_node')
            grid(2,i) = (l7 - 1)/2**(grid(3,i)-1) + 1
         end if 

         if (loopx.gt.l7) then 
            call warn (92,nopt(1),loopx,'x_node')
            grid(1,i) = (l7 - 1)/2**(grid(3,i)-1) + 1
         end if 

         if (grid(4,i).gt.l7) then 
            call warn (92,nopt(1),loopx,'1dpath')
            grid(4,i) = l7 - 1
         end if  

         if (grid(5,i).lt.1) then 
            call warn (113,rid(1,i),grid(5,i),'VARIAN')
            grid(5,i) = 1
         end if 

         if (rid(1,i).lt.1d-2) then 
            call warn (114,rid(1,i),i,'INPUT1')
         end if 

      end do
c                                 --------------------------------------
c                                 program/computation specific settings
c                                 set autorefine factor
      if (icopt.eq.1) then
         nopt(17) = nopt(19)
      else if (icopt.le.3) then 
         nopt(17) = nopt(18)
      end if 
c                                 compute resolution/number of iterations 
      nopt(24) = 2d0*nopt(21)/(1d0 + nopt(21))
c                                 effective initial resolution
      if (nopt(13).eq.0d0) then 
c                                 user wants to use solution model values, set
c                                 initial resolution to a representative value
          rid(3,1) = 0.1d0
      else 
          rid(3,1) = nopt(13)
      end if 

      rid(3,2) = rid(3,1)/nopt(17)

      do i = 1, 2

         if (rid(2,i).gt.rid(3,i)) then
c                                  requested final resolution > requested initial resolution
c                                  grid(6,i) is the iteration counter
            grid(6,i) = 0
c                                  speciation tolerance, later set to nopt(5)
            rid(5,i) = rid(3,i)/nopt(5)

         else
          
            grid(6,i) = 1

            do 

               res0 = rid(3,i)*nopt(24)/nopt(21)**grid(6,i)
c                                 actual final resolution
               rid(4,i) = res0

               if (res0.lt.rid(2,i)) exit
c                                  grid(6,i) is the iteration counter
               grid(6,i) = grid(6,i) + 1

            end do
c                                  real final resolution is res0
c                                  speciation tolerance, later set to nopt(5)
            rid(5,i) = res0/nopt(5)

         end if

      end do 
c                                 --------------------------------------
c                                 output
      if (((iam.eq.1.or.iam.eq.15).and.output).or.iam.eq.2
     *                                        .or.iam.eq.3
     *                                        .or.iam.eq.5) then
c                                 console
         call outopt (6)

         if (lopt(12).and.iam.ne.5) then 
c                                 file version, create the file name 
            if (iam.eq.1) then           
               call mertxt (tfname,prject,'_VERTEX_options.txt',0)
            else if (iam.eq.2) then 
               call mertxt (tfname,prject,'_MEEMUM_options.txt',0)
            else if (iam.eq.3) then 
               call mertxt (tfname,prject,'_WERAMI_options.txt',0)
            else if (iam.eq.15) then 
               call mertxt (tfname,prject,'_CONVEX_options.txt',0)
            end if 

            open (n8, file = tfname)
            call outopt (n8)
            close (n8) 

            write (*,1000) tfname

         end if 

      end if 
c                                 convert speciation factor to speciation tolerance
c                                 for fluids and frendly assuming final resolution 
c                                 value, for meemum and vertex this will be over-ridden
c                                 according to whether the program is in the exploratory
c                                 or auto-refine stage in setau1. 
      nopt(5) = rid(5,2)
c                                 -------------------------------------
c                                 recalculate parameters:
c                                 proportionality constant for shear modulus
      nopt(1) = nopt(16)
      nopt(16) = 1.5d0*(1d0-2d0*nopt(16))/(1d0+nopt(16))
c                                 -------------------------------------
c                                 check dynamic memory allocation:
      if (iam.lt.3) then

c        write (*,*) k21, k19
c                                 recommended value for k21, to be absolutely 
c                                 safe multiply by 2.
         ik21 = iopt(52)*icp

         if (k21.le.ik21) then 
c                                 consequent value for k1
            ik1 = (memory - (1 + k32)*ik21)/(2*k31+k32+2)

            write (*,'(a,i8,2a,i7,a/,a,i7,a,/)') 
     *            '**warning ver072** parameter k21 (',k21,') is below',
     *            ' its recommended value (',ik21,')','decrease k1 to ',
     *            ik1,' or reduce keep_max to avoid potential problems'

         end if 

         if (k21.le.100) then

            ik1 = (memory - (1 + k32)*50000)/(2*k31+k32+2)

            write (*,'(a,i7,2a/,a,i7,a,/)') 
     *            '**error ver072** k21 (',k21,') is too far below the',
     *            ' minimum safe value (~10000)','decrease k1 to < ',
     *            ik1,' and recompile Perple_X'

            call errpau

         end if

      end if

1000  format ('Context specific options are echoed in: ',a,/)
1010  format (/,'ERROR: reading option file: ',a50,/,
     *          'most likely you are using a pre-6.8.6 option',
     *          ' file with only one',/,'initial_resolution value. Corr'
     *         ,'ect the error by adding a second value,',/,'typically',
     *          ' 1/3 the first (exploratory stage) value.')
1040  format (/,'Warning: iteration keyword value must be ',
     *         ' > 1',/,'iteration will be',
     *         ' assigned its default value [2].',/)
1050  format (/,'Warning: initial_resolution values must be ',
     *         '< 1',/,'the keyword will be',
     *         ' assigned its default values.',/)
1060  format (/,'Warning: stretch_factor cannot be < 0 ',
     *          'the keyword will be assigned its default value',/)
1070  format (/,'Warning: initial_resolution values must be ',
     *         '< 1',/,'the keyword will be',
     *         ' assigned its default value (',i2,').',/)
1090  format (/,'Warning: refinement_points must be ',
     *         ' >1 and <k5+2',/,'refinement_points will be',
     *         ' assigned its default value [',i2,'].',/)
1120  format (/,'Warning: the Perple_X option file: ',a,/,
     *       'was not found, default option values will be used.',/) 
1130  format (/,'Reading Perple_X options from: ',a)
1140  format ('Writing pseudocompound glossary to: ',a)
1150  format ('Writing auto refine summary to: ',a)
1160  format ('Writing Perple_X option summary to: ',a)
1170  format ('Writing complete reaction list to: ',a)
1180  format (/,'Error: value ',a,' is invalid for Perple_X option ',
     *        'keyword ',a,/,'see www.perplex.ch/perplex_options.html ',
     *        'for a list of valid values',/)
      end 

      subroutine outopt (n)
c----------------------------------------------------------------------
c outopt - for program IAM, outputs context specific options to LUN N,
c called by redop1  
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, siz, n

      character nval1*12, text(14)*1, numb*5, nval2*12

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      logical oned
      common/ cst82 /oned

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 version
      if (n.ne.6) call vrsion (n)
c                                 generic blurb
      if (iam.eq.1) then 
         write (n,1000) 'VERTEX'
      else if (iam.eq.2) then 
         write (n,1000) 'MEEMUM'
      else if (iam.eq.3) then 
         write (n,1000) 'WERAMI'
      else if (iam.eq.5) then 
         write (n,1000) 'FRENDLY'
      else if (iam.eq.15) then 
         write (n,1000) 'CONVEX'
      end if 

      if (iam.le.2.or.iam.eq.15) then 
c                                 VERTEX/MEEMUM:
c                                 solvus tolerance text
         if (lopt(9)) then 

           nval1 = 'aut    '

         else 

           call numtxt (nopt(8),text,siz)
           write (nval1,'(14a)') (text(i),i=1,siz)

         end if

         if (lopt(13)) then 

           nval2 = 'aut    '

         else 

           call numtxt (nopt(25),text,siz)
           write (nval2,'(14a)') (text(i),i=1,siz)

         end if 

         if (iam.eq.1.or.iam.eq.15) write (n,1015) valu(6)
c                                 context specific parameters:
         if (icopt.le.3.and.(iam.eq.1.or.iam.eq.15)) then 
c                                 non-adaptive calculations
c                                 reaction format and lists
            if (icopt.gt.0) write (n,1160) grid(5,1),grid(5,2),rid(1,1),
     *                    rid(1,2),isec,valu(7),valu(9),valu(8),valu(10)
         else 
c                                 adaptive optimization
            write (n,1180) rid(2,1),rid(2,2),int(nopt(21)),
     *                     iopt(31),k5,lopt(49),lopt(52),iopt(52),
     *                     k21/10,nval2,
     *                     int(nopt(23)),valu(20),nopt(9)
c                                 gridding parameters
            if (iam.eq.1.and.icopt.eq.5.and.oned) then
c                                 1d multilevel grid
               write (n,1190) grid(2,1),grid(2,2),l7,
     *                  (grid(2,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(2,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(3,1),grid(3,2),l8

            else if (iam.eq.1.and.icopt.eq.5) then
c                                 2d multilevel grid
               write (n,1200) grid(1,1),grid(1,2),l7,
     *                  (grid(1,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(1,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(2,1),grid(2,2),l7,
     *                  (grid(2,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(2,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(3,1),grid(3,2),l8,valu(18)

            else if (iam.eq.1.and.icopt.eq.7) then 
c                                 1d fractionation grid
                write (n,1210) grid(4,1),grid(4,2),l7

            end if 
c                                 closed or open composition space
            if (iam.eq.1.and.icont.gt.1) write (n,1220) lopt(1)

         end if 
c                                 generic subdivision parameters:
         if (icopt.eq.0) then
            numb = '1/160'
         else
            numb = '1/48'
         end if 

         write (n,1010) nopt(13),nopt(13)/nopt(17),numb,nopt(14),
     *                  lopt(38),valu(13),valu(16),lopt(39)

         if (iam.eq.15)  write (n,1011) nopt(15)
c                                 generic thermo parameters:
         write (n,1012) nval1,nopt(12),nopt(20),valu(17),
     *                  lopt(8),lopt(4),nopt(5),iopt(21),
     *                  iopt(25),iopt(26),iopt(27),valu(5),
     *                  lopt(32),lopt(44),lopt(36),lopt(46),nopt(34)
c                                 for meemum add fd stuff
         if (iam.eq.2) write (n,1017) nopt(31),nopt(26),nopt(27)

         if (iam.eq.1.or.iam.eq.15) then 
c                                 vertex output options, dependent potentials
c                                 pause_on_error
            write (n,1013) valu(11),lopt(19)
c                                 auto_exclude
            write (n,1234) lopt(5)
c                                 logarithmic_p, bad_number, interim_results
            if (iam.eq.1) write (n,1014) lopt(14),nopt(7),valu(34)

         end if 

      end if

      if (iam.eq.3) then 
c                                 WERAMI input/output options
         write (n,1230) lopt(25),iopt(32),l9,valu(26),valu(27),
     *                  lopt(15),lopt(14),nopt(7),lopt(22),valu(2),
     *                  valu(21),valu(3),lopt(41),lopt(42),lopt(45),
     *                  valu(4),lopt(6),valu(22),lopt(51),lopt(21),
     *                  lopt(24),
     *                  valu(14),lopt(19),lopt(20),valu(34),lopt(48)
         write (n,1234) lopt(5)
c                                 WERAMI info file options
         write (n,1241) lopt(12)       
c                                 WERAMI thermodynamic options
         write (n,1016) lopt(8),lopt(4),iopt(25),iopt(26),iopt(27)
         write (n,1017) nopt(31),nopt(26),nopt(27)

      else if (iam.eq.2) then 
c                                 MEEMUM input/output options
         write (n,1231) lopt(25),iopt(32),l9,valu(26),valu(27),
     *                  lopt(14),nopt(7),lopt(22),valu(2),
     *                  valu(21),valu(3),lopt(6),valu(22),lopt(51),
     *                  lopt(21),
     *                  lopt(24),valu(14),lopt(19),lopt(20)
         write (n,1234) lopt(5)

      else if (iam.eq.5) then 
c                                 FRENDLY input/output options
         write (n,1232) lopt(15),lopt(14),nopt(7),lopt(6),lopt(19),
     *                  .false.

      end if 
c                                 seismic property options
      if (iam.eq.2.or.iam.eq.3.or.iam.eq.5) write (n,1233) lopt(50),
     *         valu(19),
     *         nopt(6),lopt(17),valu(15),nopt(16),valu(14),lopt(20),
     *         .false.

      if (iam.eq.5) then 
c                                 FRENDLY thermo options
         write (n,1016) lopt(8),lopt(4),iopt(25),iopt(26),iopt(27)
         write (n,1017) nopt(31),nopt(26),nopt(27)

      end if 

      if (iam.le.2) then 
c                                 info file options
         write (n,1240) lopt(12),lopt(10)
         if (iam.eq.1.or.iam.eq.15) write (n,1250) lopt(11)
         if (iam.eq.1) write (n,'(4x,a,l1,10x,a)') 
     *                     'seismic_data_file      ',lopt(50),'[F] T;'//
     *                     ' echo seismic wavespeed options'
      end if 
c                                 resolution blurb
      if ((iam.le.2.or.iam.eq.15).and.nopt(13).gt.0d0) then

         write (n,1090) rid(4,1),rid(4,2)
         write (n,1100) grid(6,1),grid(6,2)

      end if 

      write (n,1020) 

1000  format (/,'Perple_X computational option settings for ',a,':',//,
     *      '    Keyword:               Value:     Permitted values ',
     *          '[default]:')

1010  format (/,2x,'Solution subdivision options:',//,
     *        4x,'initial_resolution:    ',/,
     *        4x,'  exploratory stage    ',f6.4,5x,
     *           '0->1 [1/16], 0 => off',/,
     *        4x,'  auto-refine stage    ',f6.4,5x,
     *           '0->1 [',a,'], 0 => off',/,
     *        4x,'stretch_factor         ',f6.4,5x,'>0 [2d-3]',/,
     *        4x,'non_linear_switch      ',l1,10x,'[F] T',/,
     *        4x,'subdivision_override   ',a3,8x,'[off] lin str',/,
     *        4x,'hard_limits            ',a3,8x,'[off] on',/,
     *        4x,'refine_endmembers      ',l1,10x,'[F] T')
1011  format (4x,'pc_perturbation        ',f6.4,5x,'[5d-3]')
c                                 generic thermo options
1012  format (/,2x,'Thermodynamic options:',//,
     *        4x,'solvus_tolerance       ',a7,4x,          
     *           '[aut] or 0->1; aut = automatic, 0 => ',
     *           'p=c pseudocompounds, 1 => homogenize',/,
     *        4x,'T_stop (K)             ',f6.1,5x,'[0]',/,
     *        4x,'T_melt (K)             ',f6.1,5x,'[873]',/,
     *        4x,'order_check            ',a3,8x,'off [on]',/,
     *        4x,'approx_alpha           ',l1,10x,'[T] F',/,
     *        4x,'Anderson-Gruneisen     ',l1,10x,'[F] T',/,
     *     4x,'speciation_factor      ',f6.0,5x,'>10 [100] speciation ',
     *           'precision = final resolution/speciation_factor',/,
     *        4x,'speciation_max_it      ',i4,7x,'[100]',/,
     *        4x,'hybrid_EoS_H2O         ',i4,7x,'[4] 0-2, 4-7',/,
     *        4x,'hybrid_EoS_CO2         ',i4,7x,'[4] 0-4, 7',/,
     *        4x,'hybrid_EoS_CH4         ',i4,7x,'[0] 0-1, 7',/,
     *        4x,'aq_bad_results         ',a3,8x,'[err] 101, 102, 103,',
     *                                           ' ignore',/,
     *        4x,'aq_lagged_speciation   ',l1,10x,'[F] T',/,
     *        4x,'aq_ion_H+              ',l1,10x,'[T] F => use OH-',/,
     *        4x,'aq_oxide_components    ',l1,10x,'[F] T',/,
     *        4x,'aq_solvent_solvus      ',l1,10x,'[T] F',/,
     *        4x,'aq_vapor_epsilon       ',f3.1,8x,'[1.]')
1013  format (/,2x,'Input/Output options:',//,
     *        4x,'dependent_potentials   ',a3,8x,'off [on]',/,
     *        4x,'pause_on_error         ',l1,10x,'[T] F')
1014  format (4x,'logarithmic_p          ',l1,10x,'[F] T',/,
     *        4x,'bad_number          ',f7.1,7x,'[NaN]',/,
     *        4x,'interim_results        ',a3,8x,'[auto] off manual')
1015  format (/,2x,'Auto-refine options:',//,
     *        4x,'auto_refine            ',a3,8x,'off manual [auto]')
c                                 thermo options for frendly
1016  format (/,2x,'Thermodynamic options:',//,
     *        4x,'approx_alpha           ',l1,10x,'[T] F',/,
     *        4x,'Anderson-Gruneisen     ',l1,10x,'[F] T',/,
     *        4x,'hybrid_EoS_H2O         ',i4,7x,'[4] 0-2, 4-7',/,
     *        4x,'hybrid_EoS_CO2         ',i4,7x,'[4] 0-4, 7',/,
     *        4x,'hybrid_EoS_CH4         ',i4,7x,'[0] 0-1, 7')
1017  format (4x,'fd_expansion_factor    ',f3.1,8x,'>0 [2]',/,
     *        4x,'finite_difference_p    ',d7.1,4x,'>0 [1d4]; ',
     *           'fraction = ',d7.1,3x,'[1d-2]')
1020  format (/,'To change these options see: ',
     *        'www.perplex.ethz.ch/perplex_options.html',/)
1090  format (/,2x,
     *        'Worst-case Cartesian compositional resolution (mol)',
     *        ': ',//,4x,'Exploratory stage: ',g11.3E1,/,
     *                4x,'Auto-refine stage: ',g11.3E1)
1100  format (/,2x,'Adapative minimization will be done with: ',
     *        //,3x,i2,' iterations in the exploratory stage',/,
     *           3x,i2,' iterations in the autorefine stage')
1160  format (/,2x,'Schreinemakers and Mixed-variable diagram ',
     *           'options:',//,
     *        4x,'variance               ',i2,' /',i2,5x,
     *           '[1/99], >0, maximum true variance',/,
     *        4x,'increment           ',f5.3,'/',f5.3,3x,
     *           '[0.1/0.025], ',
     *           'default search/trace variable increment',/,
     *        4x,'efficiency               ',i1,8x,'[3] >0 < 6',/,      
     *        4x,'reaction_format        ',a3,8x,'[min] ',
     *           'full stoichiometry S+V everything',/,
     *        4x,'reaction_list          ',a3,8x,'[off] on',/,
     *        4x,'console_messages       ',a3,8x,'[on] off',/,
     *        4x,'short_print_file       ',a3,8x,'[on] off')
1180  format (/,2x,'Free energy minimization options:',//,
     *        4x,'final_resolution:      ',/,
     *        4x,'  exploratory stage    ',g7.1E1,4x,
     *           '[1e-2], target value, see actual values below',/,
     *        4x,'  auto-refine stage    ',g7.1E1,4x,
     *           '[1e-3], target value, see actual values below',/,
     *        4x,'resolution_factor      ',i2,9x,'>= 2 [2]',/,
     *        4x,'refinement_points       ',i2,8x,'[aut] or 1->',i2,
     *           '; aut = automatic',/,
     *        4x,'refinement_switch       ',l1,9x,'[T] F',/,
     *        4x,'keep_auto               ',l1,9x,'[T] F',/,
     *        4x,'keep_max              ',i7,5x,
     *           '[20000], ~100 < keep_max < ~k21/10 =',i7,/,
     *        4x,'solvus_tolerance_II     ',a7,3x,'0->1 [0.2]',/,
     *        4x,'global_reach_increment ',i2,9x,'>= 0 [0]',/,
     *        4x,'reach_increment_switch  ',a3,7x,'[on] off all',/,
     *        4x,'zero_mode               ',e7.1E2,3x,
     *           '0->1 [1e-6]; < 0 => off')
1190  format (/,2x,'1D grid options:',//,
     *        4x,'y_nodes               ',i3,' /',i3,4x,'[40/40], >0, '
     *          ,'<',i4,'; effective y-resolution ',i4,' /',i4,
     *           ' nodes',/
     *        4x,'grid_levels             ',i1,' /',i2,5x,'[1/4], >0, '
     *          ,'<',i2,/)
1200  format (/,2x,'2D grid options:',//,
     *        4x,'x_nodes               ',i3,' /',i3,4x,'[40/40], >0, '
     *          ,'<',i4,'; effective x-resolution ',i4,' /',i4
     *          ,' nodes',/
     *        4x,'y_nodes               ',i3,' /',i3,4x,'[40/40], >0, '
     *          ,'<',i4,'; effective y-resolution ',i4,' /',i4,
     *           ' nodes',/
     *        4x,'grid_levels             ',i1,' /',i2,5x,'[1/4], >0, '
     *          ,'<',i2,/,
     *        4x,'linear_model             ',a3,6x,'off [on]')
1210  format (/,2x,'Fractionation path options:',//,
     *        4x,'1d_path               ',i3,' /',i3,4x,
     *           '[20/150], >0, <',i4)
1220  format (/,2x,'Composition options:',//,
     *        4x,'closed_c_space         ',l1,10x,'F [T]')
1230  format (/,2x,'Input/Output options:',//,
     *        4x,'aqueous_output         ',l1,10x,'[F] T',/
     *        4x,'aqeuous_species        ',i3,8x,'[20] 0-',i3,/,
     *        4x,'aq_solvent_composition ',a3,8x,
     *        '[y] m: y => mol fraction, m => molality',/,
     *        4x,'aq_solute_composition  ',a3,8x,
     *        'y [m]: y => mol fraction, m => molality',/,
     *        4x,'spreadsheet            ',l1,10x,'[F] T',/,
     *        4x,'logarithmic_p          ',l1,10x,'[F] T',/,
     *        4x,'bad_number         ',f7.1,8x,'[NaN]',/,
     *        4x,'composition_constant   ',l1,10x,'[F] T',/,
     *        4x,'composition_phase      ',a3,8x,'[mol] wt',/,
     *        4x,'composition_system     ',a3,8x,'[wt] mol',/,
     *        4x,'proportions            ',a3,8x,'[vol] wt mol',/,
     *        4x,'absolute               ',l1,10x,'[F] T',/,
     *        4x,'cumulative             ',l1,10x,'[F] T',/,
     *        4x,'fancy_cumulative_modes ',l1,10x,'[F] T',/,
     *        4x,'interpolation          ',a3,8x,'[on] off ',/,
     *        4x,'melt_is_fluid          ',l1,10x,'[F] T',/,
     *        4x,'solution_names         ',a3,8x,'[model] abbreviation',
     *                                           ' full',/,
     *        4x,'structural_formulae    ',l1,10x,'[T] F',/,
     *        4x,'species_output         ',l1,10x,'[T] F',/,
     *        4x,'species_Gibbs_energies ',l1,10x,'[F] T',/,
     *        4x,'seismic_output         ',a3,8x,'[some] none all',/,
     *        4x,'pause_on_error         ',l1,10x,'[T] F',/,
     *        4x,'poisson_test           ',l1,10x,'[F] T',/,
     *        4x,'interim_results        ',a3,8x,'[auto] off manual',/,
     *        4x,'sample_on_grid         ',l1,10x,'[T] F')
1231  format (/,2x,'Input/Output options:',//,
     *        4x,'aq_output              ',l1,10x,'[T] F',/
     *        4x,'aq_species             ',i3,8x,'[20] 0-',i3,/,
     *        4x,'aq_solvent_composition ',a3,8x,
     *        '[y] m: y => mol fraction, m => molality',/,
     *        4x,'aq_solute_composition  ',a3,8x,
     *        'y [m]: y => mol fraction, m => molality',/,
     *        4x,'logarithmic_p          ',l1,10x,'[F] T',/,
     *        4x,'bad_number         ',f7.1,8x,'[NaN]',/,
     *        4x,'composition_constant   ',l1,10x,'[F] T',/,
     *        4x,'composition_phase      ',a3,8x,'[mol] wt',/,
     *        4x,'composition_system     ',a3,8x,'[wt] mol',/,
     *        4x,'proportions            ',a3,8x,'[vol] wt mol',/,
     *        4x,'melt_is_fluid          ',l1,10x,'[F] T',/,
     *        4x,'solution_names         ',a3,8x,'[mod] abb ful',/,
     *        4x,'structural_formulae    ',l1,10x,'[T] F',/,
     *        4x,'species_output         ',l1,10x,'[T] F',/,
     *        4x,'endmember_Gs           ',l1,10x,'[F] T',/,
     *        4x,'seismic_output         ',a3,8x,'[some] none all',/,
     *        4x,'pause_on_error         ',l1,10x,'[T] F',/,
     *        4x,'poisson_test           ',l1,10x,'[F] T')
1232  format (/,2x,'Input/Output options:',//,
     *        4x,'spreadsheet            ',l1,10x,'[F] T',/,
     *        4x,'logarithmic_p          ',l1,10x,'[F] T',/,
     *        4x,'bad_number          ',f7.1,7x,'[NaN]',/,
     *        4x,'melt_is_fluid          ',l1,10x,'[F] T',/,
     *        4x,'pause_on_error         ',l1,10x,'[T] F',/,
     *        4x,'Tisza_test             ',l1,10x,'[F] T')
1233  format (/,2x,'Seismic wavespeed computational options:',//,
     *        4x,'seismic_data_file      ',l1,10x,'[F] T',/,
     *        4x,'bounds                 ',a3,8x,'[VRH] HS',/,
     *        4x,'vrh/hs_weighting       ',f3.1,8x,'[0.5] 0->1',/,
     *        4x,'explicit_bulk_modulus  ',l1,10x,'[T] F',/,
     *        4x,'poisson_ratio          ',a3,8x,'[on] all off; ',
     *        'Poisson ratio = ',f4.2,/,
     *        4x,'seismic_output         ',a3,8x,'[some] none all',/,
     *        4x,'poisson_test           ',l1,10x,'[F] T',/,
     *        4x,'Tisza_test             ',l1,10x,'[F] T')
1234  format (4x,'auto_exclude           ',l1,10x,'[T] F')
1240  format (/,2x,'Information file output options:',//,
     *        4x,'option_list_files      ',l1,10x,'[F] T; ',
     *           'echo computational options',/,
     *        4x,'pseudocompound_file    ',l1,10x,'[F] T; ',
     *           'echo static pseudocompound compositions')
1241  format (/,2x,'Information file output options:',//,
     *        4x,'option_list_files      ',l1,10x,'[F] T; ',
     *           'echo computational options')
1250  format (4x,'auto_refine_file       ',l1,10x,'[T] F; ',
     *           'echo auto-refine compositions')

      end

      subroutine redcd1 (lun,ier,key,val,nval1,nval2,nval3,strg,strg1)
c----------------------------------------------------------------------
c this routine seeks a card containing a keyword and as many as 
c six values (char variables nval, nval1...), the first 3 letters of nval
c is also returned in val, if the second value is longer than 12 characters
c it is also saved in the character variable strg*80
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer lun, ier, iscan, iscnlt, ibeg, iend, ist, lend

      character card*(lchar), key*22, val*3,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
 
      ier = 0 
      key = ' '

      do 

         read (lun,'(a)',iostat=ier) card
         if (ier.ne.0) return

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

            if (ibeg.ge.com) cycle
c                                 for programs (actcor,ctransf) 
c                                 that echo data read the full card
            length = iscnlt (lchar,1,' ')

            exit 

         end if 

      end do 
c                                 find end of keyword 
      iend = ibeg + 1
      iend = iscan (iend,lchar,' ') - 1

c      if (iend-ibeg.gt.21) then
c         call warn (99,0d0,ier,'invalid keyword in '
c     *   //'REDCD1, keywords must be < 23 characters.')
c         ier = 1
c         return 
c      end if 

      if (iend-ibeg.gt.21) then
         lend = ibeg + 21
      else 
         lend = iend
      end if 
c                                 load chars into key
      write (key,'(22a)') chars(ibeg:lend)

      iend = iend + 1
c                                 now locate the value:
      ibeg = iscnlt (iend,com,' ')
c                                 now find trailing blank
      iend = iscan (ibeg,lchar,' ') 
c                                 return if just a keyword
      if (iend.gt.lchar) return
c                                 look if it contains a comment character
      ist = iscan (ibeg,iend,'|') 
      if (ist.lt.iend) iend = ist - 1
c                                 save longer versions (only on first value)
c                                 this is done in case it's long text or 
c                                 several numbers on certain options. 
      strg = ' '
      strg1 = ' '
      nval1 = '0'
      nval2 = '0'
      nval3 = '0'
      
      if (iend-ibeg.gt.39) iend = ibeg+39
      write (strg,'(40a)') chars(ibeg:iend)
      write (strg1,'(40a)') chars(ibeg:ibeg+39)
c                                 read value:
      if (ibeg+2.gt.iend) then 
         ist = iend
      else 
         ist = ibeg + 2
      end if 

      write (val,'(3a)') chars(ibeg:ist)
c                                 look for a second value
      ist = iscan (ibeg,lchar,' ')
      if (ist.gt.com) return

      ibeg = iscnlt (ist,com,' ')
      if (ibeg.gt.com) return 

      iend = iscan (ibeg,com,' ')
      if (iend-ibeg.gt.11) iend = ibeg + 11 
      write (nval1,'(12a)') chars(ibeg:iend)
c                                 look for a third value
      ist = iscan (ibeg,lchar,' ')
      if (ist.gt.com) return 

      ibeg = iscnlt (ist,com,' ')
      if (ibeg.gt.com) return 

      iend = iscan (ibeg,com,' ')
      if (iend-ibeg.gt.11) iend = ibeg + 11 
      write (nval2,'(12a)') chars(ibeg:iend) 
c                                 look for a fourth value
      ist = iscan (ibeg,lchar,' ')
      if (ist.gt.com) return

      ibeg = iscnlt (ist,com,' ')
      if (ibeg.gt.com) return

      iend = iscan (ibeg,com,' ')
      if (iend-ibeg.gt.11) iend = ibeg + 11 
      write (nval3,'(12a)') chars(ibeg:iend)

      end

      subroutine rdstrg (lun,nstrg,string,eof)
c----------------------------------------------------------------------
c rdstrg - read 240 column card images from unit lun until a non-blank
c (i.e., with data other than comments) record, then read up to three
c strings from the record. on output nstrg is the number of strings read
c from the record. 
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer lun, iscan, iscnlt, ibeg, iend, ier, nstrg, imax

      logical eof

      character card*(lchar), string(3)*8

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      eof = .false.

      do 
c                                 read cards till a non-blank
c                                 card.
         read (lun,'(a)',iostat=ier) card

         if (ier.ne.0) then
c                                 error on read = eof
            eof = .true.

            return

         else if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1

            if (com.eq.0) cycle 
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

            exit 

         end if 

      end do 
c 
c                                 we have a non-blank card
      nstrg = 1

      do 
c                                 find the end of the string
         iend = iscan (ibeg,lchar,' ') - 1

         if (iend-ibeg.gt.7) then

            imax = ibeg + 7 

         else

            imax = iend
 
         end if 
c                                 load chars into string
         write (string(nstrg),'(8a)') chars(ibeg:imax)
c                                 find the next string
         ibeg = iscnlt (iend+1,com,' ')

         if (ibeg.gt.com.or.nstrg.eq.3) return
 
         nstrg = nstrg + 1

      end do 

      end

      subroutine rdnumb (numb,def,inumb,idef,reel)
c----------------------------------------------------------------------
c rdnumb - reads a line from terminal input for numeric input, if blank
c assigns default (def, idef); if non numeric prompts for new value.
c----------------------------------------------------------------------    
      implicit none

      integer inumb, idef, ier

      double precision numb, def

      logical defalt, reel

      character card*80
c----------------------------------------------------------------------
      defalt = .true.

      do 
c                                 read input
         read (*,'(a)',iostat=ier) card
c                                 user enters a blank
         if (ier.ne.0.or.card.eq.' ') exit
c                                 read data from card
         if (reel) then 
            read (card,*,iostat=ier) numb
         else 
            read (card,*,iostat=ier) inumb
         end if 

         if (ier.ne.0) then 
            call rerr 
         else
            defalt = .false.
            exit 
         end if 

      end do 

      if (defalt) then 
         if (reel) then 
            numb = def
         else 
            inumb = idef
         end if 
      end if 

      end

      subroutine eohead (n)
c----------------------------------------------------------------------
c eohead reads cards from n until an 'END ' or 'end ' is found in
c the first 4 columns
c----------------------------------------------------------------------
      implicit none

      integer n, ier

      character tag*4

      rewind n

      do 
         read (n,'(a)',iostat=ier) tag
         if (ier.ne.0) call error (37,1d0,n,'EOHEAD')
         if (tag.eq.'end'.or.tag.eq.'END') exit
      end do 

      end

      subroutine rerror (ier,*)
c---------------------------------------------------------------------
c rerror - routine to check for errors during list directed i/o
 
      implicit none

      integer ier
c---------------------------------------------------------------------
 
      if (ier.eq.0) then
         return
      else
         write (*,1000)
         ier = 0
         return 1
      end if
 
1000  format (/,'Your input is incorrect, probably you have specified ',
     *        'an invalid numerical value',/,'or you are using ',
     *        'a character where you should be using a number ',
     *        'or vice versa.',/,'try again...',/)
 
      end

      subroutine rerr 
c---------------------------------------------------------------------
c rerror - routine to write bad input message for interactive i/o
 
      implicit none
c---------------------------------------------------------------------

      write (*,1000)
 
1000  format (/,'Your input is incorrect, probably you are using ',
     *        'a character where',/,'you should be using a number ',
     *        'or vice versa, try again...',/)
 
      end

      subroutine errpau 
c---------------------------------------------------------------------
c pause on error option
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character a*1
      
      if (lopt(19)) then 
         write (*,'(/,a,/)') 'Press Enter to quit...' 
         read (*,'(a)') a
      end if 
      
      stop
      
      end 

      subroutine error (ier,realv,int,char)
c---------------------------------------------------------------------
c write error messages and terminate execution
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, int
 
      character char*(*), tag*3

      double precision realv

      integer iam
      common/ cst4 /iam

      character tname*10
      logical refine, resub
      common/ cxt26 /refine,resub,tname

      if (refine) then
         tag = '2nd'
      else 
         tag = '1st'
      end if

      if (ier.eq.1.or.ier.eq.2) then 
         write (*,1) char,int
      else if (ier.eq.3) then 
         write (*,3) char
      else if (ier.eq.4) then
         write (*,4) char 
      else if (ier.eq.5) then 
         write (*,5) int,char,j3
      else if (ier.eq.6) then 
         write (*,6) char
      else if (ier.eq.7) then 
         write (*,7) char
      else if (ier.eq.8) then 
         write (*,8) 
      else if (ier.eq.9) then 
         write (*,9) char
      else if (ier.eq.10) then 
         write (*,10) char, int
      else if (ier.eq.11) then 
         write (*,11) char,int
      else if (ier.eq.12) then 
         write (*,12) char
      else if (ier.eq.13) then
         write (*,13) h8
      else if (ier.eq.14) then
         write (*,14) char
      else if (ier.eq.15) then
         write (*,15) char
      else if (ier.eq.16) then
         write (*,16) h5
      else if (ier.eq.17) then
         write (*,17) int
      else if (ier.eq.18) then
         write (*,18) char
      else if (ier.eq.19) then
         write (*,19) char
      else if (ier.eq.20) then
         write (*,20) int, char
      else if (ier.eq.21) then
         write (*,21) char 
      else if (ier.eq.22) then
         write (*,22) int, char
      else if (ier.eq.23) then
         write (*,23) char
      else if (ier.eq.24) then
         write (*,24) int
      else if (ier.eq.25) then
         write (*,25) h9
      else if (ier.eq.26) then
         write (*,26) int, char
      else if (ier.eq.27) then
         write (*,27) char
      else if (ier.eq.28) then 
         write (*,28) int, char
      else if (ier.eq.29) then 
         write (*,29) int, char
      else if (ier.eq.32) then 
         write (*,32)
      else if (ier.eq.33) then 
         write (*,33) char, int
      else if (ier.eq.34) then
         write (*,34)
      else if (ier.eq.35) then
         write (*,35)
      else if (ier.eq.36) then
         write (*,36)
      else if (ier.eq.37) then
         write (*,37) int
      else if (ier.eq.38) then 
         write (*,38) 
      else if (ier.eq.39) then
         write (*,39) int
      else if (ier.eq.40) then
         write (*,40) int, char
      else if (ier.eq.41) then

         write (*,41)

         if (int.eq.0) then 
            write (*,410)
         else if (int.eq.1) then 
            write (*,411)
         else if (int.eq.2) then 
            write (*,*) ' are you kidding me? please report'
         end if

         write (*,414) char
         write (*,413)

      else if (ier.eq.42) then 
         write (*,42) char
      else if (ier.eq.43) then
         write (*,43) char
      else if (ier.eq.44) then 
         write (*,44) 
      else if (ier.eq.45) then 
         write (*,45) 
      else if (ier.eq.46) then 
         write (*,46) iopt(16), int, char
      else if (ier.eq.47) then 
         write (*,47) char
      else if (ier.eq.48) then 
         write (*,48) char,int
      else if (ier.eq.49) then 
         write (*,49) char,int
      else if (ier.eq.50) then 
         write (*,50) realv, char, int
      else if (ier.eq.51) then
         write (*,51) char
      else if (ier.eq.52) then 
         write (*,52) h9
      else if (ier.eq.53) then 
         write (*,53) 
      else if (ier.eq.54) then 
         write (*,54)
      else if (ier.eq.55) then 
         write (*,55) k16
      else if (ier.eq.56) then 
         write (*,56) k17
      else if (ier.eq.57) then 
         write (*,57) char, char, char
      else if (ier.eq.60) then 
         write (*,60) k22, char
      else if (ier.eq.61) then 
         write (*,61) k18, char
      else if (ier.eq.62) then
         write (*,62) char, int, realv
      else if (ier.eq.63) then 
         write (*,63) 
      else if (ier.eq.64) then
         write (*,64) char
      else if (ier.eq.65) then
         write (*,65) char
      else if (ier.eq.66) then
         write (*,66) 
      else if (ier.eq.67) then
         write (*,67) char
      else if (ier.eq.68) then
         write (*,68) char
      else if (ier.eq.69) then 
         write (*,69) char
      else if (ier.eq.72) then 
         write (*,72) char
      else if (ier.eq.73) then 
         write (*,73) char
      else if (ier.eq.74) then 
         write (*,74) int
      else if (ier.eq.75) then 
         write (*,75) char
      else if (ier.eq.76) then 
         write (*,76) char, char
      else if (ier.eq.78) then 
         write (*,78) char,char
      else if (ier.eq.89) then
         write (*,89) 
      else if (ier.eq.90) then
         write (*,90) l6
      else if (ier.eq.106) then
         write (*,106) char
      else if (ier.eq.107) then
         write (*,107) int
      else if (ier.eq.108) then
         write (*,108) int
      else if (ier.eq.109) then
         write (*,109) int
      else if (ier.eq.110) then
         write (*,110)
      else if (ier.eq.111) then
         write (*,111)
      else if (ier.eq.112) then
         write (*,112) char
      else if (ier.eq.116) then
         write (*,116)
      else if (ier.eq.117) then
         write (*,117)
      else if (ier.eq.118) then
         write (*,118)
      else if (ier.eq.120) then
         write (*,120) char
      else if (ier.eq.125) then 
         write (*,125) char, int
      else if (ier.eq.169) then
         write (*,169) int
      else if (ier.eq.180) then

         if (resub) then
c                                 dynamic 
            write (*,180) k13,int,char,tag
            write (*,1804)

         else
c                                 static
            write (*,180) k13,int,tname,tag
            write (*,1803)
            write (*,1804)

         end if

         if (realv.gt.0d0) then 
c                                 non-linear schemes are in use, instruct
c                                 accordingly:
            if (.not.lopt(38).and.iam.ne.15) then
               write (*,1805) tname,'T (true).'
            else if (lopt(38).and.iam.eq.15) then
               write (*,1805) tname,'F (false).'
            end if

         end if

      else if (ier.eq.181) then
         write (*,181) int
      else if (ier.eq.182) then
         write (*,182) k2
      else if (ier.eq.183) then
         write (*,183) k2,char
      else if (ier.eq.197) then
         write (*,197) int, k5, char
      else if (ier.eq.200) then
         write (*,200)
      else if (ier.eq.204) then
         write (*,204) int
      else if (ier.eq.206) then
         write (*,206) int
      else if (ier.eq.207) then
         write (*,207) realv,char
      else if (ier.eq.208) then
         write (*,208) char
      else if (ier.eq.227) then
         write (*,227) char, int
      else
         write (*,999) ier, realv, int, char
      end if
 
      call errpau

1     format (/,'**error ver001** increase parameter ',a,' to ',i7,' in'
     *       ,' perplex_parameters.h and recompile Perple_X',/)
2     format (/,'**error ver002** solution model file versions ',
     *         '6.8.2-6.8.7 are inconsistent with',/,
     *         'this version of Perple_X. Update the solution ',
     *         'file.',/)
3     format (/,'**error ver003** the solution model file format (',a,
     *         ') is inconsistent with',/,
     *         'this version of Perple_X. Update the file and/or ',
     *         'Perple_X',/)
4     format (/,'**error ver004** you must use ',a,' to analyze this '
     *         ,'type of calculation.',/)
5     format (/,'**error ver005** too many ordered species (',i2,') in',
     *        ' solution model ',a,/,'increase dimension j3 (',i2,')',/)
6     format (/,'**error ver006** fractionation path coordinate file: '
     *          ,a,/,'does not exist.',/)
7     format (/,'**error ver007** reference phase ',a,' is not in the ',
     *          'thermodynamic data file.',/)
8     format (/,'**error ver008** the thermodynamic data file ',
     *          'is out of date, download',/,'the current ',
     *          'version from: www.perplex.ethz.ch',/)
9     format (/,'**error ver009** invalid tag (',a,') in the ',
     *          'thermodynamic data file.',/)
10    format (/,'**error ver010** the text string:',/,a,/,'is too long',
     *        ' to be processed.',/,'The maximum allowed length is ',i3,
     *        ' characters.',/)
11    format (/,'**error ver011** invalid ',a,' choice (',i3,')',a,/)
12    format (/,'**error ver012** file:',a,/,'is missing or formatted ',
     *        'incorrectly, run paralyzer or create/edit it manually',/)
13    format ('**error ver013** too many excluded phases, ',
     *        'increase dimension h8 (',i3,')')
14    format ('**error ver014** programming error, routine ',a)
15    format (/,'**error ver015** missing composant for: ',a,/)
16    format (/,'**error ver016** too many saturated components, ',
     *        'increase dimension h5 (',i2,')')
17    format (/,'**error ver017** too many composants for a saturation',
     *        ' constraint increase dimension h6 (',i3,')')
18    format (/,'**error ver018** ',a,' is defined as a saturated ',
     *        'phase component in the thermodynamic data file.')
19    format (/,'**error ver019** probable cause missing composant,',
     *        ' executing routine ',a)
20    format (/,'**error ver020** error ',i2,' reading solution model',
     *        ' file.',/,'   Reading model: ',a,' Check format.',/)
21    format (/,'**error ver021**error reading ',
     *        'header section of',/,'the thermodynamic data ',
     *        'file, last data read:',/,a,/,'Check formatting',/)
22    format (/,'**error ver022** too many divariant assemblages, ',
     *        'increase dimension j9 (',i8,') routine: ',a)
23    format (/,'**error ver023**error reading',
     *        ' thermodynamic data file.',/,'Last data read:',
     *        /,a,/,'Check formatting.',/)
24    format (/,'**error ver024** too many solution models in',
     *        ' solution model file',/,' increase parameter i9 (',
     *        i3,')')
25    format (/,'**error ver025** too many solution models ',
     *          'increase parameter h9 (',i3,')')
26    format (/,'**error ver026** the number of fixed components (',
     *        i2,') in ',a,/,' is >= the number of components ',/)
27    format (/,'**error ver027** Error reading the problem',
     *        ' definition file:',//,a,//,
     *        'Probable cause: You are using an ',
     *        'input file created by an out-of-date',/,
     *        '                version of BUILD, or you have',
     *        ' incorrectly edited the',/'                input file',
     *        ' created by BUILD',/)
28    format (/,'**error ver028** invalid buffer choice (',i3,') in',
     *          ' routine: ',a,/)
29    format (/,'**error ver029** unknown term type ',i6,' for',
     *          ' solution model: ',a,/)
32    format (/,'**error ver032** stability field calculations (',
     *          'option 2) are disabled in this version of PERPLEX',/)
33    format (/,'**error ver033** expression with too many terms in ',a
     *       ,/,'increase m0 or j6 to',i2,'.',/)
34    format (/,'**error ver034** vmax is lt vmin, check input.')
35    format (/,'**error ver035** dv is lt 0, check input.')
36    format (/,'**error ver036** missing composant for the saturated',
     *        ' phase,',/,'you have probably excluded either H2O or',
     *        ' CO2,',/,'or a composant is duplicated in the',
     *        ' thermodynamic data file',/)
37    format (/,'**error ver037** no end marker in header',/,
     *        'section of thermodynamic data file unit ',i2,/)
38    format (/,'**error ver038** you have configured a ',
     *       'problem with only one independent variable.',/,
     *       'This case cannot be handled by constrained minimization',
     *       ' use the unconstrained computational mode.'/)
39    format (/,'**error ver039** too many end-members, ',
     *        'increase dimension k12 (',i2,') Routine: ',a)
40    format (/,'**error ver040** too many compositional coordinates, ',
     *        'increase dimension k13 (',i7,')  Routine: ',a)
41    format (/,'**error ver041** too many static compositions this ',
     *          'error can be eliminated',/,'by one of the ',
     *          'following actions (best listed first):',/)
410   format (2x,'- increase exploratory stage initial_resolution',
     *           ' in perplex_option.dat',/,
     *        2x,'- restrict subdivision ranges of complex solutio',
     *           'ns in solution_model.dat')
411   format (2x,'- increase auto-refine stage initial_resolution ',
     *           'in perplex_option.dat')
413   format (2x,'- simplify the calculation, e.g., eliminate ',
     *           'components and/or simplify solution models')
414   format (2x,'- increase parameter ',a,' and recompile ',
     *           'Perple_X')
42    format (/,'**error ver042** cannot open file:',a,/,'check that it'
     *       ,' is not being used by another program',/)
43    format (/,'**error ver043** you cannot simultaneously treat: ',
     *          a,/,'as a thermodynamic solution and as a saturated',
     *          ' phase.',/)
44    format (/,'**error ver044** too many saturated phase components.'
     *        /)
45    format (/,'**error ver045** too many mobile components.'/)
46    format (/,'**error ver046** the first value of the iteration ',
     *          'keyword exceeds (',i2,') the value',/,'of MRES (',i3,
     *          ') specified in routine ',a,'. Either reduce the',/,
     *          'iteration keyword value or increase MRES.',/) 
47    format (/,'**error ver047** solution model ',a,' is incorrectly ',
     *        'formatted (van Laar).',/)
48    format (/,'**error ver048** too many terms in solution model ',a,
     *        ' increase parameter m1 (',i2,').',/)
49    format (/,'**error ver049** the order of solution model ',a,
     *        ' is too high, increase parameter m2 (',i2,').',/)
50    format (/,'**error ver050** requested resolution ',
     *          '(',f6.0,') for a component in solution:',a,/,
     *          'exceeds 1/MRES (MRES=',i5,') ',
     *          'reduce requested resolution or inrease',/,
     *          'MRES in routine CARTES',/)
51    format (/,'**error ver051** DUMMY1 could not find the auxilliary'
     *         ,' input file:',/,a,/,'required for open system model ',
     *          'computations (ICOPT=9).',/)
52    format (/,'**error ver052** too many solution models in your'
     *         ,' calculation',/,'reduce the number of models or ',
     *          'increase parameter H9 (',i2,').',/)
53    format (/,'**error ver053** phase fractionation calculations '
     *         ,'require >1 thermodynamic component.',/)
54    format (/,'**error ver054** unanticipated condition, probable ',
     *          'cause is incorrect ordering of',/,'endmembers in the',
     *          ' solution model, which leads to inconsistent site ',
     *          'occupancies',/)
55    format (/,'**error ver055** too many make definitions, delete '
     *         ,'unused definitions from the',/
     *         ,'thermodynamic data file or '
     *         ,'increase parameter K16 (',i3,') and recompile.',/)
56    format (/,'**error ver056** too many phases in a make definition'
     *         ,', increase parameter K17 (',i2,') and recompile.',/)
57    format (/,'**error ver057** failed on an accepted make definition'
     *         ,' for ',a,/,'routine INPUT2. Exclude ',a,' and restart',
     *          ' the calculation.',/,'If ',a,/,' is legitimate. please'
     *         ,' report this error.',/)
60    format (/,'**error ver060** too many coordinates generated by ',
     *        'refinement, increase dimension k22 (',i8,') routine: ',a)
61    format (/,'**error ver061** too many solution coordinates, ',
     *        'increase dimension k18 (',i8,') routine: ',a)
62    format (/,'**error ver062** solution model ',a,' specifies non-',
     *          'Cartesian subdivision (',i1,')',/,' and must be refor',
     *          'mulated for adapative minimization, but VERTEX cannot',
     *        /,' do the reformulation because the initial_reolution ',
     *          'keyword specified in',/,' perplex_option.dat (',f5.2,
     *          ') is invalid',/)
63    format (/,'**error ver063** inconsistent auto-refine data.',
     *        ' Suppress or reinitialize auto-refinement.',/) 
64    format (/,'**error ver064** PSVDRAW plots only ',
     *          'binary mixed-variable and ',/,
     *          'ternary composition diagrams (',a,').',/)
65    format (/,'**error ver065** dimensioning error (',a,').',/)
66    format (/,'**error ver066** invalid format, most probably this',
     *          ' result should be plotted with PSSECT.',/)
67    format (/,'**error ver067** file ',a,/,
     *        'is not formatted correctly for PSVDRAW.',/)
68    format (/,'**error ver068** solution model: ',a,
     *          ' is in a format that is no longer supported',/,
     *          'Use a more recent solution model file, e.g., copy ',
     *          'the current version from: ',//,
     *          'www.perplex.ethz.ch/datafiles/solution_model.dat',/)
69    format (/,'**error ver069** too many points (',a,'), increase ',
     *          'parameter L5',/)
72    format (/,'**error ver072** ',a,/)
73    format (/,'**error ver073** the thermodynamic data file has ', 
     *          'more than one entity named: ',a,/,'delete or rename ',
     *          'the entities, this error is often caused by make ',
     *          'definitions.')
74    format (/,'**error ver074** unrecognized EoS pointer (',i3, 
     *          ') in routine SETINS',/)
75    format (/,'**error ver075** more than one solution model is ',
     *          'named ',a,/,'delete or rename the replicate models in',
     *          ' the solution model file.',/)
76    format (/,'**error ver076** the ',a,' solution model was not ',
     *          'reformulated correctly.',//,
     *          'If this error was preceeded by warning ver114 for a ',
     *          'dependent endmember, then',/,'it is possible that ',
     *          'the endmember definition is missing from the model.',
     *           //,
     *          'Otherwise excluding either more or fewer ',a,
     *          ' endmembers may bypass/rectify this error.',/)
78    format (/,'**error ver078** ',a,' has dependent endmembers with ',
     *        'invalid site populations',/,'it cannot be used unless ',
     *        'it is corrected or the site_check_override keyword is',/,
     *        'specified at the end of the ',a,' model.',/)
89    format (/,'**error ver089** SMPLX programming error. Change ',
     *        'minimnization method.',/)
90    format (/,'**error ver090** SMPLX failed to converge within ', 
     *        i6,' iterations.',/,'Probable cause: the possible ',
     *        'phases do not span the systems composition',/,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/,'Alternatively, although less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h',/,
     *        'and recompiling VERTEX permit SMPLEX to converge.',/)
106   format (/,'**error ver106** programming error in ',a)
107   format (/,'**error ver107** the assemblage you input is ',
     *        'metastable (ier=',i3,').')
108   format (/,'**error ver108** the assemblage you input is ',
     *        'degenerate (ier=',i3,').')
109   format (/,'**error ver109** the assemblage you input does not '
     *       ,'span the specified bulk composition (ier=',i3,').')
110   format (/,'**error ver110** you have requested a calculation ',
     *        'with the composition ',/,'of a saturated phase as a ',
     *        'variable, but you have not defined its composition')
111   format (/,'**error ver111** you have requested a calculation ',
     *        'with the composition ',/,'of a saturated phase as a ',
     *        'variable, but the phase has only one component')
112   format (/,'**error ver112** the maximum value of an independent '
     *       ,'variable',/,'is ',a,' to the minimum value')
116   format (/,'**error ver116** an independent variable, or at least'
     *       ,' its name, is undefined')
117   format (/,'**error ver117** vmax(iv(3)) ne vmin(iv(3) but no ',
     *        'sectioning variable v(iv(3)) is defined')
118   format (/,'**error ver118** the default increment of the ',
     *        'sectioning variable will result ',/,
     *        'in the generation of more ',
     *        'than 10 sections, to avoid this',/,' error increase ',
     *        'the increment or modify this test')
120   format (/,'**error ver120** file:',/,a,/,
     *        'could not be opened, check that it exists or that it is',
     *        ' not in use by another program.',/) 
125   format (/,'**error ver125** a site fraction is out of range for ',
     *          'solution: ',a,' endmember ',i2,/,'The configurational',
     *          ' entropy model is probably incorrect.',/)
169   format (/,'**error ver169** cart, imod=',i2,' is an invalid ',
     *          'request')
180   format (/,'**error ver180** too many pseudocompounds ',i8,' gener'
     *       ,'ated subdividing site ',i1,' of solution: ',a,/,
     *        'this error can usually be eliminated by one of the ',
     *        'following actions (best listed first):',//,
     *        2x,'- increase the ',a,' value of the initial_resolution '
     *          ,'keyword in perplex_option.dat')
1803  format (2x,'- restrict the subdivision range for the solution')
1804  format (2x,'- if nonlinear subdivision is specified for the solu'
     *          ,'tion then increase the',/,4x,'resolution parameters '
     *          ,'or change to linear subdivision',/,
     *        2x,'- increase parameter k13 and recompile',/)
1805  format ('NOTE!! ',a,' specifies nonlinear subdivision sch',
     *        'emes, the number of pseudo-',/,'compounds generated by ',
     *        'nonlinear subdivision is independent of initial_resolut',
     *        'ion',/,'unless the non_linear_switch option is ',a,
     *        ' For additional information on',/,'nonlinear subdivi',
     *        'sion refer to the commentary in the header of the ',
     *        'solution model',/,'file. If you do not understand nonli',
     *        'near subdivision use linear subdivision.',/)
181   format (/,'**error ver181** too many reactions,',
     *          ' increase dimension k2 (',i6,')')
182   format (/,'**error ver182** too many invariant points,',
     *           ' increase parameter k2 (',i6,')')
183   format (/,'**error ver183** too many assemblages; increase ',
     *        ' parameter k2 (',i6,'), routine ',a)
197   format (/,'**error ver197** to many components (',i2,'), increase'
     *         ,' parameter k5 (',i2,'), routine:',a,/)
200   format (/,'**error ver200** you are trying to use a fluid ',
     *        'equation of state for an invalid component',/)
204   format ('**error ver204** too many stable assemblages, i',
     *        'ncrease dimension j9 (',i8,')',/)
206   format ('**error ver206** too many univariant assemblages ',
     *        'intersect the edges of the diagram, i',
     *        'ncrease dimension k2 (',i6,')',/)
207   format (/,'**error ver207** the value of the stretching ',
     *          ' parameter (',g13.6,')',/,'for solution ',a,
     *          ' is invalid (<1) for transform subdivision,',/,
     *          'check section 4 of PERPLEX documentation.',/)
208   format (/,'**error ver208** too many phases on one side of a',/
     *        ' reaction.',/,'Do not use the full reaction',
     *        ' equation option (',a,').')
227   format (/,'**error ver227** in solution model ',a,' a DQF ',
     *          'correction is specified for endmember: ',i2,/,
     *          'DQF corrections can only be made on the idependent ',
     *          'endmembers of a solution model',/)
999   format (/,'**error vertex** unspecified error ier=',i3,/,
     *        ' real=',g15.7,/,' i=',i12,/,' char=',a)
      end

      subroutine warn (ier,realv,int,char)
c---------------------------------------------------------------------
c write warning message and continue execution

c generic warnings:

c 49  - future instances of warning int will not be repeated.
c 99  - just dump char
c 100 - use int (i3) as error number and dump char.
c 999 - unspecified real, int, char dump.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      integer ier,int

      double precision realv
 
      character char*(*)

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)
c----------------------------------------------------------------------
      if (ier.eq.1) then 
         write (*,1) 
      else if (ier.eq.2) then 
         write (*,2) realv
      else if (ier.eq.3) then 
         write (*,3)
      else if (ier.eq.4) then 
         write (*,4) char
      else if (ier.eq.5) then
         write (*,5) 
      else if (ier.eq.6) then
         write (*,6) 
      else if (ier.eq.7) then
         write (*,7) 
      else if (ier.eq.8) then
         write (*,8) h8
      else if (ier.eq.9) then
         write (*,9) char
      else if (ier.eq.10) then
         write (*,10) int, realv, char
      else if (ier.eq.11) then
         write (*,11) char
      else if (ier.eq.12) then
         write (*,12) char
      else if (ier.eq.13) then
         write (*,13) char, char
      else if (ier.eq.14) then
         write (*,14) char
      else if (ier.eq.15) then
         write (*,15)
      else if (ier.eq.16) then
         write (*,16) char
      else if (ier.eq.18) then
         write (*,18) realv
      else if (ier.eq.19) then
         write (*,19) 
      else if (ier.eq.20) then
         write (*,20)
      else if (ier.eq.21) then
         write (*,21) realv, char
      else if (ier.eq.22) then
         write (*,22) realv, char
      else if (ier.eq.23) then
         write (*,23) char     
      else if (ier.eq.24) then
         write (*,24) realv
      else if (ier.eq.25) then 
         write (*,25) int, char
      else if (ier.eq.26) then 
         write (*,26) char
      else if (ier.eq.27) then 
         write (*,27) int
      else if (ier.eq.28) then
         write (*,28)
      else if (ier.eq.29) then
         write (*,29) char
      else if (ier.eq.30) then
         write (*,30) char
      else if (ier.eq.31) then 
         write (*,31)
      else if (ier.eq.32) then
         write (*,32) char
      else if (ier.eq.33) then
         write (*,33) char
      else if (ier.eq.34) then
         write (*,34) char
      else if (ier.eq.35) then
         write (*,35) char,realv
      else if (ier.eq.36) then 
         write (*,36) realv, char 
      else if (ier.eq.37) then
         write (*,37)  
      else if (ier.eq.38) then
         write (*,38) 
      else if (ier.eq.39) then
         write (*,39) 
      else if (ier.eq.40) then
         write (*,40) 
      else if (ier.eq.42) then
         write (*,42)     
      else if (ier.eq.45) then
         write (*,45) char
      else if (ier.eq.46) then 
         write (*,46) realv, char, char
      else if (ier.eq.47) then
         write (*,47) int, realv
      else if (ier.eq.48) then 
         write (*,48) 
      else if (ier.eq.49) then 
         write (*,49) int, char
      else if (ier.eq.50) then
         write (*,50) char
      else if (ier.eq.51) then 
         write (*,51) char
      else if (ier.eq.52) then 
         write (*,52) char
      else if (ier.eq.53) then 
         write (*,53) realv
      else if (ier.eq.54) then 
         write (*,54)
      else if (ier.eq.55) then 
         write (*,55) char
      else if (ier.eq.56) then 
         write (*,56) char
      else if (ier.eq.58) then

         write (*,58)
         write (*,582)
         if (nopt(21).gt.2d0) write (*,587)
         write (*,583)
         if (lopt(49)) write (*,584)
         write (*,585)
         if (lopt(32)) write (*,586)
         write (*,581)
         write (*,413)
         write (*,580) char

58    format (/,'**warning ver058** exhausted memory during adaptive o',
     *       'ptimization, the result will',/,'be reported as a failed',
     *       ' optimization, if such failures are excessive the ',/,
     *       'problem may be eliminated by the ',
     *       'following actions (best listed first):',/)
580   format (2x,'- increase parameter ',a,' and recompile ',
     *           'Perple_X')
581   format (2x,'- set the low_reach flag for high-dimension ',
     *           'solution models in solution_model.dat')
582   format (2x,'- set keep_auto to T or default ',
     *           'in perplex_option.dat')
583   format (2x,'- reduce any reach_increments specified ',
     *           'in solution_model.dat',/,
     *        2x,'  or set reach_increment_switch to off ',
     *           'in perplex_option.dat')
584   format (2x,'- set refinement_switch to F ',
     *           'in perplex_option.dat')
585   format (2x,'- reduce refinement_points (< c+2, > 0) ',
     *           'in perplex_option.dat')
586   format (2x,'- set aq_solvent_solvus to F ',
     *           'in perplex_option.dat')
587   format (2x,'- reduce resolution_factor to 2 (default) ',
     *           'in perplex_option.dat')
413   format (2x,'- simplify the calculation, e.g., eliminate ',
     *           'components and/or simplify solution models')

      else if (ier.eq.589) then
         write (*,589)
      else if (ier.eq.59) then
         write (*,59) char
      else if (ier.eq.60) then
         write (*,60) char
         if (int.eq.1) then 
            write (*,601) char
         else 
            write (*,602)
         end if
      else if (ier.eq.61) then
         write (*,61) char
      else if (ier.eq.63) then
         write (*,63)
      else if (ier.eq.68) then
         write (*,68)
      else if (ier.eq.73) then
         write (*,73) char, realv, int
      else if (ier.eq.74) then
         write (*,74)
      else if (ier.eq.79) then
         write (*,79) char
      else if (ier.eq.87) then
         write (*,87)
      else if (ier.eq.88) then
         write (*,88)
      else if (ier.eq.89) then
         write (*,89)
      else if (ier.eq.90) then
         write (*,90) 
      else if (ier.eq.91) then
         write (*,91)
      else if (ier.eq.92) then 
         write (*,92) int, l7, char, (l7 - 1)/2**(grid(3,2)-1) + 1
      else if (ier.eq.99) then
         write (*,99) char
      else if (ier.eq.100) then
         write (*,100) int, char
      else if (ier.eq.106) then
         write (*,106) char
      else if (ier.eq.108) then
         write (*,108)
      else if (ier.eq.109) then
         write (*,109)
      else if (ier.eq.113) then
         write (*,113) int
      else if (ier.eq.114) then
         write (*,114)
      else if (ier.eq.172) then
         write (*,172) 
      else if (ier.eq.173) then
         write (*,173) 
      else if (ier.eq.175) then
         write (*,175) char,ier,realv
      else if (ier.eq.176) then
         write (*,176) char, iopt(21)
      else if (ier.eq.177) then
         write (*,177) nopt(5)
      else if (ier.eq.190) then
         write (*,190) int
      else if (ier.eq.205) then
         write (*,205) int
         write (*,900)
      else if (ier.eq.228) then 
         write (*,228) char, realv, int, char
      else
         write (*,999) ier, char, realv, int
      end if

1     format (/,'**warning ver001** the amount of a saturated phase is'
     *       ,' < 0, this indicates that',/,'the specified amount of a '
     *       ,'saturated component is inadequate to saturate the system'
     *       ,/)
2     format (/,'**warning ver002** the amount of a phase is <',g12.3,
     *        ' (-zero_mode) this may be',/,'indicative of numeric ',
     *        'instability',/)
3     format (/,'**warning ver003** the solution model file is ',
     *         ' inconsistent with this',/,
     *         'this version of Perple_X. Update the file and/or ',
     *         'Perple_X',/)
4     format (/,'**warning ver004** the data includes ',a,' values, '
     *      ,'probably because bad_number',/,'in perplex_option.dat = '
     *      ,'NaN, these values will be replaced by zeros. To avoid ',/
     *      ,'this behavior set bad_number to a numeric value or use a'
     *      ,' plotting program capable',/,'of handling NaNs, e.g., ',
     *       'MatLab or PYWERAMI',/)
5     format (/,'**warning ver005** fluid components are specified',
     *        ' as thermodynamic AND as either',/,'saturated phase',   
     *      ' or saturated components; almost certainly a BAD idea.',/)
6     format (/,'**warning ver006** fluid components are specified',
     *        ' as both thermodynamic AND',/,'saturated ',    
     *        'components; almost certainly a BAD idea.',/)
7     format (/,'**warning ver007** fluid components are specified as'
     *       ,' a saturated phase component',/,'AND as a thermodynamic', 
     *        'or saturated component; almost certainly a BAD idea.',/)
8     format ('**warning ver08** too exclude more phases ',
     *        'increase paramter h8 (',i3,')')
9     format ('**warning ver009** unable to deconstruct transition,'
     *       ,/,'data for ',a,' will not be output.')
10    format (/,'**warning ver010** not able to traverse ',
     *          'the entire  extent of equilibrium',/,'(',i6,')',
     *          ' at v(3)=',g12.6,/,
     *          'this error can usually be avoided by increasing the ',
     *          'finite difference',/,'increment delt(iv(1)) or delv',
     *          '(iv(2)), as defined on card 6 of',/,'the file on n2.',
     *          ' In routine:',a,/)
11    format (/,'**warning ver011** ',a,' has > 1',
     *          ' transition with dp/dT ne 0 and may not be treated ',/,
     *          ' correctly')
12    format (/,'**warning ver012** ',a,' has a transition ',
     *          ' with dp/dT < 0 and may not be treated ',/,
     *          ' correctly')
13    format (/,'**warning ver013** the total amount of the thermodynam'
     *         ,'ic components in: ',a,'is < 0,',/,'it will be rejected'
     *         ,' from the thermodynamic composition space. If non-elem'
     *         ,'ental (e.g., ',/,'oxide) ',
     *          'components are in use this rejection criterion may be '
     *         ,'incorrect. When this is',/,'the case redefine the ',
     *          'data base components so that the total amount of the'/,
     *          'thermodynamic components in ',a,' is > 0.')
14    format (/,'**warning ver014** You can not redefine the ',
     *          'saturated phase component:',a,/,'To circumvent this ',
     *          'restriction use CTRANSF to make a data base with the',/
     *         ,'the desired component transformations',/)
15    format (/,'**warning ver015** if you select > 1 saturated ',
     *          'component, then the order you',/,'enter the ',
     *          'components determines the saturation heirarchy and may'
     *          ,' effect your',/,'results (see Connolly 1990).',/)
16    format (/,'**warning ver016** ',a,' has been rejected because it',
     *       ' has no associated volumetric EoS.',/,'To override this ',
     *        'behavior set auto_exclude to false or add an ',
     *        'association.',/)
18    format (/,'**warning ver018** the value of the default dependen',
     *         't variable (',g14.6,') for the following',/,
     *         'equilibrium was inconsistent with the an earlier ',
     *         'determination of the invariant condition',/,
     *         'and will be reset. This may cause the curve to ',
     *         'kink near the invariant point',/)
19    format ('**warning ver019** you must specify at least ',
     *        'one thermodynamic component, try again',/)
20    format ('**warning ver020** sfol2')
21    format ('**warning ver021** xmax (',g12.6,') > 1 for '
     *         ,' solution model ',a,/,' xmax will be reset to 1',
     *        /,' see documentation, section 4.0.')
22    format ('**warning ver022** xmin (',g12.6,') < 0 for '
     *         ,' solution model ',a,/,' xmin will be reset to 1',
     *        /,' see documentation, section 4.0.')
23    format ('**warning ver023** xmin > xmax for solution ',a,/,
     *        'xmin will be set to xmax NO PSEUDOCOMPOUNDS WILL BE',
     *        ' GENERATED.',/,'see documentation, section 4.0',/)
24    format (/,'**warning ver024** wway, increment refined out of',
     *          ' range (',g8.1,')',/,'before the stable',
     *          ' extension of the equilibria was located')
25    format ('**warning ver025** ',i1,' endmembers for ',a,
     *          ' The solution will not be considered.')
26    format ('**warning ver026** only one endmember for ',a,
     *          ' The solution will not be considered.')
27    format (/,'**warning ver027** only ',i2,' user defined ',
     *      'compositions permitted.',/,'do multiple runs with WERAMI',
     *      'or redimension common block comps.',/)
28    format (/,'**warning ver028** minimization failed, ill-',
     *        'conditioned?',/)
29    format ('**warning ver029** programming error, routine ',a,/)
30    format (/,'**warning ver030** Because of missing endmembers, ',
     *        'or that the',/,'subdivision',
     *        ' scheme specified for solution model ',a,/,'is too',
     *        ' restrictive, there are no valid compositions for', 
     *        ' this model.',/)
31    format (/,'**warning ver031** this choice is disabled because ',
     *        'the dependent_potentials',/,'keyword is missing or off',
     *        ' in perplex_option.dat, to use this choice set the',/,
     *        'keyword to on and re-run VERTEX.',/)
32    format ('**warning ver032** fixed activity option requested',
     *          ' for ',a,/,'This option is disabled, the',
     *          ' solution will not be considered.')
33    format ('**warning ver033** missing endmembers for ',a,/,
     *        'The model may be recast in > 1 way for',
     *        ' the endmember subset.',/,'To control this choice',
     *        ' eliminate undesired endmembers.')
34    format ('**warning ver034** ',a,' could not be recast as',
     *          ' a simpler model.',/,'The solution will not be',
     *          ' considered. Add the missing endmembers or eliminate'
     *          ,/,'additional endmembers to allow this model.',/)
35    format (/,'**warning ver035** ',a,' is only for pure fluids',
     *        /,' XCO2 will be reset to: ',f4.2,/)
36    format ('**warning ver021** xinc (',g12.6,') < 0 for'
     *         ,' solution model ',a,/,'xinc will be reset to 1.'
     *         ,' see documentation, section 4.0. ',/)
37    format (/,'**warning ver37** you will not be able to plot the ',
     *       'results of this',/,'calculation with PSVDRAW. PSVDRAW ',
     *       'only plots ternary composition diagrams.',/)
38    format (/,'**warning ver38** you will not be able to plot the ',
     *       'results of this',/,'calculation with PSVDRAW. PSVDRAW ',
     *       'only plots mixed-variable diagrams for',/,'a binary ',
     *       'system with one independent potential variable.',/)
39    format (/,'**warning ver39** PSVDRAW will plot the results of ',
     *       'this calculation as a',/,'projected section, such plots ',
     *       'may be difficult to interpret. To plot',/,
     *       'pseudosections as a an explicit function of a systems ',
     *       'composition use',/, 'gridded minimization.',/)
40    format (/,'**warning ver040** you have configured a ',
     *       'problem with only one independent variable.',/)
42    format (/,'**warning ver042** an optimization failed due ',
     *          'to numerical instability',/,
     *          'or because the phases of the system do not span ',
     *          'its bulk composition.',//,
     *          4x,'In the 1st case:',/,
     *          8x,'increase (sic) final_resolution and/or',/,
     *          8x,'increase resolution_factor and/or',/,
     *          8x,'increase reach_increment and/or',/,
     *          8x,'increase speciation_factor and/or',/,
     *          8x,'increase speciation_max_it.',/,
     *          4x,'see: www.perplex.ch/perplex_options.html for ',
     *          'explanation.',//,
     *          4x,'In the 2nd case: ',
     *          'change the bulk composition or add phases.',/)
45    format (/,'**warning ver045** ',a,' involves a nonlinear EoS.',/,
     *          ' Output properties are stoichiometric averages.',/)
46    format (/,'**warning ver046** temperature (',g12.6,' K) is out',
     *        ' of range for endmember/phase: ',a,/,a,
     *        ' will be destabilized at this condition. In some cases',
     *        ' this problem can be corrected by',/,'setting ',
     *        'Anderson_Gruneisen to TRUE in the Perple_X option file.')
47    format (/,'**warning ver047** univariant field ',i6,' terminates',
     *        ' at an invariant field',/,'that could not be located ',
     *         'within the tolerance specified in the thermodynamic',/,
     *         'data file (PTOL= ',g12.6,').',/)
48    format (/,'**warning ver048** fluid phase pseudocompound data ',
     *         'does not include',/,' volumetric properties (SWASH).',/)
49    format (/,'**warning ver049** warning ',i3,' will not be repeated'
     *         ,' for future instances of this problem.',/,
     *          'currently in routine: ',a,//)
50    format (/,'**warning ver050** reformulating polytopic ',
     *          'solution: ',a,' because of missing endmembers. ',
     *        /,'(reformulation can be controlled explicitly ',
     *          'by excluding additional endmembers).',/)
51    format (/,'**warning ver051** cannot make ',a,' because of ',
     *          'missing data or an'
     *       ,/,'invalid definition in the thermodynamic data file.',/)
52    format (/,'**warning ver052** rejecting ',a,'; excluded or '
     *       ,'invalid composition.',/)
53    format (/,'**warning ver053** the failure rate during speciation',
     *          ' (o/d) calculations is ',f5.1,'%.',/,
     *          'A high failure rate may cause failed optimizations. ',
     *          'Usually the failure rate can be',/,'reduced by ',
     *          'increasing speciation_max_it in ',
     *          'perplex_option.dat',/)
54    format (/,'**warning ver054** property choices 25, 36, and 38 are'
     *         ,' not allowed in combination',/,'with other property '
     *         ,'choices',/)
55    format (/,'**warning ver055** possible compositions of ',a,'lie ',
     *         'entirely within the saturated',/,'component composition'
     *         ,' space. the compositions will not be considered.',/,
     *         'If this is problematic, then eliminate the component ',
     *         'saturation constraints',/,'or use convex.',/)
56    format (/,'**warning ver056** the EoS specified for ',a,' by the',
     *        ' hybrid_EoS option will be',/,'overridden by the EoS sp',
     *        'ecified in the problem definition file. To prevent this',
     *      /,'behavior delete the special_component section from the ',
     *        'header of the',/,'thermodynamic data file.',/)
589   format (/,'**warning ver589** wway, the equilibrium of the '
     *         ,'following reaction',/,'is inconsistent with the ',
     *          'invariant equilibrium.',/)
59    format (/,'**warning ver059** endmember ',a,
     *        ' has invalid site populations.',/)
60    format (/,'**warning ver060** Tait EoS conditions out of range ',
     *          'for endmember/phase: ',a,'(at T=',g12.6,' K)')
601   format ('base v1 is negative, endmember/phase ',a,' will be ',
     *        'destabilized.',/)
602   format ('base v2 will be zeroed on the assumption that it is ',
     *        'insignificant.',/)
61    format (/,'**warning ver061** the data includes NaN values, '
     *      ,'probably because bad_number',/,'in perplex_option.dat = '
     *      ,'NaN, these values will be replaced by zeros. To avoid ',/,
     *       'this behavior set bad_number to a numeric value or use a',
     *       ' plotting program capable',/,'of handling NaNs, e.g., ',
     *       'MatLab or PYWERAMI.',//,'program/routine: ',a,/)
63    format (/,'**warning ver063** wway, invariant point on an edge?',
     *        /)
68    format (/,'**warning ver068** degenerate initial assemblage in ',
     *          'COFACE, this should never occur',/,'if you see this ',
     *          'message please report the problem',/)
73    format (/,'**warning ver073** an invariant point has been ',
     *          'skipped, routine: ',a,/,
     *          'decreasing DTOL (',g9.3,') in the thermodynamic ', 
     *          'data file for variable ',i1,/,
     *          'may eliminate this problem',/)
74    format (/,'**warning ver074** no new equilibria identified,',
     *          ' if degenerate segments have',/,' been skipped',
     *          ' increase the computational reliability level.',/)
79    format (/,'**warning ver079** univeq failed on an edge for ',
     *          'the following equilibrium.',/,' Probable cause is ',
     *          'extreme independent variable limits (e.g., xco2=0)',/
     *          ' or poor convergence criteria ',
     *          'in the thermodynamic data file. In routine:',a,/)
87    format (/,'**warning ver087** wway-univeq did not converge ',
     *          'when div was refined',/)
88    format (/,'**warning ver088** SMPLX converged to a non-unique ',
     *        'solution.',/,3x,'Probable cause: system composition ',
     *        'coincides with that of ',//,3x,'a compound or a ',
     *        'tieline between compounds.',//,3x,'This may lead to ',
     *        'erratic results.',//,3x,'To avoid this problem ',
     *        'perturb the systems composition.',/)
89    format (//,'**warning ver089** BUILD you did not request',
     *        'plot file output.',/,' You will not be able to process',
     *        ' the results of the requested calculation.',//)
90    format (/,'**warning ver090** optimization failed. ',
     *        'Most probably, the possible ',
     *        'phases do not span',/,'the systems composition.',
     *        'In this case, add phases or modify the bulk ',
     *        'composition.',/,'Less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h may permit convergence.',/)
91    format (/,'**warning ver091** optimization failed. Change ',
     *        'minimnization method',/)
92    format (/,'**warning ver092** you have requested ',i4,
     *        ' grid points. Current',/,'dimensioning is for ',
     *        i4,' points. To obtain the requested resolution',/,
     *        'increase parameter L7 and recompile; or reduce the ',
     *        'required resolution via',/,'the ',a,' keyword in ',
     *        'perplex_option.dat',/,'The program will continue ',
     *        'with an effective grid resolution of ',i4,
     *        ' points.')
99    format (/,'**warning ver099** ',a,/)
100   format (/,'**warning ver',i3,'** ',a,/)
106   format ('**warning ver106** programming error in ',a)
108   format (/,'**warning ver108** wway, a phase field with the '
     *         ,'following',/,' reaction is stable on both ',
     *          'sides of an invariant point',/,' this error can ',
     *          'usually be avoided by increasing the finite ',/,
     *          ' difference increment delt(iv(1)) or delv',
     *          '(iv(2)), defined',
     *           /,' on card 6 of the thermodynamic data file',/)
109   format (/,'**warning ver109** you may ',
     *        'have assigned a mobile component as an independent ',/,
     *        ' variable without defining the component',/)
113   format (/,'**warning ver113** maximum variance for equilibrium',
     *        ' tracing [the variance keyword in ',/,
     *        'perplex_option.dat] must be > 0, but is ',i2,
     *        '. Set to 1 for the current calculation',/)
114   format (/,'**warning ver114** the default increment of an ',
     *        'independent variable is <',/,'1 percent of ',
     *        'its range, this is usually inefficient',/)
172   format (/,'**warning ver172** you cannot use this equation of',
     *          ' state with Y(CO2)',/,' as an indepedent variable, ',
     *          ' pick another one:',/)
173   format (/,'**warning ver173** invalid buffer choice ',/)
175   format (/,'**warning ver175** speciation routine ',a,' did',
     *          ' not converge ',/,' possibly due to graphite super-',
     *          'saturation. ier = ',i1,' real = ',g16.5,/)
176   format (/,'**warning ver176** fluid equation of state routine ',
     *        a,' did not converge.',/,'If execution continues this ',
     *        'may lead to incorrect results. To avoid this ',/,
     *        'problem increase speciation_max_it (',i4,') in ',
     *        'perplex_option.dat or choose',/,
     *        'a different equation of state.',//,
     *        'NOTE: at compositional extremes fluid speciation ',
     *        'calculations may require ',/,
     *        'thousands of iterations',/)
177   format (/,'**warning ver177** Invalid fluid speciation. ',
     *          'Reducing speciation tolerance (',g14.6,') in ',
     *          'perplex_option.dat',/,'may resolve this problem',/)
190   format (/,'**warning ver190** SMPLX failed to converge within ',
     *        i6,' iterations.',/,3x,'Probable cause: the possible ',
     *        'phases do not span the systems composition',/,3x,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/,3x,'Alternatively, although less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h',/,3x,
     *        'and recompiling VERTEX permit SMPLEX to converge.',/)
205   format (/,'**error ver205** too many new phase assemblages, ',
     *        'found by routine newhld',/,'increase dimension j9 (',
     *        i8,')',/)
228   format (/,'**warning ver228** in solution model ',a,' negative ',
     *          'composition (',g12.6,') for component ',i2,/,
     *          'this may indicate an incorrect stoichiometric ',
     *          'dependent endmember definition.',//,
     *          'This warning is issued only for the first negative ',
     *          'composition of ',a,/)
900   format ('the calculation may be incomplete !!!!',/)
999   format (/,'**warning unspecified** ier =',i3,' routine ',a6
     *         ,' r = ',g12.6,' int = ',i9,/)
      end

      subroutine rmakes (jopt)
c----------------------------------------------------------------------
c rmakes is called by topn2 to read make definitions of thermodynamic
c entities, these entities are defined as a linear combination of 
c exisiting entities as defined in the thermodynamic file, with an 
c optional pressure/temperature dependent DQF correction. The format
c assumes data on one line of less than 240 characters, the expected format
c is

c name = num_1 * name_1 + num_2 * name_2 ....+ num_int * name_int
c dnum_1 dnum_2 dnum_3

c where i num_j is a number or fraction (i.e., two numbers separated by a 
c '/') and name_j is the name of the int existing entities. 
c and the dqf correction to the entity 'name' is
c Gdqf(J/mol) = dnum_1 + T[K]*dnum_2 + P[bar]*dnum_3

c end_of_data is either a "|" or the end of the record.

c make definitions are preceded by the keyword:

c begin_makes 

c and truncated by the keyword:

c end_makes

c if jopt > 3, data is echoed to LUN n8 (for ctransf/actcor).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, ier, iscan, i, nreact, jopt

      double precision rnum, nums(m3)

      character tname*8, name*8, rec*(lchar), tag*3

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      call readcd (n2,ier,.true.)
      if (ier.ne.0) goto 90 
c                                 echo data for ctransf/actcor
      if (jopt.gt.3) write (n8,'(400a)') chars(1:com)

      nmak = 0 

      write (rec,'(400a)') chars
      read (rec,'(a3)') tag

      do while (tag.ne.'end')   

         nmak = nmak + 1
         if (nmak.gt.k16) call error (55,mkcoef(1,1),nmak,'RMAKES')
c                                 get first name
         ibeg = 1
         call readnm (ibeg,iend,com,ier,tname)
         if (ier.ne.0) goto 90
c                                 find start of data marker '='
         ibeg = iscan (1,com,'=') + 1
c                                 the rest of the data should
c                                 consist of coefficients followed
c                                 by names
         nreact = 0 

         do while (ibeg.lt.com) 
c                                 find the number
            call readfr (rnum,ibeg,iend,com,ier)
            if (ier.eq.2) then 
c                                 ier = 2 = a read error
               goto 90
            else if (ier.eq.1) then 
c                                 ier = 1, end-of-definition
               exit 
            end if 
c                                 find the name
            call readnm (ibeg,iend,com,ier,name)
            if (ier.ne.0) goto 90

            nreact = nreact + 1
            if (nreact.gt.k17) call error (56,mkcoef(1,1),nmak,'RMAKES')

            mkcoef(nmak,nreact) = rnum 
            mknam(nmak,nreact) = name
           
         end do

         if (nreact+1.gt.k17) call error (56,mkcoef(1,1),nmak,'RMAKES')
         mknam(nmak,nreact+1) = tname
         mknum(nmak) = nreact
c                                 now the dqf
         call readcd (n2,ier,.true.)
         if (ier.ne.0) goto 90
c                                 echo data for ctransf/actcor 
         if (jopt.gt.3) write (n8,'(400a)') chars(1:com)
c                                 read the DQF coefficients
         ibeg = 1
         call redlpt (nums,ibeg,iend,ier) 
         if (ier.ne.0) goto 90

         do i = 1, m3 
            mdqf(nmak,i) = nums(i)
         end do 
c                                 start next make definition
         call readcd (n2,ier,.true.)
         write (rec,'(400a)') chars
         read (rec,'(a3)') tag
c                                 echo data for ctransf/actcor
         if (jopt.gt.3) write (n8,'(400a)') chars(1:com)

c                                 reject excluded makes
         do i = 1, ixct
            if (tname.eq.exname(i)) then 
               nreact = nreact - 1
               exit 
            end if
         end do 

      end do 

      goto 99

90    write (*,1000) chars(1:com)
      stop
      
1000  format (/,'**error ver200** READMK bad make definition in the',
     *        ' thermodynamic data file',/,'currently reading: ',/
     *        ,400a)

99    end 

      subroutine readnm (ibeg,iend,siz,ier,name)
c----------------------------------------------------------------------
c readnm looks for the first word in a record chars, ibeg is the index
c of the 1st letter, iend is the index of the last letter.

c input 
c         ibeg - starting index for search
c         siz  - end index for search
c output
c         ibeg - starting index of word
c         iend - end index of word
c         ier  - error code
c         name - word
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, iscan, iscnlt, ier, siz

      external iscan, iscnlt

      character name*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 
c                                 find start of name
      ibeg = iscnlt (ibeg,siz,' ') 
c                                 find next blank
      iend = iscan (ibeg,siz,' ') - 1
c                                 initialize to be safe:
      name = ' '

      if ( iend - ibeg.le.7) then

         write (name,'(20a)') chars(ibeg:iend)

      else 
c                                 can't be a valid name, save it
c                                 anyway in case it's a tag
         write (name,'(20a)') chars(ibeg:ibeg+7)
         ier = 4

      end if 

      ibeg = iend + 1

      end 

      subroutine readcd (nloc,ier,strip)
c----------------------------------------------------------------------
c readcd - read 240 column card image from unit 9, strip out unwanted
c characters if strip. ier = 1 no card found.
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      logical strip

      integer ier, iscan, ict, i, iscnlt, ibeg, nloc

      external iscan, iscnlt

      character card*(lchar)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 

      ibeg = 0
  
      com = 0 

      card = ' '

      do while (ibeg.ge.com) 

         read (nloc,'(a)',end=90) card

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1
c                                 '|' in first column
            if (com.eq.0) cycle
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

         end if 

      end do 
c                                 there is a non-blank data character
      if (strip) then 

         ict = 1

         do i = 2, com
c                                 strip out '+' and '*' chars
            if (chars(i).eq.'+'.or.chars(i).eq.'*') chars(i) = ' '
c                                 eliminate blanks after '/' and '-'
c                                 and double blanks
            if ((chars(ict).eq.'/'.and.chars(i  ).ne.' ') .or. 
     *          (chars(ict).eq.'-'.and.chars(i  ).ne.' ') .or.
     *          (chars(ict).eq.' '.and.chars(i  ).ne.' ') .or.
     *          (chars(ict).ne.'-'.and.chars(ict).ne.'/'.and.
     *           chars(ict).ne.' ') ) then
                ict = ict + 1
                chars(ict) = chars(i)
            end if

         end do 

         com = ict

      else
c                                 scan backwards to the last non-blank
         com = iscnlt (com,1,' ')

      end if

      goto 99

90    ier = 3

99    end

      subroutine readfr (rnum,ibeg,iend,siz,ier)
c----------------------------------------------------------------------
c readfr looks for a number or two numbers separated by a backslash / in
c array elements chars(iend:ibeg), the latter case is interpreted as a ratio. 
c the result is returned as num
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision rnum, rnum1 

      integer ibeg, iend, iback, ier, iscan, iscnlt, siz

      external iscan, iscnlt

      character num*30

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0 
c                                 now find start of a number
      ibeg = iscnlt (ibeg,siz,' ')  
c                                 find backslash
      iback = iscan (ibeg,siz,'/') - 1
c                                 find next blank
      iend = iscan (ibeg,siz,' ') - 1
c                                 three cases:
      if (iend.ge.com) then

         ier = 1
         goto 99 

      else if (iback.gt.iend) then
c                                 no fraction
         if (iend-ibeg+1.gt.30) goto 90
c                                 first constant
         write (num,'(30a)') chars(ibeg:iend)
         read (num,*,err=90) rnum

      else 
c                                 fraction write numerator
         if (iback+1-ibeg.gt.30) goto 90
c                                 first number
         write (num,'(30a)') chars(ibeg:iback)       
         read (num,*,err=90) rnum
c                                 second number 

         if (iend-iback-1.gt.30) goto 90
         write (num,'(30a)') chars(iback+2:iend)
         read (num,*,err=90) rnum1

         rnum = rnum/rnum1

      end if 

      ibeg = iend + 1

      goto 99

90    ier = 2

99    end 

      subroutine redfr0 (rnum,ibeg,iend,ier)
c----------------------------------------------------------------------
c redfr0 looks for a number or two numbers separated by a backslash / in
c that array chars(iend:ibeg), the latter case is interpreted as a ratio. 
c the result is returned as rnum. differs from readfr in that redfr0
c expects ibeg/iend are known on input.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision rnum, rnum1 

      integer ibeg, iend, iback, ier, iscan

      external iscan

      character num*30

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0 

c                                 find backslash
      iback = iscan (ibeg,iend,'/') - 1
c                                 two cases:
      if (iback.ge.iend) then
       
         iback = iscan(ibeg,iend,' ') - 1
c                                 no fraction
         if (iback-ibeg+1.gt.30) goto 90
c                                 simple number
         write (num,'(30a)') chars(ibeg:iback)
         read (num,*,err=90) rnum

      else 
c                                 fraction write numerator
         if (iback+1-ibeg.gt.30) goto 90
c                                 first number
         write (num,'(30a)') chars(ibeg:iback)       
         read (num,*,err=90) rnum
c                                 second number 
         if (iend-iback-1.gt.30) goto 90
         write (num,'(30a)') chars(iback+2:iend)      
         read (num,*,err=90) rnum1

         rnum = rnum/rnum1

      end if 

      return

90    ier = 2

      end

      subroutine getnam (name,ids)
c----------------------------------------------------------------------
c subroutine to retrieve phase name corresponding to index ids
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids

      character name*14

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c-----------------------------------------------------------------------
      if (ids.lt.0) then
c                                 simple compound:
         name = names(-ids)

      else  
c                                 solution phases:
         if (iopt(24).eq.0.or.lname(ids).eq.'unclassified') then
c                                 use model name
            name = fname(ids)

         else if (iopt(24).eq.1) then
c                                 use phase abbreviation
            name = aname(ids) 

         else 
c                                 use full name
            name = lname(ids)

         end if 

      end if 
      
      end 

      subroutine sopen 
c-----------------------------------------------------------------------
c simple file open for data echo programs (e.g., actcor, rewrit, ctransf)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character*100 n2name

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
c                                 first the thermo data file
      do 

         call fopen2 (2,n2name)
 
         if (iam.eq.6) then 
            write (*,1070) 'ctransf.dat'
            open (n8,file='ctransf.dat')
         else if (iam.eq.9) then 
            write (*,1070) 'actcor.dat'
            open (n8,file='actcor.dat')
         else if (iam.eq.10) then 
            write (*,1070) 'new_'//n2name
            open (n8,file='new_'//n2name)
         end if 

         exit

      end do 
 
1070  format (/,'Output will be written to file: ',a,/)
 
      end

      block data
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision ptx
      integer ipt2
      common/ cst32 /ptx(l5),ipt2

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)

      character*80 com
      common/delet/com 

      integer hs2p
      double precision hsb
      common/ cst84 /hsb(i8,4),hs2p(6)

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      double precision vrt
      integer irt
      logical sroot
      common/ rkroot /vrt,irt,sroot

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)
c-----------------------------------------------------------------------
      data hs2p/4, 5, 18, 19, 20, 21/

      data iff,ipt2,goodc,badc/3*0,6*0d0/
c
      data mu, uf/ k8*0d0, 2*0d0/

      data r/8.3144126d0/

      data gflu, sroot/ 2*.false./

      data com/' '/
c                                 tags for thermo data i/o
      data strgs/'G0 ','S0 ','V0 ','c1 ','c2 ','c3 ','c4 ','c5 ','c6 ',
     *           'c7 ','b1 ','b2 ','b3 ','b4 ','b5 ','b6 ','b7 ','b8 ',
     *           'b9 ','b10','c8 ','c9 ','c10','c11',
     *           'Tc ','B  ','p  ','v  ','cs1','cs2','cs3','cs4'/
      data mstrg/'m0','m1','m2','k0','k1','k2'/
      data dstrg/'d1','d2','d3','d4','d5','d6','d7','d8','d9'/
      data tstrg/'t1 ','t2 ','t3 ','t4 ','t5 ','t6 ','t7 ','t8 ','t9 ',
     *           't10','t11','t12','t13','t14','t15'/
      data e16st/'G0 ','S0 ','V0 ','Cp0','w ','q ','a1 ','a2 ','a3 ',
     *           'a4 ','c1 ','c2 ','HOH'/
c     data estrg/'eG0','eS0','eV0','ec1','ec2','ec3','ec4','ec5','ec6',
c    *           'ec7','eb1','eb2','eb3','eb4','eb5','eb6','eb7','eb8'/
c                                 tags for interaction coefficients (Redlich-Kister polynomial)
      data wstrg/'w0 ','wT ','wP ','wP1','wP2','wP0'/
c                                 fluid eos species
      data specie /
     *      'H2O ','CO2 ','CO  ','CH4 ','H2  ','H2S ','O2  ',
     *      'SO2 ','COS ','N2  ','NH3 ','O   ','SiO ','SiO2',
     *      'Si  ','C2H6','HF  '/

      data times,btime,etime/90*0d0/

      end

      subroutine getphi (name,aq,eof)
c----------------------------------------------------------------------
c read phase data from the thermodynamic data file from lun N2, assumes
c topn2 has read the header section of the file.

c on input aq is a flag which determines if solute species data is
c accepted (ieos = 15-16).
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i, it, j, ier

      double precision ct

      logical eof, aq

      character key*22, val*3, name*8,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      eof = .false.

      do 

         call redcd1 (n2,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (ier.lt.0) then
 
            eof = .true.
            exit

         else if (ier.gt.0) then
                 
            call error (23,ct,i,name)

         end if 
c                                 name          
         read (key,'(a)',iostat=ier) name
         if (ier.ne.0) exit
c                                 on the off chance of a loose end keyword
         if (key.eq.'end') cycle
c                                 EoS
         read (nval2,*,iostat=ier) ieos
         if (ier.ne.0) exit    
c                                 look for comments
c        write (commnt,'(80a)') chars(com:com+79)
c                                 composition
         call formul (n2)
c                                 thermodynamic data
         call indata (n2)
c                                 do component transformation if
c                                 itrans is not zero
         if (itrans.gt.0) then
 
            do i = 1, itrans
               it = ictr(i)
               if (comp(it).ne.0d0.and.ctrans(it,i).ne.0d0) then
c                                 ct is how much of the new
c                                 component is in the phase.
                  ct =  comp(it) / ctrans(it,i)
 
                  do j = 1, icmpn
                     comp(j) = comp(j) - ct * ctrans(j,i)
                  end do 
 
                  comp(it) = ct
               end if 
            end do
         end if

         if (.not.aq.and.(ieos.eq.15.or.ieos.eq.16)) cycle

         if (ieos.gt.0.and.ieos.lt.5.and.thermo(3,k10).eq.0d0) then 
c                                 standard form with no volumetric EoS, 
c                                 reset ieos internally:
            ieos = 0
         end if 
 
         exit 

      end do

      end

      subroutine indata (lun)
c----------------------------------------------------------------------
c called by getphi to decompile thermodynamic data cards read from lun.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lun, ier, iscan, iscnlt, i, j, ibeg, iend, ic2p(k4), jkind

      character key*22, values*80, strg*80

      double precision var

      logical ok

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      double precision emodu
      common/ cst318 /emodu(k15)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(k4),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer ic
      common/ cst42 /ic(k0)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer iam
      common/ cst4 /iam

      save ic2p
      data ic2p/31,32,22,1,2,3,4,5,6,7,12,13,14,15,16,17,18,19,20,21,8,
     *          9,10,11,23,24,25,26,27,28,29,30/
c-----------------------------------------------------------------------
c                                 initialize data
c                                 flag for t-dependent disorder
      idiso = 0
c                                 flag for mock-lambda transitions
      ilam = 0 
c                                 counter of mock-lambda transitions
      jlam = 0 
c                                 flag for bulk + shear moduli
      ikind = 0
c                                 flag for just shear
      jkind = 0
c                                 hsc conversion 
      hsc(k10) = .false.
c                                 standard thermo parameters
      do i = 1, k4
         thermo(i,k10) = 0d0
      end do 
c                                 shear modulus
      do i = 1, k15
         emodu(i) = 0d0
      end do
c                                 lamda transitions
      do j = 1, m6
         do i = 1, m7
            tm(i,j) = 0d0
         end do
      end do
c                                 t-dependent disorder
      do i = 1, m8
         td(i) = 0d0
      end do 

      do 
c                                 find a data card
         call redcd0 (lun,ier,key,values,strg)
         if (ier.ne.0) call error (23,tot,ier,strg) 

         ibeg = 1

         if (key.eq.'end') then 

            exit 

         else if (key.eq.'transition') then 

            ibeg = iscan (ibeg,com,'=') + 1
            ibeg = iscnlt (ibeg,com,' ')
            iend = iscan (ibeg+1,com,'=') + 1
c                                 write ilam data to values
            write (values,'(80a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg)
c                                 ilam as read is the counter, code
c                                 currently assumes the data is entered
c                                 sequentially, therefore this isn't necessary.
            read (values,*,iostat=ier) ilam
            if (ier.ne.0) call error (23,tot,ier,strg)
c                                 next get the type flag jlam.
            ibeg = iend
            iend = iscnlt (ibeg,com,'9')

            write (values,'(80a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg)
            read (values,*,iostat=ier) jlam
            if (ier.ne.0) call error (23,tot,ier,strg) 
c                                 position for next keyword
            ibeg = iend

         end if
c                                 read remaining keywords and values
c                                 from card
         do 

            key = ''
c                                 locate end of keyword
            if (ibeg.ge.com) exit 
            iend = iscan (ibeg,com,'=') - 1
            if (iend.ge.com) exit
c                                 write keyword
            write (key,'(22a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg) 
c                                 locate data
            ibeg = iscnlt (iend+2,com,' ')
            iend = iscan (ibeg,com,' ')
c                                 write data 
            write (values,'(80a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg) 
c                                 shift pointer to next key
            ibeg = iscnlt(iend,com,' ')
c                                 assign data
            ok = .false.
c                                 =====================================
c                                 thermo data 
            if (ieos.eq.12.or.ieos.eq.14.or.ieos.eq.17) then
c                                 calphad format
               do i = 1, k4
                  if (key.eq.strgs(i)) then 
                     read (values,*,iostat=ier) thermo(ic2p(i),k10)
                     if (ier.ne.0) call error (23,tot,ier,key) 
                     ok = .true.
                     exit 
                  end if 
               end do

            else if (ieos.eq.16) then 
c                                 DEW/HKF aqueous data
               do i = 1, 13
                  if (key.eq.e16st(i)) then 
                     read (values,*,iostat=ier) thermo(i,k10)
                     if (ier.ne.0) call error (23,tot,ier,key) 
                     ok = .true.
                     exit 
                  end if 
               end do

            else 
c                                 generic thermo data 
               do i = 1, 21

                  if (key.eq.strgs(i)) then 

                     read (values,*,iostat=ier) thermo(i,k10)
                     if (ier.ne.0) call error (23,tot,ier,strg) 
                     ok = .true.
                     exit

                  else if (key.eq.'GH') then

                     read (values,*,iostat=ier) thermo(1,k10)
                     if (ier.ne.0) call error (23,tot,ier,strg)
                     hsc(k10) = .true.

                     if (hscon.and.iam.ne.5) then 
c                                 convert HSC G0 to SUP G0
                        do j = 1, icomp
                           thermo(1,k10) = thermo(1,k10) 
     *                                   + tr*comp(ic(j))*sel(j)
                        end do

                     else if (hscon) then 

                        do j = 1, icmpn
                           thermo(1,k10) = thermo(1,k10) 
     *                                   + tr*comp(j)*sel(j)
                        end do

                     end if 
 
                     ok = .true.
                     exit

                  end if 

               end do 

            end if 

            if (ok) cycle
c                                 =====================================
c                                 shear mod data 
            do i = 1, 6

               if (key.eq.mstrg(i)) then 
c                                 set shear/bulk mod flag
                  if (i.lt.4) then
c                                 shear
                     ikind = 1

                  else if (i.gt.3) then 
c                                 bulk
                     jkind = 1

                  end if 
                  
                  read (values,*,iostat=ier) emodu(i)
                  if (ier.ne.0) call error (23,tot,ier,strg) 
                  ok = .true.

                  exit

               end if

            end do

            if (ok) cycle
c                                 =====================================
c                                 explicit temperature-dependent disorder data 
            do i = 1, m8
               if (key.eq.dstrg(i)) then 
c                                 set disorder flag
                  idiso = 1
                  read (values,*,iostat=ier) td(i)
                  if (ier.ne.0) call error (23,tot,ier,strg) 
                  ok = .true.
                  exit 
               end if 
            end do 

            if (ok) cycle
c                                 =====================================
c                                 mock-lambda transition data 
            do i = 1, m7
               if (key.eq.tstrg(i)) then 
                  read (values,*,iostat=ier) tm(i,ilam)
                  if (ier.ne.0) call error (23,tot,ier,strg) 
                  ok = .true.
                  exit 
               end if 
            end do 

            if (ok) cycle

            call error (9,var,i,key)

         end do

      end do
c                                 reset ikind flag:
c                                 ikind = 0, no explicit moduli
c                                 ikind = 1 and jkind = 0 set ikind = 1 => just shear
c                                 ikind = 1 and jkind = 1 set ikind = 2 => both
c                                 ikind = 0 and jkind = 1 set ikind = 3 => just bulk
      ikind = ikind + jkind
      if (ikind.eq.1.and.jkind.eq.1) ikind = 3

      end

      subroutine getkey (n,ier,key,values,strg)
c----------------------------------------------------------------------
c getkey calls redcd0 and outputs error message 21 on error.
c----------------------------------------------------------------------
      implicit none
 
      integer n, ier 

      character key*22, values*80, strg*80
c----------------------------------------------------------------------

      call redcd0 (n,ier,key,values,strg)

      if (ier.ne.0) call error (21,0d0,n,strg)

      end 

      subroutine redcd0 (lun,ier,key,values,strg)
c----------------------------------------------------------------------
c this routine seeks a non-blank card that contains data, i.e., something
c other than a comment (text preceded by a "|") character. 
c the first word of the data is saved as key, the remaining words are
c are saved in values, and the complete data is saved as strg.
c the full record (including comments) is saved in chars.

c the card is also loaded into chars with:

c  length - position of last non-blank character
c  com    - position of the comment character
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer lun, ier, iscan, iscnlt, ibeg, iend, iblank

      character card*(lchar), key*(*), values*(*), strg*(*)

      external iscan, iscnlt

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 
      key = ' '

      do 

         read (lun,'(a)',iostat=ier) card

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')
c                                 find the next blank
            iblank = iscan (ibeg,com,' ')
c                                 len < ibeg => only comments
            if (ibeg.ge.com) cycle
c                                 full record length
            length = iscnlt (lchar,1,' ')

            exit 

         else if (ier.ne.0) then 

            exit 

         end if 

      end do 

      if (ier.eq.0) then 
c                                 find end of keyword 
         iend = ibeg + 1
         iend = iscan (iend,lchar,' ') - 1
         if (iend.gt.22) iend = 22
c                                 load chars into key
         write (key,'(22a)') chars(ibeg:iend)
c                                 now the values
         ibeg = iscnlt (iend+1,lchar,' ') 

         if (ibeg.lt.com) then 

            iend = iscnlt (com,ibeg,' ')
            if (iend-ibeg.gt.79) iend = ibeg + 79
c                                 load chars into value
            write (values,'(80a)') chars(ibeg:iend)
c                                 load chars into strg
            if (iend.gt.80) iend = 80
            write (strg,'(80a)') chars(1:iend)
        
         else
c                                 no values
            strg = key

         end if 

      end if 

      end 

      subroutine formul (lun)
c----------------------------------------------------------------------
c formul reads a text formula and coverts the result to the composition
c array comp. allowed component names are read from xcmpnt array to 
c avoid problems associated with component transformations. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lun, len0, len1, ier, iscan, i, ibeg, iend

      character key*22, values*80, strg*80, ctemp*5

      logical ok

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)
c-----------------------------------------------------------------------
      do i = 1, icmpn
         comp(i) = 0d0
      end do

      call getkey (lun,ier,key,values,strg)

      if (ier.ne.0) call error (23,0d0,i,strg)

      ibeg = 1
      iend = iscan (ibeg,lchar,' ') - 1

      do 
c                                 find the "(" and ")"
         len0 = iscan (ibeg,iend,'(') 
         len1 = iscan (len0,iend,')')
c                                 write the name and number
         write (ctemp,'(5a)')   chars(ibeg:len0-1)
c                                 identify the component
         ok = .false.

         do i = 1, icmpn

            if (xcmpnt(i).eq.ctemp) then
               call redfr0 (comp(i),len0+1,len1-1,ier)
               if (ier.eq.0) ok = .true.
               exit 
            end if 

         end do      

         if (.not.ok) call error (23,0d0,i,strg)

         if (len1.eq.iend) exit

         ibeg = len1 + 1

      end do                  

      end 

      subroutine outdat (lun,id,option)
c----------------------------------------------------------------------
c lun    - output LUN
c id     - phase pointer
c option - source of compositional data 

c called by frendly, vertex, ctransf, actcor, rewrit. for processed data 
c requires unver and unlam to recover original data; unlam puts the 
c transition data into the local array tm; other data is the primary 
c arrays (cp or comp[see option below], thermo, therdi)

c if option = 0  then formula for entity id is created from the 
c    composition array comp and name array cmpnt.
c if option = 1 then formula for entity id created from the 
c    composition array cp and name array cmpnt via pointer ic.
c if option = 2  then formula for entity id created from the 
c    composition array cp0 and name array cmpnt.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lun, i, j, ibeg, iend, id, option, jcomp, siz

      character text(14)*1

      double precision var, dg

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer ic
      common/ cst42 /ic(k0)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer eos
      common/ cst303 /eos(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      character*80 commnt
      common/delet/commnt

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)
c-----------------------------------------------------------------------
c                                 =====================================
c                                 name & EoS
      write (lun,*) 
      read (names(id),'(8a)') chars(1:8)
      ibeg = 9
      var = eos(id)
      call outthr (var,' EoS',4,ibeg) 

      if (commnt.ne.' ') then 
         chars(ibeg) = '|'
         read (commnt,'(80a)') chars(ibeg+1:ibeg+80)
         ibeg = ibeg + 80
      end if 

      write (lun,'(400a)') chars(1:ibeg)
c                                 =====================================
c                                 formula
      ibeg = 1
      iend = 0

      if (option.eq.1) then 
         jcomp = icomp
      else 
          jcomp = icmpn
      end if 

      dg = 0d0 

      do i = 1, jcomp

         if (option.eq.0) then
            var = comp(i)
         else if (option.eq.1) then 
            var = cp(i,id)
         else if (option.eq.2) then 
            var = cp0(i,id)
         end if 

         if (var.ne.0) then 
c                                 load text name
            iend = ibeg + cl(ic(i)) - 1

            read (cmpnt(ic(i)),'(5a)') chars(ibeg:iend)
c                                 left parenthesis
            chars(iend + 1) = '('
c                                 get number
            call numtxt (var,text,siz)
c                                 load number into chars
            ibeg = iend + 2
            iend = ibeg + siz - 1

            do j = ibeg, iend
               chars(j) = text(j-ibeg+1)
            end do
c                                get the delta g HSC correction
            dg = dg + var*sel(i)

            chars(j) = ')'
 
            ibeg = iend + 2

         end if 

      end do 
c                                 write the formula
      write (lun,'(400a)') chars(1:iend+1)
c                                 =====================================
c                                 thermo data
      if (eos(id).eq.16) then 
c                                 HKF aqueous electrolyte data (13 values)
         ibeg = 1
 
         do i = 1, 5
            call outthr (thermo(i,id),e16st(i),3,ibeg)
         end do

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

         ibeg = 1
 
         do i = 6, 10
            call outthr (thermo(i,id),e16st(i),3,ibeg)
         end do

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

         ibeg = 1
 
         do i = 11, 13
            call outthr (thermo(i,id),e16st(i),3,ibeg)
         end do

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      else if (eos(id).eq.12.or.eos(id).eq.14.or.eos(id).eq.17) then 

         call error (72,r,id,'routine OUTDAT is not programmed to '//
     *              'output CALPHAD data tags.')

      else 

         ibeg = 1

         if (hscon.and.hsc(id)) then
c                                 convert back to HSC apparent G
            call outthr (thermo(1,id) - tr*dg,'GH',2,ibeg)

         else if (hsc(id)) then 
c                                 direct output of HSC apparent G
            call outthr (thermo(1,id),'GH',2,ibeg)

         else

            call outthr (thermo(1,id),strgs(1),2,ibeg)

         end if
 
         do i = 2, 3
            call outthr (thermo(i,id),strgs(i),2,ibeg)
         end do
c                                 write G,S,V
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 c1->c7 of thermo data
         ibeg = 1
  
         do i = 4, 10
            call outthr (thermo(i,id),strgs(i),2,ibeg)
         end do
c                                 write c1->c7
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 b1->b10 of thermo data
         ibeg = 1

         do i = 11, 20
            call outthr (thermo(i,id),strgs(i),2,ibeg)
         end do
c                                 write b1->b8, c8
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      end if 
c                                 =====================================
c                                 shear/bulk modulus
      ibeg = 1

      do i = 1, 6
         call outthr (emod(i,id),mstrg(i),2,ibeg)
      end do

      if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 =====================================
c                                 disorder parameters
      if (idis(id).ne.0) then

         ibeg = 1

         do i = 1, 8
            call outthr (td(i),dstrg(i),2,ibeg)
         end do 

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      end if 
c                                 =====================================
c                                 transition parameters
      do i = 1, lct(id)

         ibeg = 1
         var = i
         call outthr (var,'transition',10,ibeg)
         var = ltyp(id)
         call outthr (var,'type',4,ibeg) 

         do j = 1, m7
            call outthr (tm(j,i),tstrg(j),3,ibeg)
         end do 

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      end do 

      write (lun,1000)

1000  format ('end',/)

      end

      subroutine outthr (num,strg,siz,ibeg)
c----------------------------------------------------------------------
c output prettified data.

c    num - the numeric data 
c    strg - a text tag for the data
c    siz  - length of strg
c    ibeg - pointer to the location for the data in the output array (chars)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision num

      character strg*(*), text(14)*1

      integer i, ibeg, iend, siz, len0, jend

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      if (num.ne.0d0) then 
c                                 pad with 2 blanks, if not at line begining
         if (ibeg.gt.1) then
            chars(ibeg) = ' '
            ibeg = ibeg + 1
         end if

         iend = ibeg + siz - 1

         read (strg,'(14a)') chars(ibeg:iend)
c                                 trim out trailing blanks
         jend = ibeg

         do i = ibeg + 1, iend
            if (chars(i).eq.' ') cycle
            jend = jend + 1
         end do

         iend = jend

         chars(iend+1) = ' '
         chars(iend+2) = '='
         chars(iend+3) = ' '

         call numtxt (num,text,len0)

         do i = 1, len0
            chars(iend+3+i) = text(i)
         end do

         chars(iend+3+i) = ' '
         chars(iend+4+i) = ' '

         ibeg = iend + 4 + i

      end if 

      end 

      subroutine numtxt (num,text,siz)
c----------------------------------------------------------------------
c convert a g14.7e2 number to simplest possible text
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision num

      character text(*)*1, strg*14

      logical dec

      integer i, siz, inum, ier, ibeg, iend, jscnlt, jscan

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      inum = int(num)

      siz = 14

      if (dabs(num-inum).lt.zero) then 
c                                 the number can be represented as 
c                                 an integer
         write (strg,'(i14)',iostat=ier) inum

      else 

         write (strg,'(g14.7E2)',iostat=ier) num

      end if

      read (strg,'(14a)') text(1:siz)

      ibeg = jscnlt (1,siz,' ',text)
      iend = jscan (ibeg,siz,' ',text) - 1
c                                 shift text left
      siz = 0 

      dec = .true.

      do i = ibeg, iend

         siz = siz + 1
         text(siz) = text(i)

         if (text(siz).gt.'A') dec = .false. 

      end do 
c                                 pruning:
      if (text(1).eq.'0') then
c                                 cut leading zero/+
         do i = 1, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1 

      else if (text(1).eq.'-'.and.text(2).eq.'0') then
c                                 cut leading zero
         do i = 2, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1

      end if

      if (dec) then 
c                                decimal number
         iend = jscan (1,siz,'.',text)
c                                reduce len to cut trailing zeroes
         if (iend.lt.siz) siz = jscnlt (siz,iend,'0',text)

      else if (num-inum.ne.0d0) then 
c                                 find the E char
         iend = jscnlt (1,siz,'A',text)
         ibeg = jscnlt (iend-1,1,'0',text) + 1
         inum = iend - ibeg
c             
         do i = ibeg, siz - inum
            text(i) = text(i + inum)
         end do   

         siz = siz - inum
c                                 the E character is now at
         ibeg = iend - inum   

         if (text(ibeg+1).eq.'+') then

            inum = 1 
            if (text(ibeg+2).eq.'0') inum = 2
c                                 delete superfluous + and 0
            do i = ibeg+1, siz - inum
               text(i) = text(i + inum)
            end do
          
            siz = siz - inum

         else if (text(ibeg+1).eq.'-') then
c                                 delete superfluous 0
            if (text(ibeg+2).eq.'0') then
               do i = ibeg+2, siz - 1
                  text(i) = text(i + 1)
               end do
           
               siz = siz - 1

            end if 

         end if 
       
      end if 

      end 

      subroutine znmtxt (num,text,siz)
c----------------------------------------------------------------------
c convert a f7.3 number to simplest possible text, 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision num

      character text(*)*1, strg*7

      integer i, siz, inum, ier, ibeg, iend, jscnlt, jscan

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      inum = int(num)

      siz = 7

      if (dabs(num-inum).lt.zero) then 
c                                 the number can be represented as 
c                                 an integer
         write (strg,'(i7)',iostat=ier) inum

      else 

         write (strg,'(f7.4)',iostat=ier) num

      end if

      read (strg,'(7a)') text(1:siz)

      ibeg = jscnlt (1,siz,' ',text)
      iend = jscan (ibeg,siz,' ',text) - 1
c                                 shift text left
      siz = 0 

      do i = ibeg, iend

         siz = siz + 1
         text(siz) = text(i)

      end do 
c                                 pruning:
      if (text(1).eq.'0') then
c                                 cut leading zero/+
         do i = 1, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1 

      else if (text(1).eq.'-'.and.text(2).eq.'0') then
c                                 cut leading zero
         do i = 2, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1

      end if

      do i = siz + 1, 7
         text(i) = ' '
      end do

      iend = jscan (1,siz,'.',text)
c                                reduce len to cut trailing zeroes
c     if (iend.lt.siz) siz = jscnlt (siz,iend,'0',text)

      end

      subroutine fopen1 
c-----------------------------------------------------------------------
c fopen1 gets the project name and opens the problem definition file
c    n1name = project_name.dat
c iam is a flag indicating the calling program:
c    4 - build
c    1 - vertex
c    2 - meemum
c   13 - unsplt, global
c   14 - unplst, local
c------------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character y*1,n1name*100

      integer ierr

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iam
      common/ cst4 /iam

      integer jx, jy, lev, xn, yn
      common/ cst58 /jx, jy, lev, xn, yn
c-----------------------------------------------------------------------
      do 
c                                 get the root for all output files
c                                 except if unsplt-local
         if (iam.ne.14) then 

            if (iam.eq.4) then 
c                                 BUILD
               write (*,1040)
c                                 readrt loads the root into prject
               call readrt

            else  
c                                 VERTEX, MEEMUM, and plotting programs
               write (*,1030) 
               call readrt

            end if 

         end if 
c                                 make the problem definition file name
         call mertxt (n1name,prject,'.dat',0)

         if (iam.eq.4) then 

            write (*,1070) n1name
c                                 BUILD
            open (n1,file=n1name,iostat=ierr,status='new')

            if (ierr.ne.0) then
c                                 name exists
               write (*,1050) n1name
               read (*,'(a)') y

               if (y.eq.'Y'.or.y.eq.'y') then 
c                                 overwrite it
                  open (n1,file=n1name)

               else
c                                 try again 
                  cycle 

               end if 

            end if
         
         else 
c                                 VERTEX, MEEMUM, UNSPLT
            open (n1,file=n1name,iostat=ierr,status='old')

            if (ierr.ne.0) then
c                                 name does not exist
               write (*,1080) n1name
               read (*,'(a)') y
               if (y.eq.'Y'.or.y.eq.'y') then 
c                                 try again
                  cycle 

               else 
c                                 quit
                  stop 
 
               end if 

            end if

            if (iam.eq.13) then
c                                 unsplt, read my_project.spt
               call mertxt (tfname,prject,'.spt',0)
               open (n8,file=tfname,iostat=ierr,status='old')
               if (ierr.ne.0) then
c                                 file does not exist
                  call error (12,0d0,ierr,tfname)
               end if 

               read (n8,*,iostat=ierr) jx
               if (ierr.ne.0) call error (12,0d0,ierr,tfname)
               read (n8,*,iostat=ierr) jy
               if (ierr.ne.0) call error (12,0d0,ierr,tfname)

            end if
         end if

         exit 

      end do 

1030  format (/,'Enter the project name (the name assigned ',
     *        'in BUILD) [default = my_project]:')
1040  format (/,'Enter a name for this project (the name',
     *        ' will be used as the',/,'root for all output file names)'
     *       ,' [default = my_project]:')
1050  format (/,'The file: ',a,/,'exists, overwrite it (y/n)?')
1070  format (/,'The problem definition file will be named: ',a) 
1080  format (/,'**warning ver191** no problem definition file named: ',
     *       a,/,'Run BUILD to create the file or change project names.'
     *       ,//,'Enter a different project name (y/n)?')

      end 

      subroutine readrt
c----------------------------------------------------------------------
c readrt - read file name root from terminal
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'
 
      integer kscan, iscnlt, ierr, siz

      character*100 prject, tfname
      common/ cst228 /prject,tfname

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      do 

         read (*,'(a)') prject

         if (prject.ne.' ') then 

            read (prject,'(100a)') chars(1:100)
c                                 find end of name ' '
            length = iscnlt (100,1,' ') 
c                                 check length
            if (length.gt.90) then 
               write (*,1010) 
               cycle 
            end if 
c                                 look for path characters / or \
            siz = kscan (100,1,'/')
            if (siz.eq.0) siz = kscan (100,1,'\')

            if (siz.eq.length) then 
               write (*,1030)
               cycle
            end if 
c                                 check if directory is valid
            if (siz.ne.0) then

               write (tfname,'(100a)') chars(1:siz)
               call mertxt (tfname,tfname,'delete_me',0)

               open (n1,file=tfname,iostat = ierr)
               close (n1,status='delete')

               if (ierr.ne.0) then 
                  write (*,1040)
                  cycle 
               end if  
c                                 mertxt uses chars, so re-read chars
               read (prject,'(100a)') chars(1:100)

            end if
c                                 look for illegal "." character
            if (kscan(siz+1,length,'.').lt.length) then 
               write (*,1000)
               cycle 
            end if 
c                                 look for illegal " " character
            if (kscan(siz+1,length,' ').lt.length) then 
               write (*,1020)
               cycle
            end if 

         else

            prject = 'my_project'

         end if 

         exit

      end do 

1000  format (/,'file/project names cannot include . characters, ',
     *          'try again',/)
1010  format (/,'file/project names must be < 91 characters, '
     *         ,'try again',/)
1020  format (/,'file/project names cannot include blanks, ',
     *          'try again',/)
1030  format (/,'file/project names cannot end with a / or \ character',
     *        ', try again',/)
1040  format (/,'the path specified in your project name is invalid,',
     *          ' check that all the ',/,
     *          'directories in the path exist, try again.',/)
      end

      subroutine getrt
c----------------------------------------------------------------------
c getrt - extracts root file name from the full file name in tfname
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'
 
      integer kscan, siz

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      read (tfname,'(100a)') chars(1:100)
c                                 find end of name ' '
      length = kscan (1,100,' ') - 1
c                                 look for dot character
      siz = kscan (length,1,'.') - 1

      if (siz.le.0) siz = length

      write (prject,'(100a)') chars(1:siz)

      end

      subroutine fopen2 (jam,name)
c-----------------------------------------------------------------------
c fopen2 - choose and open a thermodynamic data file on unit n2, jam 
c indicates behavior required by the calling program:
c  0 - name passed as argument, error if not found.
c  1 - BUILD, queries for name and writes it to N1
c  2 - queries for name
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*100 name, y*1, ddata*14, text*140

      integer ierr, jam

      data ddata/'hp02ver.dat   '/
c-----------------------------------------------------------------------

      do 

         if (jam.ne.0) then 
            write (*,1000)
            read (*,'(a)') name
            if (name.eq.' ') name = ddata
         end if 

         open (n2,file=name,iostat=ierr,status='old')

         if (ierr.ne.0) then
c                                 system could not find the file
            if (jam.eq.0) call error (120,0d0,n2,name)
c                                 if not vertex allow another try
            write (*,1010) name
            read (*,'(a)') y

            if (y.ne.'Y'.and.y.ne.'y') then
               write (*,1060)
               stop
            end if             
c                                 try again
            cycle

         end if
 
         if (jam.ne.1) exit 
c                                 BUILD, echo name to n1: 
         call mertxt (text,name,'thermodynamic data file',5)
         write (n1,'(a)') text 

         exit

      end do 
 
1000  format (/,'Enter thermodynamic data file name',
     *          ' [default = hp02ver.dat]:')
1010  format (/,'**warning ver191** FOPEN2 cannot find file:',/,a
     *         ,//,'try again (y/n)?')
1060  format (/,'O.K., I quit too.')

      end

      subroutine mertxt (text,text1,text2,nblank)
c----------------------------------------------------------------------- 
c mertxt - merge text two strings cutting trailing/leading blanks, with 
c          nblank characters between the two strings if text1 eq 
c          blank then pad by an initial 40 blanks.

c     input  - text1, text2 - character strings
c              nblank - number of blanks between the strings in text
c     output - text - character string 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer nchar1, nchar2, nblank

      character text*(*), text1*(*), text2*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      chars(1:lchar) = ' '
c                                 strip leading blanks in text1 and
c                                 get pointer to end of string
      call leblnk (text1,1,nchar1)

      if (nchar1.le.0) then 
c                                 text1 is blank
         nchar1 = 40

      else 
c                                 put nblank blanks between 1st and 
c                                 2nd strings, this is necessary despite
c                                 initialization because leblnk may 
c                                 shift strings left.
         chars(nchar1+1:nchar1 + nblank) = ' '

      end if
c                                 nchar1 points to the first empty char
      nchar1 = nchar1 + nblank + 1 
c                                 strip leading blanks from string in text2 and
c                                 get pointer to end of string in chars
      call leblnk (text2,nchar1,nchar2)

      text = ' '

      if (nchar2.gt.len(text)) call error (10,0d0,len(text),text2)

      write (text,'(400a)') chars(1:nchar2)

      end

      subroutine gettrn (jopt)
c----------------------------------------------------------------------
c jopt = 3 -> build
c jopt = 5 -> ctransf
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,jopt,ict,ier,jscan
 
      character*5 pname, rname, y*1

      double precision sum, ssum

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec
c-----------------------------------------------------------------------
c                                 recombine components:
      do 
         write (*,1030)
         write (*,1040) (cmpnt(i), i = 1, icmpn)
         write (*,1050)
         read (*,'(a)') y
         if (y.ne.'Y'.and.y.ne.'y') exit

         write (*,1060)
         read (*,'(a)') pname
         if (pname.eq.' ') exit
c                                 get the identity of the real comp
c                                 to be replaced.
50       write (*,1070) pname
         read (*,'(a)') rname

         do i = 1, icmpn
            if (rname.eq.cmpnt(i)) then
c                                 matches a name, check if it's a 
c                                 special component
               do j = 1, ispec

                  if (i.ne.idspe(j)) cycle 
c                                 matches special component id(j) 
                  if (jopt.eq.3) then
c                                 don't allow build users 
c                                 to transform saturated
c                                 phase components
                     call warn (14,atwt(1),i,cmpnt(i))
                     goto 60
                  else 
c                                 ctransf, ask the user if the 
c                                 new component will be a special 
c                                 component
                     write (*,1010) cmpnt(i),pname
                     read (*,'(a)') y
                     if (y.eq.'y'.or.y.eq.'Y') cycle
                     idspe(j) = 0 

                  end if 
               end do 

               icout(1) = i
               goto 70

            end if
         end do 
 
60       write (*,1080) 
         write (*,1040) (cmpnt(i), i = 1, icmpn)
         goto 50
c                                 get the identities of the other 
c                                 components in the new component:
70       itrans = itrans + 1
         ict = 1
         if (itrans.gt.k0) call error (999,atwt(1),ict,'GETTRN')
       
         write (*,4050) k5-1,pname
30       read (*,'(a)') rname
         if (rname.eq.'     ') goto 80

         do i = 1, icmpn
            if (rname.eq.cmpnt(i)) then 
               ict = ict + 1
               icout(ict) = i
               goto 30
            end if
         end do 
c                                 no match, try again message
         write (*,2300)
         goto 30
c                                 get the component stoichiometries:
80       write (*,4030) (cmpnt(icout(i)),i=1,ict)
         write (*,4040) pname

         do 
            read (*,*,iostat=ier) (ctrans(icout(i),itrans), i= 1, ict)
            if (ier.eq.0) exit
            call rerr
         end do 
 
         write (*,1100) pname,(ctrans(icout(i),itrans),
     *                      cmpnt(icout(i)), i = 1, ict)
         write (*,1110)
         read (*,'(a)') y
 
         if (y.eq.'y'.or.y.eq.'Y') then
            sum = 0d0
            ssum = 0d0             
            do i = 1, ict
               sum = sum + ctrans(icout(i),itrans) * atwt(icout(i))
               ssum = ssum + ctrans(icout(i),itrans) * sel(icout(i))
            end do 
            atwt(icout(1)) = sum
            sel(icout(1)) = ssum
            cmpnt(icout(1)) = pname
            cl(icout(1)) = jscan(1,5,' ',pname) - 1
            tcname(itrans) = pname
            ictr(itrans) = icout(1)
         else
            itrans = itrans - 1
            write (*,1000) 
         end if

      end do 

1000  format ('Try again.')
1010  format (/,a,' is a possible saturated phase component. Is ',
     *        'the new component ',a,/,'also a possible saturated ',
     *        'phase component (Y/N)?')
1030  format (/,'The current data base components are:')
1040  format (12(1x,a))
1050  format ('Transform them (Y/N)? ')
1060  format ('Enter new component name, < 6 characters,',
     *          ' left justified: ')
1070  format ('Enter old component to be replaced',
     *          ' with ',a,': ')
1080  format ('Select the component from the set: ')
1100  format (1x,a,' = ',6(f6.2,1x,a),/,9x,6(f6.2,1x,a))
1110  format ('Is this correct (Y/N)? ')
2300  format (/,'You made a mistake, try again.',/
     *          'Check spelling and upper/lower case matches.',/)
4030  format ('Enter stoichiometric coefficients of:',/,
     *        2x,12(a,1x))
4040  format ('in ',a,' (in above order): ')
4050  format ('Enter other components (< ',i2,') in ',a,' 1 per',
     *        ' line, <enter> to finish:')

      end

      subroutine topn2 (option)
c----------------------------------------------------------------------
c topn2 reads the header of the thermodynamic data file, if option > 3 
c echoes data to n8
c
c     option  calling program
c       0     vertex, meemum (all calls from clib.f)
c       1     frendly 
c       2     build (2nd, 3rd calls)
c       3     build (first call) 
c       4     actcor
c       5     ctransf
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character tag*4, string*140, key*22, values*80, strg*80

      integer option, i, j, ier, iscan

      double precision sum, ssum

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
  
      character vname*8, xname*8
      common/ csta2 /xname(k5),vname(l2)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct
c-----------------------------------------------------------------------
      rewind n2
c                                 frendly or actcor is reading
c                                 the header, no transformations
      if (option.eq.1.or.option.eq.4) itrans = 0
c                                 test for old (pre 4/2010) data file format:
      read (n2,*,iostat=ier) i
      if (ier.eq.0) call error (8,r,i,dname)

      rewind n2
c                                 database name
      call getkey (n2,ier,key,values,strg)

      dname = strg
c                                 extrinsic variable names & reference values
c                                 read "begin"
      call getkey (n2,ier,key,values,strg)

      do i = 1, l2

         call getkey (n2,ier,key,values,strg)
c                                 over ride default mobile component naming to 
c                                 allow build to use fugacities or activities.
         if (option.gt.3.or.i.lt.4) read (key,'(a8)') vname(i)
         read (values,*) v(i), delt(i)

      end do
c                                 reset the chemical potential tolerances if the
c                                 potential is to be replaced by activity/fugacity
      do i = 1, ipot
         if (jv(i).gt.3) then
            if (imaf(jv(i)-3).ne.1) delt(jv(i)) = delt(jv(i))/1e4
         end if
      end do 
c                                 set log p variable name.
      if (icopt.gt.4.and.lopt(14)) vname(1) = 'log[P,b]'
c                                 read end key
      call getkey (n2,ier,key,values,strg)
c                                  set reference conditions
      pr = v(1)
      tr = v(2)
c                                  this block of code probably never gets executed?
      if (option.lt.4) then 
         if ((ifug.ge.10.and.ifug.le.12).or.
     *       ifug.eq.15.or.ifug.eq.17.or.ifug.eq.18) then 
            vname(3) = ' X(O) ' 
         else if (ifug.eq.25) then 
            vname(3) = 'Y(CO2)*'
         else if (ifug.eq.13) then 
            vname(3) = 'X(H2)'
         end if 
      end if  
c                                 read tolerance dtol
      call getkey (n2,ier,key,values,strg)

      read (values,*) dtol
      dtol = -dabs(dtol)
c                                 set utol and ptol, the 
c                                 tolerances (in energy units) for determination of
c                                 the stability of divariant, univariant
c                                 and invariant equilibria or reactions.
c                                 utol must be smaller than -utol
c                                 ptol must be > 2*-dtol
      utol = -dtol/1d1
      ptol = -dtol*3d0 

      hscon = .false.
      oxchg = .false.

      do i = 1, k0
         sel(i) = 0d0
         cox(i) = 0d0
      end do

      do
c                                 component names, formula weights, and, optionally, 3rd law s_elements
         call getkey (n2,ier,key,values,strg)
c                                 look for optional HSC_conversion key
         if (key.eq.'HSC_conversion') then 

            hscon = .true.

         else if (key.eq.'reference_oxidation_st') then

            oxchg = .true.

         else if (key.eq.'begin_components') then 

            exit 

         else

            call error (72,utol,i,'invalid thermodynamic data file '//
     *                             'keyword '//key)

         end if

      end do 


      icmpn = 0 

      do

         call getkey (n2,ier,key,values,strg)
         
         if (key.eq.'end_components') exit

         icmpn = icmpn + 1

         read (key,'(a5)') cmpnt(icmpn)
c                                 get component string length
         cl(icmpn) = iscan(1,length,' ') - 1

         if (hscon.and.oxchg) then
            read (values,*) atwt(icmpn), sel(icmpn), cox(icmpn)
         else if (hscon) then 
            read (values,*) atwt(icmpn), sel(icmpn)
         else 
            read (values,*) atwt(icmpn)
         end if 

      end do
c                                 save old names for component transformations
c                                 in ctransf or build.
      do i = 1, k0
         xcmpnt(i) = cmpnt(i)
      end do
c                                 for programs that read the composition 
c                                 vector, check that the cp-dimension (k5)
c                                 is adequate
c     if (option.ne.0.and.option.ne.2.and.icmpn.gt.k5) 
c     *    call error (197,r,icmpn,'TOPN2')
c                                 read special components.
      lopt(7) = .false.
  
      call getkey (n2,ier,key,values,strg)

      if (key.eq.'begin_special_componen') then

         ispec = 0 

         do 
 
            call getkey (n2,ier,key,values,strg)
         
            if (key.eq.'end_special_components') exit

            do j = 1, icmpn
               if (key.eq.cmpnt(j)) then 
                  ispec = ispec + 1
                  idspe(ispec) = j
                  lopt(7) = .true.
                  exit
               end if 
            end do 
           
         end do 

      else
c                                 no saturated phase constraint possible
c                                 set lopt(7) to false. 
         backspace(n2)

      end if 

      if (option.eq.3.or.option.eq.5) then 
c                                 get transformations if build or ctransf
         call gettrn (option)
c                                 check special components
         if (lopt(7)) then 
            j = 0
            do i = 1, ispec
               if (idspe(i).ne.0) then
                  j = j + 1
                  idspe(j) = idspe(i)
               end if
            end do 
            ispec = j
            if (ispec.eq.0) lopt(7) = .false.
         end if 

      else if (option.ne.2) then 
c                                 substitute transformed component names
c                                 and compute the new formula wieghts, this
c                                 is done in gettrns for ioption 3/5.
         do i = 1, itrans

            cmpnt(ictr(i)) = tcname(i)

            sum = 0d0
            ssum = 0d0 

            do j = 1, icmpn
               sum =  sum  + ctrans(j,i) * atwt(j)
               ssum = ssum + ctrans(j,i) * sel(j)
            end do
 
            atwt(ictr(i)) = sum
            sel(ictr(i))  = ssum

         end do 

      end if 

      if (option.gt.3) then 
c                                 echo formatted header data for ctransf/actcor:
         write (n8,1000) 
         write (n8,'(a,a,/)') dname,' |<= data base title'

         write (n8,'(a,a)') 'begin_standard_variables |<= name (<9 ',
     *                      'characters), reference value, tolerance'
         do i = 1, l2 
            write (n8,'(a8,1x,f7.2,3x,g6.1E1)') vname(i),v(i),delt(i)
         end do

         write (n8,'(a,/)') 'end_standard_variables'

         write (n8,'(a,g6.1E1,a,/)') 'tolerance  ',dtol,
     *         '  |<= DTOL for unconstrained minimization, energy units'

         if (hscon) then
            write (n8,'(a,//,a)') 'HSC_conversion |<= tag enabling HSC '
     *                          //'to SUP apparent energy conversion, '
     *                          //'requires elemental entropies in the '
     *                          //'component list below',
     *                         'begin_components | < 6 chars, '//
     *                         'molar weight (g), elemental entropy (R)'
            write (n8,'(a5,2x,f9.4,3x,f9.4)') (cmpnt(i),atwt(i),sel(i),
     *                                        i = 1, icmpn)

         else 

            write (n8,'(a)') 'begin_components | < 6 chars, '//
     *                       'molar weight (g)'
            write (n8,'(a5,1x,f9.4)') (cmpnt(i),atwt(i), i = 1, icmpn)

         end if 

         write (n8,'(a,/)') 'end_components'

         if (lopt(7)) then 
            write (n8,'(a)') 'begin_special_components'
            do i = 1, ispec
               write (n8,'(a)') cmpnt(idspe(i))
            end do 
            write (n8,'(a,/)') 'end_special_components'
         end if  

      end if 
c                                 read and echo unformatted comments and make data                            
      do 

         read (n2,'(a)',iostat=ier) string
         if (ier.ne.0) call error (21,r,i,dname)
         read (string,'(a)') tag

         if (option.gt.3) then
            call mytrim (string)
            write (n8,'(400a)') chars(1:length)
         end if 

         if (string.eq.'begin_makes'.and.option.lt.4) then
 
            call rmakes (option)

            cycle 

         else if (tag.ne.'end') then     

            cycle

         else

            exit

         end if

      end do  

1000  format (/,' | comments are indicated by the | character.',/,
     *     ' | check for warnings at the end of the header section.',/) 

      end 

      double precision function dnan()
c----------------------------------------------------------------------
c george helffrich's function to make a safe NaN
c----------------------------------------------------------------------
      implicit none 

      integer i4

      character c8(8)

      double precision val

      equivalence (c8,i4,val)
c----------------------------------------------------------------------
      i4 = 1

      if (ichar(c8(1)).eq.1) then
C                                 Little-endian
         c8(8)=char(127)
         c8(7)=char(248)
         c8(6)=char(0)
         c8(5)=char(0)
         c8(4)=char(0)
         c8(3)=char(0)
         c8(2)=char(0)
         c8(1)=char(0)

      else
C                                 Big-endian
         c8(1)=char(127)
         c8(2)=char(248)
         c8(3)=char(0)
         c8(4)=char(0)
         c8(5)=char(0)
         c8(6)=char(0)
         c8(7)=char(0)
         c8(8)=char(0)

      endif

      dnan = val

      end

      subroutine plblrb (typ)
c----------------------------------------------------------------------
c write a blurb on plotting options for the calculation
c
c type 1 - tab, ctr, table format
c type 2 - plt, psvdraw format
c type 3 - pts, pspts format
c type 4 - phm, phemgp table format
c----------------------------------------------------------------------
      implicit none

      integer typ
c----------------------------------------------------------------------

      if (typ.eq.1) then
c                                 2d - tab format
         write (*,1000) 
         write (*,1010)
      else if (typ.eq.2) then 
c                                 plt format
         write (*,1020)
      else if (typ.eq.3) then
c                                 pts format
         write (*,1030)
      else if (typ.eq.4) then 
c                                 1d - tab format
         write (*,1000) 
         write (*,1040)
      end if 

1000  format (/,'The tabulated data from this calculation can be ',
     *          'plotted with:',/)
1010  format (5x,'PSTABLE - a Perple_X plotting program',
     *      /,5x,'PYWERAMI - github.com/ondrolexa/pywerami',
     *      /,5x,'PERPLE_X_PLOT - a MATLAB plotting script',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *           'for details of the table format refer to:',/,
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *           '.txt',/)
1020  format (/,'The output from this calculation can be plotted with ',
     *          'PSVDRAW',/)
1030  format (/,'The output from this calculation can be plotted with ',
     *          'PSPTS or converted to',/,'table/plot format with ',
     *          'PT2CURV',/)
1040  format (5x,'pstable - a Perple_X plotting program',
     *      /,5x,'perple_x_plot - a Matlab plotting script',
     *      /,5x,'spread-sheet programs, e.g., Excel',//,
     *          'for details of the table format refer to:',/,
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *           '.txt',/)

      end

c-------------------------------------------------------------------
c text manipulation routines
c
c evntually all such rountines will be in one block or file, this is
c not the case now. 6/20/2011. 
c-------------------------------------------------------------------


      subroutine unblnk (text)
c------------------------------------------------------------------- 
c unblnk - subroutine to remove blanks from text
 
c     text - character string 
c     length - length of unblanked character string, 0 on input 
c             if unknown.
c-------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar

      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c-------------------------------------------------------------------
      nchar = len(text)
 
      read (text,'(400a)') (chars(i), i = 1, nchar)
c                                 scan for blanks:
      length = 0

      do i = 1, nchar
         if (chars(i).eq.' ') cycle 
         length = length + 1
         chars(length) = chars(i)
      end do 

      write (text,'(400a)') (chars(i), i = 1, length)

      end

      subroutine enblnk (text)
c---------------------------------------------------------------------- 
c enblnk - scan text to first blank and cut remaining text.
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ict,nchar
 
      character text*(*), bitsy(lchar)*1 

      nchar = len(text) 
      read (text,'(400a)') (bitsy(i), i = 1, nchar)
c                                 scan for blanks:
      ict = 0

      do i = 1, nchar
         if (bitsy(i).eq.' ') exit
         ict = ict + 1
      end do 

      text = ' '
      write (text,'(400a)') (bitsy(i), i = 1, ict)

      end

      subroutine reblnk (text)
c------------------------------------------------------------------- 
c reblnk - subroutine to replace blanks followed by a character
c          with an underscore.
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c----------------------------------------------------------------------
      implicit none 

      integer i,ict
 
      character*8 text, bitsy(8)*1 
 
      read (text,'(400a)') bitsy
c                                 scan for blanks:
      ict = 0

      do i = 1, 7

         if (i.eq.1.and.bitsy(i).eq.' ') cycle

         if (bitsy(i).eq.' '.and.bitsy(i+1).ne.' ') then
            ict = ict + 1
            bitsy(ict) = '_'
         else if (bitsy(i).eq.' ') then 
            cycle
         else
            ict = ict + 1
            bitsy(ict) = bitsy(i)
         end if

      end do 
 
      bitsy(ict+1) = bitsy(8)

      write (text,'(400a)') (bitsy(i), i = 1, ict + 1)

      end

      subroutine ftext (ist,iend)

c subprogram to filter blanks from text 

      implicit none

      include 'perplex_parameters.h'

      integer ist,iend,i,itic,igot,jend

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      itic = ist - 1
      igot = 0
      jend = iend 

      do i = ist, iend-1

         if (chars(i).eq.' '.and.chars(i + 1).eq.' '.or. 
     *       chars(i).eq.' '.and.chars(i + 1).eq.')'.or. 
     *       chars(i).eq.' '.and.chars(i + 1).eq.'('.or.
     *       igot.eq.0.and.chars(i).eq.' ') cycle
         if (i.gt.ist) then
            if (chars(i-1).eq.'-'.and.chars(i).eq.' ') cycle 
         end if 

         itic = itic + 1
         igot = 1
         chars(itic) = chars(i)

      end do 

      if (chars(iend).ne.' ') then
         itic = itic + 1
         chars(itic) = chars(iend)
      end if 

      iend = itic + 1

      do i = iend, jend 
         chars(i) = ' '
      end do 
 
      end

      subroutine leblnk (text,ibeg,nchar)
c----------------------------------------------------------------------
c leblnk - scan text and strip leading blanks load into chars from
c position ibeg, nchar is the last non blank in char
c non-blank character in stripped string
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ibeg, nchar, ist
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      nchar = len(text) + ibeg -1 
      if (nchar.gt.lchar) nchar = lchar

      read (text,'(400a)') chars(ibeg:nchar)
c                                 find last non-blank   
      
      do i = ibeg, nchar
         if (chars(i).gt.' ') exit
      end do  

      ist = i 

      if (ist.gt.nchar) then 
c                                 all blank
         nchar = 0

      else 

         if (ist.gt.ibeg) then 
c                                 shift chars ist-1 left
            do i = ist, nchar
               chars(i+ibeg-ist) = chars(i)
            end do
         end if  
c                                 find last non-blank
         nchar = nchar + ibeg - ist

         do i = nchar, ibeg, -1
            if (chars(i).gt.' ') exit
         end do 

         nchar = i 

      end if 

      end

      integer function iscan (ibeg,iend,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of char in chars(ibeg..iend), if not
c found iscan is ?; kscan does the forward/inverse scans and sets value
c for bad ranges
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character char*1

      integer ibeg, iend

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      do iscan = ibeg, iend

         if (chars(iscan).eq.char) exit

      end do 

      end 

      integer function kscan (ibeg,iend,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of char in chars(ibeg..iend), if not
c found iscan is iend + inc
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character char*1

      integer ibeg, iend, inc

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      if (ibeg.gt.iend) then 
         inc = -1
      else 
         inc = 1
      end if 

      do kscan = ibeg, iend, inc 

         if (chars(kscan).eq.char) exit

      end do 
c                                 normal loop exit kscan = iend + inc
      end 

      integer function iscnlt (ibeg,iend,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of a character in chars(ibeg..iend) that
c is greater than char. assuming ascii collating sequence +/- < 0 < a
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character char*1

      integer ibeg, iend, inc

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      if (ibeg.le.iend) then 
         inc = 1
      else 
         inc = -1
      end if 

      do iscnlt = ibeg, iend, inc

         if (chars(iscnlt).gt.char) exit

      end do 

      end 

      integer function jscnlt (ibeg,iend,char,chars)
c----------------------------------------------------------------------
c iscan finds the first occurence of a character in chars(ibeg..iend) that
c is greater than char. assuming ascii collating sequence +/- < 0 < a
c----------------------------------------------------------------------
      implicit none

      character char*1, chars(*)*1

      integer ibeg, iend, inc
c----------------------------------------------------------------------

      if (ibeg.le.iend) then 
         inc = 1
      else 
         inc = -1
      end if 

      do jscnlt = ibeg, iend, inc

         if (chars(jscnlt).gt.char) exit

      end do 

      end 

      integer function jscan (ibeg,iend,char,chars)
c----------------------------------------------------------------------
c jscan finds the first occurence of char in chars(ibeg..iend)
c----------------------------------------------------------------------
      implicit none

      character char*1, chars(*)*1

      integer ibeg, iend
c----------------------------------------------------------------------
      do jscan = ibeg, iend

         if (chars(jscan).eq.char) exit

      end do 

      end 

      subroutine mytrim (text)
c----------------------------------------------------------------------
c mytrim - scan text and delete trailing blank characters.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c---------------------------------------------------------------------- 
      nchar = len(text) 

      read (text,'(400a)') chars(1:nchar)
c                                find last non-blank
      length = 1 
      
      do i = 1, nchar
         if (chars(i).gt.' ') length = i
      end do

      end 

      subroutine deblnk (text)
c----------------------------------------------------------------------
c deblnk - scan text and delete multiple blank characters, strip
c out sequential + - or - + operators, trailing -/+ operators, 
c leading blanks, and leading + operators. also deletes blanks
c preceding punctuation (',' and ':' and ';').
 
c     text - character string 

c no test is made for a blank string or a string of "+" signs.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar

      logical strip
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c---------------------------------------------------------------------- 
      nchar = len(text) 

      read (text,1000) chars(1:nchar)
c                                find last non-blank
      length = 1 
      
      do i = 1, nchar
         if (chars(i).gt.' ') length = i
      end do

      nchar = length
c                                 kill any trailing +/- or ','
      if (chars(nchar).eq.'+'.or.chars(nchar).eq.'-'.or.
     *    chars(nchar).eq.',') nchar = nchar - 1
         
c                                 scan for first non blank/+ character:
      length = 0 
      
      do i = 1, nchar
         if (chars(i).eq.' '.or.chars(i).eq.'+') cycle
         length = i
         exit 
      end do 
c                                 shift everything right
      if (length.gt.1) then 

         length = length - 1
         
         do i = length + 1, nchar
            chars(i-length) = chars(i)
         end do 

         nchar = nchar - length

      end if 

      length = 1
      
      do i = 2, nchar
c                                 strip out double blanks
         if ((chars(i).eq.' '.and.chars(i+1).eq.' ').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.':').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.';').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.',').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.')')) cycle

         length = length + 1
         chars(length) = chars(i)

      end do

      nchar = length

      if (nchar.eq.1) return
c                                 strip put + - and - + strings
      strip = .false.

      do i = 1, nchar - 2

         if (chars(i).eq.'+'.and.chars(i+1).eq.'-'.or.
     *       chars(i).eq.'-'.and.chars(i+1).eq.'+') then

             chars(i) = '-'
             chars(i+1) = ' '
             strip = .true.

         else if (chars(i).eq.'+'.and.chars(i+2).eq.'-'.or.
     *            chars(i).eq.'-'.and.chars(i+2).eq.'+') then
c                                allow +/- or -/+
             if (chars(i+1).eq.'/') cycle 

             chars(i) = '-'
             chars(i+2) = ' '
             strip = .true.

         end if 

      end do 
c                                 special cases:
      if (chars(nchar).eq.'*'.and.chars(nchar-1).eq.' '.and.
     *    chars(nchar-2).eq.',') then
          chars(nchar-2) = '*'
          chars(nchar) = ' '
      end if 

      if (strip) then 
c                                 strip out new double blanks
         length = 1

         do i = 2, nchar

            if (chars(i).eq.' '.and.chars(i-1).eq.' ') cycle

            length = length + 1
            chars(length) = chars(i)

         end do 

      end if 

      write (text,1000) chars(1:length),(' ',i = length+1, len(text))
 
1000  format (400a)

      end

      subroutine getstg (text)
c---------------------------------------------------------------------- 
c getstg - subroutine to get first non-blank string  
c          
c     text - character string on input, first non-blank strg on output
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar, ist
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      nchar = len(text) 
      if (nchar.gt.lchar) nchar = lchar

      read (text,'(400a)') chars(1:nchar)
c                                 scan for blanks:
      ist = 1

      do i = 1, nchar
         if (chars(i).eq.' ') cycle 
         ist = i
         exit 
      end do 

      do i = ist, nchar
         if (chars(i).ne.' ') cycle 
         nchar = i-1
         exit 
      end do 

      text = ' '

      write (text,'(400a)') chars(ist:nchar)

      end

      subroutine redlpt (coeffs,ibeg,iend,ier)
c----------------------------------------------------------------------
c redlpt - read coefficients of a linear p-t function:

c                f = c0 + c1 T + c2 P

c from chars array. 

c on input

c    ibeg - the first possible location of the data

c assumes one of two formats:

c    c0 c1 c2

c or

c    c0 +/- ci Tag ....

c where the first character of a < 9 character Tag is used to identify P or T coefficients
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ibeg, iend, ier, iscan, iscnlt, itag

      double precision coeffs(3)

      external iscan, iscnlt

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      do i = 2, 3
         coeffs(i) = 0d0
      end do 
c                                 scan for an equals sign from ibeg
      iend = iscan (ibeg,com,'=') + 1
      if (iend.lt.com) ibeg = iend
c                                 get the first number
      ibeg = iscnlt (ibeg,com,' ') 

      call readfr (coeffs(1),ibeg,iend,com,ier)
      if (ier.ne.0.or.iend+1.ge.com) return

      ibeg = iend + 2
      itag = ibeg
c                                 try reading as though no tags 
c                                 are present (pre-6.7.3)
      do i = 2, 3

         call readfr (coeffs(i),ibeg,iend,com,ier)
         if (ier.ne.0) exit

      end do 

      if (ier.eq.0) return

      do i = 2, 3
         coeffs(i) = 0d0
      end do 
c                                 if an error, numbs/tags must be present
c                                 locate the number
      ibeg = itag
      iend = iscan (ibeg,com,' ') 
c                                 locate the first character of the tag
      itag = iend + 1

      if (chars(itag).eq.'T'.or.chars(itag).eq.'t') then
         i = 2
      else if (chars(itag).eq.'P'.or.chars(itag).eq.'p') then
c                                 must be c2, but check for invalid tag
         i = 3
      else 
         ier = 1
         return
      end if 
c                                 read the number
      call readfr (coeffs(i),ibeg,iend,com,ier)
c                                 the next number, if present begins at
      ibeg = iscan (itag,com,' ') + 1
      iend = iscan (ibeg,com,' ')

      if (ier.ne.0.or.iend.ge.com) return 
c                                 swap indexes
      if (i.eq.2) then 
          i = 3
       else 
          i = 2
      end if 
c                                 read the second tag
      call readfr (coeffs(i),ibeg,iend,com,ier)

      end 

      subroutine nanchk (x,y,text)
c----------------------------------------------------------------------
c nanchk - check x-y coordinate pair for bad_number and resets to 0.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical warn1

      character*(*) text

      double precision x, y

      save warn1
      data warn1/.true./
c----------------------------------------------------------------------
      if (warn1.and.(isnan(x).or.isnan(y))) then 
c                                 
         call warn (61,x,1,text)
         warn1 = .false.

      end if 

      if (isnan(x)) x = 0d0
      if (isnan(y)) y = 0d0

      end 

      subroutine begtim (itime)
c----------------------------------------------------------------------
c begin timer itime
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer itime
c----------------------------------------------------------------------
      call CPU_TIME(btime(itime))

      end 

      subroutine endtim (itime,output,chars)
c----------------------------------------------------------------------
c begin timer itime
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical output

      character chars*(*)

      integer itime
c----------------------------------------------------------------------

      call CPU_TIME(etime(itime))
      times(itime) = times(itime) + (etime(itime)-btime(itime))

      if (output) then 

      write (*,'(/,a,3(2x,g14.7))') chars,
     *                              times(itime), 
     *                              etime(itime)-btime(itime)
      write (666,'(/,a,3(2x,g14.7))') chars,
     *                              times(itime), 
     *                              etime(itime)-btime(itime)

      end if 

      end 

      subroutine outsei
c-----------------------------------------------------------------------
c output details of modulus computation for endmembers and solutions
c if lopt(50).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character stag*9, ptag*8, made*12

      logical stx, notstx, lstx, lmake

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer eos
      common/ cst303 /eos(k10)

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      integer make
      common / cst335 /make(k10)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character prject*100,tfname*100
      common/ cst228 /prject,tfname
c-----------------------------------------------------------------------
      call mertxt (tfname,prject,'_seismic_data.txt',0)
      open (n8, file = tfname)

      stx = .false.
      notstx = .false.
      lmake = .false.

      write (n8,1233) valu(19),nopt(6),lopt(17),valu(15),
     *                nopt(1),valu(14),lopt(20),lopt(4),.false.

      write (n8,1030)

      write (n8,'(/,a)') 'Endmembers and stoichiometric compounds:'
      write (n8,1000)

      do i = istct, ipoint

         if (eos(i).eq.5.or.eos(i).eq.6) then 
            stx = .true.
            lstx = .true.
         else
            lstx = .false.
            notstx = .true.
         end if

         if (iemod(i).eq.3) then

            ptag = 'explicit'
            stag = 'missing'

         else if (iemod(i).eq.0) then

            ptag = 'implicit'
            stag = 'missing'

         else if (iemod(i).eq.1) then

            ptag = 'implicit'
            stag = 'explicit'

         else if (iemod(i).eq.2) then

            ptag = 'explicit'
            stag = 'explicit'

         end if

         if (.not.lopt(17)) ptag = 'implicit'

         if (iopt(16).eq.1) then

            if (stag.eq.'missing') stag = 'Poisson'

         else if (iopt(16).eq.2) then 

            stag = 'Poisson'

         end if

         if (lstx) then
            if (iemod(i).ge.1) stag = 'implicit'
            ptag = 'implicit'
         end if

         made = ' '

         if (make(i).ne.0) then 
            made = 'made entity*'
            lmake = .true.
         end if 

         if (ifp(i).ne.0) stag = 'fluid'

         write (n8,1050) names(i),ptag,stag,made

      end do

      if (lmake) write (n8,1040)

      if (isoct.gt.0) then

         if (stx.and.notstx) write (n8,1010) 

         write (n8,'(/,a)') 'Solutions:'
         write (n8,1000)

         do i = 1, isoct

            if (.not.pmod(i).and..not.smod(i)) then

               ptag = 'implicit'
               stag = 'missing'

            else if (pmod(i).and..not.smod(i)) then

               ptag = 'explicit'
               stag = 'missing'

            else if (.not.pmod(i).and.smod(i)) then

               ptag = 'implicit'
               stag = 'explicit'

            else if (pmod(i).and.smod(i)) then

               ptag = 'explicit'
               stag = 'explicit'

            end if

            if (lopt(17).and.pmod(i)) ptag = 'explicit'

            if (iopt(16).gt.0.and.(.not.smod(i).or.iopt(16).eq.2))
     *                                stag = 'Poisson'

            if (stx) then
               ptag = 'implicit'
               if (stag.eq.'explicit') stag = 'implicit*'
            end if

            if (fp(i)) stag = 'fluid'

            if (lname(i).eq.'liquid') stag = 'liquid'

            write (n8,1050) fname(i),ptag,stag

         end do

         if (stx) write (n8,1020)

      end if 

      close (n8)

1000  format (/,20x,'  Bulk Mod    Shear Mod ',/,
     *          20x,'  ---------   ---------')
1010  format (/,'**warning ver119** this computation mixes inconsistent'
     *         ,' thermodynamic data',/,'the following table may not be'
     *         ,' reliable.',/)
1020  format (/,'*computed as the Reuss average of the implicit endmemb'
     *       ,'er shear moduli.')
1030  format ('In the tables below: implicit moduli are calculated rigo'
     *       ,'rously from the EoS,',/,'explicit moduli are computed '
     *       ,'from empirical functions provided in the',/
     *       ,'thermodynamic data file.',/)
1040  format (/,'*explicit moduli of made endmembers are computed as a '
     *       ,'linear combination of ',/,'the real endmembers specified'
     *       ,' in the corresponding make definition.',/)
1050  format (6x,a10,6x,a8,4x,a9,4x,a)
1233  format (/,'Seismic wavespeed computational options:',//,
     *        4x,'bounds                 ',a3,8x,'[VRH] HS',/,
     *        4x,'vrh/hs_weighting       ',f3.1,8x,'[0.5] 0->1',/,
     *        4x,'explicit_bulk_modulus  ',l1,10x,'[T] F',/,
     *        4x,'poisson_ratio          ',a3,8x,'[on] all off; ',
     *        'Poisson ratio = ',f4.2,/,
     *        4x,'seismic_output         ',a3,8x,'[some] none all',/,
     *        4x,'poisson_test           ',l1,10x,'[F] T',/,
     *        4x,'Anderson-Gruneisen     ',l1,10x,'[F] T',/,
     *        4x,'Tisza_test             ',l1,10x,'[F] T',/)

      end
