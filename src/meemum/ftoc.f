      module mod_meemum_ftoc

        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env

        implicit none

        ! global variables
        include "perplex_parameters.h"

      contains

        ! part 1 of wrapper for meemm (meemum.f)
        ! read the input file (only needs to be called once)
        subroutine init(filename) bind(c)
          character(c_char), dimension(*), intent(in) :: filename

          integer :: i

          character*100 prject,tfname
          common/ cst228 /prject,tfname

          integer iam
          common/ cst4 /iam

          ! ----- PREP WORK -----

          ! prevent input1 crashing by setting project name
          do i = 1, 100
          if (filename(i) == c_null_char) then
            exit
          end if
          prject(i:i) = filename(i)
          end do

          iam = 2
          call vrsion (6)
          call iniprp_wrapper
        end subroutine

        !> part 2 of wrapper for meemm (meemum.f)
        ! called at each iteration
        ! pretty much imitates meemm
        subroutine minimize() bind(c)
          integer i

          logical bad

          double precision num

          integer iwt
          common/ cst209 /iwt

          double precision atwt
          common/ cst45 /atwt(k0) 

          integer io3,io4,io9
          common / cst41 /io3,io4,io9

          ! -----------------------------------------

          ! convert to moles if needed
          if (iwt.eq.1) then 
            do i = 1, jbulk
            cblk(i) = cblk(i)/atwt(i)
            end do 
          end if

          call meemum (bad)

          if (.not.bad) then
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)
         end if 
       end subroutine

       subroutine set_pressure(pressure) bind(c)
         real(c_double), intent(in), value :: pressure
 
         ! source: ?
         double precision v, tr, pr, r, ps
         common / cst5  / v(l2), tr, pr, r, ps
 
         v(1) = pressure
       end subroutine
 
       subroutine set_temperature(temperature) bind(c)
         real(c_double), intent(in), value :: temperature
 
         ! source: ?
         double precision v, tr, pr, r, ps
         common / cst5  / v(l2), tr, pr, r, ps
 
         v(2) = temperature
       end subroutine

       subroutine set_composition_component(id, amount) bind(c)
         integer(c_size_t), intent(in), value :: id
         real(c_double), intent(in), value :: amount

         cblk(id) = amount
       end subroutine

        !> @return n_soln_models Number of solution models
        function get_n_soln_models() bind(c) result(n_soln_models)
          integer(c_size_t) :: n_soln_models

          integer isoct
          common / cst79 / isoct

          n_soln_models = isoct
        end function

        !> @param  soln_id        Solution model index
        !! @return abbr_soln_name Abbreviated solution model name
        function get_abbr_soln_name(soln_id) bind(c)
     >      result(abbr_soln_name)
          integer(c_size_t), intent(in), value :: soln_id
          type(c_ptr) :: abbr_soln_name

          character fname*10, aname*6, lname*22
          common / csta7 / fname(h9), aname(h9), lname(h9)

          abbr_soln_name = alloc_c_str(aname(soln_id))
        end function

        !> @param  soln_id        Solution model index
        !! @return full_soln_name Full solution model name
        function get_full_soln_name(soln_id) bind(c)
     >      result(full_soln_name)
          integer(c_size_t), intent(in), value :: soln_id
          type(c_ptr) :: full_soln_name

          character fname*10, aname*6, lname*22
          common / csta7 / fname(h9), aname(h9), lname(h9)

          full_soln_name = alloc_c_str(fname(soln_id))
        end function

        !> @return n_phases Number of phases
        function get_n_phases() bind(c) result(n_phases)
          integer(c_size_t) :: n_phases

          integer kkp, np, ncpd, ntot
          double precision cp3, amt
          common / cxt15 / cp3(k0,k19), amt(k19), kkp(k19), 
     >    np, ncpd, ntot

          n_phases = ntot
        end function

        !> @param  phase_id   Phase index
        !! @return phase_name Phase name
        function get_phase_name(phase_id) bind(c) result(phase_name)
          integer(c_size_t), intent(in), value :: phase_id
          type(c_ptr) :: phase_name

          character pname*14
          common / cxt21a / pname(k5)

          phase_name = alloc_c_str(pname(phase_id))
        end function

        !> @param  phase_id          Phase index
        !! @return phase_weight_frac Fractional phase weight
        function get_phase_weight_frac(phase_id) bind(c)
     >      result(phase_weight_frac)
          integer(c_size_t), intent(in), value :: phase_id
          real(c_double) :: phase_weight_frac

        double precision props, psys, psys1, pgeo, pgeo1
        common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >  pgeo(i8),pgeo1(i8)
          ! source: olib.f
          phase_weight_frac = props(17, phase_id) * props(16, phase_id) 
     >    / psys(17)
        end function

        !> @param  phase_id       Phase index
        !! @return phase_vol_frac Fractional volume of a phase
        function get_phase_vol_frac(phase_id) bind(c)
     >      result(phase_vol_frac)
          integer(c_size_t), intent(in), value :: phase_id
          real(c_double) :: phase_vol_frac

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          ! source: olib.f
          phase_vol_frac = props(1, phase_id) * props(16, phase_id) 
     >    / psys(1)
        end function

        !> @param  phase_id       Phase index
        !! @return phase_mol_frac Fractional number of moles of a phase
        function get_phase_mol_frac(phase_id) bind(c)
     >      result(phase_mol_frac)
          integer(c_size_t), intent(in), value :: phase_id
          real(c_double) :: phase_mol_frac

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          ! source: olib.f
          phase_mol_frac = props(16, phase_id) / psys(16)
        end function

        !> @param  phase_id  Phase index
        !! @return phase_mol Number of moles of a phase
        function get_phase_mol(phase_id) bind(c)
     >      result(phase_mol)
          integer(c_size_t), intent(in), value :: phase_id
          real(c_double) :: phase_mol

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          ! source: olib.f
          phase_mol = props(16, phase_id)
        end function

        !> @return sys_density The system density (kg/m3)
        function get_sys_density() bind(c) result(sys_density)
          real(c_double) :: sys_density

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)
          
          ! source: olib.f
          sys_density = psys(10)
        end function

        !> @return sys_expansivity The system expansivity
        function get_sys_expansivity() bind(c) 
     >      result(sys_expansivity)
          real(c_double) :: sys_expansivity

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          ! source: olib.f
          sys_expansivity = psys(13)
        end function

        !> @return sys_mol_entropy The system molar entropy
        function get_sys_mol_entropy() bind(c)
     >      result(sys_mol_entropy)
          real(c_double) :: sys_mol_entropy

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          ! source: olib.f
          sys_mol_entropy = psys(15)
        end function

        !> @return sys_mol_heat_capacity The system molar heat capacity
        function get_sys_mol_heat_capacity() bind(c)
     >      result(sys_mol_heat_capacity)
          real(c_double) :: sys_mol_heat_capacity

          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          ! source: olib.f
          sys_mol_heat_capacity = psys(12)
        end function

        subroutine iniprp_wrapper()
          ! wrapper for iniprp (resub.f)

          ! ----- VARIABLES -----
          ! resub.f : subroutine iniprp
          logical first, err 


          ! ----- WRAPPER -----
          first = .true.

          ! *** prevent prompt for project name ***
          outprt = .true.

          call input1 (first,err)
          call input2 (first)
          call setau1 
          call input9 (first)
          if (lopt(50)) call outsei
          call initlp
        end subroutine

        !> Convert a Fortran string to a C string.
        !! @param f_str The Fortran string to be converted.
        !! @return ptr The location of the allocated C string in memory.
        function alloc_c_str(f_str) result(ptr)
          character(len=*), intent(in) :: f_str
          type(c_ptr) :: ptr

          character(c_char), dimension(:), pointer :: c_str

          ! allocate memory for C string
          allocate(c_str(len_trim(f_str)+1))

          ! copy across Fortran string
          c_str = transfer(trim(f_str) // c_null_char, c_str)

          ! return pointer to C string
          ptr = c_loc(c_str)
        end function
      end module
       
