      module mod_meemum_c_interface

        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env

        implicit none

        ! global variables
        include "perplex_parameters.h"

      contains

        ! ------------------------------------------------------------
        ! --------------------- SOLVER FUNCTIONS ---------------------
        ! ------------------------------------------------------------

        subroutine solver_init(filename) bind(c)
          ! part 1 of wrapper for meemm (meemum.f)
          ! read the input file (only needs to be called once)
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
        subroutine solver_minimize() bind(c)
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

        subroutine solver_set_pressure(pressure) bind(c)
          real(c_double), intent(in), value :: pressure

          ! source: ?
          double precision v, tr, pr, r, ps
          common / cst5  / v(l2), tr, pr, r, ps

          v(1) = pressure
        end subroutine

        subroutine solver_set_temperature(temperature) bind(c)
          real(c_double), intent(in), value :: temperature

          ! source: ?
          double precision v, tr, pr, r, ps
          common / cst5  / v(l2), tr, pr, r, ps

          v(2) = temperature
        end subroutine

        ! ------------------------------------------------------------
        ! ------------------ COMPOSITION PROPERTIES ------------------
        ! ------------------------------------------------------------

        function composition_props_get_n_components() bind(c) 
     >      result(res)
          integer(c_size_t) :: res

          ! source: olib.f
          integer icomp,istct,iphct,icp
          common / cst6  / icomp, istct, iphct, icp

          res = icomp
        end function

        function composition_props_get_name(component_idx) bind(c) result(res)
          integer(c_size_t), intent(in), value :: component_idx
          type(c_ptr) :: res

          ! source: olib.f
          character cname*5
          common / csta4  / cname(k5) 

          res = alloc_c_str(cname(component_idx+1))
        end function

        ! -----------------------------------------------------------
        ! --------------------- BULK PROPERTIES ---------------------
        ! -----------------------------------------------------------

        function bulk_props_get_composition(component_idx) bind(c) result(res)
          integer(c_size_t), intent(in), value :: component_idx
          real(c_double) :: res

          res = cblk(component_idx+1)
        end function

        subroutine bulk_props_set_composition(component_idx, amount) bind(c)
          integer(c_size_t), intent(in), value :: component_idx
          real(c_double), intent(in), value :: amount

          cblk(component_idx+1) = amount
        end subroutine

        ! ---------------------------------------------------------
        ! --------------- SOLUTION PHASE PROPERTIES ---------------
        ! ---------------------------------------------------------

        function soln_phase_props_get_n() bind(c) result(res)
          integer(c_size_t) :: res

          ! source: ???
          integer isoct
          common / cst79 / isoct

          res = isoct
        end function

        function soln_phase_props_get_abbr_name(soln_phase_idx) bind(c)
     >      result(res)
          integer(c_size_t), intent(in), value :: soln_phase_idx
          type(c_ptr) :: res

          ! source: ???
          character fname*10, aname*6, lname*22
          common / csta7 / fname(h9), aname(h9), lname(h9)

          res = alloc_c_str(aname(soln_phase_idx+1))
        end function

        function soln_phase_props_get_full_name(soln_phase_idx) bind(c)
     >      result(res)
          integer(c_size_t), intent(in), value :: soln_phase_idx
          type(c_ptr) :: res

          ! source: ???
          character fname*10, aname*6, lname*22
          common / csta7 / fname(h9), aname(h9), lname(h9)

          res = alloc_c_str(fname(soln_phase_idx+1))
        end function

        ! ----------------------------------------------------------
        ! -------------------- PHASE PROPERTIES --------------------
        ! ----------------------------------------------------------

        function res_phase_props_get_n() bind(c) result(res)
          integer(c_size_t) :: res

          ! source: ???
          integer kkp, np, ncpd, ntot
          double precision cp3, amt
          common / cxt15 / cp3(k0,k19), amt(k19), kkp(k19), 
     >    np, ncpd, ntot

          res = ntot
        end function

        function res_phase_props_get_name(res_phase_idx) bind(c) 
     >      result(res)
          integer(c_size_t), intent(in), value :: res_phase_idx
          type(c_ptr) :: res

          ! source: ???
          character pname*14
          common / cxt21a / pname(k5)

          res = alloc_c_str(pname(res_phase_idx+1))
        end function

        function res_phase_props_get_weight_frac(res_phase_idx) bind(c)
     >      result(res)
          integer(c_size_t), intent(in), value :: res_phase_idx
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          res = props(17, res_phase_idx+1) * props(16, res_phase_idx+1) 
     >    / psys(17)
        end function

        function res_phase_props_get_vol_frac(res_phase_idx) bind(c)
     >      result(res)
          integer(c_size_t), intent(in), value :: res_phase_idx
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          res = props(1, res_phase_idx+1) * props(16, res_phase_idx+1) 
     >    / psys(1)
        end function

        function res_phase_props_get_mol_frac(res_phase_idx) bind(c)
     >      result(res)
          integer(c_size_t), intent(in), value :: res_phase_idx
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8), pgeo1(i8)

          res = props(16,res_phase_idx+1) / psys(16)
        end function

        function res_phase_props_get_mol(res_phase_idx) bind(c)
     >      result(res)
          integer(c_size_t), intent(in), value :: res_phase_idx
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8), pgeo1(i8)

          res = props(16, res_phase_idx+1)
        end function

        function res_phase_props_get_composition(res_phase_idx, 
     >      component_idx) bind(c) result(res)
          integer(c_size_t), intent(in), value :: res_phase_idx
          integer(c_size_t), intent(in), value :: component_idx
          real(c_double) :: res

          ! source: olib.f
          double precision pcomp
          common / cst324 / pcomp(k0,k5)

          res = pcomp(component_idx+1,res_phase_idx+1)
        end function

        ! -----------------------------------------------------------
        ! -------------------- SYSTEM PROPERTIES --------------------
        ! -----------------------------------------------------------
        
        function sys_props_get_density() bind(c) result(res)
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8), pgeo1(i8)

          res = psys(10)
        end function

        function sys_props_get_expansivity() bind(c) result(res)
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8), pgeo1(i8)

          res = psys(13)
        end function

        function sys_props_get_mol_entropy() bind(c) result(res)
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8), pgeo1(i8)

          res = psys(15)
        end function

        function sys_props_get_mol_heat_capacity() bind(c) result(res)
          real(c_double) :: res

          ! source: olib.f
          double precision props, psys, psys1, pgeo, pgeo1
          common / cxt22 / props(i8,k5), psys(i8), psys1(i8),
     >    pgeo(i8), pgeo1(i8)

          res = psys(12)
        end function

        ! -----------------------------------------------------------
        ! -------------------- PRIVATE FUNCTIONS --------------------
        ! -----------------------------------------------------------

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

