      module props

        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env

        implicit none

        include "perplex_parameters.h"

      contains

        function abbr_soln_name(soln_id) bind(c) result(ptr)
          integer(c_size_t), intent(in) ::soln_id
          type(c_ptr) :: ptr

          character(c_char), dimension(:), pointer :: soln_name

          ! source: rlib.f
          character fname*10, aname*6, lname*22
          common / csta7 / fname(h9), aname(h9), lname(h9)

          integer :: i, strlen

          strlen = len_trim(aname(soln_id))
      
          allocate(soln_name(1:strlen+1))

          soln_name = transfer(aname(soln_id), " ", strlen)
          soln_name(strlen+1) = c_null_char

          ptr = c_loc(soln_name)
        end function

      ! subroutine abbr_soln_name(soln_id, soln_name) bind(c)
      !   integer(c_size_t), intent(in) ::soln_id
      !   character(c_char), dimension(*), intent(out) :: soln_name

      !   ! source: rlib.f
      !   character fname*10, aname*6, lname*22
      !   common/ csta7 /fname(h9),aname(h9),lname(h9)

      !   integer :: i, strlen

      !   strlen = len(trim(aname(soln_id)))

      !   do i = 1, strlen
      !   soln_name(i) = aname(soln_id)(i:i)
      !   end do
      !   soln_name(strlen+1) = c_null_char
      ! end subroutine

        function n_soln_models() bind(c) result(n)
          integer(c_size_t) :: n

          ! source: perplex_parameters.h
          integer isoct
          common / cst79 / isoct

          n = isoct
        end function

       subroutine load_full_soln_name(soln_id, soln_name)
     >     bind(c)
          integer(c_size_t), intent(in) ::soln_id
          character(c_char), dimension(*), intent(out) :: soln_name

         ! source: rlib.f
         character fname*10, aname*6, lname*22
         common/ csta7 /fname(h9),aname(h9),lname(h9)
         integer :: i, strlen

         strlen = len(trim(lname(soln_id)))

         do i = 1, strlen
           soln_name(i) = lname(soln_id)(i:i)
         end do
         soln_name(strlen+1) = c_null_char

        end subroutine

        subroutine load_phase_name(phase_id, phase_name)
     >      bind(c)
          integer(c_size_t), intent(in) :: phase_id
          character(c_char), dimension(*), intent(out) :: phase_name

          ! meemum_trimmed_subprogram.f
          character pname*14
          common/ cxt21a /pname(k5)

          integer :: i, strlen

          strlen = len(trim(pname(phase_id)))

          do i = 1, strlen
            phase_name(i:i) = pname(phase_id)(i:i)
          end do

          ! null-terminate C string
          phase_name(strlen+1) = c_null_char
        end subroutine

        function phase_weight_frac(phase_id) result(frac)
     >      bind(c)
          integer(c_size_t), intent(in) :: phase_id
          real(c_double) :: frac

          ! olib.f
          double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

          frac = props(17, phase_id) * props(16, phase_id) / psys(17)
        end function

        function phase_vol_frac(phase_id) result(frac)
     >      bind(c)
          integer(c_size_t), intent(in) :: phase_id
          real(c_double) :: frac

          ! olib.f
          double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

          frac = props(1, phase_id) * props(16, phase_id ) / psys(1)
        end function

        function phase_mol_frac(phase_id) result(frac)
     >      bind(c)
          integer(c_size_t), intent(in) :: phase_id
          real(c_double) :: frac

          ! olib.f
          double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

          frac = props(16, phase_id) / psys(16)
        end function

        function phase_mol(phase_id) result(frac)
     >      bind(c)
          integer(c_size_t), intent(in) :: phase_id
          real(c_double) :: frac

          ! olib.f
          double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

          frac = props(16, phase_id)
        end function

        function density() bind(c) result(rho)
          real(c_double) :: rho

          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)
          
          rho = psys(10)
        end function

        function n_phases() bind(c) result(n)
          integer(c_int) :: n

          ! meemum_trimmed_subprogram.f
          integer kkp,np,ncpd,ntot
          double precision cp3,amt
          common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

          n = ntot
        end function

        function get_phase_amount(phase_id) bind(c) result(amount)
          integer(c_int), intent(in) :: phase_id
          real(c_double) :: amount

          ! meemum_trimmed_subprogram
          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)

          amount = props(16, phase_id) 
        end function

        function expansivity() bind(c) 
     >      result(alpha)
          real(c_double) :: alpha

          ! meemum_trimmed_subprogram
          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)

          alpha = psys(13)
        end function

        function heat_capacity() result(Cp)
     >      bind(c)
          real(c_double) :: Cp

          ! meemum_trimmed_subprogram
          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),
     >    pgeo(i8),pgeo1(i8)

          Cp = psys(12) / psys(1) * 1d5 / psys(10)
        end function

        function get_melt_frac() bind(c) result(melt_frac)
          real(c_double) :: melt_frac
          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)

          ! melt_frac = props(17,imelt)*props(16,imelt) / psys(17)
          melt_frac = 0.1
        end function

        function entropy() bind(c) result(S)
          real(c_double) :: S

          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)

          S = psys(15)
        end function

        function has_melt() bind(c)
          logical(c_bool) :: has_melt

          logical gflu,aflu,fluid,shear,lflu,volume,rxn
          common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

          has_melt = aflu
        end function

        function is_melt(phase_id) bind(c)
          integer(c_int), intent(in) :: phase_id
          logical(c_bool) :: is_melt

          logical gflu,aflu,fluid,shear,lflu,volume,rxn
          common / cxt20 / gflu, aflu, fluid(k5), shear, lflu, volume,
     >    rxn

          is_melt = fluid(phase_id)
        end function

        function get_n_components() bind(c) result(n)
          integer(c_int) :: n

          n = k5
        end function

        subroutine get_component_name(comp_id, comp_name) bind(c)
          integer(c_int), intent(in) :: comp_id
          character(c_char), dimension(*), intent(out) :: comp_name

          ! component names (from where?)
          character*5 cname
          common / csta4 / cname(k5)

          integer :: i

          do i = 1, 5
          comp_name(i:i) = cname(comp_id)(i:i)
          end do
        end subroutine

        subroutine get_component_amount(comp_id, comp_amount) bind(c)
          integer(c_int), intent(in) :: comp_id
          real(c_double), intent(out) :: comp_amount

          integer :: i

          ! where is this from? file reference
          double precision a,b,c
          common / cst313 / a(k5,k1), b(k5), c(k1)

          comp_amount = b(comp_id)
        end subroutine
       
      end module
