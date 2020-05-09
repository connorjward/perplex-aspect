      module meemum_wrapper_mod

        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env

        implicit none

        include "perplex_parameters.h"

      contains

        subroutine init(filename) bind(c, name="meemum_init")
          character(c_char), dimension(*), intent(in) :: filename

          integer :: i

          ! rlib.f : subroutine input1
          character*100 prject,tfname
          common/ cst228 /prject,tfname

          ! part 1 of wrapper for meemm (meemum.f)
          integer iam
          common/ cst4 /iam

         ! ----- PREP WORK -----

         ! prevent input1 crashing
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

        subroutine minimize(temperature, pressure)
     >      bind(c, name="meemum_minimize")
          ! arguments
          real(c_double) :: temperature, pressure

          ! part 2 of wrapper for meemm (meemum.f)
          integer i, ier

          logical bulk, bad

          character amount*6, yes*1

          double precision num

          integer iwt
          common/ cst209 /iwt

          integer npt,jdv
          logical fulrnk
          double precision cptot,ctotal
          common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

          double precision atwt
          common/ cst45 /atwt(k0) 

          double precision v,tr,pr,r,ps
          common/ cst5  /v(l2),tr,pr,r,ps

          integer ipot,jv,iv
          common / cst24 /ipot,jv(l2),iv(l2)

          character*8 vname,xname
          common/ csta2  /xname(k5),vname(l2)

          character*5 cname
          common/ csta4 /cname(k5)

          double precision a,b,c
          common/ cst313 /a(k5,k1),b(k5),c(k1)

          integer icomp,istct,iphct,icp
          common/ cst6  /icomp,istct,iphct,icp

          integer io3,io4,io9
          common / cst41 /io3,io4,io9

          logical gflu,aflu,fluid,shear,lflu,volume,rxn
          common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

          double precision goodc, badc
          common/ cst20 /goodc(3),badc(3)

          integer iam
          common/ cst4 /iam

c----------------------------------------------------------------------- 

          ! set potential (P, T) values
          v(2) = temperature
          v(1) = pressure

          ! convert to moles if needed
          if (iwt.eq.1) then 
            do i = 1, jbulk
            cblk(i) = cblk(i)/atwt(i)
          end do 
        end if

c                                 meemum does the minimization and outputs
c                                 the results to the print file.
        call meemum_wrapper

         if (.not.bad) then
c                                 print summary to LUN 6
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)

         end if 

         if (goodc(1)+badc(1).gt.0d0) then

            num = badc(1)/(badc(1)+goodc(1))*1d2
            if (num.gt.1d-1) call warn (53,num,i,'MEEMUM')

         end if 
       end subroutine


        function density() bind(c, name="meemum_density") result(rho)
          real(c_double) :: rho

          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)
          
          rho = psys(10)
        end function

        function n_phases() bind(c, name="meemum_n_phases") result(n)
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

        function expansivity() bind(c, name="meemum_expansivity") 
     >      result(alpha)
          real(c_double) :: alpha

          ! meemum_trimmed_subprogram
          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)

          alpha = psys(13)
        end function

        function heat_capacity() result(Cp)
     >      bind(c, name="meemum_heat_capacity") 
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

        function entropy() bind(c, name="meemum_entropy") result(S)
          real(c_double) :: S

          double precision props,psys,psys1,pgeo,pgeo1
          common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),
     >    pgeo1(i8)

          S = psys(15)/psys(1)*1d5/psys(10)
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

        ! PRIVATE FUNCTIONS
       
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

       subroutine meemum_wrapper()
         ! wrapper for meemum (resub.f)

         ! ----- VARIABLES -----

      integer i, idead

      logical nodata, bad

      integer itri(4),jtri(4),ijpt

      double precision wt(3), cum

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5)

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
         ! ----- PREP WORK -----

         ! ----- WRAPPER -----
      ! normalize vector
      ctotal = 0d0
      do i = 1, icp
          ctotal = ctotal + cblk(i)
      end do 

      do i = 1, icp
         b(i) = cblk(i)/ctotal
      end do
      call incdp0

      call lpopt0 (idead)

      if (idead.eq.0) then
c                                 compute derivative properties
         call getloc (itri,jtri,ijpt,wt,nodata)

         bad = .false.

      else 

         bad = .true.

      end if
       end subroutine

       
       end module
