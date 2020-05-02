      module meemum_wrapper_mod

        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env

        implicit none

        include "perplex_parameters.h"

      contains

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

        subroutine init() bind(c)
          ! part 1 of wrapper for meemm (meemum.f)
          integer iam
          common/ cst4 /iam

          iam = 2
          call vrsion (6)
          call iniprp_wrapper
        end subroutine

        subroutine minimize(temperature, pressure, component_amounts)
     *      bind(c)
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

          ! arguments
          double precision :: pressure, temperature
          double precision, dimension(k5) :: component_amounts

c----------------------------------------------------------------------- 

          !
          ! set potential (P, T) values
          pressure = 100
          temperature = 120

          v(1) = temperature
          v(2) = pressure


          ! load composition
          do i = 1, jbulk
          cblk(i) = component_amounts(i)
          end do

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
       
       subroutine iniprp_wrapper()
         ! wrapper for iniprp (resub.f)

         ! ----- VARIABLES -----
         ! resub.f : subroutine iniprp
         logical first, err 

         ! rlib.f : subroutine input1
         character*100 prject,tfname
         common/ cst228 /prject,tfname

         ! ----- PREP WORK -----

         ! prevent input1 crashing
         prject = "my_project"

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
