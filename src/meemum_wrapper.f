! requires molar masses

       module meemum_wrapper

       use iso_c_binding, only: c_double
       implicit none

       contains

       subroutine c_init() bind(c)
         call init()
       end subroutine

       subroutine c_minimize() bind(c)
         call minimize()
       end subroutine
       
       subroutine init()
         include "perplex_parameters.h"

         ! meemum.f : program meemum
         integer iam
         common/ cst4 /iam

         ! resub.f : subroutine iniptp
         logical first, err 

         ! rlib.f : subroutine input1
         character*100 prject,tfname
         common/ cst228 /prject,tfname

         ! initial work

         ! prevent input1 crashing
         prject = "my_project"

         ! now start the wrapper
         ! meemum.f : program meemm

         iam = 2
         call vrsion (6)

         ! resub.f : subroutine iniptp

         first = .true.

         ! prevent prompt for project name
         outprt = .true.

         call input1 (first,err)
         call input2 (first)
         call setau1 
         call input9 (first)
         if (lopt(50)) call outsei
         call initlp
       end subroutine
       
       subroutine minimize()
         include "../perplex/perplex_parameters.h"

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
         ! avoid bulk usage prompt
         bulk = .true.

c                                 iwt is set by input, it is only used below to determine
c                                 whether to convert weight compositions to molar. the 
c                                 computations are done solely in molar units. 
      amount = 'molar '

      if (iwt.eq.1) amount = 'weight'

      if (lopt(28)) open (666,file='times.txt')
c                                 computational loop
      do 
c                                 read potential variable values    
c                                 v(1) is P(bar), v(2) is T(K) the pointer jv used 
c                                 for general problems but can be eliminated for calculations 
c                                 simply as a f(P,T)       
         write (*,1070) (vname(jv(i)), i = 1, ipot)
         read (*,*,iostat=ier) (v(jv(i)), i = 1, ipot)
         if (ier.ne.0) cycle
         if (v(jv(1)).eq.0d0) exit 
          
         if (bulk) then 
c                                 load the composition into b, the component names are  
c                                 in cname, if iwt = 1 the composition is in mass fractions
c                                 otherwise in molar units. 
            do 
               write (*,1060) amount
               write (*,'(12(a,1x))') (cname(i),i=1,jbulk)
               read (*,*,iostat=ier) (cblk(i),i=2,jbulk)
               if (ier.eq.0) exit
            end do  
         
            if (iwt.eq.1) then 
c                                 convert mass to molar 
               do i = 1, jbulk
                  cblk(i) = cblk(i)/atwt(i)
               end do 

            end if

         end if 
c                                 meemum does the minimization and outputs
c                                 the results to the print file.
         call meemum (bad)

         if (.not.bad) then
c                                 print summary to LUN 6
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)

         end if 

         if (goodc(1)+badc(1).gt.0d0) then

            num = badc(1)/(badc(1)+goodc(1))*1d2
            if (num.gt.1d-1) call warn (53,num,i,'MEEMUM')

         end if 

      end do

1000  format (/,'Interactively enter bulk compositions (y/n)?',/,
     *          'If you answer no, MEEMUM uses the bulk composition',
     *         ' specified in the input file.',/)
1060  format (/,'Enter ',a,' amounts of the components:')
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))
       end subroutine
       
       subroutine print_stuff
         call calpr0(6)
       end subroutine
       
       end module
