PROGRAM SquareBlock_CPFE
!-----------------------------------------------------------------------
! Anton Shterenlikht (University of Bristol)
! Luis Cebamanos (EPCC, University of Edinburgh)
!
! Program Rect_block - linking ParaFEM with CGPACK, specifically
! modifying p121 from 5th edition to link with the cgca module.
!
! The model is a 3D cylinder under tension.
!
! Reproducible RND seed routine (cgca_ins) is used here.
!
! 12.1 is a three dimensional analysis of an elastic solid
! using 20-node brick elements, preconditioned conjugate gradient
! solver; diagonal preconditioner diag_precon; parallel version
! loaded_nodes only
!
! This version must conform to F2008 standard.
! This mainly means none of the routines using Cray extensions
! can be used.
!
! The CA CS is aligned with the FE CS, but not the same origin.
! Anyway, all cells are within the FE model, so mapping is easy.
!-----------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation

!*** CGPACK part *****************************************************72
! The CGPACK module must be used
use casup
!*** end CGPACK part *************************************************72

USE precision
USE global_variables
USE mp_interface
USE input
USE output
USE loading
USE timing
USE maths
USE gather_scatter
USE steering
USE new_library
USE cgca_m3pfem
USE read_umat_params__genmod
USE init_umat_params__genmod
USE umat__genmod

IMPLICIT NONE
!include 'module_interface.h'

! neq, ntot are now global variables - must not be declared

INTEGER, PARAMETER :: nodof = 3, ndim = 3, nst = 6
INTEGER :: loaded_nodes, iel, i, j, k, iters, limit, nn, nr, nip, nod, &
   nels, ndof, npes_pp, node_end, node_start, nodes_pp, meshgen,       &
   partitioner, nlen, nelements, ios

REAL( iwp ), PARAMETER :: zero=0.0_iwp
REAL( iwp ) :: e, v, det, tol, up, alpha, beta, q

LOGICAL :: converged=.false.

CHARACTER( LEN=50 ) :: argv
CHARACTER( LEN=15 ) :: element
CHARACTER( LEN=6 ) :: ch

!---------------------------- dynamic arrays -------------------------72
REAL( iwp ), ALLOCATABLE :: points(:,:), dee(:,:), weights(:),         &
   val(:,:),                                                           &
   disp_pp(:), g_coord_pp(:,:,:), jac(:,:), der(:,:), deriv(:,:),      &
   bee(:,:), storkm_pp(:,:,:), eps(:), sigma(:), diag_precon_pp(:),    &
   p_pp(:), r_pp(:), x_pp(:), xnew_pp(:), u_pp(:), pmul_pp(:,:),       &
   utemp_pp(:,:), d_pp(:), timest(:), diag_precon_tmp(:,:),            &
   eld_pp(:,:), temp(:), tot_r_pp(:)
INTEGER, ALLOCATABLE :: rest(:,:), g_num_pp(:,:), g_g_pp(:,:), node(:)
REAL, ALLOCATABLE    :: centr_ca(:,:)

! --------------------------Additional parameters----------------------------
REAL(iwp)              :: stot_upfc, smean_upfc,a_upfc,b_upfc,         &
                          TLoad_upfc, pres_upfc, dtot_upfc,            &
                          dmean_upfc
REAL(iwp)              :: EigenVec(3)=(/1,1,1/),EigenVal
INTEGER, ALLOCATABLE   :: elem_id_upfc(:),etype(:),nf(:,:),node_id_upfc(:)
INTEGER                :: fixed_freedoms,inod,iseed,incr1,incr2,incr3
INTEGER                :: steps,nel,rel,el,rimg,nintp,maxindx,maxelemN
!INTEGER                :: node_id_upfc(4)=(/1,25,225,249/)
!INTEGER                :: node_id_upfc(1)=(/25/)
INTEGER, PARAMETER     :: ncomp=6
REAL                   :: stressTensor(3,3),maxeval
REAL( kind=rdef ), ALLOCATABLE :: str_tmp( : , : ),stot(:,:)
!INTEGER                :: seed(1) = (/100/),n,iseedx,iseedy,iseedz,errstat
INTEGER                :: n,iseedx,iseedy,iseedz,errstat
REAL                   :: num
INTEGER, ALLOCATABLE   :: icgca_cx(:),icgca_cy(:),icgca_cz(:)
CHARACTER(LEN=15)      :: elemtype
REAL(iwp), ALLOCATABLE :: g_coord(:,:)
INTEGER( kind=iarr )   :: cgca_test(10,20,5)[*]
!-------------------------------------------------------------------------
!------------------------CPFE UMAT declaration----------------------------
REAL*8, ALLOCATABLE    :: props(:) 
character(len=50)      :: fname
integer                :: nprops
real*8                 :: sse,spd,scd,rpl,drpldt,celent,dtime
real*8                 :: temp_umat,dtemp,predef,dpred,pnewdt
real*8                 :: kspt,kstep,kinc
real*8, allocatable    :: drot(:,:),coords(:),dfgrd0(:,:)
real*8, allocatable    :: stress(:),statev(:),dfgrd1(:,:)
real*8, allocatable    :: ddsddt(:),ddsdde(:,:),drplde(:)
real*8, allocatable    :: stran(:),dstran(:),time(:)
character(len=8)       :: cmname
integer                :: ndi,nshr,nstatv,noel,layer
!-------------------------------------------------------------------------
interface
  subroutine read_elem_upfc(elem_id_upfc)
    integer, allocatable, intent(inout) :: elem_id_upfc(:)
  end subroutine read_elem_upfc
end interface
interface
  subroutine read_node_upfc(node_id_upfc)
    integer, allocatable, intent(inout) :: node_id_upfc(:)
  end subroutine read_node_upfc
end interface
interface
  subroutine FINDEIGEN(Matrix, n, x, EigenVal, steps)
    INTEGER, INTENT(IN) :: n, steps  !n = order of matrix, steps = number of iterations
    REAL, INTENT(IN), DIMENSION(n,n) :: Matrix(n,n)  !Input Matrix
    REAL, INTENT(INOUT), DIMENSION(n) :: x !Eigenvector
    REAL, INTENT(INOUT) :: EigenVal !Eigenvalue
  end subroutine FINDEIGEN
end interface
interface
  subroutine MULMATRIX(a, b, n)
    INTEGER, INTENT(IN) :: n !matrix size
    REAL, INTENT(IN), DIMENSION(n,n) :: a  !Matrix of order > 1
    REAL, INTENT(INOUT), DIMENSION(n) :: b !1x1 matrix 
  end subroutine
end interface

interface
  subroutine FINDLARGEST(x, n, l)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN), DIMENSION(n) :: x
    REAL, INTENT(INOUT) :: l !Largest value
  end subroutine
end interface

!*** CGPACK part *****************************************************72
! CGPACK parameters
!integer, parameter :: cgca_linum=100 ! number of loading iterations
integer, parameter :: cgca_linum=20 ! number of loading iterations
logical( kind=ldef ), parameter :: cgca_yesdebug = .true.,             &
 cgca_nodebug = .false.
real( kind=rdef ), parameter :: cgca_zero = 0.0_rdef,                  &
 cgca_one = 1.0_rdef,                                                  &
 ! cleavage stress on 100, 110, 111 planes for BCC,
 ! see the manual for derivation, GPa.
!! cgca_scrit(3) = (/ 110000.0_rdef, 118200.0_rdef, 125500.0_rdef /) !Load=3000MPa
cgca_scrit(3) = (/ 1.05e1_rdef, 1.25e1_rdef, 4.90e1_rdef /)

! CGPACK variables
integer( kind=idef ) ::                                                &
   cgca_ir(3),            & ! coarray codimensions
   cgca_img,              &
   cgca_nimgs,            &
   cgca_ng,               & ! number of grains in the whole model
   cgca_clvg_iter,        & ! number of cleavage iterations
   cgca_liter               ! load iteration number
integer( kind=iarr ) :: cgca_c(3) ! coarray dimensions
integer( kind=iarr ), allocatable :: cgca_space(:,:,:,:) [:,:,:]

real( kind=rdef ) ::                                                   &
  cgca_qual,             & ! quality
  cgca_bsz(3),           & ! the given and the updated "box" size
  cgca_origin(3),        & ! origin of the "box" cs, in FE cloads
  cgca_rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
  cgca_dm,               & ! mean grain size, linear dim, phys units
  cgca_res,              & ! resolutions, cells per grain
  cgca_bcol(3),          & ! lower phys. coords of the coarray on image
  cgca_bcou(3),          & ! upper phys. coords of the coarray on image
  cgca_stress(3,3),      & ! stress tensor
  cgca_length,           & ! fracture length scale
  cgca_time_inc,         & ! time increment
  cgca_lres,             & ! linear resolution, cells per unit of length
  cgca_charlen,          & ! characteristic element length
  cgca_fracvol             ! volume (number) of fractured cells per img
real( kind=rdef ), allocatable :: cgca_grt(:,:,:)[:,:,:]
logical( kind=ldef )           :: cgca_solid
character( len=6 )             :: cgca_citer
real(kind=rdef)                :: cgca_cxyz(3)
!*** end CGPACK part *************************************************72
! Read element numbers of the upper surface
call read_elem_upfc(elem_id_upfc)

! Read node numbers of the upper surface
call read_node_upfc(node_id_upfc)
write(*,*)'node_id_upfc',node_id_upfc

!------------------------ input and initialisation ---------------------
ALLOCATE( timest( 20 ) )
timest = zero
timest( 1 ) = elap_time()

!*    Get the rank of the processes and the total number of processes
!* intent( out ):
!*          numpe - integer, process number (rank)
!*           npes - integer, total number of processes (size)
CALL find_pe_procs( numpe, npes )
!*           argv - character(*), data file base name
!*           nlen - integer, number of characters in data file base name
CALL getname( argv, nlen )

!*    Master processor reads the general data for the problem
!     and broadcasts it to the slave processors.
!  in:
!        argv - character, file name to read from
!       numpe - MPI rank
!           e - Young's modulus
!     element - character, element type
!       limit - max number of iterations
!loaded_nodes - number of nodes with applied forces
!fixed_freedoms - number of nodes with fixed degrees of freedom
!     meshgen - mesh numbering scheme
!        nels - total number of elements
!         nip - number of Gauss points
!          nn - total number of nodes in the mesh
!         nod - number of nodes per element
!         tol - tolerance
!           v - Poisson's ratio
!https://code.google.com/p/parafem/source/browse/trunk/parafem/src/modules/mpi/input.f90
CALL read_p121f( argv, numpe, e, element, limit, loaded_nodes, fixed_freedoms, meshgen, &
                nels, nip, nn, nod, nr, partitioner, tol, v )

! Calculates the number of elements, nels_pp, assigned to each
! processor.
! It is effectively a very naive method of mesh partitioning.
! The subroutine also computes, iel_start, the first element number
! on each processor. iel_start and nels_pp are the external variables
! modified by this subroutine. Note that they are global variables
! that are not passed through the list of arguments.
!
! nels_pp is indeed a global var, defined in
! http://sourceforge.net/p/parafem/code/HEAD/tree/trunk/parafem/src/modules/shared/global_variables.f90
CALL calc_nels_pp( argv, nels, npes, numpe, partitioner, nels_pp )

write(*,*)'nels_pp',nels_pp

!   nod - number of nodes per element
! nodof = 3 (see the beginning of this file), probably
!           degrees of freedom per node.
ndof = nod * nodof

! ntot - a global variable,
! the total number of degrees of freedom per element
ntot = ndof

! g_num_pp(nod,nels_pp) - integer, elements connectivity
! g_coord_pp - global coord?
ALLOCATE( g_num_pp( nod, nels_pp ) )
allocate( g_coord_pp( nod, ndim, nels_pp ) )
allocate( rest( nr,nodof+1) )
allocate( nf(nodof+1,nn) )
allocate(etype(nels_pp))
allocate(g_coord(ndim,nn))
g_num_pp = 0
g_coord_pp = zero
rest = 0

! Read element connectivity
CALL read_g_num_pp( argv, iel_start, nn, npes, numpe, g_num_pp )

IF ( meshgen == 2 ) CALL abaqus2sg( element, g_num_pp )
CALL read_g_coord_pp( argv, g_num_pp, nn, npes, numpe, g_coord_pp )

CALL read_g_coord(argv,nn,g_coord)

CALL read_rest( argv, numpe, rest )

do i = 1,nn
  nf(1,i)=i
  do j = 2,nodof+1
    nf(j,i)=1
  enddo
enddo

do i = 1,nn
  do j = 1,ubound(rest,1)
    if (i.eq.rest(j,1)) then
      do k = 2,nodof+1
        nf(k,i)=rest(j,k)
      enddo      
    endif
  enddo
enddo

do i = 1,nels_pp
  etype(i)=1
enddo

timest(2) = elap_time()

ALLOCATE( points(nip,ndim) )
allocate( dee(nst,nst) )
allocate( jac(ndim,ndim) )
allocate( der(ndim,nod) )
allocate( deriv( ndim, nod ) )
allocate( bee( nst, ntot ) )
allocate( weights( nip ) )
allocate( eps( nst ) )
allocate( sigma(nst) )
allocate( storkm_pp( ntot, ntot, nels_pp ) )
allocate( pmul_pp( ntot, nels_pp ) )
allocate( utemp_pp(ntot,nels_pp) )
allocate( g_g_pp(ntot,nels_pp) )

IF ( numpe==1 ) THEN
   open(  11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='write' )
   write( 11,'(A,I7,A)') "This job ran on ",npes," processes"
   write( 11,'(A,3(I12,A))') "There are ",nn," nodes", nr, &
        " restrained and ",neq," equations"
   write( 11,'(A,F10.4)') "Time to read input is:", timest(2)-timest(1)
   write( 11,'(A,F10.4)') "Time after setup is:", elap_time()-timest(1)
END IF
!*** end of ParaFEM input and initialisation *************************72

100 continue

!*** CGPACK part *****************************************************72
! *** CGPACK first executable statement ***
! In this test set the number of images via the env var,
! or simply as an argument to aprun.
! the code must be able to cope with any value >= 1.
  cgca_img = this_image()
cgca_nimgs = num_images()

! Need to init RND before cgca_pfem_cenc, where there order of comms
! is chosen at *random* to even the remote access pattern.

! Initialise random number seed. Choose either a routine with
! reproducible seeds (cgca_ins), or random seeds (cgca_irs).
! Multiple runs with cgca_ins *on the same number of cores (images)*
! on the same platform should produce reproducible results.
!
! Argument:
! .false. - no debug output
!  .true. - with debug output
call cgca_ins( .true. )

! dump CGPACK parameters and some ParaFEM settings
if ( cgca_img .eq. 1 ) then
  call cgca_pdmp
  write (*,*) "Young's mod:", e, "Poisson's ratio", v
end if

! Try to separate stdout output
sync all

! The disk is of 10mm diameter, 5mm thick.
! The bottom surface nodes are fixed and the corresponding top nodes are applied
! vertical loads.

! The physical dimensions of the box, must be the same
! units as in the ParaFEM. 
! 
! The disk diameter is 10mm, so the box fully encloses
! the diameter.
cgca_bsz = (/ 62.5, 60.0, 25.0 /)

! Origin of the box in the finite element cs, in the same units.
cgca_origin = (/ 0.0, 0.0, 0.0 /)

! Rotation tensor *from* FE cs *to* CA cs.
! The box cs is aligned with the box.
cgca_rot         = cgca_zero
cgca_rot( 1, 1 ) = 1.0
cgca_rot( 2, 2 ) = 1.0
cgca_rot( 3, 3 ) = 1.0

! mean grain size, also mm
cgca_dm = 10.0e0_rdef

! Scale the grain size in all 3 crystalline directions
cgca_cxyz(1)=1.0
cgca_cxyz(2)=1.0
cgca_cxyz(3)=1.0

! resolution, cells per grain
cgca_res = 1.0e5_rdef

! cgpack length scale, also in mm
! Equivalent to crack propagation distance per unit of time,
! i.e. per second. Let's say 1 km/s = 1.0e3 m/s = 1.0e6 mm/s. 
cgca_length = 1.0e6_rdef

! In p121_medium, each element is 0.25 x 0.25 x 0.25 mm, so
! the charlen must be bigger than that.
cgca_charlen = 0.4

! each image calculates the coarray grid dimensions
call cgca_gdim( cgca_nimgs, cgca_ir, cgca_qual )

! calculate the resolution and the actual phys dimensions of the box
! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
call cgca_cadim( cgca_bsz, cgca_res, cgca_dm, cgca_ir, cgca_c,         &
                 cgca_lres, cgca_ng, cgca_cxyz(:),argv(1:nlen)//".gid" )

! dump some stats from img 1
if (cgca_img .eq. 1 ) then
  write ( *, "(9(a,i0),tr1,g0,tr1,g0,3(a,g0),a)" )                     &
    "img: ", cgca_img  , " nimgs: ", cgca_nimgs,                       &
     " ("  , cgca_c (1), ","       , cgca_c (2), ",", cgca_c (3),      &
     ")["  , cgca_ir(1), ","       , cgca_ir(2), ",", cgca_ir(3),      &
     "] "  , cgca_ng   ,                                               &
    cgca_qual, cgca_lres,                                              &
         " (", cgca_bsz(1), ",", cgca_bsz(2), ",", cgca_bsz(3), ")"
  write (*,*) "dataset sizes for ParaView", cgca_c*cgca_ir
  write (*,"(a, es10.2, a, i0)") "Total cells in the model (real): ",  &
    product( real(cgca_c) * real(cgca_ir) ), " (int): ",               &
    product( int(cgca_c, kind=ilrg) * int(cgca_ir, kind=ilrg) )
  write(*,*)'Linear resolution:',cgca_lres
  write(*,*)'No. of grains:',cgca_ng
end if

! Allocate space coarray with 2 layers, implicit SYNC ALL inside.
!subroutine cgca_as( l1, u1, l2, u2, l3, u3,                           &
!             col1, cou1, col2, cou2, col3, props, coarray )
call cgca_as( 1, cgca_c(1),  1, cgca_c(2),  1, cgca_c(3),              &
              1, cgca_ir(1), 1, cgca_ir(2), 1, 2, cgca_space )

! Calculate the phys. dim. of the coarray on each image
!subroutine cgca_imco( space, lres, bcol, bcou )
call cgca_imco( cgca_space, cgca_lres, cgca_bcol, cgca_bcou )

! dump box lower and upper corners from every image
write ( *,"(a,i0,2(a,3(es9.2,tr1)))" ) "img ", cgca_img,               &
       " bcol: ", cgca_bcol, "bcou: ", cgca_bcou

! and now in FE cs:
write ( *,"(a,i0,2(a,3(es9.2,tr1)),a)" ) "img: ", cgca_img,            &
   " FE bcol: (",                                                      &
    matmul( transpose( cgca_rot ),cgca_bcol ) + cgca_origin,           &
  ") FE bcou: (",                                                      &
    matmul( transpose( cgca_rot ),cgca_bcou ) + cgca_origin, ")"

! confirm that image number .eq. MPI process number
write (*,*) "img",cgca_img," <-> MPI proc", numpe

! Allocate the tmp centroids array: cgca_pfem_centroid_tmp%r ,
! an allocatable array component of a coarray variable of derived type.
call cgca_pfem_ctalloc( ndim, nels_pp )

! Set the centroids array component on this image, no remote comms.
! first dim - coord, 1,2,3
! second dim - element number, always starting from 1
! g_coord_pp is allocated as g_coord_pp( nod, ndim, nels_pp )
cgca_pfem_centroid_tmp%r = sum( g_coord_pp(:,:,:), dim=1 ) / nod

         ! set cgca_pfem_centroid_tmp[*]%r
sync all ! must add execution segment
         ! use cgca_pfem_centroid_tmp[*]%r
! Set lcentr private arrays on every image. Choose one of the two
! routines that do this:
! - cgca_pfem_cenc - uses all-to-all algorithm.
! - cgca_pfem_map  - uses CO_SUM, CO_MAX and *large* tmp arrays
! Both routines have identical sets of input arguments.
call cgca_pfem_cenc( cgca_origin, cgca_rot, cgca_bcol, cgca_bcou )

! Allocate cgca_pfem_integrity%i(:), array component of a coarray of
! derived type. Allocating *local* array.
! i is set to 1.0 on allocation.
call cgca_pfem_integalloc( nels_pp )
 
! Allocate the Young's modulus 2D array
call cgca_pfem_ealloc( nip, nels_pp )
  
! initially set the Young's modulus to "e" everywhere
cgca_pfem_enew = e

! Generate microstructure
! Allocate rotation tensors, implicit SYNC ALL inside.
call cgca_art( 1, cgca_ng, 1, cgca_ir(1), 1, cgca_ir(2), 1, cgca_grt )

! Set initial values to all layers of the space array.
cgca_space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
cgca_space( :, :, :, cgca_state_type_frac  ) = cgca_intact_state

! Make sure all images set their space arrays, before calling
! the grain nucleation routine, which will update the state of
! the space array.
sync all

! Set grain nuclei, SYNC ALL inside.
! last argument:
! .false. - no debug output
!  .true. - with debug output
call cgca_nr( cgca_space, cgca_ng, .false. )

! Assign rotation tensors, SYNC ALL inside
call cgca_rt( cgca_grt,cgca_cxyz,argv(1:nlen)//".tra" )

!
!! Solidify, SYNC ALL inside.
!subroutine cgca_sld( coarray, periodicbc, iter, heartbeat, solid )
! second argument:
!  .true. - periodic BC
! .false. - no periodic BC
call cgca_sld( cgca_space, .false., 0, 10, cgca_solid )

! initiate grain boundaries
call cgca_igb( cgca_space )
sync all
call cgca_hxi( cgca_space )
sync all

! Smoothen the GB, several iterations. cgca_gbs has not remote comms.
! cgca_hxi has remote comms, so need to sync before and after it.
call cgca_gbs( cgca_space )
sync all
call cgca_hxi( cgca_space )
sync all
call cgca_gbs( cgca_space )
sync all
call cgca_hxi( cgca_space )
sync all

! update grain connectivity, local routine, no sync needed
call cgca_gcu( cgca_space )

! Set crack nuclei in the centre of the x1=x2=0 edge
if ( cgca_img .eq. 1 ) then
!  cgca_space( 1, 1, cgca_c(3)/2, cgca_state_type_frac )                &
!!              [ 1, 1, cgca_ir(3)/2 ] = cgca_clvg_state_100_edge
!                [ 1, 1, 1 ] = cgca_clvg_state_100_edge

!!  cgca_space( 10, 15, 10, cgca_state_type_frac )                &
!!              [ 1, 1, 1 ] = cgca_clvg_state_100_edge
!!   iseed=1
!--------------------------------------
! Uniform seeding of CA array
!--------------------------------------
incr1=100
incr2=100
incr3=50
iseed=0
do i = 1,cgca_c(1),incr1
  do j = 1,cgca_c(2),incr2
    do k = 1,cgca_c(3),incr3
        iseed=iseed+1
!!        cgca_space(i,j,k,cgca_state_type_frac)[1,1,1] = cgca_clvg_state_100_edge
    enddo
  enddo
enddo
write(*,*)'iseed=',iseed
end if
!----------------------------------------------
! Random seeding of CA array
!----------------------------------------------
! Compute the no. of seeds in each direction for the given increments
if ( cgca_img .eq. 1 ) then
  iseedx=0
  do i = 1,cgca_c(1),incr1
    iseedx=iseedx+1
  enddo

  iseedy=0
  do i = 1,cgca_c(2),incr2
    iseedy=iseedy+1
  enddo

  iseedz=0
  do i = 1,cgca_c(3),incr3
    iseedz=iseedz+1
  enddo

  allocate(icgca_cx(iseedx),source=0,stat=errstat)
  allocate(icgca_cy(iseedy),source=0,stat=errstat)
  allocate(icgca_cz(iseedz),source=0,stat=errstat)

! Generate random numbers in each direction
! Factoring by cgca_c is necessary, as the numbers generated are restricted
! within 0...1
! Restart seeding and generate new set of numbers at each call of RANDOM_NUMBER
  call random_seed()                    
  do i =1,iseedx
    CALL RANDOM_NUMBER(num)
    icgca_cx(i)=floor(num*(cgca_c(1)+1))
    if (icgca_cx(i).eq.0) icgca_cx(i)=1 ! Exclude zero numbers
  enddo

  call random_seed()
  do i =1,iseedy
    CALL RANDOM_NUMBER(num)
    icgca_cy(i)=floor(num*(cgca_c(2)+1))
    if (icgca_cy(i).eq.0) icgca_cy(i)=1
  enddo

  call random_seed()
  do i =1,iseedz
    CALL RANDOM_NUMBER(num)
    icgca_cz(i)=floor(num*(cgca_c(3)+1))
    if (icgca_cz(i).eq.0) icgca_cz(i)=1
  enddo

  write(*,*)'iseedx,iseedy,iseedz',iseedx,iseedy,iseedz
  write(*,*)'icgca_cx',icgca_cx
  write(*,*)'icgca_cy',icgca_cy
  write(*,*)'icgca_cz',icgca_cz

  iseed=0
  do i = 1,iseedx
    do j = 1,iseedy
      do k = 1,iseedz
          iseed=iseed+1
          cgca_space(icgca_cx(i),icgca_cy(j),icgca_cz(k),                 &
          cgca_state_type_frac)[1,1,1] = cgca_clvg_state_100_edge
      enddo
    enddo
  enddo
  write(*,*)'iseed',iseed
endif

         ! cgca_space changed locally on every image
sync all !
         ! cgca_space used

! Now can deallocate the temp array cgca_pfem_centroid_tmp%r.
! Could've done this earlier, but best to wait until sync all is
! required, to avoid extra sync.
call cgca_pfem_ctdalloc

! Img 1 dumps space arrays to files
! Remote comms, no sync inside, so most likely want to sync afterwards
if ( cgca_img .eq. 1 ) write (*,*) "dumping model to files"
call cgca_pswci( cgca_space, cgca_state_type_grain, "zg0.raw" )
call cgca_pswci( cgca_space, cgca_state_type_frac,  "zf0.raw" )
if ( cgca_img .eq. 1 ) write (*,*) "finished dumping model to files"

! Allocate the stress array component of cgca_pfem_stress coarray
! subroutine cgca_pfem_salloc( nels_pp, intp, ncomp )
call cgca_pfem_salloc( nels_pp, nip, nst )

! Need a sync after file write, because not sure what is coming next...
sync all

!*** end CGPACK part *************************************************72

! **************************Start parafem part************************
  !----------  find the steering array and equations per process ---------
  CALL rearrange(rest)
  g_g_pp = 0
  neq = 0
  ! map equations to elements
  elements_0: DO iel=1,nels_pp
     CALL find_g3( g_num_pp(:,iel), g_g_pp(:,iel), rest )
  END DO elements_0

  neq = MAXVAL(g_g_pp)
  neq = max_p(neq)
  CALL calc_neq_pp
  CALL calc_npes_pp( npes, npes_pp )
  CALL make_ggl( npes_pp, npes, g_g_pp )

  ALLOCATE( p_pp(neq_pp), source=zero )
  ! loads in the current increment, a fraction of tot_r_pp
  allocate( r_pp(neq_pp), source=zero )
  ! total loads
  allocate( tot_r_pp(neq_pp), source=zero )
  allocate( x_pp(neq_pp), source=zero )
  allocate( xnew_pp(neq_pp), source=zero )
  allocate( u_pp(neq_pp), source=zero )
  allocate( d_pp(neq_pp), source=zero )
  allocate( diag_precon_pp(neq_pp), source=zero )

  ! do this only once per program
  IF ( loaded_nodes > 0 ) THEN
     ALLOCATE( node(loaded_nodes), source=0 )
     allocate( val(ndim, loaded_nodes), source=zero )
     CALL read_loads( argv, numpe, node, val )
     CALL load( g_g_pp, g_num_pp, node, val, tot_r_pp(1:) )
     DEALLOCATE( node )
!     deallocate( val )
  end if
  DEALLOCATE(g_g_pp)
  ! allocate arrays outside of the loop
  ALLOCATE( eld_pp(ntot,nels_pp), source=zero )
  fname='CPFE_UMAT_PARAM.txt'
  call read_umat_params(fname,props,nprops)
  call init_umat_params(stress,statev,ddsdde,sse,spd,scd,rpl,       &
                        ddsddt,drplde,drpldt,stran,dstran,          &
                        time,dtime,temp_umat,dtemp,predef,dpred,    &
                        cmname, ndi,nshr,nst,nstatv,props,          &
                        nprops,coords,drot,pnewdt,celent,           &
                        dfgrd0,dfgrd1,noel,nip,layer,kspt,          &
                        kstep,kinc)
  cgca_pfem_enew=props(1)
  v=props(2)
  ! Since the Young's modulus can be updated,
  ! all below must be inside the loading iterations!
  ! Make sure not to allocate/deallocate arrays within the loop
  load_iter: do cgca_liter=1, cgca_linum
     !------ element stiffness integration and build the preconditioner ---72

     ! deemat needs to go into a loop elements_1
     ! need an element array of E

     dee = zero

     CALL sample( element, points, weights )
     storkm_pp = zero

!     pause
     elements_1: DO iel=1,nels_pp
        gauss_pts_1: DO i=1,nip
!           CALL deemat( dee, cgca_pfem_enew(i,iel), v )
           props(1)=cgca_pfem_enew(i,iel)
           call umat(stress,statev,ddsdde,sse,spd,scd,rpl,       &
                     ddsddt,drplde,drpldt,stran,dstran,          &
                     time,dtime,temp_umat,dtemp,predef,dpred,    &
                     cmname, ndi,nshr,nst,nstatv,props,          &
                     nprops,coords,drot,pnewdt,celent,           &
                     dfgrd0,dfgrd1,noel,nip,layer,kspt,          &
                     kstep,kinc)

           CALL shape_der( der, points, i )
           jac = MATMUL( der, g_coord_pp(:,:,iel) )
           det = determinant(jac)
           CALL invert(jac)
           deriv = MATMUL( jac, der )
           CALL beemat( bee, deriv )
           storkm_pp(:,:,iel) = storkm_pp(:,:,iel) +                             &
!                MATMUL( MATMUL( TRANSPOSE(bee), dee ), bee ) * det * weights(i)
                MATMUL( MATMUL( TRANSPOSE(bee), ddsdde ), bee ) * det * weights(i)
        END DO gauss_pts_1
     END DO elements_1

     ALLOCATE( diag_precon_tmp( ntot, nels_pp ), source=zero )

     elements_2: DO iel=1,nels_pp
        DO i=1,ndof
           diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkm_pp(i,i,iel)
        END DO
     END DO elements_2

     CALL scatter( diag_precon_pp, diag_precon_tmp )
     DEALLOCATE( diag_precon_tmp )
     
     !----------------------------- get starting r --------------------------

     ! loading increases from 1/cgca_linum of tot_r_pp
     ! to tot_r_pp
     r_pp(1:) = tot_r_pp(1:) * cgca_liter/cgca_linum

     IF ( loaded_nodes > 0 ) THEN
        q = SUM_P( r_pp(1:) )
        write(*,*)'Total loads q-->',q
        IF ( numpe==1 ) then
           write (11, *) "Load iter:", cgca_liter
           write (11,'(A,E12.4)') "The total load is:", q
        end if
     END IF

     diag_precon_pp = 1._iwp/diag_precon_pp
     d_pp = diag_precon_pp*r_pp
     p_pp = d_pp
     x_pp = zero

     !--------------------- preconditioned cg iterations --------------------

     write(*,*)'PCG iterations start'
     iters=0
     timest(3) = elap_time()

     iterations: DO
        iters = iters + 1
        u_pp = zero
        pmul_pp = zero
        utemp_pp = zero
        CALL gather(p_pp,pmul_pp)
        elements_3: DO iel=1,nels_pp
           utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
        END DO elements_3

        CALL scatter(u_pp,utemp_pp)

        !-------------------------- pcg equation solution ----------------------
        up = DOT_PRODUCT_P(r_pp,d_pp)
        alpha = up/DOT_PRODUCT_P(p_pp,u_pp)
        xnew_pp = x_pp+p_pp*alpha

        r_pp = r_pp-u_pp*alpha
        d_pp = diag_precon_pp*r_pp
        beta = DOT_PRODUCT_P(r_pp,d_pp)/up
        p_pp = d_pp+p_pp*beta
        CALL checon_par(xnew_pp,tol,converged,x_pp)
        IF ( converged .OR. iters == limit ) EXIT
     END DO iterations

     IF ( numpe==1 ) THEN
        write(11,'(A,I6)')"The number of iterations to convergence was ",iters
        write(11,'(A,F10.4)')"Time to solve equations was  :",                &
             elap_time()-timest(3)
        write(11,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(1)
     END IF


!---------------   Recover stresses at Gauss points   ----------------72

     CALL gather(xnew_pp(1:),eld_pp)

     elmnts: DO iel = 1, nels_pp
        intpts: DO i = 1, nip
!           call deemat(dee, cgca_pfem_enew( i, iel ), v)
           props(1)=cgca_pfem_enew(i,iel)
           call umat(stress,statev,ddsdde,sse,spd,scd,rpl,       &
                     ddsddt,drplde,drpldt,stran,dstran,          &
                     time,dtime,temp_umat,dtemp,predef,dpred,    &
                     cmname, ndi,nshr,nst,nstatv,props,          &
                     nprops,coords,drot,pnewdt,celent,           &
                     dfgrd0,dfgrd1,noel,nip,layer,kspt,          &
                     kstep,kinc)
           ! Compute the derivatives of the shape functions at a Gauss point.
           CALL shape_der(der,points,i)
           jac = MATMUL(der,g_coord_pp(:,:,iel))
           CALL invert(jac)
           deriv = MATMUL(jac,der)
           CALL beemat(bee,deriv)
           eps = MATMUL(bee,eld_pp(:,iel))
!           sigma = MATMUL(dee,eps)
           sigma = MATMUL(ddsdde,eps)

           ! set the stress array on this image
           ! from cgca_m3pfem.f90:
           cgca_pfem_stress%stress( iel, i, : ) = sigma
        END DO intpts
     end do elmnts

!**************** end ParaFEM part ***********************************72
! Calculate mean stresses for the upper surface in the vertical direction
     stot_upfc=0
     a_upfc=62.5d0
     b_upfc=25.0d0
     do j=1,size(elem_id_upfc)
       iel=elem_id_upfc(j)
       do i = 1,nip
         stot_upfc=stot_upfc+cgca_pfem_stress%stress( iel, i, 2 )
       enddo
     enddo
     smean_upfc=stot_upfc/(size(elem_id_upfc)*nip)

!*** CGPACK part *****************************************************72
! debug: dump stresses to stdout
!call cgca_pfem_sdmp
! dump stresses from last image for element 1
if ( cgca_img .eq. cgca_nimgs ) then
  do i = 1, nip
    write (*,"(2(a,i0),a,6es10.2)") "img: ", cgca_nimgs,               &
          " FE 1 int. p. ", i, " stress: ",                            &
          cgca_pfem_stress % stress( 1 , i , : )
  end do
end if

! all images sync here
sync all

! Calculate the mean stress tensor per image
!   subroutine cgca_pfem_simg( simg )
!   real( kind=rdef ), intent(out) :: simg(3,3)
call cgca_pfem_simg( cgca_stress )
write (*,*) "img:", cgca_img, " mean s tensor:", cgca_stress

nel=size(lcentr)
nintp=size(cgca_pfem_stress%stress, dim=2)
allocate(str_tmp(nintp,ncomp), source=0.0_rdef)
allocate(stot(nel,ncomp),source=0.0_rdef)

maxeval=0
! Loop through the elements to obtain the max. principal stress
do el=1,nel
  rimg = lcentr(el) % image
  rel = lcentr(el) % elnum
  str_tmp=cgca_pfem_stress[ rimg ] % stress( rel, : , : )

  ! Mean of the stresses for all gauss pts in the element
  do j = 1,nintp
    stot(el,:)=stot(el,:)+str_tmp(j,:) 
  enddo
  do i=1,ncomp
    stot(el,i)=stot(el,i)/real(nintp,kind=rdef)
  enddo

  ! Construct a symmetric stress tensor
  stressTensor(1,1)=stot(el,1)
  stressTensor(2,2)=stot(el,2)
  stressTensor(3,3)=stot(el,3)
  stressTensor(1,2)=stot(el,4)
  stressTensor(1,3)=stot(el,5)
  stressTensor(2,3)=stot(el,6)
  stressTensor(2,1)=stressTensor(1,2)
  stressTensor(3,1)=stressTensor(1,3)
  stressTensor(3,2)=stressTensor(2,3)

  EigenVal=1
  steps=10

  ! The eigen values corresponds to the magnitude of the principal stress
  ! components. Here the maximum principal stress arising out of the elements
  ! are chosen and are applied to simulate the cleavege fracture
  call FINDEIGEN(stressTensor, 3, EigenVec, EigenVal, steps)
    if(EigenVal .gt. maxeval) then
      maxeval=EigenVal
      maxindx=el
      maxelemN=rel
    endif  
enddo
write(20,*)'maxelemN',maxelemN

! Stress tensor of the element corresponding to maximum principal stress
stressTensor(1,1)=stot(maxindx,1)
stressTensor(2,2)=stot(maxindx,2)
stressTensor(3,3)=stot(maxindx,3)
stressTensor(1,2)=stot(maxindx,4)
stressTensor(1,3)=stot(maxindx,5)
stressTensor(2,3)=stot(maxindx,6)
stressTensor(2,1)=stressTensor(1,2)
stressTensor(3,1)=stressTensor(1,3)
stressTensor(3,2)=stressTensor(2,3)

! all images wait for each other, to make sure the stress arrays
! are not modified until all images calculate their mean values
sync all

! no real time increments in this problem
! I use the inverse of the length scale,
! which gives 1mm of crack propagation per increment maximum.
! I then can multiply it by a factor, e.g. a factor of 3 will mean
! that I get 3mm max ( 1/3 of the model ) per load increment.
cgca_time_inc = 2 * 1.0_rdef / cgca_length

! run cleavage for a correct number of iterations, which is a function
! of the characteristic length and the time increment
cgca_clvg_iter = nint( cgca_length * cgca_lres * cgca_time_inc )
if ( cgca_img .eq. 1 ) write (*,*) "load inc:", cgca_liter,            &
                                   "clvg iter:", cgca_clvg_iter

! ===>>> sync all inside <<<===
! lower the crit stresses by a factor of 100.
! On Intel 16: no support for CO_SUM yet, so use _nocosum.
! On Intel 16: subroutine cgca_clvgp_nocosum( coarray, rt, t, scrit, &
!                     sub, gcus, periodicbc, iter, heartbeat, debug )
cgca_stress=stressTensor
write(*,*)'cgca_stress',stressTensor

call cgca_clvgp_nocosum( cgca_space, cgca_grt, cgca_stress,            &
     0.01_rdef * cgca_scrit, cgca_clvgsd, cgca_gcupdn, .false. ,       &
     cgca_clvg_iter, 10, cgca_yesdebug)

! dump the model out, no sync inside
if ( cgca_img .eq. 1 ) write (*,*) "dumping model to file"
write ( cgca_citer, "(i0)" ) cgca_liter
if (mod(cgca_liter,5).eq.0.0d0) then
  call cgca_pswci( cgca_space, cgca_state_type_frac,                   &
                 "zf"//trim( cgca_citer )//".raw" )
endif
if ( cgca_img .eq. 1 ) write (*,*) "finished dumping model to file"

sync all
     
! Calculate number (volume) of fractured cells on each image.
! cgca_fracvol is a local, non-coarray, array, so no sync needed.
call cgca_fv( cgca_space, cgca_fracvol )

! calculate integrity, remote write
!!call cgca_pfem_intcalc1( cgca_c, cgca_fracvol )
! Calculate integrity due to damage occuring on a single element alone
call cgca_pfem_intcalc2( iseed, maxelemN, cgca_c, cgca_fracvol )

! dump min integrity for all FE stored on this image
write (*,*) "img:", cgca_img, "min. integrity:",                       &
  minval( cgca_pfem_integrity % i )

! wait for integrity on all images to be calculated
sync all

! Young's modulus need to be updated on each image, local arrays only.
! The original Young's modulus value must be given as input.
! For each FE stored on this image.
call cgca_pfem_uym( e, nels_pp )

! dump updated Young's modulus
write (*,*) "img:", cgca_img, "*min* Young's mod:",                    &
            minval( cgca_pfem_enew )

sync all

!*** end CGPACK part *************************************************72


!*** ParaFEM part ****************************************************72

! end loading iterations
!!end do load_iter

!------------------------ write out displacements ----------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)

  IF ( numpe==1 ) THEN
     write(ch,'(I6.6)') numpe
     open( 12, file=argv(1:nlen)//".ensi.DISPL-"//ch,status='unknown',    &
          action='write')
     write(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
     write(12,'(A/A/A)') "part", "     1","coordinates"
  END IF
  
  ALLOCATE( disp_pp(nodes_pp*ndim), source=zero )
  allocate( temp(nodes_pp), source=zero )
  CALL scatter_nodes( npes, nn, nels_pp, g_num_pp, nod, ndim, nodes_pp,  &
       node_start, node_end, eld_pp, disp_pp, 1 )

  DO i=1,ndim
     temp=zero
     DO j=1,nodes_pp
        k = i+(ndim*(j-1))
        temp(j) = disp_pp(k)
     END DO
     write(*,*)'Writing displacements for step #:',cgca_liter
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
  END DO

! Obtain mean displacements on the upper surface in the vertical direction
  dtot_upfc=0
  DO i=1,size(node_id_upfc)
    inod=node_id_upfc(i)
     k=2+(ndim*(inod-1))
!     k=1+(ndim*(inod-1))
     dtot_upfc=dtot_upfc+disp_pp(k)
  ENDDO
  dmean_upfc=dtot_upfc/size(node_id_upfc)

! Output applied loads and displacements on the upper surface
     pres_upfc = val(2,1)
     pres_upfc = pres_upfc * cgca_liter/cgca_linum
     TLoad_upfc=smean_upfc*a_upfc*b_upfc
     open(UNIT=OUTPUT_UNIT, FILE='LDST.txt', ACTION='write',      &
          STATUS='unknown', ACCESS='append', iOSTAT=ios)
     if (cgca_liter.eq.1) write(OUTPUT_UNIT,'(5A15)')'Pressure(MPa)','Loads(N)','Displacements(mm)','Stiffness(L/D)', 'Y.modulus(GPa)'
     write(OUTPUT_UNIT,'(5F15.4)')pres_upfc, TLoad_upfc, dmean_upfc,TLoad_upfc/dmean_upfc, minval( cgca_pfem_enew )

DEALLOCATE( disp_pp )
DEALLOCATE( temp )
DEALLOCATE( str_tmp )
DEALLOCATE( stot )

close( 12 )

! end loading iterations
end do load_iter

! Generate paraview resultant files 
elemtype='hexahedron'
CALL mesh_ensi(argv,nlen,g_coord,g_num_pp,elemtype,etype,nf,    &
              tot_r_pp(1:),1,1,1.0d0,.true.)

! deallocate all arrays, moved from inside the loop
DEALLOCATE( p_pp )
deallocate( r_pp )
deallocate( x_pp )
deallocate( u_pp )
deallocate( d_pp )
deallocate( diag_precon_pp )
deallocate( storkm_pp )
deallocate( pmul_pp )
DEALLOCATE( xnew_pp )
DEALLOCATE( g_coord_pp )

  IF ( numpe==1 ) then
     write( 11, '(A,F10.4)') "This analysis took: ", elap_time()-timest(1)
     close( 11 )
!     close( 12 )
  end if
!*** end ParaFEM part ************************************************72


!*** CGPACK part *****************************************************72
! deallocate all CGPACK arrays
call cgca_ds(  cgca_space )
call cgca_drt( cgca_grt )
call cgca_pfem_sdalloc
call cgca_pfem_edalloc
call cgca_pfem_integdalloc
!*** end CGPACK part *************************************************72

write (*,*) "dataset sizes for ParaView", cgca_c*cgca_ir

110 format(F7.2,A1,F7.2,A1,F7.2,A1,I5)

END PROGRAM SquareBlock_CPFE

!-------------------------------------------------------------------------
subroutine read_elem_upfc(elem_id_upfc)
!--------------------------------------------------------------------------
integer :: funit,ios,i,elemid,istat,errstat
integer :: ielem
character(50) :: keywd
integer, allocatable, intent(inout) :: elem_id_upfc(:)
!--------------------------------------------------------------------------
funit=50
open(UNIT=funit, FILE='abaqus_Elem_Reac.rpt', ACTION='read', IOSTAT=ios)

! Compute size for holding the element array 
do while (keywd(1:3) .ne. 'Loc')
  read(funit,*)keywd
enddo
do i =1,5
  read(funit,*)
enddo
read(funit,*)keywd,elemid
backspace funit

ielem=0
do while (keywd(1:1) .eq. 'P')
  ielem=ielem+1
  read(funit,*,iostat=istat)keywd,elemid
enddo

allocate(elem_id_upfc(ielem-1), source=0, stat=errstat)

rewind funit

! Read actual element list
do while (keywd(1:3) .ne. 'Loc')
  read(funit,*)keywd
end do
do i =1,5
  read(funit,*)
enddo

read(funit,*)keywd,elemid
backspace funit

ielem=0
do while (keywd(1:1) .eq. 'C')
  ielem=ielem+1
  read(funit,*,iostat=istat)keywd,elemid
  elem_id_upfc(ielem)=elemid
enddo
close(UNIT=funit)

end subroutine read_elem_upfc


!-------------------------------------------------------------------------
subroutine read_node_upfc(node_id_upfc)
!--------------------------------------------------------------------------
integer :: funit,ios,i,nodeid,istat,errstat
integer :: inode
character(50) :: keywd
integer, allocatable, intent(inout) :: node_id_upfc(:)
!--------------------------------------------------------------------------
open(UNIT=funit, FILE='abaqus_Nodal_Disp.rpt', ACTION='read', IOSTAT=ios)

! Compute size for holding the element array 
do while (keywd(1:3) .ne. 'Var')
  read(funit,*)keywd
end do
do i =1,6
  read(funit,*)
enddo

read(funit,*)keywd,nodeid
backspace funit

inode=0
do while (keywd(1:1) .eq. 'P')
!  inode=inode+1
  read(funit,*,iostat=istat)keywd,nodeid
  if (nodeid.ne.0) inode=inode+1
enddo

allocate(node_id_upfc(inode), source=0, stat=errstat)
rewind funit

! Read actual node list
do while (keywd(1:3) .ne. 'Var')
  read(funit,*)keywd
end do
do i =1,6
  read(funit,*)
enddo

read(funit,*)keywd,nodeid
backspace funit

inode=0
do while (keywd(1:1) .eq. 'P')
  inode=inode+1
  read(funit,*,iostat=istat)keywd,nodeid
  if (nodeid.ne.0) node_id_upfc(inode)=nodeid
enddo
end subroutine read_node_upfc

!--------------------------------------------------------------------------------------------------------
!Author : Louisda16th a.k.a Ashwith J. Rego
!These set of subroutines find the largest eigenvalue and eigenmatrix of the
!matrix.
!The algorithm is based on Rayleigh's power method
!Please note that the subroutine used to multiply the two matrices is not
!general.
!Also note that the number of iterations must be specified
!--------------------------------------------------------------------------------------------------------
SUBROUTINE FINDEIGEN(Matrix, n, x, EigenVal, steps)
    use cgca_m1co, only : rdef
    INTEGER, INTENT(IN) :: n, steps  !n = order of matrix, steps = number of iterations
    REAL, INTENT(IN), DIMENSION(n,n) :: Matrix(n,n)  !Input Matrix
    REAL, INTENT(INOUT), DIMENSION(n) :: x !Eigenvector
    REAL, INTENT(INOUT) :: EigenVal !Eigenvalue
    INTEGER :: i, j

    x  = 1 !Initialize eigen vector to any value.
    DO i = 1, steps
        CALL MULMATRIX(Matrix, x, n)       !Multiply input matrix by eigenvector
        CALL FINDLARGEST(x, n, EigenVal)   !Find eigenvalue
        IF(EigenVal == 0) EXIT         
        DO j = 1, n                        !Find eigenvector
            x(j) = x(j)/EigenVal
        END DO 
    END DO
END SUBROUTINE FINDEIGEN

SUBROUTINE MULMATRIX(a, b, n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n !matrix size
    REAL, INTENT(IN), DIMENSION(n,n) :: a  !Matrix of order > 1
    REAL, INTENT(INOUT), DIMENSION(n) :: b !1x1 matrix  
    INTEGER i, j
    REAL, DIMENSION(n) :: temp !temporary matrix
 
    temp = 0
    !These two loops to the multiplication
    DO i = 1, n
        DO j = 1, n
            temp(i) = temp(i) + a(i,j)*b(j)
        END DO
    END DO
    b = temp

END SUBROUTINE MULMATRIX
 
SUBROUTINE FINDLARGEST(x, n, l)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN), DIMENSION(n) :: x
    REAL, INTENT(INOUT) :: l !Largest value
     
    INTEGER :: i
    !Algorithm is easy
    !Let the largest number be the first one.
    !If you find a number larger than it, store this number and then continue

    l = ABS(x(1))
    DO i = 2, n
        IF (ABS(x(i)) > l) l = ABS(x(i))
    END DO
         
END SUBROUTINE FINDLARGEST

