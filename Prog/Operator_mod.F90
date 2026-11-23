!  Copyright (C) 2016 - 2023 The ALF project
! 
!  This file is part of the ALF project.
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Defines the Operator type, and provides a number of operations on this type. 
!
!--------------------------------------------------------------------
Module Operator_mod

  Use mpi_shared_memory
  Use mat_subroutines
  Use MyMats
  Use Fields_mod
  Use runtime_error_mod
  Use iso_fortran_env, only: error_unit
  Use Natural_Constants, only: Eps_machine, Eps_small
  
  Implicit none
  
  
  
  Type Operator
     Integer          :: N, N_non_zero                      !> dimension of Operator (P and O), number of non-zero eigenvalues
     Integer, private :: win_M_exp, win_U                   !> MPI_windows which can be used for fences (memory synch.) and dealloc.
     logical          :: diag                               !> encodes if Operator is diagonal
     logical, private :: U_alloc, M_exp_alloc, g_t_alloc    !> logical to track if memory is allocated
     complex (Kind=Kind(0.d0)), pointer :: O(:,:), U (:,:)  !> Storage for operator matrix O and its eigenvectors U
     complex (Kind=Kind(0.d0)), pointer, private :: M_exp(:,:,:), E_exp(:,:)  !> Internal storage for exp(O) and exp(E)
     Real    (Kind=Kind(0.d0)), pointer :: E(:)             !> Eigenvalues of O
     Integer, pointer :: P(:)                               !> Projector P encoding DoFs that contribute in Operator
     complex (Kind=Kind(0.d0)) :: g                         !> coupling constant
     complex (Kind=Kind(0.d0)), allocatable :: g_t(:)       !> time dependent  coupling constant
     complex (Kind=Kind(0.d0)) :: alpha                     !> operator shift
     Integer          :: Type                               !> Type of the operator: 1=Ising; 2=discrete HS; 3=continuous scalar HS, 4=  Complex for  three  bondy term.
     Integer          :: Flip_protocol=1                    !> Flip protocol  for  local  updates.  Only  relevant  for  type =3  fields. 
     ! P is an N X Ndim matrix such that  P.T*O*P*  =  A  
     ! P has only one non-zero entry per column which is specified by P
     ! All in all.   g * Phi(s,type) * ( c^{dagger} A c  + alpha )
     ! The variable Type allows you to define the type of HS. 
     ! The first N_non_zero elements of diagonal matrix E are non-zero. The rest vanish.

     ! !!!!! M_exp and E_exp  are for storage   !!!!!
     ! If Type =1   then the Ising field  takes the values  s = +/- 1 
     !              and M_exp   has dimensions  M_exp(N,N,3)   last index=1 rep field -1 and index=3 rep field 1 
     ! If Type =2   then the Ising field  takes the values  s = +/- 1,  +/- 2
     !              and M_exp   has dimensions  M_exp(N,N,5)
     ! M_exp(:,:,s) =  e^{g * Phi(s,type) *  O(:,:) }
     !
     ! 
     ! E_exp(:,s) = e^{g * Phi(s,type) *E(:) } and has dimensions E_exp(N,-1:1)  for Type = 1
     !              and dimensions E_exp(N,-2:2)   for Type = 2
     !
     ! !!! If Type .neq. 1,2  then  E_exp  and  M_exp  are  not allocated !!!

     ! Low-rank path (Op_set decides; see op_lowrank_wanted): the mmults apply
     ! e^{g phi O} = 1 + V (e^{g phi E} - 1) V^dag with V = U(:,1:N_non_zero), and the
     ! wraps rotate by the product Q of N_non_zero Householder reflectors (QR of V)
     ! instead of U. Q's first N_non_zero columns are U's up to a phase, which the
     ! update in the rotated basis does not see; only the null-space basis differs.
     ! U and O stay whole (Op%N x Op%N): Op_exp and other readers use every column.
     logical          :: lowrank = .false.
     complex (Kind=Kind(0.d0)), allocatable :: HY(:,:), HT(:,:)  !> Q = I - HY HT HY^dag (compact WY)
     Integer          :: P0 = 0                               !> P(1) if P is a contiguous ascending run, else 0
   contains
     procedure  :: get_g_t_alloc => operator_get_g_t_alloc     !>    External access to   private  logical  variable g_t_alloc
  end type Operator

  !> -1: read ALF_LOWRANK_VERTEX at the first Op_set (unset or nonzero: on); 0 off; 1 on.
  Integer :: Op_lowrank_mode = -1


Contains

!!$  Real (Kind=Kind(0.d0)) :: Phi(-2:2,2),  Gaml(-2:2,2)
!!$  Integer ::  NFLIPL(-2:2,3)
!!$      Subroutine Op_setHS
!!$
!!$        Implicit none
!!$        !Local
!!$        Integer :: n
!!$        
!!$
!!$        Phi = 0.d0
!!$        do n = -2,2
!!$           Phi(n,1) = real(n,Kind=Kind(0.d0))
!!$        enddo
!!$        Phi(-2,2) = - SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
!!$        Phi(-1,2) = - SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
!!$        Phi( 1,2) =   SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
!!$        Phi( 2,2) =   SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
!!$        
!!$        Do n = -2,2
!!$           gaml(n,1) = 1.d0
!!$        Enddo
!!$        GAML(-2,2) = 1.D0 - SQRT(6.D0)/3.D0
!!$        GAML( 2,2) = 1.D0 - SQRT(6.D0)/3.D0
!!$        GAML(-1,2) = 1.D0 + SQRT(6.D0)/3.D0
!!$        GAML( 1,2) = 1.D0 + SQRT(6.D0)/3.D0
!!$        
!!$        NFLIPL(-2,1) = -1
!!$        NFLIPL(-2,2) =  1
!!$        NFLIPL(-2,3) =  2
!!$        
!!$        NFLIPL(-1,1) =  1
!!$        NFLIPL(-1,2) =  2
!!$        NFLIPL(-1,3) = -2
!!$        
!!$        NFLIPL( 1,1) =  2
!!$        NFLIPL( 1,2) = -2
!!$        NFLIPL( 1,3) = -1
!!$        
!!$        NFLIPL( 2,1) = -2
!!$        NFLIPL( 2,2) = -1
!!$        NFLIPL( 2,3) =  1 
!!$
!!$      end Subroutine Op_setHS

  

!--------------------------------------------------------------------
!> @author
!> 
!
!> @brief 
!> calculate the phase of a given set of operators and HS fields.
!
!> @param[inout] Phase  Complex
!> * On entry: phase of \f$ \prod_{f} \det M_f(C) \f$, f is the flavor index.
!> * On exit:  phase of  \f$ W(C)  =   \left[ \left( \prod_{n,\tau,f}  \exp \left[ g_f(n) \alpha_f(n) \phi(\sigma(n,\tau)) \right] \right) \det(M_f(C))\right]^{N_{SUN}} \prod_{n,\tau }\gamma(\sigma(n,\tau)) \f$ 
!> @param[in] Op_V  Dimension(:,:)  Type(Operator)
!> * List of interaction operators. OP_V(n,f) has no tau index.  
!> @param[in] Nsigma 
!> Type(Fields)
!> * Fields
!> @param[in] nf
!> Integer
!> * flavor index
!--------------------------------------------------------------------
  Subroutine  Op_phase(Phase,OP_V,Nsigma,nf) 
    Implicit none

    Complex  (Kind=Kind(0.d0)), Intent(Inout) :: Phase
    Integer,           Intent(IN)    :: nf
    Type  (Fields),    Intent(IN)    :: Nsigma
    Type (Operator),   dimension(:,:), Intent(In) :: Op_V
    Real  (Kind=Kind(0.d0))                       :: angle
    
    Integer :: n, nt
    complex (kind=kind(0.d0)) :: g_loc
    
    do n = 1,size(Op_V,1)
       do nt = 1,size(nsigma%f,2)
          g_loc = Op_V(n,nf)%g
          if (op_v(n,nf)%g_t_alloc) g_loc = Op_V(n,nf)%g_t(nt)
          angle = Aimag( g_loc * Op_V(n,nf)%alpha * nsigma%Phi(n,nt) )
          Phase = Phase*CMPLX(cos(angle),sin(angle), Kind(0.D0))
       enddo
    enddo
    
  end Subroutine Op_phase
  
!--------------------------------------------------------------------
!> @author
!> 
!> @brief 
!> Set up _core_ data of the operator required for type=0 operators
!
!> @param[inout] Op
!> @param[in] N 
!--------------------------------------------------------------------

  Pure subroutine Op_make(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Allocate (Op%O(N,N), Op%P(N) )
    ! F.F.A  Op%M_exp and Op%E_exp are allocated  in Op_set once the type is available.
    
    Op%O = cmplx(0.d0, 0.d0, kind(0.D0))
    Op%P = 0
    Op%N = N
    Op%N_non_zero = N
    Op%g     = cmplx(0.d0,0.d0, kind(0.D0))
    Op%alpha = cmplx(0.d0,0.d0, kind(0.D0))
    Op%diag  = .false.
    Op%type = 0
    OP%flip_protocol = 1
    Op%U_alloc = .false.
    Op%M_exp_alloc = .false.
    Op%g_t_alloc = .false.
    Op%lowrank = .false.
  end subroutine Op_make

!--------------------------------------------------------------------
!> @author
!> 
!> @brief 
!> Deallocate the operator
!
!> @param[inout] Op
!> @param[in] N 
!--------------------------------------------------------------------

  subroutine Op_clear(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N

!     If ( associated(OP%O) )   deallocate(OP%O)
!     If ( associated(OP%P) )   deallocate(OP%P)
!     If ( associated(OP%E_exp) )   deallocate(OP%E_exp)
!     If ( associated(OP%E) )   deallocate(OP%E)
    if (use_mpi_shm) then
      If ( Op%U_alloc ) then
         deallocate(OP%O, OP%P, OP%E)
         !call deallocate_shared_memory(OP%win_U)
      endif
      If ( Op%M_exp_alloc ) then
         deallocate(OP%E_exp)
         !call deallocate_shared_memory(OP%win_M_exp)
      endif
    else
      If ( Op%M_exp_alloc )   deallocate(OP%M_exp, OP%E_exp)
      If ( Op%U_alloc )   deallocate(OP%U, OP%O, OP%P, OP%E)
      if ( Op%g_t_alloc )  deallocate(Op%g_t)
    endif
    if (allocated(Op%HY)) deallocate(Op%HY, Op%HT)
    Op%lowrank = .false.

  end subroutine Op_clear

!--------------------------------------------------------------------
!> @author
!> 
!> @brief 
!> Setup storage for type=1,2 and 3 vertices.  Setup the exponential of the operator for
!> type 1 and 2
!
!> @param[inout] Op  Type (Operator)
!--------------------------------------------------------------------
  subroutine Op_set(Op)
    Implicit none
    Type (Operator), intent(INOUT) :: Op

    Complex (Kind=Kind(0.d0)), allocatable :: U(:,:), TMP(:, :)
    Real    (Kind=Kind(0.d0)), allocatable :: E(:)
    Real    (Kind=Kind(0.d0)), parameter :: Zero = Eps_small
    Real    (Kind=Kind(0.d0)) :: herm_tol = Eps_machine
    Real    (Kind=Kind(0.d0)) :: herm_dev
    Integer :: N, I, J, np,nz, noderank, arrayshape2d(2), arrayshape(3)
#ifdef MPI
    Integer :: ierr
#endif
    Complex (Kind=Kind(0.d0)) :: Z
    Type  (Fields)   :: nsigma_single
    
    ! --- Validate Op%type ---
    ! Valid types: 0 (hopping/no HS field), 1 (Ising), 2 (discrete HS), 3 (continuous scalar HS), 4 (complex three-body)
    if (Op%type < 0 .or. Op%type > 4) then
       Write(error_unit, '(A,I0)') ' Op_set: Invalid operator type: ', Op%type
       Write(error_unit, '(A)')    '   Valid types are 0 (hopping), 1 (Ising), 2 (discrete HS), 3 (continuous scalar HS), 4 (complex three-body).'
       Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
    endif

    ! --- Validate projector Op%P ---
    do I = 1, Op%N
       if (Op%P(I) < 1) then
          Write(error_unit, '(A,I0,A,I0)') ' Op_set: Projector index P(', I, ') out of bounds: ', Op%P(I)
          Write(error_unit, '(A)')          '   All projector entries must be >= 1.'
          Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
       endif
    enddo

    ! --- Validate Hermiticity of Op%O and detect if diagonal (single pass) ---
    N = Op%N
    Op%diag = .true.
    if (N > 1) then
       do I = 1, N
          do J = I+1, N
             herm_dev = abs(Op%O(I,J) - conjg(Op%O(J,I)))
             if (herm_dev > herm_tol * max(abs(Op%O(I,J)), abs(Op%O(J,I)), 1.D-30)) then
                Write(error_unit, '(A)')        ' Op_set: Operator matrix Op%O is not Hermitian.'
                Write(error_unit, '(A,I0,A,I0,A,2ES15.7,A)') &
                     '   O(', I, ',', J, ') = (', real(Op%O(I,J)), aimag(Op%O(I,J)), ')'
                Write(error_unit, '(A,I0,A,I0,A,2ES15.7,A)') &
                     '   O(', J, ',', I, ') = (', real(Op%O(J,I)), aimag(Op%O(J,I)), ')'
                Write(error_unit, '(A,ES15.7)') '   |O(i,j) - conjg(O(j,i))| = ', herm_dev
                Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
             endif
             ! Binary comparison is OK here as Op%O was initialized to zero during Op_make.
             ! Checking only the upper-triangle element suffices: after the Hermiticity check above,
             ! Op%O(I,J)==0 implies Op%O(J,I)==conjg(0)==0.
             if (Op%diag) then
                if (Op%O(I,J) .ne. cmplx(0.d0,0.d0, kind(0.D0))) Op%diag = .false.
             endif
          enddo
       enddo
       ! Check that diagonal elements are real
       do I = 1, N
          if (abs(aimag(Op%O(I,I))) > herm_tol * max(abs(Op%O(I,I)), 1.D-30)) then
             Write(error_unit, '(A)')        ' Op_set: Operator matrix Op%O has complex diagonal elements (not Hermitian).'
             Write(error_unit, '(A,I0,A,I0,A,2ES15.7,A)') &
                  '   O(', I, ',', I, ') = (', real(Op%O(I,I)), aimag(Op%O(I,I)), ')'
             Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
          endif
       enddo
    else
       ! N = 1: diagonal element must be real for Hermiticity
       if (abs(aimag(Op%O(1,1))) > herm_tol * max(abs(Op%O(1,1)), 1.D-30)) then
          Write(error_unit, '(A)')        ' Op_set: 1x1 operator matrix has non-zero imaginary part (not Hermitian).'
          Write(error_unit, '(A,2ES15.7,A)') '   O(1,1) = (', real(Op%O(1,1)), aimag(Op%O(1,1)), ')'
          Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
       endif
    endif

    if (allocated(OP%g_t)) Op%g_t_alloc = .true.
    
    Call nsigma_single%make(1,1)
    noderank=0

    N = OP%N
    Allocate ( Op%E(N) )
    if (use_mpi_shm) then
       arrayshape2d=(/ Op%N,Op%N /)
       call allocate_shared_memory(Op%U,Op%win_U,noderank,arrayshape2d)
       if (noderank == 0) Op%U = cmplx(0.d0, 0.d0, kind(0.D0))
    else
       Allocate ( Op%U(N,N) )
       Op%U = cmplx(0.d0, 0.d0, kind(0.D0))
    endif
    Op%U_alloc = .true.
    Op%E = 0.d0
    
    If (Op%N > 1) then
       if (Op%diag) then
          do I=1,N
             Op%E(I)=DBLE(Op%O(I,I))
             if (noderank == 0) Op%U(I,I)=cmplx(1.d0,0.d0, kind(0.D0))
          enddo
          Op%N_non_zero = N
          ! FFA Why do we assume that Op%N_non_zero = N for a diagonal operator? 
       else
          Allocate (U(N,N), E(N), TMP(N, N))
          Call Diag(Op%O,U, E)  
          Np = 0
          Nz = 0
          do I = 1,N
             if ( abs(E(I)) > Zero ) then
                np = np + 1
                if (noderank == 0) Op%U(:, np) = U(:, i)
                Op%E(np)   = E(I)
             else
                if (noderank == 0) Op%U(:, N-nz) = U(:, i)
                Op%E(N-nz)   = E(I)
                nz = nz + 1
             endif
          enddo
          Op%N_non_zero = np
          ! Write(6,*) "Op_set", np,N
          if (noderank == 0) then
             TMP = Op%U ! that way we have the changes to the determinant due to the permutation
             Z = Det_C(TMP, N)
             if (abs(Z) < 1.D-30) then
                Write(error_unit, '(A,ES15.7)') ' Op_set: Eigenvector matrix has near-zero determinant: |det| = ', abs(Z)
                Write(error_unit, '(A)')        '   Cannot normalize to SU(N). Check operator matrix.'
                Call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
             endif
             ! Scale Op%U to be in SU(N) 
             DO I = 1, N
                Op%U(I,1) = Op%U(I, 1)/Z 
             ENDDO
          endif
          deallocate (U, E, TMP)
          ! Op%U,Op%E)
          ! Write(6,*) 'Calling diag 1'
       endif
    else
       Op%E(1)   = REAL(Op%O(1,1), kind(0.D0))
       if (noderank == 0) Op%U(1,1) = cmplx(1.d0, 0.d0 , kind(0.D0))
       Op%N_non_zero = 1
       Op%diag = .true.
    endif
#ifdef MPI
    if (use_mpi_shm) call MPI_WIN_FENCE(0, Op%win_U, ierr)
#endif
    if (allocated(Op%HY)) deallocate(Op%HY, Op%HT)
    Op%lowrank = op_lowrank_wanted(Op)
    if (Op%lowrank) call op_lowrank_setup(Op)

    select case(OP%type)
    case(1)
       if (use_mpi_shm) then
          Allocate(Op%E_exp(Op%N, -Op%type : Op%type))
          arrayshape=(/ Op%N,Op%N,3 /)
          call allocate_shared_memory(Op%M_exp,Op%win_M_exp,noderank,arrayshape)
       else
          ! The low-rank path never reads M_exp, which would cost Op%N**2 * 3 per operator.
          if (Op%lowrank) then
             Allocate(Op%E_exp(Op%N, -Op%type : Op%type), Op%M_exp(0, 0, 3))
          else
             Allocate(Op%E_exp(Op%N, -Op%type : Op%type), Op%M_exp(Op%N, Op%N, 3))
          endif
       endif
       Op%M_exp_alloc = .true.
       nsigma_single%t(1) = 1
       Do I=1,Op%type
          nsigma_single%f(1,1) = real(I,kind=kind(0.d0))
          do n = 1, Op%N
             Op%E_exp(n,I)  = cmplx(1.d0, 0.d0, kind(0.D0))
             Op%E_exp(n,-I) = cmplx(1.d0, 0.d0, kind(0.D0))
             if ( n <= Op%N_non_Zero) then
                !Op%E_exp(n,I) = exp(Op%g*Op%E(n)*Phi_st(I,1))
                Op%E_exp(n,I) = exp(Op%g*Op%E(n)*nsigma_single%Phi(1,1))
                Op%E_exp(n,-I) = 1.D0/Op%E_exp(n,I)
             endif
          enddo
          !call Op_exp(Op%g*Phi_st( I,1),Op,Op%M_exp(:,:,I))
          !call Op_exp(Op%g*Phi_st(-I,1),Op,Op%M_exp(:,:,-I))
          if (noderank == 0 .and. size(Op%M_exp, 1) > 0) then
             call Op_exp( Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:, I+2))
             call Op_exp(-Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:,-I+2))
          endif
#ifdef MPI
          if (use_mpi_shm) call MPI_WIN_FENCE(0, Op%win_M_exp, ierr)
#endif
       enddo
    case(2)
       if (use_mpi_shm) then
          Allocate(Op%E_exp(Op%N, -Op%type : Op%type))
          arrayshape=(/ Op%N,Op%N,5 /)
          call allocate_shared_memory(Op%M_exp,Op%win_M_exp,noderank,arrayshape)
       else
          ! The low-rank path never reads M_exp, which would cost Op%N**2 * 5 per operator.
          if (Op%lowrank) then
             Allocate(Op%E_exp(Op%N, -Op%type : Op%type), Op%M_exp(0, 0, 5))
          else
             Allocate(Op%E_exp(Op%N, -Op%type : Op%type), Op%M_exp(Op%N, Op%N, 5))
          endif
       endif
       Op%M_exp_alloc = .true.
       nsigma_single%t(1) = 2
       Do I=1,Op%type  ! Here Op%type = 2    
          nsigma_single%f(1,1) = real(I,kind=kind(0.d0))
          do n = 1, Op%N
             Op%E_exp(n ,I) = cmplx(1.d0, 0.d0, kind(0.D0))
             Op%E_exp(n,-I) = cmplx(1.d0, 0.d0, kind(0.D0))
             if ( n <= Op%N_non_Zero) then
                !Op%E_exp(n,I)  = exp(Op%g*Op%E(n)*Phi_st(I,2))
                Op%E_exp(n, I) = exp(Op%g*Op%E(n)*nsigma_single%Phi(1,1))
                Op%E_exp(n,-I) = 1.D0/Op%E_exp(n,I)
             endif
          enddo
          if (noderank == 0 .and. size(Op%M_exp, 1) > 0) then
             call Op_exp( Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,: ,I+3))
             call Op_exp(-Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:,-I+3))
          endif
#ifdef MPI
          if (use_mpi_shm) call MPI_WIN_FENCE(0, Op%win_M_exp, ierr)
#endif
          !call Op_exp(Op%g*Phi_st( I,2),Op,Op%M_exp(:,:,I))
          !call Op_exp(Op%g*Phi_st(-I,2),Op,Op%M_exp(:,:,-I))
       enddo
    case default
    end select

    Call nsigma_single%clear() 

  end subroutine Op_set

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> Calculate the exponentiated operator and returns the full matrix.
!>
!> @details
!> It is assumed that  the U and E arrays of the
!> operator are present.
!> @param[in]  g  Complex
!> @param[in]  Op Type(Operator)
!> @param[out]
!> Mat   Complex, Dimension(:,:)
!> * On output  \f$ M = e^{g O} \f$
!
!--------------------------------------------------------------------
  Subroutine Op_exp(g,Op,Mat)

    Implicit none 
    
    Type (Operator), Intent(IN)  :: Op
    Complex (Kind=Kind(0.d0)), Dimension(:,:), INTENT(OUT) :: Mat
    Complex (Kind=Kind(0.d0)), INTENT(IN) :: g
    Complex (Kind=Kind(0.d0)) :: Z, Z1, y, t
    Complex (Kind=Kind(0.d0)), allocatable, dimension(:,:) :: c
    
    Integer :: n, j, I, iters
    
    iters = Op%N
    Mat = cmplx(0.d0, 0.d0, kind(0.D0))
    if (Op%diag) then
       Do n = 1, iters
          Mat(n,n)=exp(g*Op%E(n))
       enddo
    else
       Allocate (c(iters, iters))
       c = 0.D0
       Do n = 1, iters
          Z = exp(g*Op%E(n))
          do J = 1, iters
             Z1 = Z*conjg(Op%U(J,n))
             do I = 1, iters
                ! This performs Kahan summation so as to improve precision.
                y = Z1 * Op%U(I, n) - c(I, J)
                t = Mat(I, J) + y
                c(I, J) = (t - Mat(I,J)) - y
                Mat(I, J) = t
                !  Mat(I, J) = Mat(I, J) + Z1 * Op%U(I, n)
             enddo
             ! Mat(1:iters, J) = Mat(1:iters, J) + Z1 * Op%U(1:iters, n)
          enddo
       enddo
       Deallocate(c)
    endif
  end subroutine Op_exp


!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief 
!> Out   With phi_loc = Phi_of_field(HS_Field, Op%type) the routine computes
!>       Op%Mat = Mat* Op ( exp( sign * phi_loc*g* P^T O T) )
!>       For  Op%type = 1,2  the  exponential is  stored. Otherwise  it is
!>       computed  on the fly
!>
!> @param[inout] Mat Complex Dimension(:,:)
!> * On exit Mat = Mat*Op ( exp(phi_loc*g* P^T O T) )
!> @param[in] Op Type(Operator)
!> * The Operator containing g and the sparse matrix P^T O P 
!> @param[in]  HS_Field Complex
!> * The field
!> @param[in] cop  Character
!> * cop = N,  Op = None
!> * cop = T,  Op = Transposed
!> * cop = C,  Op = Transposed + Complex conjugation
!> 
  
!--------------------------------------------------------------------
  subroutine Op_mmultL(Mat,Op,HS_Field,cop, nt, sign)
    Implicit none 
    Type (Operator)          , INTENT(IN)    :: Op
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (:,:)
    Integer                  , INTENT(IN)    :: sign
    Complex (Kind=Kind(0.d0)), INTENT(IN)    :: Hs_Field
    Character                , Intent(IN)    :: cop
    integer                  , intent(in)    :: nt
    
    ! Local 
    Integer :: I, N1, N2, sp
    Complex (Kind=Kind(0.d0)) :: ExpMat (Op%n,Op%n), g_loc
    Complex (Kind=Kind(0.d0)) :: phi_loc
    
    phi_loc = Phi_of_field(HS_field, op%type)

    N1=size(Mat,1)
    N2=size(Mat,2)

    g_loc = OP%g
    if (op%g_t_alloc)  g_loc = Op%g_t(nt)
    if ( abs(g_loc) < Eps_machine ) return
    if ( Op%lowrank ) then
       call op_lowrank_mmult(Mat, N1, N2, Op, HS_Field, cop, nt, sign, 'R')
       return
    endif

    if ( op%type < 3 ) then
       sp = sign*nint(Real(HS_Field))
       if ( Op%diag ) then
        if (Op%g_t_alloc) then
            do I=1,Op%N
                if ( cop == 'c' .or. cop =='C' ) then
                   call ZSCAL(N1,conjg(exp(sign*phi_loc*Op%g_t(nt)*Op%E(I))),Mat(1,Op%P(I)),1)
                else
                   call ZSCAL(N1,exp(sign*phi_loc*Op%g_t(nt)*Op%E(I)),Mat(1,Op%P(I)),1)
                endif
            enddo
        else
            do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N1,conjg(Op%E_exp(I,sp)),Mat(1,Op%P(I)),1)
             else
                call ZSCAL(N1,Op%E_exp(I,sp),Mat(1,Op%P(I)),1)
             endif
          enddo
        endif
       else
          if (Op%g_t_alloc) then
            call Op_exp(sign*phi_loc*Op%g_t(nt),Op,expmat)
            call ZSLGEMM('r',cop,Op%N,N1,N2,expmat,Op%P,Mat)
          else
            call ZSLGEMM('r',cop,Op%N,N1,N2,Op%M_exp(:,:,sp+op%type+1),Op%P,Mat)
          endif
       endif
    else
       if ( Op%diag ) then
          do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N1,conjg(exp(sign*phi_loc*g_loc*Op%E(I))),Mat(1,Op%P(I)),1)
             else
                call ZSCAL(N1,exp(sign*phi_loc*g_loc*Op%E(I)),Mat(1,Op%P(I)),1)
             endif
          enddo
       else
          call Op_exp(sign*g_loc*phi_loc,Op,expmat)
          call ZSLGEMM('r',cop,Op%N,N1,N2,expmat,Op%P,Mat)
       endif
    endif
    
  end subroutine Op_mmultL

  
!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief 
!> Out   With phi_loc = Phi_of_field(HS_Field, Op%type) the routine computes
!>       Op%Mat =  Op ( exp( phi_loc*g* P^T O T) ) * Mat
!>       For  Op%type = 1,2  the  exponential is  stored. Otherwise  it is
!>       computed  on the fly
!
!> @param[inout] Mat Complex Dimension(:,:)
!> * On exit Mat = Op ( exp(phi_loc*g* P^T O T) )* Mat
!> @param[in] Op Type(Operator)
!> * The Operator containing g and the sparse matrix P^T O P 
!> @param[in] Hs_Field Complex
!> * The field
!> @param[in] cop  Character
!> * cop = N,  Op = None
!> * cop = T,  Op = Transposed
!> * cop = C,  Op = Transposed + Complex conjugation
!> 
!--------------------------------------------------------------------
  subroutine Op_mmultR(Mat,Op,HS_Field,cop,nt)
    Implicit none
    Type (Operator)          , INTENT(IN )   :: Op
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (:,:)
    Complex (Kind=Kind(0.d0)), INTENT(IN )   :: Hs_Field
    Character                , Intent(IN)    :: cop
    integer                  , intent(in)    :: nt
    
    ! Local 
    Integer :: I, N1, N2, sp
    Complex (Kind=Kind(0.d0)) :: ExpMat (Op%n,Op%n), g_loc
    Complex (Kind=Kind(0.d0)) :: phi_loc
    
    phi_loc = Phi_of_field(Hs_Field, op%type)

    N1=size(Mat,1)
    N2=size(Mat,2)
    
    
    ! quick return if possible
    g_loc = Op%g
    if (op%g_t_alloc) g_loc = Op%g_t(nt)
    if ( abs(g_loc) < Eps_machine ) return
    if ( Op%lowrank ) then
       call op_lowrank_mmult(Mat, N1, N2, Op, HS_Field, cop, nt, 1, 'L')
       return
    endif

    if ( op%type < 3 ) then
       sp = nint(Real(Hs_Field))
       if ( Op%diag ) then
        if (op%g_t_alloc) then
            do I=1,Op%N
                if ( cop == 'c' .or. cop =='C' ) then
                    call ZSCAL(N2,conjg(exp(phi_loc*Op%g_t(nt)*Op%E(I))),Mat(Op%P(I),1),N1)
                else
                    call ZSCAL(N2,exp(phi_loc*Op%g_t(nt)*Op%E(I)),Mat(Op%P(I),1),N1)
                endif
            enddo
        else
            do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N2,conjg(Op%E_exp(I,sp)),Mat(Op%P(I),1),N1)
             else
                call ZSCAL(N2,Op%E_exp(I,sp),Mat(Op%P(I),1),N1)
             endif
            
            enddo
        endif
       else
        if (op%g_t_alloc) then
            call Op_exp(Op%g_t(nt)*phi_loc,Op,expmat)
            call ZSLGEMM('L',cop,Op%N,N1,N2,expmat,Op%P,Mat)
        else
          call ZSLGEMM('L',cop,Op%N,N1,N2,Op%M_exp(:,:,sp+op%type+1),Op%P,Mat)
        endif
       endif
    else
       if ( Op%diag ) then
          do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N2,conjg(exp(phi_loc*g_loc*Op%E(I))),Mat(Op%P(I),1),N1)
             else
                call ZSCAL(N2,exp(phi_loc*g_loc*Op%E(I)),Mat(Op%P(I),1),N1)
             endif
          enddo
       else
          call Op_exp(g_loc*phi_loc,Op,expmat)
          call ZSLGEMM('L',cop,Op%N,N1,N2,expmat,Op%P,Mat)
       endif
    endif
  end subroutine Op_mmultR

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief 
!> Wrapup the Green function
!>
!> @param[inout] Mat(Ndim,Ndim)  Complex
!> \verbatim
!>  N_type = 1, Mat = exp(Op%g*nsigma%phi(1,1)*Op%E)*(Op%U^{dagger}) * Mat * Op%U*exp(-Op%g*nsigma%phi(1,1)*Op%E)
!>  N_type = 2, Mat = Op%U * Mat * (Op%U^{dagger})
!> \endverbatim
!>  **Do not** mix up  \p N_type  and    \p OP\%type 
!> @param[in] Op Type(Operator)
!> \verbatim
!>  The operator containing g, U, P
!> \endverbatim
!> @param[in] HS_Field Complex
!> \verbatim
!>  The field
!> \endverbatim
!> @param[in] Ndim Integer
!> 
!--------------------------------------------------------------------
  Subroutine Op_Wrapup(Mat,Op,HS_Field,Ndim,N_Type, nt)

    Implicit none 

    Integer :: Ndim
    Type (Operator)           , INTENT(IN )   :: Op
    Complex (Kind=Kind(0.d0)) , INTENT(INOUT) :: Mat (Ndim,Ndim)
    Complex (Kind=Kind(0.d0)) , INTENT(IN )   :: HS_Field
    Integer                   , INTENT(IN)    :: N_Type
    Integer                   , INTENT(IN)    :: nt

    ! Local 
    Complex (Kind=Kind(0.d0)) :: VH1(Op%N,Op%N)
    Integer :: I,sp
    Complex (kind=kind(0.d0)) :: g_loc
    Complex (Kind=Kind(0.d0)) :: phi_loc

    if ( Op%lowrank ) then
       call op_lowrank_wrap(Mat, Ndim, Op, HS_Field, N_Type, nt, .true.)
       return
    endif
    phi_loc = Phi_of_field(HS_Field, op%type)

    if ( op%type < 3 ) then
       sp = nint(Real(HS_Field))
       If (N_type == 1) then
          if(Op%diag) then
            if (op%g_t_alloc) then
                do I=1,Op%N
                    call ZSCAL(Ndim,exp( phi_loc*Op%g_t(nt)*Op%E(I)),Mat(Op%P(I),1),Ndim)
                enddo
                do I=1,Op%N
                    call ZSCAL(Ndim,exp(-phi_loc*Op%g_t(nt)*Op%E(I)),Mat(1,Op%P(I)),1)
                enddo
            else
             do I=1,Op%N
                call ZSCAL(Ndim,Op%E_Exp(I,sp),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,Op%E_Exp(I,-sp),Mat(1,Op%P(I)),1)
             enddo
            endif
          else
            if (op%g_t_alloc) then
                Do i = 1,Op%N
                   VH1(:,i)=Op%U(:,i)*exp(-phi_loc*Op%g_t(nt)*Op%E(I))
                Enddo
                call ZSLGEMM('r','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
                Do i = 1,Op%N
                   VH1(:,i)=exp(phi_loc*Op%g_t(nt)*Op%E(I))*conjg(Op%U(:,i))
                Enddo
                call ZSLGEMM('l','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)

            else
                 Do i = 1,Op%N
                    VH1(:,i)=Op%U(:,i)*Op%E_Exp(I,-sp)
                 Enddo
                 call ZSLGEMM('r','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
                 Do i = 1,Op%N
                    VH1(:,i)=Op%E_Exp(I, sp)*conjg(Op%U(:,i))
                 Enddo
                 call ZSLGEMM('l','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
            endif
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('l','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('r','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    else
       If (N_type == 1) then
        g_loc = Op%g
        if (op%g_t_alloc) g_loc = Op%g_t(nt)
          if(Op%diag) then
             do I=1,Op%N
                call ZSCAL(Ndim,exp( phi_loc*g_loc*Op%E(I)),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,exp(-phi_loc*g_loc*Op%E(I)),Mat(1,Op%P(I)),1)
             enddo
          else
             Do i = 1,Op%N
                VH1(:,i)=Op%U(:,i)*exp(-phi_loc*g_loc*Op%E(I))
             Enddo
             call ZSLGEMM('r','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             Do i = 1,Op%N
                VH1(:,i)=exp(phi_loc*g_loc*Op%E(I))*conjg(Op%U(:,i))
             Enddo
             call ZSLGEMM('l','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('l','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('r','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    endif
  end Subroutine Op_Wrapup

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief 
!> Wrapup the Green function
!>
!> @param[inout] Mat(Ndim,Ndim)  Complex
!> \verbatim
!>  N_type = 1, Mat = Op%U*exp(-Op%g*nsigma%phi(1,1)*Op%E)*Mat*exp(Op%g*nsigma%phi(1,1)*Op%E)*(Op%U^{dagger})
!>  N_type = 2, Mat = (Op%U^{dagger}) * Mat * Op%U
!> \endverbatim
!>  **Do not** mix up  \p N_type  and    \p OP\%type 
!> @param[in] Op Type(Operator)
!> \verbatim
!>  The operator containing g, U, P
!> \endverbatim
!> @param[in] HS_Field Complex
!> \verbatim
!>  The field
!> \endverbatim
!> @param[in] Ndim Integer
!> 
!--------------------------------------------------------------------

  Subroutine Op_Wrapdo(Mat,Op,HS_Field,Ndim,N_Type,nt)
    Implicit none 
    
    Integer :: Ndim
    Type (Operator) , INTENT(IN)   :: Op
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Complex   (Kind=Kind(0.d0)), INTENT(IN )   :: HS_Field
    Integer, INTENT(IN) :: N_Type, nt

    ! Local 
    Integer :: n, i, sp
    Complex (Kind = Kind(0.D0)) :: VH1(Op%N,OP%N)
    Complex (Kind = Kind(0.D0)) :: g_loc
    Complex (Kind=Kind(0.d0)) :: phi_loc

    if ( Op%lowrank ) then
       call op_lowrank_wrap(Mat, Ndim, Op, HS_Field, N_Type, nt, .false.)
       return
    endif
    phi_loc = Phi_of_field(HS_Field, op%type)

    if ( op%type < 3 ) then
       sp = nint(Real(HS_Field))
       If (N_type == 1) then
          if(Op%diag) then
             if (op%g_t_alloc) then
                do I=1,Op%N
                   call ZSCAL(Ndim,exp(-phi_loc*Op%g_t(nt)*Op%E(I)),Mat(Op%P(I),1),Ndim)
                enddo
                do I=1,Op%N
                   call ZSCAL(Ndim,exp(phi_loc*Op%g_t(nt)*Op%E(I)),Mat(1,Op%P(I)),1)
                enddo
             else
                do I=1,Op%N
                   call ZSCAL(Ndim,Op%E_Exp(I,-sp),Mat(Op%P(I),1),Ndim)
                enddo
                do I=1,Op%N
                   call ZSCAL(Ndim,Op%E_Exp(I, sp),Mat(1,Op%P(I)),1)
                enddo
             endif
          else
             if (op%g_t_alloc) then
                Do n = 1,Op%N
                   VH1(:,n)=Op%U(:,n)*exp(-phi_loc*Op%g_t(nt)*Op%E(n))
                Enddo
                call ZSLGEMM('l','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
                Do n = 1,Op%N
                   VH1(:,n)=exp(phi_loc*Op%g_t(nt)*Op%E(n))*conjg(Op%U(:,n))
                Enddo
                call ZSLGEMM('r','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             else
                Do n = 1,Op%N
                   VH1(:,n)=Op%U(:,n)*Op%E_Exp(n,-sp)
                Enddo
                call ZSLGEMM('l','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
                Do n = 1,Op%N
                   VH1(:,n)=Op%E_Exp(n, sp)*conjg(Op%U(:,n))
                Enddo
                call ZSLGEMM('r','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             endif
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('r','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('l','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    else
       If (N_type == 1) then
          g_loc = Op%g
          if (op%g_t_alloc) g_loc = Op%g_t(nt)
          if(Op%diag) then
             do I=1,Op%N
                call ZSCAL(Ndim,exp(-phi_loc*g_loc*Op%E(I)),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,exp( phi_loc*g_loc*Op%E(I)),Mat(1,Op%P(I)),1)
             enddo
          else
             Do n = 1,Op%N
                VH1(:,n)=Op%U(:,n)*exp(-phi_loc*g_loc*Op%E(n))
             Enddo
             call ZSLGEMM('l','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             Do n = 1,Op%N
                VH1(:,n)=exp(phi_loc*g_loc*Op%E(n))*conjg(Op%U(:,n))
             Enddo
             call ZSLGEMM('r','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('r','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('l','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    endif
  end Subroutine Op_Wrapdo

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Apply the same conjugation Op_Wrapup/Op_Wrapdo apply to a Green's function,
!> but to the two panels of a factored one.
!>
!> @details
!> A delayed (rank-k) update holds the Green's function as `G = G_stale + X*Y^T`.
!> Op_Wrapup and Op_Wrapdo act as `G -> L*G*R` with L and R the identity outside
!> `Op%P(1:Op%N)`, hence
!>
!>     L*(G_stale + X*Y^T)*R = (L*G_stale*R) + (L*X)*(R^T*Y)^T.
!>
!> This routine applies L to X and R^T to Y, costing O(Op%N**2*ncols).
!>
!> Converting the wrap's right-multiply into a left-multiply on Y:
!>
!>     wrap does ZSLGEMM('r','n',A)  ->  here ZSLGEMM('l','T',A)
!>     wrap does ZSLGEMM('r','T',A)  ->  here ZSLGEMM('l','N',A)
!>     wrap does ZSLGEMM('r','c',A)  ->  here ZSLGEMM('l','N',conjg(A))
!>
!> Note: ZSLGEMM has no conjugate-without-transpose and therefore A must be
!> conjugated into a temporary first. For most typical models Op%N is small,
!> and hence the cost of this operation is negligible.
!>
!> This subroutine mirrors the sixteen branches of Op_Wrapup/Op_Wrapdo and hence
!> care must be taken to keep the logic synchronised between them. Within
!> the testsuite, test 37-delayed-wrap forms an explicit check to catch any
!> such drift.
!>
!> @param[inout] Xpan(Ndim,ncols), Ypan(Ndim,ncols) The panels.
!> @param[in] updo 'u' to mirror Op_Wrapup, 'd' to mirror Op_Wrapdo.
!--------------------------------------------------------------------
  Subroutine Op_Wrap_panels(Xpan, Ypan, Op, HS_Field, Ndim, ncols, N_Type, nt, updo)
    use iso_fortran_env, only: error_unit
    Implicit none

    Integer                   , INTENT(IN)    :: Ndim, ncols, N_Type, nt
    Type (Operator)           , INTENT(IN)    :: Op
    Complex (Kind=Kind(0.d0)) , INTENT(INOUT) :: Xpan(Ndim,ncols), Ypan(Ndim,ncols)
    Complex (Kind=Kind(0.d0)) , INTENT(IN)    :: HS_Field
    Character (Len=1)         , INTENT(IN)    :: updo

    ! Local
    Complex (Kind=Kind(0.d0)) :: VH1(Op%N,Op%N), VH2(Op%N,Op%N)
    Complex (Kind=Kind(0.d0)) :: g_loc, phi_loc
    Integer :: I, sp
    Logical :: up

    if (ncols <= 0) return

    if (updo == 'u' .or. updo == 'U') then
       up = .true.
    elseif (updo == 'd' .or. updo == 'D') then
       up = .false.
    else
       write(error_unit,*) 'Op_Wrap_panels: updo must be one of [u,U,d,D], got ', updo
       Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
       return
    endif

    if ( Op%lowrank ) then
       call op_lowrank_panels(Xpan, Ypan, Ndim, ncols, Op, HS_Field, N_Type, nt, up)
       return
    endif
    phi_loc = Phi_of_field(HS_Field, op%type)

    g_loc = Op%g
    if (op%g_t_alloc) g_loc = Op%g_t(nt)

    if ( op%type < 3 ) then
       sp = nint(Real(HS_Field))
       If (N_type == 1) then
          if (Op%diag) then
             ! Row scale on X, the reciprocal row scale on Y. The table entries
             ! are read rather than recomputed, so the panel product tracks the
             ! matrix to the last bit.
             if (op%g_t_alloc) then
                do I = 1,Op%N
                   if (up) then
                      call ZSCAL(ncols, exp( phi_loc*Op%g_t(nt)*Op%E(I)), Xpan(Op%P(I),1), Ndim)
                      call ZSCAL(ncols, exp(-phi_loc*Op%g_t(nt)*Op%E(I)), Ypan(Op%P(I),1), Ndim)
                   else
                      call ZSCAL(ncols, exp(-phi_loc*Op%g_t(nt)*Op%E(I)), Xpan(Op%P(I),1), Ndim)
                      call ZSCAL(ncols, exp( phi_loc*Op%g_t(nt)*Op%E(I)), Ypan(Op%P(I),1), Ndim)
                   endif
                enddo
             else
                do I = 1,Op%N
                   if (up) then
                      call ZSCAL(ncols, Op%E_Exp(I, sp), Xpan(Op%P(I),1), Ndim)
                      call ZSCAL(ncols, Op%E_Exp(I,-sp), Ypan(Op%P(I),1), Ndim)
                   else
                      call ZSCAL(ncols, Op%E_Exp(I,-sp), Xpan(Op%P(I),1), Ndim)
                      call ZSCAL(ncols, Op%E_Exp(I, sp), Ypan(Op%P(I),1), Ndim)
                   endif
                enddo
             endif
          else
             if (op%g_t_alloc) then
                Do I = 1,Op%N
                   VH1(:,I) = Op%U(:,I)*exp(-phi_loc*Op%g_t(nt)*Op%E(I))
                Enddo
                Do I = 1,Op%N
                   VH2(:,I) = exp(phi_loc*Op%g_t(nt)*Op%E(I))*conjg(Op%U(:,I))
                Enddo
             else
                Do I = 1,Op%N
                   VH1(:,I) = Op%U(:,I)*Op%E_Exp(I,-sp)
                Enddo
                Do I = 1,Op%N
                   VH2(:,I) = Op%E_Exp(I, sp)*conjg(Op%U(:,I))
                Enddo
             endif
             if (up) then
                ! Wrapup: R from ('r','n',VH1), L from ('l','T',VH2).
                call ZSLGEMM('l','T',Op%n,Ndim,ncols,VH1,Op%P,Ypan)
                call ZSLGEMM('l','T',Op%n,Ndim,ncols,VH2,Op%P,Xpan)
             else
                ! Wrapdo: L from ('l','n',VH1), R from ('r','T',VH2).
                call ZSLGEMM('l','n',Op%n,Ndim,ncols,VH1,Op%P,Xpan)
                call ZSLGEMM('l','N',Op%n,Ndim,ncols,VH2,Op%P,Ypan)
             endif
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          VH1 = conjg(Op%U)
          if (up) then
             ! Wrapup: L from ('l','n',U), R from ('r','c',U).
             call ZSLGEMM('l','n',Op%n,Ndim,ncols,Op%U,Op%P,Xpan)
             call ZSLGEMM('l','N',Op%n,Ndim,ncols,VH1,Op%P,Ypan)
          else
             ! Wrapdo: R from ('r','n',U), L from ('l','c',U).
             call ZSLGEMM('l','c',Op%n,Ndim,ncols,Op%U,Op%P,Xpan)
             call ZSLGEMM('l','T',Op%n,Ndim,ncols,Op%U,Op%P,Ypan)
          endif
       endif
    else
       If (N_type == 1) then
          if (Op%diag) then
             do I = 1,Op%N
                if (up) then
                   call ZSCAL(ncols, exp( phi_loc*g_loc*Op%E(I)), Xpan(Op%P(I),1), Ndim)
                   call ZSCAL(ncols, exp(-phi_loc*g_loc*Op%E(I)), Ypan(Op%P(I),1), Ndim)
                else
                   call ZSCAL(ncols, exp(-phi_loc*g_loc*Op%E(I)), Xpan(Op%P(I),1), Ndim)
                   call ZSCAL(ncols, exp( phi_loc*g_loc*Op%E(I)), Ypan(Op%P(I),1), Ndim)
                endif
             enddo
          else
             Do I = 1,Op%N
                VH1(:,I) = Op%U(:,I)*exp(-phi_loc*g_loc*Op%E(I))
             Enddo
             Do I = 1,Op%N
                VH2(:,I) = exp(phi_loc*g_loc*Op%E(I))*conjg(Op%U(:,I))
             Enddo
             if (up) then
                call ZSLGEMM('l','T',Op%n,Ndim,ncols,VH1,Op%P,Ypan)
                call ZSLGEMM('l','T',Op%n,Ndim,ncols,VH2,Op%P,Xpan)
             else
                call ZSLGEMM('l','n',Op%n,Ndim,ncols,VH1,Op%P,Xpan)
                call ZSLGEMM('l','N',Op%n,Ndim,ncols,VH2,Op%P,Ypan)
             endif
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          VH1 = conjg(Op%U)
          if (up) then
             call ZSLGEMM('l','n',Op%n,Ndim,ncols,Op%U,Op%P,Xpan)
             call ZSLGEMM('l','N',Op%n,Ndim,ncols,VH1,Op%P,Ypan)
          else
             call ZSLGEMM('l','c',Op%n,Ndim,ncols,Op%U,Op%P,Xpan)
             call ZSLGEMM('l','T',Op%n,Ndim,ncols,Op%U,Op%P,Ypan)
          endif
       endif
    endif

  end Subroutine Op_Wrap_panels

!--------------------------------------------------------------------
!> @brief
!> Whether Op_set gives the operator the low-rank path: a non-diagonal HS vertex
!> (type >= 1) of rank under half its width, where it is cheaper than the dense
!> Op%N x Op%N products. ALF_LOWRANK_VERTEX (op_lowrank_env_mode) can turn it off.
!--------------------------------------------------------------------
  logical function op_lowrank_wanted(Op)
    Implicit none
    Type (Operator), intent(in) :: Op
    Character (len=32) :: text
    Integer :: length, status

    if (Op_lowrank_mode < 0) then
       Op_lowrank_mode = 1
       call get_environment_variable("ALF_LOWRANK_VERTEX", text, length, status)
       if (status /= 1 .and. length > 0) then
          if (status == 0) Op_lowrank_mode = op_lowrank_env_mode(text(1:length))
          if (status /= 0 .or. Op_lowrank_mode < 0) then
             Write(error_unit, '(A)') ' ALF_LOWRANK_VERTEX must be one of 0/off/false/no or 1/on/true/yes, got "' &
                  & // trim(text) // '"'
             Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
          endif
       endif
    endif
    op_lowrank_wanted = Op_lowrank_mode == 1 .and. Op%type >= 1 .and. .not. Op%diag &
         &              .and. Op%N > 1 .and. Op%N_non_zero >= 1 .and. 2*Op%N_non_zero < Op%N
  end function op_lowrank_wanted

!> ALF_LOWRANK_VERTEX's value as a mode: 0 off (0, off, false, no), 1 on (1, on, true,
!> yes), -1 for anything else; case and surrounding blanks are ignored.
  integer function op_lowrank_env_mode(text)
    Implicit none
    Character (len=*), intent(in) :: text
    Character (len=len(text)) :: low
    Integer :: i, c

    low = adjustl(text)
    do i = 1, len(low)
       c = iachar(low(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) low(i:i) = achar(c + 32)
    enddo
    select case (trim(low))
    case ('0', 'off', 'false', 'no')
       op_lowrank_env_mode = 0
    case ('1', 'on', 'true', 'yes')
       op_lowrank_env_mode = 1
    case default
       op_lowrank_env_mode = -1
    end select
  end function op_lowrank_env_mode

!> Householder QR of V = Op%U(:,1:N_non_zero) in compact WY form, Q = I - HY HT HY^dag,
!> and the start of P when it is a contiguous ascending run (0 otherwise). U is left
!> whole.
  subroutine op_lowrank_setup(Op)
    Implicit none
    Type (Operator), intent(inout) :: Op
    Complex (Kind=Kind(0.d0)), allocatable :: work(:), tau(:)
    Complex (Kind=Kind(0.d0)) :: wq(1)
    Integer :: n, nz, info, i

    n  = Op%N
    nz = Op%N_non_zero
    allocate(Op%HY(n, nz), Op%HT(nz, nz), tau(nz))
    Op%HY = Op%U(:, 1:nz)
    call ZGEQRF(n, nz, Op%HY, n, tau, wq, -1, info)
    allocate(work(max(1, nint(real(wq(1))))))
    call ZGEQRF(n, nz, Op%HY, n, tau, work, size(work), info)
    if (info /= 0) then
       Write(error_unit, '(A,I0)') ' op_lowrank_setup: ZGEQRF failed, info = ', info
       Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
    endif
    Op%HT = cmplx(0.d0, 0.d0, kind(0.d0))
    call ZLARFT('F', 'C', n, nz, Op%HY, n, tau, Op%HT, nz)
    ! ZLARFT reads only the reflectors below the diagonal; make HY the explicit
    ! unit lower trapezoidal Y, as the products below use it whole.
    do i = 1, nz
       Op%HY(1:i-1, i) = cmplx(0.d0, 0.d0, kind(0.d0))
       Op%HY(i, i)     = cmplx(1.d0, 0.d0, kind(0.d0))
    enddo
    Op%P0 = Op%P(1)
    do i = 2, n
       if (Op%P(i) /= Op%P(1) + i - 1) Op%P0 = 0
    enddo
  end subroutine op_lowrank_setup

!> d(n) = e^{sgn g phi E(n)}, n <= N_non_zero, read as the dense path reads it.
  subroutine op_lowrank_factors(Op, HS_Field, nt, sgn, d)
    Implicit none
    Type (Operator), intent(in) :: Op
    Complex (Kind=Kind(0.d0)), intent(in) :: HS_Field
    Integer, intent(in) :: nt, sgn
    Complex (Kind=Kind(0.d0)), intent(out) :: d(:)
    Complex (Kind=Kind(0.d0)) :: g_loc, phi_loc
    Integer :: n, sp

    if (Op%type < 3 .and. .not. Op%g_t_alloc) then
       sp = sgn*nint(real(HS_Field))
       do n = 1, Op%N_non_zero
          d(n) = Op%E_exp(n, sp)
       enddo
    else
       phi_loc = Phi_of_field(HS_Field, Op%type)
       g_loc = Op%g
       if (Op%g_t_alloc) g_loc = Op%g_t(nt)
       do n = 1, Op%N_non_zero
          d(n) = exp(sgn*phi_loc*g_loc*Op%E(n))
       enddo
    endif
  end subroutine op_lowrank_factors

!--------------------------------------------------------------------
!> @brief
!> op(Q) = I + B S B^dag for which = 'N' (Q), 'C' (Q^dag), 'T' (Q^T), 'J' (conjg(Q)):
!> B = HY or conjg(HY), S = -HT, -HT^dag, -HT^T or -conjg(HT).
!--------------------------------------------------------------------
  subroutine op_lowrank_qform(Op, which, B, S)
    Implicit none
    Type (Operator), intent(in) :: Op
    Character, intent(in) :: which
    Complex (Kind=Kind(0.d0)), intent(out) :: B(:,:), S(:,:)

    select case (which)
    case ('N')
       B = Op%HY
       S = -Op%HT
    case ('C')
       B = Op%HY
       S = -conjg(transpose(Op%HT))
    case ('T')
       B = conjg(Op%HY)
       S = -transpose(Op%HT)
    case ('J')
       B = conjg(Op%HY)
       S = -conjg(Op%HT)
    end select
  end subroutine op_lowrank_qform

!--------------------------------------------------------------------
!> @brief
!> Rows P of Mat(M1,M2): X <- X + B S (B^dag X), X = Mat(P,:), B (N x nz), S (nz x nz).
!> In place when P is a contiguous run; otherwise X is gathered and scattered.
!--------------------------------------------------------------------
  subroutine op_lowrank_rows(Op, Mat, M1, M2, B, S)
    Implicit none
    Type (Operator), intent(in) :: Op
    Integer, intent(in) :: M1, M2
    Complex (Kind=Kind(0.d0)), intent(inout) :: Mat(M1, M2)
    Complex (Kind=Kind(0.d0)), intent(in) :: B(:,:), S(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: X(:,:)
    Integer :: i

    if (Op%P0 > 0) then
       call op_lowrank_rows_kernel(Op%N, M2, Op%N_non_zero, B, S, Mat(Op%P0, 1), M1)
    else
       allocate(X(Op%N, M2))
       do i = 1, Op%N
          X(i, :) = Mat(Op%P(i), :)
       enddo
       call op_lowrank_rows_kernel(Op%N, M2, Op%N_non_zero, B, S, X, Op%N)
       do i = 1, Op%N
          Mat(Op%P(i), :) = X(i, :)
       enddo
    endif
  end subroutine op_lowrank_rows

!> Columns P of Mat(M1,M2): X <- X + (X B) S B^dag, X = Mat(:,P); in place as above.
  subroutine op_lowrank_cols(Op, Mat, M1, M2, B, S)
    Implicit none
    Type (Operator), intent(in) :: Op
    Integer, intent(in) :: M1, M2
    Complex (Kind=Kind(0.d0)), intent(inout) :: Mat(M1, M2)
    Complex (Kind=Kind(0.d0)), intent(in) :: B(:,:), S(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: X(:,:)
    Integer :: i

    if (Op%P0 > 0) then
       call op_lowrank_cols_kernel(M1, Op%N, Op%N_non_zero, B, S, Mat(1, Op%P0), M1)
    else
       allocate(X(M1, Op%N))
       do i = 1, Op%N
          X(:, i) = Mat(:, Op%P(i))
       enddo
       call op_lowrank_cols_kernel(M1, Op%N, Op%N_non_zero, B, S, X, M1)
       do i = 1, Op%N
          Mat(:, Op%P(i)) = X(:, i)
       enddo
    endif
  end subroutine op_lowrank_cols

!> X(1:n,j) += B S B^dag X(1:n,j) for each column j <= m2 of X(ldx,*). Loops, not
!> ZGEMM: with nz of one or two the products are thin, and a BLAS call's overhead
!> outweighs its arithmetic; columns are unit stride.
  subroutine op_lowrank_rows_kernel(n, m2, nz, B, S, X, ldx)
    Implicit none
    Integer, intent(in) :: n, m2, nz, ldx
    Complex (Kind=Kind(0.d0)), intent(in) :: B(n, nz), S(nz, nz)
    Complex (Kind=Kind(0.d0)), intent(inout) :: X(ldx, *)
    Complex (Kind=Kind(0.d0)) :: w(nz), w2(nz), z
    Integer :: i, j, k

    do j = 1, m2
       do k = 1, nz
          z = cmplx(0.d0, 0.d0, kind(0.d0))
          do i = 1, n
             z = z + conjg(B(i, k))*X(i, j)
          enddo
          w(k) = z
       enddo
       do k = 1, nz
          w2(k) = sum(S(k, :)*w)
       enddo
       do k = 1, nz
          z = w2(k)
          do i = 1, n
             X(i, j) = X(i, j) + B(i, k)*z
          enddo
       enddo
    enddo
  end subroutine op_lowrank_rows_kernel

!> X(1:m1,1:n) += X B S B^dag on X(ldx,*), as op_lowrank_rows_kernel but from the right.
  subroutine op_lowrank_cols_kernel(m1, n, nz, B, S, X, ldx)
    Implicit none
    Integer, intent(in) :: m1, n, nz, ldx
    Complex (Kind=Kind(0.d0)), intent(in) :: B(n, nz), S(nz, nz)
    Complex (Kind=Kind(0.d0)), intent(inout) :: X(ldx, *)
    Complex (Kind=Kind(0.d0)) :: W(m1, nz), W2(m1, nz), z
    Integer :: j, k, l

    W = cmplx(0.d0, 0.d0, kind(0.d0))
    do j = 1, n
       do k = 1, nz
          z = B(j, k)
          W(:, k) = W(:, k) + z*X(1:m1, j)
       enddo
    enddo
    W2 = cmplx(0.d0, 0.d0, kind(0.d0))
    do k = 1, nz
       do l = 1, nz
          W2(:, k) = W2(:, k) + W(:, l)*S(l, k)
       enddo
    enddo
    do j = 1, n
       do k = 1, nz
          z = conjg(B(j, k))
          X(1:m1, j) = X(1:m1, j) + z*W2(:, k)
       enddo
    enddo
  end subroutine op_lowrank_cols_kernel

!> Scales rows (rows = .true.) or columns P(1:N_non_zero) of Mat(M1,M2) by d.
  subroutine op_lowrank_scale(Op, Mat, M1, M2, d, rows)
    Implicit none
    Type (Operator), intent(in) :: Op
    Integer, intent(in) :: M1, M2
    Complex (Kind=Kind(0.d0)), intent(inout) :: Mat(M1, M2)
    Complex (Kind=Kind(0.d0)), intent(in) :: d(:)
    Logical, intent(in) :: rows
    Integer :: i

    do i = 1, Op%N_non_zero
       if (rows) then
          call ZSCAL(M2, d(i), Mat(Op%P(i), 1), M1)
       else
          call ZSCAL(M1, d(i), Mat(1, Op%P(i)), 1)
       endif
    enddo
  end subroutine op_lowrank_scale

!--------------------------------------------------------------------
!> @brief
!> Mat(P,:) = op(M) Mat(P,:) (side 'L') or Mat(:,P) = Mat(:,P) op(M) (side 'R'), with
!> M = e^{sgn g phi O} = 1 + V diag(c) V^dag, c = e^{sgn g phi E} - 1: O(Op%N N_non_zero)
!> per row or column instead of O(Op%N**2). op(M) - 1 = A diag(c) A^dag with A = V
!> (cop N, C; c conjugated for C) or conjg(V) (cop T).
!--------------------------------------------------------------------
  subroutine op_lowrank_mmult(Mat, M1, M2, Op, HS_Field, cop, nt, sgn, side)
    Implicit none
    Integer, intent(in) :: M1, M2
    Complex (Kind=Kind(0.d0)), intent(inout) :: Mat(M1, M2)
    Type (Operator), intent(in) :: Op
    Complex (Kind=Kind(0.d0)), intent(in) :: HS_Field
    Character, intent(in) :: cop, side
    Integer, intent(in) :: nt, sgn
    Complex (Kind=Kind(0.d0)) :: A(Op%N, Op%N_non_zero), S(Op%N_non_zero, Op%N_non_zero)
    Complex (Kind=Kind(0.d0)) :: c(Op%N_non_zero)
    Integer :: i

    call op_lowrank_factors(Op, HS_Field, nt, sgn, c)
    c = c - cmplx(1.d0, 0.d0, kind(0.d0))
    if (cop == 'C' .or. cop == 'c') c = conjg(c)
    A = Op%U(:, 1:Op%N_non_zero)
    if (cop == 'T' .or. cop == 't') A = conjg(A)
    S = cmplx(0.d0, 0.d0, kind(0.d0))
    do i = 1, Op%N_non_zero
       S(i, i) = c(i)
    enddo
    if (side == 'L') then
       call op_lowrank_rows(Op, Mat, M1, M2, A, S)
    else
       call op_lowrank_cols(Op, Mat, M1, M2, A, S)
    endif
  end subroutine op_lowrank_mmult

!--------------------------------------------------------------------
!> @brief
!> Op_Wrapup (up) / Op_Wrapdo with Q in place of U, D = diag(e^{g phi E}):
!>   up,   N_type 1: D Q^dag G Q D^-1        up,   N_type 2: Q G Q^dag
!>   down, N_type 1: Q D^-1 G D Q^dag        down, N_type 2: Q^dag G Q
!--------------------------------------------------------------------
  subroutine op_lowrank_wrap(Mat, Ndim, Op, HS_Field, N_Type, nt, up)
    Implicit none
    Integer, intent(in) :: Ndim
    Complex (Kind=Kind(0.d0)), intent(inout) :: Mat(Ndim, Ndim)
    Type (Operator), intent(in) :: Op
    Complex (Kind=Kind(0.d0)), intent(in) :: HS_Field
    Integer, intent(in) :: N_Type, nt
    Logical, intent(in) :: up
    Complex (Kind=Kind(0.d0)) :: dp(Op%N_non_zero), dm(Op%N_non_zero)
    Complex (Kind=Kind(0.d0)) :: B(Op%N, Op%N_non_zero), S(Op%N_non_zero, Op%N_non_zero)

    if (N_Type == 1) then
       call op_lowrank_factors(Op, HS_Field, nt,  1, dp)
       call op_lowrank_factors(Op, HS_Field, nt, -1, dm)
       if (up) then
          call op_lowrank_qform(Op, 'C', B, S)
          call op_lowrank_rows(Op, Mat, Ndim, Ndim, B, S)
          call op_lowrank_scale(Op, Mat, Ndim, Ndim, dp, .true.)
          call op_lowrank_qform(Op, 'N', B, S)
          call op_lowrank_cols(Op, Mat, Ndim, Ndim, B, S)
          call op_lowrank_scale(Op, Mat, Ndim, Ndim, dm, .false.)
       else
          call op_lowrank_scale(Op, Mat, Ndim, Ndim, dm, .true.)
          call op_lowrank_qform(Op, 'N', B, S)
          call op_lowrank_rows(Op, Mat, Ndim, Ndim, B, S)
          call op_lowrank_scale(Op, Mat, Ndim, Ndim, dp, .false.)
          call op_lowrank_qform(Op, 'C', B, S)
          call op_lowrank_cols(Op, Mat, Ndim, Ndim, B, S)
       endif
    elseif (N_Type == 2) then
       if (up) then
          call op_lowrank_qform(Op, 'N', B, S)
          call op_lowrank_rows(Op, Mat, Ndim, Ndim, B, S)
          call op_lowrank_qform(Op, 'C', B, S)
          call op_lowrank_cols(Op, Mat, Ndim, Ndim, B, S)
       else
          call op_lowrank_qform(Op, 'C', B, S)
          call op_lowrank_rows(Op, Mat, Ndim, Ndim, B, S)
          call op_lowrank_qform(Op, 'N', B, S)
          call op_lowrank_cols(Op, Mat, Ndim, Ndim, B, S)
       endif
    endif
  end subroutine op_lowrank_wrap

!--------------------------------------------------------------------
!> @brief
!> Op_Wrap_panels for op_lowrank_wrap's G -> L G R: X -> L X as G's rows, Y -> R^T Y:
!>   up,   N_type 1: R^T = D^-1 Q^T          up,   N_type 2: R^T = conjg(Q)
!>   down, N_type 1: R^T = conjg(Q) D        down, N_type 2: R^T = Q^T
!--------------------------------------------------------------------
  subroutine op_lowrank_panels(Xpan, Ypan, Ndim, ncols, Op, HS_Field, N_Type, nt, up)
    Implicit none
    Integer, intent(in) :: Ndim, ncols
    Complex (Kind=Kind(0.d0)), intent(inout) :: Xpan(Ndim, ncols), Ypan(Ndim, ncols)
    Type (Operator), intent(in) :: Op
    Complex (Kind=Kind(0.d0)), intent(in) :: HS_Field
    Integer, intent(in) :: N_Type, nt
    Logical, intent(in) :: up
    Complex (Kind=Kind(0.d0)) :: dp(Op%N_non_zero), dm(Op%N_non_zero)
    Complex (Kind=Kind(0.d0)) :: B(Op%N, Op%N_non_zero), S(Op%N_non_zero, Op%N_non_zero)

    if (N_Type == 1) then
       call op_lowrank_factors(Op, HS_Field, nt,  1, dp)
       call op_lowrank_factors(Op, HS_Field, nt, -1, dm)
       if (up) then
          call op_lowrank_qform(Op, 'C', B, S)
          call op_lowrank_rows(Op, Xpan, Ndim, ncols, B, S)
          call op_lowrank_scale(Op, Xpan, Ndim, ncols, dp, .true.)
          call op_lowrank_qform(Op, 'T', B, S)
          call op_lowrank_rows(Op, Ypan, Ndim, ncols, B, S)
          call op_lowrank_scale(Op, Ypan, Ndim, ncols, dm, .true.)
       else
          call op_lowrank_scale(Op, Xpan, Ndim, ncols, dm, .true.)
          call op_lowrank_qform(Op, 'N', B, S)
          call op_lowrank_rows(Op, Xpan, Ndim, ncols, B, S)
          call op_lowrank_scale(Op, Ypan, Ndim, ncols, dp, .true.)
          call op_lowrank_qform(Op, 'J', B, S)
          call op_lowrank_rows(Op, Ypan, Ndim, ncols, B, S)
       endif
    elseif (N_Type == 2) then
       if (up) then
          call op_lowrank_qform(Op, 'N', B, S)
          call op_lowrank_rows(Op, Xpan, Ndim, ncols, B, S)
          call op_lowrank_qform(Op, 'J', B, S)
          call op_lowrank_rows(Op, Ypan, Ndim, ncols, B, S)
       else
          call op_lowrank_qform(Op, 'C', B, S)
          call op_lowrank_rows(Op, Xpan, Ndim, ncols, B, S)
          call op_lowrank_qform(Op, 'T', B, S)
          call op_lowrank_rows(Op, Ypan, Ndim, ncols, B, S)
       endif
    endif
  end subroutine op_lowrank_panels

!--------------------------------------------------------------------
!> @brief
!> Compare two Operator instances on their defining properties.
!> Returns .true. if operators are equal (within tolerance for
!> floating-point fields), .false. otherwise.
!> Only compares N, type, P, O, g, alpha, and g_t.
!> Derived fields (E, U, M_exp) are skipped.
!--------------------------------------------------------------------
  function Op_equal(Op_a, Op_b) result(equal)
    Implicit None
    Type (Operator), INTENT(IN) :: Op_a, Op_b
    Logical :: equal
    Real (Kind=Kind(0.d0)), parameter :: eq_tol = 1.D-10
    Real (Kind=Kind(0.d0)) :: tol, scale
    Integer :: i, j

    equal = .false.

    ! Integer fields: exact match
    if (Op_a%N /= Op_b%N) return
    if (Op_a%type /= Op_b%type) return

    ! Projector indices: exact match
    do i = 1, Op_a%N
       if (Op_a%P(i) /= Op_b%P(i)) return
    enddo

    ! Operator matrix: tolerance-based
    scale = max(1.d0, maxval(abs(Op_a%O(1:Op_a%N, 1:Op_a%N))), &
                      maxval(abs(Op_b%O(1:Op_b%N, 1:Op_b%N))))
    tol = eq_tol * scale
    do j = 1, Op_a%N
       do i = 1, Op_a%N
          if (abs(Op_a%O(i,j) - Op_b%O(i,j)) > tol) return
       enddo
    enddo

    ! Coupling constant g
    scale = max(1.d0, abs(Op_a%g), abs(Op_b%g))
    if (abs(Op_a%g - Op_b%g) > eq_tol * scale) return

    ! Operator shift alpha
    scale = max(1.d0, abs(Op_a%alpha), abs(Op_b%alpha))
    if (abs(Op_a%alpha - Op_b%alpha) > eq_tol * scale) return

    ! Time-dependent coupling g_t
    if (allocated(Op_a%g_t) .neqv. allocated(Op_b%g_t)) return
    if (allocated(Op_a%g_t)) then
       if (size(Op_a%g_t) /= size(Op_b%g_t)) return
       scale = max(1.d0, maxval(abs(Op_a%g_t)), maxval(abs(Op_b%g_t)))
       tol = eq_tol * scale
       do i = 1, size(Op_a%g_t)
          if (abs(Op_a%g_t(i) - Op_b%g_t(i)) > tol) return
       enddo
    endif

    equal = .true.
  end function Op_equal

!--------------------------------------------------------------------
!> @brief
!> Check if the Op_V operator sequence is symmetric, i.e., the product
!> of (I + g_n * V_n) in forward order equals the product in reverse order,
!> where V_n is operator O_n embedded in the Ndim space via projector P_n.
!> This is satisfied when operators are arranged palindromically,
!> or when all operators commute (e.g., on-site operators in the Hubbard model).
!>
!> @param[in] Op_V_arr  Operator array of shape (N_op, N_FL)
!> @param[in] N_FL_in   Number of flavors
!> @param[in] Ndim_in   Dimension of the full single-particle Hilbert space
!> @return    .true. if symmetric, .false. otherwise
!--------------------------------------------------------------------
  logical function Op_V_is_symmetric(Op_V_arr, N_FL_in, Ndim_in)
    Implicit None
    Type(Operator), intent(in) :: Op_V_arr(:,:)
    Integer, intent(in) :: N_FL_in, Ndim_in

    Integer :: n_op, n, nf, i, j, NN
    Complex(Kind=Kind(0.d0)), allocatable :: M_fwd(:,:), M_rev(:,:), M_tmp(:,:), M_op(:,:)
    Real(Kind=Kind(0.d0)) :: diff, scale
    Real(Kind=Kind(0.d0)), parameter :: sym_tol = 1.D-10

    Op_V_is_symmetric = .true.
    n_op = size(Op_V_arr, 1)
    NN = Ndim_in

    if (n_op <= 1) return

    allocate(M_fwd(NN, NN), M_rev(NN, NN), M_tmp(NN, NN), M_op(NN, NN))

    do nf = 1, N_FL_in
       ! Initialize M_fwd = Identity
       M_fwd = cmplx(0.d0, 0.d0, kind(0.d0))
       do i = 1, NN
          M_fwd(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
       enddo

       ! Forward product: (I + g_1*V_1) * ... * (I + g_N*V_N)
       do n = 1, n_op
          M_op = cmplx(0.d0, 0.d0, kind(0.d0))
          do i = 1, NN
             M_op(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
          enddo
          do j = 1, Op_V_arr(n, nf)%N
             do i = 1, Op_V_arr(n, nf)%N
                M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) = &
                   M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) + &
                   Op_V_arr(n, nf)%g * Op_V_arr(n, nf)%O(i, j)
             enddo
          enddo
          M_tmp = matmul(M_fwd, M_op)
          M_fwd = M_tmp
       enddo

       ! Initialize M_rev = Identity
       M_rev = cmplx(0.d0, 0.d0, kind(0.d0))
       do i = 1, NN
          M_rev(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
       enddo

       ! Reverse product: (I + g_N*V_N) * ... * (I + g_1*V_1)
       do n = n_op, 1, -1
          M_op = cmplx(0.d0, 0.d0, kind(0.d0))
          do i = 1, NN
             M_op(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
          enddo
          do j = 1, Op_V_arr(n, nf)%N
             do i = 1, Op_V_arr(n, nf)%N
                M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) = &
                   M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) + &
                   Op_V_arr(n, nf)%g * Op_V_arr(n, nf)%O(i, j)
             enddo
          enddo
          M_tmp = matmul(M_rev, M_op)
          M_rev = M_tmp
       enddo

       ! Compare forward and reverse products
       diff = maxval(abs(M_fwd - M_rev))
       scale = max(1.d0, maxval(abs(M_fwd)), maxval(abs(M_rev)))
       if (diff > sym_tol * scale) then
          Op_V_is_symmetric = .false.
          deallocate(M_fwd, M_rev, M_tmp, M_op)
          return
       endif
    enddo

    deallocate(M_fwd, M_rev, M_tmp, M_op)
  end function Op_V_is_symmetric

  function Op_is_real(Op) result(retval)
    Implicit None
    
    Type (Operator) , INTENT(IN)   :: Op
    Logical :: retval
    Real (Kind=Kind(0.d0)) :: myzero
    integer :: i,j
    
    retval = (Abs(aimag(Op%g)) < Abs(Op%g)*epsilon(1.D0))
    ! calculate a matrix scale
    myzero = maxval(abs(Op%E))*epsilon(Op%E)
    
    do i = 1, Op%N
       do j = 1, Op%N
          retval = retval .and. (Abs(aimag(Op%O(i,j))) < myzero)
       enddo
    enddo
  end function Op_is_real
  
  logical function operator_get_g_t_alloc(this)
    Implicit None
    
    class (Operator) , INTENT(IN)   :: this
    
    operator_get_g_t_alloc = this%g_t_alloc
    
  end function operator_get_g_t_alloc
end Module Operator_mod
