! The low-rank vertex path (Op%lowrank) must act as the dense one.
!
! For a vertex of rank nz < N/2 Op_set swaps the dense N x N products for
! e^{g phi O} = 1 + V (e^{g phi E} - 1) V^dag in the mmults, and for Householder
! reflectors Q in the wraps. Q's first nz columns are U's up to a phase, but its
! null-space columns differ, so a single wrap leaves G in a different basis than
! the dense wrap. What must agree is the cycle the sequential update runs: rotate
! (Wrapup N_type 1 / Wrapdo N_type 2), apply a rank-nz update in the eigenbasis,
! rotate back. The update below has the phase covariance of Upgrade2's (it reads
! only rows and columns P(1:nz) and a diagonal Delta), so the cycle is
! basis-independent exactly when the rotation's first nz columns are right.
!
! Checked per type (1..4), with and without g_t, rank 1 and 2 on a width-6 vertex
! in Ndim = 9, scattered (gathered path) and contiguous (in place):
!   * Op_mmultR / Op_mmultL for cop N, T, C (and sign +-1), low-rank vs dense;
!   * the up and down update cycles, low-rank vs dense;
!   * Op_Wrap_panels vs the explicit wrap on the low-rank operator (test 37's check);
!   * that Op_set picks the path only where it should (full rank, diagonal, type 0),
!     and leaves U whole where it does, so Op_exp matches the dense operator's;
!   * how ALF_LOWRANK_VERTEX's value is read (op_lowrank_env_mode).

Program LowrankVertex

   Use Operator_mod
   Use Fields_mod

   Implicit None

   Integer, parameter :: Ndim = 9, opn = 6, ncols = 3
   ! Scattered (gathered path) and a contiguous run (in place, as a full-lattice support).
   Integer, parameter :: Playout(opn,2) = reshape([8, 2, 5, 3, 9, 6,  3, 4, 5, 6, 7, 8], [opn, 2])
   Integer :: Pidx(opn), ip
   Complex (Kind=Kind(0.D0)), parameter :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)
   Complex (Kind=Kind(0.D0)) :: g0(Ndim,Ndim), gd(Ndim,Ndim), gl(Ndim,Ndim), gb(Ndim,Ndim)
   Complex (Kind=Kind(0.D0)) :: mr0(Ndim,5), mrd(Ndim,5), mrl(Ndim,5)
   Complex (Kind=Kind(0.D0)) :: ml0(4,Ndim), mld(4,Ndim), mll(4,Ndim)
   Complex (Kind=Kind(0.D0)) :: xp(Ndim,ncols), yp(Ndim,ncols), field
   Complex (Kind=Kind(0.D0)) :: vecs(opn,2), ed(opn,opn), el(opn,opn), eye(opn,opn)
   Real    (Kind=Kind(0.D0)) :: lam(2)
   Integer :: nt, nz, igt, icop, isg, iud, N_Type, i, j, nfail
   Character (Len=1) :: cop
   Character (Len=1), parameter :: cops(3) = ['n', 't', 'c']
   Type (Operator) :: Opd, Opl
   Type (Fields)   :: nsigma_single

   nfail = 0
   Call Fields_init()
   Call nsigma_single%make(1,1)
   lam = [0.9d0, -0.6d0]

   ! Two orthonormal complex vectors on the vertex's width (Gram-Schmidt).
   Do i = 1, opn
      vecs(i,1) = cmplx(sin(dble(2*i+1)), cos(dble(i*i)), kind(0.D0))
      vecs(i,2) = cmplx(cos(dble(3*i-2)), sin(dble(i+4)), kind(0.D0))
   End Do
   vecs(:,1) = vecs(:,1)/sqrt(sum(abs(vecs(:,1))**2))
   vecs(:,2) = vecs(:,2) - dot_product(vecs(:,1), vecs(:,2))*vecs(:,1)
   vecs(:,2) = vecs(:,2)/sqrt(sum(abs(vecs(:,2))**2))

   Do j = 1, Ndim
      Do i = 1, Ndim
         g0(i,j) = cmplx(0.4d0*sin(dble(i+2*j)), 0.4d0*cos(dble(3*i-j)), kind(0.D0))
      End Do
      g0(j,j) = g0(j,j) + 0.5d0
   End Do
   Do j = 1, 5
      Do i = 1, Ndim
         mr0(i,j) = cmplx(cos(dble(i*j)), sin(dble(i-2*j)), kind(0.D0))
      End Do
   End Do
   Do j = 1, Ndim
      Do i = 1, 4
         ml0(i,j) = cmplx(sin(dble(i+j*j)), cos(dble(2*i+j)), kind(0.D0))
      End Do
   End Do

   Do ip = 1, 2
   Pidx = Playout(:, ip)
   Do nt = 1, 4
   Do nz = 1, 2
   Do igt = 1, 2

      Op_lowrank_mode = 0
      Call build(Opd)
      Op_lowrank_mode = 1
      Call build(Opl)
      If (Opd%lowrank .or. .not. Opl%lowrank .or. Opl%P0 /= merge(0, 3, ip == 1)) Then
         Write (*,*) "PATH not chosen as expected: type", nt, "rank", nz, Opd%lowrank, Opl%lowrank, Opl%P0
         nfail = nfail + 1
      End If
      If (any(shape(Opl%U) /= [opn, opn]) .or. any(shape(Opd%U) /= [opn, opn])) Then
         Write (*,*) "U not whole: type", nt, "rank", nz, shape(Opl%U), shape(Opd%U)
         nfail = nfail + 1
      End If
      ! Op_exp reads every column of U, low-rank or not.
      eye = zero
      Do i = 1, opn
         eye(i,i) = one
      End Do
      Call Op_exp(Opd%g, Opd, ed)
      Call Op_exp(Opl%g, Opl, el)
      Call check("Op_exp", ed, el, eye)

      Select Case (nt)
      Case (1)
         field = cmplx(-1.d0, 0.d0, kind(0.D0))
      Case (2)
         field = cmplx(2.d0, 0.d0, kind(0.D0))
      Case (3)
         field = cmplx(0.813d0, 0.d0, kind(0.D0))
      Case (4)
         field = cmplx(-1.d0, 0.5d0, kind(0.D0))
      End Select
      nsigma_single%t(1) = nt
      nsigma_single%f(1,1) = field

      ! mmults
      Do icop = 1, 3
         cop = cops(icop)
         mrd = mr0; mrl = mr0
         Call Op_mmultR(mrd, Opd, field, cop, 1)
         Call Op_mmultR(mrl, Opl, field, cop, 1)
         Call check("mmultR "//cop, mrd, mrl, mr0)
         Do isg = -1, 1, 2
            mld = ml0; mll = ml0
            Call Op_mmultL(mld, Opd, field, cop, 1, isg)
            Call Op_mmultL(mll, Opl, field, cop, 1, isg)
            Call check("mmultL "//cop, mld, mll, ml0)
         End Do
      End Do

      ! Update cycles: up = Wrapup(1), update, Wrapup(2); down = Wrapdo(2), update, Wrapdo(1).
      gd = g0; gl = g0
      Call Op_Wrapup(gd, Opd, field, Ndim, 1, 1)
      Call Op_Wrapup(gl, Opl, field, Ndim, 1, 1)
      Call fake_update(gd); Call fake_update(gl)
      Call Op_Wrapup(gd, Opd, field, Ndim, 2, 1)
      Call Op_Wrapup(gl, Opl, field, Ndim, 2, 1)
      Call check("up cycle", gd, gl, g0)

      gd = g0; gl = g0
      Call Op_Wrapdo(gd, Opd, field, Ndim, 2, 1)
      Call Op_Wrapdo(gl, Opl, field, Ndim, 2, 1)
      Call fake_update(gd); Call fake_update(gl)
      Call Op_Wrapdo(gd, Opd, field, Ndim, 1, 1)
      Call Op_Wrapdo(gl, Opl, field, Ndim, 1, 1)
      Call check("down cycle", gd, gl, g0)

      ! Panels on the low-rank operator against its own explicit wrap.
      Do N_Type = 1, 2
      Do iud = 1, 2
         Do j = 1, ncols
            Do i = 1, Ndim
               xp(i,j) = cmplx(0.3d0*sin(dble(i*j+1)), 0.3d0*cos(dble(i-j)), kind(0.D0))
               yp(i,j) = cmplx(0.3d0*cos(dble(2*i+j)), 0.3d0*sin(dble(i+3*j)), kind(0.D0))
            End Do
         End Do
         gd = g0
         Call ZGEMM('N','T',Ndim,Ndim,ncols,one,xp,Ndim,yp,Ndim,one,gd,Ndim)
         gb = gd
         gl = g0
         If (iud == 1) Then
            Call Op_Wrapup(gd, Opl, field, Ndim, N_Type, 1)
            Call Op_Wrapup(gl, Opl, field, Ndim, N_Type, 1)
            Call Op_Wrap_panels(xp, yp, Opl, field, Ndim, ncols, N_Type, 1, 'u')
         Else
            Call Op_Wrapdo(gd, Opl, field, Ndim, N_Type, 1)
            Call Op_Wrapdo(gl, Opl, field, Ndim, N_Type, 1)
            Call Op_Wrap_panels(xp, yp, Opl, field, Ndim, ncols, N_Type, 1, 'd')
         End If
         Call ZGEMM('N','T',Ndim,Ndim,ncols,one,xp,Ndim,yp,Ndim,one,gl,Ndim)
         Call check("panels", gd, gl, gb)
      End Do
      End Do

      Call Op_clear(Opd, opn)
      Call Op_clear(Opl, opn)

   End Do
   End Do
   End Do
   End Do

   ! Where the path must not be taken.
   Op_lowrank_mode = 1
   Call path_case("full rank", 2, 1, .false.)
   Call path_case("diagonal", 2, 2, .true.)
   Call path_case("type 0", 0, 3, .false.)

   ! ALF_LOWRANK_VERTEX: words as well as digits, anything else refused (-1).
   Call env_case("0", 0);   Call env_case(" OFF ", 0); Call env_case("False", 0); Call env_case("no", 0)
   Call env_case("1", 1);   Call env_case("on", 1);    Call env_case("TRUE", 1);  Call env_case("yes", 1)
   Call env_case("2", -1);  Call env_case("of", -1);   Call env_case("0.0", -1);  Call env_case("", -1)

   Call nsigma_single%clear()
   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If
   Write (*,*) "SUCCESS"

Contains

   ! O = sum_k lam_k v_k v_k^dag on P = Pidx, type nt, rank nz, g_t if igt == 2.
   Subroutine build(Op)
      Type (Operator), Intent(InOut) :: Op
      Integer :: a, b, k
      Call Op_make(Op, opn)
      Op%P = Pidx
      Do b = 1, opn
         Do a = 1, opn
            Op%O(a,b) = zero
            Do k = 1, nz
               Op%O(a,b) = Op%O(a,b) + lam(k)*vecs(a,k)*conjg(vecs(b,k))
            End Do
         End Do
      End Do
      Do a = 1, opn
         Op%O(a,a) = cmplx(real(Op%O(a,a)), 0.d0, kind(0.D0))
         Do b = a+1, opn
            Op%O(b,a) = conjg(Op%O(a,b))
         End Do
      End Do
      Op%type  = nt
      Op%g     = cmplx(0.7d0, 0.d0, kind(0.D0))
      Op%alpha = zero
      If (igt == 2) Then
         Allocate (Op%g_t(1))
         Op%g_t(1) = cmplx(0.41d0, 0.23d0, kind(0.D0))
      End If
      Call Op_set(Op)
   End Subroutine build

   ! G <- G - G(:,S) Delta X^-1 (1 - G)(S,:), X = 1 + Delta (1 - G)(S,S), S = P(1:nz):
   ! covariant under a phase per eigenvector and blind to the null-space basis.
   Subroutine fake_update(G)
      Complex (Kind=Kind(0.D0)), Intent(InOut) :: G(Ndim,Ndim)
      Complex (Kind=Kind(0.D0)) :: X(nz,nz), R(nz,Ndim), C(Ndim,nz), delta(2)
      Integer :: ipiv(nz), info, a, b
      delta = [cmplx(0.3d0, 0.1d0, kind(0.D0)), cmplx(-0.2d0, 0.25d0, kind(0.D0))]
      Do a = 1, nz
         R(a,:) = -G(Pidx(a),:)
         R(a,Pidx(a)) = R(a,Pidx(a)) + one
         C(:,a) = G(:,Pidx(a))
      End Do
      Do b = 1, nz
         Do a = 1, nz
            X(a,b) = delta(a)*R(a,Pidx(b))
         End Do
         X(b,b) = X(b,b) + one
      End Do
      Do a = 1, nz
         R(a,:) = delta(a)*R(a,:)
      End Do
      Call ZGESV(nz, Ndim, X, nz, ipiv, R, nz, info)
      Call ZGEMM('N','N',Ndim,Ndim,nz,-one,C,Ndim,R,nz,one,G,Ndim)
   End Subroutine fake_update

   ! Relative agreement, and that the operation moved the matrix at all.
   Subroutine check(label, ref, cand, before)
      Character (Len=*), Intent(In) :: label
      Complex (Kind=Kind(0.D0)), Intent(In) :: ref(:,:), cand(:,:), before(:,:)
      Real (Kind=Kind(0.D0)) :: scale, err, moved
      scale = maxval(abs(ref))
      err   = maxval(abs(ref - cand))/scale
      moved = maxval(abs(ref - before))/scale
      If (err > 1.d-12 .or. moved < 1.d-6) Then
         Write (*,*) "ERROR ", label, " type", nt, "rank", nz, "g_t", igt == 2, &
            &        " rel err", err, " moved", moved
         nfail = nfail + 1
      End If
   End Subroutine check

   Subroutine env_case(text, expected)
      Character (Len=*), Intent(In) :: text
      Integer, Intent(In) :: expected
      If (op_lowrank_env_mode(text) /= expected) Then
         Write (*,*) "ALF_LOWRANK_VERTEX='", text, "' read as", op_lowrank_env_mode(text), &
            &        " expected", expected
         nfail = nfail + 1
      End If
   End Subroutine env_case

   Subroutine path_case(label, nz_case, mode, diag)
      Character (Len=*), Intent(In) :: label
      Integer, Intent(In) :: nz_case, mode
      Logical, Intent(In) :: diag
      Type (Operator) :: Op
      Integer :: a
      Call Op_make(Op, 4)
      Op%P = [1, 2, 3, 4]
      If (diag) Then
         Op%O(1,1) = one
      Else
         ! Full rank: a 4-site ring hopping has eigenvalues +-2, 0, 0 -- so give it
         ! a diagonal too, making all four nonzero.
         Do a = 1, 4
            Op%O(a, mod(a,4)+1) = one
            Op%O(mod(a,4)+1, a) = one
            Op%O(a,a) = cmplx(0.3d0*dble(a), 0.d0, kind(0.D0))
         End Do
      End If
      Op%type = nz_case
      Op%g = cmplx(0.5d0, 0.d0, kind(0.D0))
      If (mode == 3) Then
         ! Rank 1 and non-diagonal, so only the type excludes it.
         Op%type = 0
         Op%O = zero
         Op%O(1:2,1:2) = one
      End If
      Call Op_set(Op)
      If (Op%lowrank) Then
         Write (*,*) "PATH taken for ", label
         nfail = nfail + 1
      End If
      Call Op_clear(Op, 4)
   End Subroutine path_case

End Program LowrankVertex
