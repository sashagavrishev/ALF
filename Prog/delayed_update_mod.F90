!  Copyright (C) 2016 - 2026 The ALF project
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
!     Under Section 7 of GPL version 3 we require you to fulfill the following
!     additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific
!       community. Being part of that community we feel that it is reasonable to
!       require you to give an attribution back to the original authors if you
!       have benefitted from this program. Guidelines for a proper citation can
!       be found on the project's homepage http://alf.physik.uni-wuerzburg.de.
!
!     - We require the preservation of the above copyright notice and this
!       license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source
!       files. To obtain the original source files please visit the homepage
!       http://alf.physik.uni-wuerzburg.de.
!
!     - If you make substantial changes to the program we require you to either
!       consider contributing to the ALF project or to mark your material in a
!       reasonable way as different from the original version.

!-------------------------------------------------------------------------------
!> @brief
!> This module enables us to hold the Green's function in factored form across
!> one imaginary time slice; an accepted update will now append a column pair
!> instead of updating the whole matrix. Thus this is a method of "delaying" the
!> Green's function updates, with the aim of reducing the burden on the memory
!> system of the execution environment. This module is designed only for the
!> sequential auxiliary field sampler.
!>
!> @details
!> When a proposed field change is accepted a rank-d update has to be applied to
!> the entire Green's function (see upgrade_mod). An "instantaneous" scheme
!> requires this rank-d BLAS operation at every such update, which is memory
!> bandwidth bound and hence can cause performance degradation. This is
!> especially true for large matrix dimensions and when the compute environment
!> is fully saturated by many other memory intensive jobs.
!>
!> We can take advantage of the fact that no step of the sequential sweep needs
!> the full Green's function, G. With P = Op%P and d = Op%N_non_zero:
!>
!>    - every proposal, accepted or not, needs only the d x d block G(P,P)
!>      to form the Metropolis ratio;
!>
!>    - an accepted proposal additionally needs the d rows G(P,:) and d
!>      columns G(:,P) of the same vertex, which form the factors of its
!>      rank-d update.
!>
!> Thus we are at liberty to split the Green's function into
!>
!>             G = G_stale + X * Y^T,   X, Y of shape (Ndim, ~k)
!>
!> and pay Ndim**2 only once every k panel columns, i.e. every ~k/d accepted
!> updates. Traffic per accepted field update is then expected to fall from
!> ~2*Ndim**2 to ~d*Ndim*k + 2*d*Ndim**2/k.
!> The matrices X and Y are referred to as "panels". The scheme, including its
!> generalisation to vertices of rank d > 1, follows F. Sun and X. Y. Xu,
!> Phys. Rev. B 109, 235140 (2024); see the "Delayed (rank-k) updates" section
!> of the ALF documentation.
!>
!> The implementation is such that the Green's function is in its "factorised"
!> form only within the sequential vertex loop of a single time slice;
!> thus stabilisation, measurement and global moves always receive the fully
!> flushed G.
!>
!> By default delayed updates are off. The environment variable ALF_DELAY_K set
!> to the appropriate value enables it; see delay_resolve for more details.
!-------------------------------------------------------------------------------

module delayed_update_mod

   use operator_mod
   use runtime_error_mod
   use iso_fortran_env, only: error_unit

   implicit none

   private
   public :: delay_alloc, delay_dealloc, delay_resolve
   public :: delay_assert_inactive, delay_open, delay_close
   public :: delay_block, delay_row, delay_col, delay_append, delay_flush
   public :: delay_wrap, delay_pending, delay_log

   ! Define the panels X and Y with their live column count, ncol, one set per
   ! flavour. The panels are only Ndim*(k+dmax) per flavour hence we can
   ! afford to allocate them once.
   complex (Kind=Kind(0.d0)), private, save, allocatable :: xp(:,:,:), yp(:,:,:)
   integer,                   private, save, allocatable :: ncol(:)

   integer, private, save :: kmax    = 0 ! Flush threshold, the k of the scheme
   integer, private, save :: panel_w = 0 ! Allocated panel width, kmax + dmax
   integer, private, save :: ndim_s  = 0 ! Ndim, recorded for the GR dummies
   integer, private, save :: nfl_s   = 0 ! N_FL, likewise

   ! Record if a factored region is open
   logical, public, protected, save :: delay_active = .false.

   ! ALF_DELAY_K as parsed by delay_resolve, kept for delay_log: 0 disables the
   ! delay and a positive value fixes the depth, so named requests are negative.
   integer, private, parameter :: K_AUTO    = -1 ! "auto": measure the depth
   integer, private, parameter :: K_FORMULA = -2 ! "formula": k ~ sqrt(2*Ndim)
   integer, private, save      :: k_request = 0

   ! Constants on the depth "auto" may resolve to: a ceiling and floor on the
   ! possible delay depths we can set.
   integer, private, parameter :: K_FLOOR   = 8
   integer, private, parameter :: K_CEILING = 256

   ! How the depth in force was arrived at, for the info file: off, fixed,
   ! probe or formula. Protected rather than behind an accessor, since
   ! Wrapgr_mod hands it straight to Control_set_delay_depth.
   character (Len=16), public, protected, save :: delay_source = 'off'

   ! Depths the delay depth probe times, ascending.
   integer, private, parameter :: K_CAND(*) = [8, 16, 32, 64, 128, 256]
   integer, private, parameter :: N_CAND = size(K_CAND)

   ! Delay depth probe timing details.
   real (Kind=Kind(0.d0)), private, parameter :: PROBE_MIN_SECONDS = 5.d-3
   integer,                private, parameter :: PROBE_MAX_REPS    = 4096

   ! Sweep the probe multiple times to average over the background of the
   ! execution environment.
   integer, private, parameter :: PROBE_SWEEPS = 3

   ! Define a baseline margin in the comparison of the results from the
   ! delay depth probe.
   real (Kind=Kind(0.d0)), private, parameter :: PROBE_MARGIN = 1.05d0

   ! Complex identity "z-one"
   complex (Kind=Kind(0.d0)), private, parameter :: ZONE = (1.d0, 0.d0)

   ! Kernel selector for the flush / panel
   integer, private, parameter :: PROBE_FLUSH = 1
   integer, private, parameter :: PROBE_PANEL = 2

   ! Details for the delay_log: the probe's cost per candidate (-1 where it
   ! never ran), its own wall clock, and ALF_DELAY_K as it was given.
   real (Kind=Kind(0.d0)), private, save :: probe_cost(N_CAND) = -1.d0
   real (Kind=Kind(0.d0)), private, save :: probe_seconds = 0.d0
   real (Kind=Kind(0.d0)), private, save :: probe_scratch_mb = 0.d0
   character (Len=32),     private, save :: k_request_text = '<unset>'

   ! Every GR dummy below is explicit shape rather than assumed shape. The
   ! arrays are handed element-first to ZGEMM and ZCOPY, which have no explicit
   ! interface, and sequence association from an assumed-shape actual is not
   ! something the standard guarantees -- a compiler may pass a descriptor or a
   ! copy. upgrade_mod declares the same array the same way for the same
   ! reason.

   ! Argument names shared by the routines below: nf is the flavour, d the rank
   ! of the vertex being applied -- its number of column pairs -- and P(d) the
   ! operator's support Op%P, the rows and columns of G that vertex touches.
   ! The local c is that flavour's live column count, ncol(nf).

contains

!-------------------------------------------------------------------------------
!> @brief
!> Resolve the delay depth for this run; 0 when the delayed update is disabled.
!>
!> @details
!> Reads ALF_DELAY_K, case-insensitively:
!>
!>    - unset, empty, "0", negative, unreadable or over 32 characters: delay
!>      off, with the reason in delay_log
!>    - a positive integer: used verbatim, not clamped
!>    - "formula": delay_formula, k ~ sqrt(2*Ndim)
!>    - "auto": delay_probe times the flush and panel costs at this Ndim and
!>      takes the largest k within PROBE_MARGIN of the cheapest, falling back
!>      to delay_formula if the measurement is impossible or inconclusive
!>
!> "formula" and "auto" are clamped to [K_FLOOR, K_CEILING], and turn the delay
!> off when Ndim < K_FLOOR. Also sets delay_source and what delay_log reports.
!>
!> Call once and pass the result to delay_alloc: "auto" is a timing, so a
!> second call may answer differently. Under MPI only one rank should call it
!> and broadcast the result, as Wrapgr_delay_alloc does: delay_probe holds an
!> Ndim**2 scratch, and ranks probing at once would contend for the very memory
!> system they are measuring.
!>
!> A measured depth may differ between runs of one chain. This is harmless:
!> different k give the same chain up to rounding, which can only flip
!> Metropolis decisions that are near-ties.
!>
!> The depth is deliberately not a simulation parameter, so that a
!> run-to-run choice does not show up as a difference between runs' parameters.
!-------------------------------------------------------------------------------

   integer function delay_resolve(Ndim)

      implicit none

      integer, intent(in) :: Ndim
      character (Len=32) :: text   ! ALF_DELAY_K as the environment gave it
      character (Len=32) :: word   ! the same, trimmed and lower-cased
      integer :: length, status    ! from get_environment_variable
      integer :: value             ! the depth, where the request was a number
      integer :: i

      k_request      = 0
      k_request_text = '<unset>'
      delay_source   = 'off'

      call get_environment_variable("ALF_DELAY_K", text, length, status)
      if (status == -1) then
         ! Longer than text holds: turn delays off, but say so.
         k_request_text = '<too long>'
      else if (status == 0 .and. length > 0) then
         word           = trim(adjustl(text(1:length)))
         k_request_text = word
         ! Capitalisation agnostic
         do i = 1, len_trim(word)
            if (word(i:i) >= 'A' .and. word(i:i) <= 'Z') &
            &  word(i:i) = achar(iachar(word(i:i)) + 32)
         enddo
         select case (word)
          case ('auto')
            k_request = K_AUTO
          case ('formula')
            k_request = K_FORMULA
          case default
            ! Unreadable or negative turns delays off, but say so.
            read (word, *, iostat=status) value
            if (status /= 0) then
               k_request_text = trim(k_request_text)//' (unreadable)'
            else if (value < 0) then
               k_request_text = trim(k_request_text)//' (negative)'
            else
               k_request = value
            endif
         end select
      endif

      select case (k_request)
       case (K_FORMULA)
         delay_resolve = delay_formula(Ndim)
         delay_source  = 'formula'
       case (K_AUTO)
         ! delay_probe sets delay_source itself, since it may have to fall back
         ! to the formula.
         delay_resolve = delay_probe(Ndim)
       case default
         delay_resolve = k_request
         delay_source  = 'fixed'
      end select
      ! Whichever path gave it, a zero depth is off.
      if (delay_resolve <= 0) delay_source = 'off'

      ! A named request only resolves to zero on a matrix below the floor.
      if (delay_resolve <= 0 .and. k_request < 0) then
         write (word, '(i0)') K_FLOOR
         k_request_text = trim(k_request_text)//' (Ndim < '//trim(word)//')'
      endif
   end function delay_resolve

!-------------------------------------------------------------------------------
!> @brief
!> Log what the delay is set to and the decision pathway.
!>
!> @details
!> Called once at setup, from main under the rank guard.
!-------------------------------------------------------------------------------

   subroutine delay_log(unit)
      implicit none
      integer, intent(in) :: unit      ! where to write; main passes 6
      integer :: i
      real (Kind=Kind(0.d0)) :: lo          ! best cost, the curve's normaliser
      character (Len=13) :: mark, tag       ! marker for the depth in force

      write (unit,'(a)')  ' Delayed update:'
      write (unit,'(2a)') '   ALF_DELAY_K            : ', trim(k_request_text)
      if (kmax <= 0) then
         write (unit,'(a)') '   status                 : off (immediate)'
         return
      endif
      write (unit,'(a,i0)') '   Ndim                   : ', ndim_s
      write (unit,'(a,i0)') '   depth k                : ', kmax
      write (unit,'(2a)')   '   chosen by              : ', trim(delay_source)
      write (unit,'(a,i0)') '   panel width (k + dmax) : ', panel_w
      write (unit,'(a,i0,a,i0,a)') '   validated range        : [', &
      & K_FLOOR, ', ', K_CEILING, ']'

      ! Write a note when one had to fall back to the formula on a refused probe
      if (trim(delay_source) == 'formula' .and. k_request == K_AUTO) &
      & write (unit,'(a)') '   WARN: the probe was refused; this is the formula'

      ! Still at its -1 default: the probe never ran, so there is no curve.
      if (probe_cost(1) < 0.d0) return

      write (unit,'(a,f6.3,a,f10.3,a)') '   probe cost             : ', &
      & probe_seconds, ' s, scratch ', probe_scratch_mb, ' MB'
      write (unit,'(a)') '        k   rel. cost   (1.00 = best)'

      lo = minval(probe_cost, mask=(probe_cost > 0.d0))

      ! Mark the depth in force; on a fallback it came from the formula.
      mark = '   <- chosen'
      if (trim(delay_source) == 'formula') mark = '   <- formula'
      do i = 1, N_CAND
         tag = ''
         if (K_CAND(i) == kmax) tag = mark
         if (K_CAND(i) > ndim_s) then
            write (unit,'(a,i5,a)') '   ', K_CAND(i), '       --  (above Ndim)'
         else if (probe_cost(i) <= 0.d0) then
            write (unit,'(a,i5,a)') '   ', K_CAND(i), '       --  (not timed)'
         else
            write (unit,'(a,i5,f12.3,a)') '   ', K_CAND(i), probe_cost(i)/lo, &
            & trim(tag)
         endif
      enddo
   end subroutine delay_log

!-------------------------------------------------------------------------------
!> @brief
!> The closed-form depth: nint(sqrt(2*Ndim)), clamped; 0 when Ndim < K_FLOOR.
!>
!> @details
!> Minimises the traffic per accepted update, d*Ndim*k for the panel
!> reconstructions (2*d ZGEMVs against a half-full panel) plus 2*d*Ndim**2/k
!> for the flush; d cancels. It assumes both terms move bytes at the same
!> cost. In practice the panel is often cache-resident and the flush partly
!> compute-bound, so the true optimum tends to lie above this estimate.
!-------------------------------------------------------------------------------

   integer function delay_formula(Ndim)
      implicit none
      integer, intent(in) :: Ndim
      ! Below the floor the matrix is too small for the delay to pay. Above it
      ! the clamp keeps k within the range the delayed path has been exercised
      ! over, which for Ndim >= K_FLOOR is also never wider than the matrix.
      if (Ndim < K_FLOOR) then
         delay_formula = 0
      else
         delay_formula = min(K_CEILING, &
         &               max(K_FLOOR, nint(sqrt(2.d0*real(Ndim, Kind(0.d0))))))
      endif
   end function delay_formula

!-------------------------------------------------------------------------------
!> @brief
!> Pick the delay depth by timing the two k-dependent costs at this Ndim.
!>
!> @details
!> Per accepted flip the delayed scheme pays a flush ZGEMM('N','T',Ndim,Ndim,k)
!> once every k/d flips, and 2*d panel ZGEMVs against a panel that is half full
!> on average:
!>
!>     cost(k) = t_gemm(k)*d/k + 2*d*t_gemv(k/2)
!>             = d * [ t_gemm(k)/k + 2*t_gemv(k/2) ]
!>
!> We notice here that d factors out; the result is operator rank independent.
!> t_gemm and t_gemv are times that are probed at runtime.
!>
!> Falls back to delay_formula when fewer than two candidates fit Ndim, the
!> scratch cannot be allocated, the clock gives no reading, or the curve is
!> flat to within PROBE_MARGIN.
!-------------------------------------------------------------------------------

   integer function delay_probe(Ndim)

      implicit none

      integer, intent(in) :: Ndim

      ! Scratch the kernels run on: g stands in for the Green's function, xs
      ! and ys for the panels, v for a panel coefficient vector and w for the
      ! ZGEMV target.
      complex (Kind=Kind(0.d0)), allocatable :: g(:,:), xs(:,:), ys(:,:)
      complex (Kind=Kind(0.d0)), allocatable :: v(:), w(:)

      ! cost(i) is the modelled per-flip cost at K_CAND(i), built from the flush
      ! time tg and the panel time tv of one reading, this. lo and hi bound the
      ! finished curve.
      real (Kind=Kind(0.d0)) :: cost(N_CAND), tg, tv, lo, hi, this

      ! k is the candidate depth, c the half occupancy the panel is timed at,
      ! kwide the widest candidate that fits Ndim and so the width allocated,
      ! best the winner, stat the allocation status.
      integer :: i, k, c, kwide, stat, sweep, best

      ! Wall clock over the whole ladder, which delay_log reports.
      integer (Kind=8) :: wall0, wall1, wall_rate

      ! The fallback
      delay_source = 'formula'
      delay_probe  = delay_formula(Ndim)

      ! No candidate may exceed Ndim
      kwide = 0
      do i = 1, N_CAND
         if (K_CAND(i) <= Ndim) kwide = K_CAND(i)
      enddo
      ! A curve needs at least two points.
      if (kwide < K_CAND(2)) return

      allocate (g(Ndim,Ndim), xs(Ndim,kwide), ys(Ndim,kwide), &
      & v(kwide), w(Ndim), stat=stat)
      if (stat /= 0) return

      ! Everything just allocated, at 16 bytes per complex entry.
      probe_scratch_mb = 16.d0*(real(Ndim, Kind(0.d0))**2 &
      &                  + 2.d0*real(Ndim*kwide, Kind(0.d0)) &
      &                  + real(kwide + Ndim, Kind(0.d0)))/1048576.d0

      call probe_fill(g,  Ndim*Ndim)
      call probe_fill(xs, Ndim*kwide)
      call probe_fill(ys, Ndim*kwide)
      call probe_fill(v,  kwide)
      call probe_fill(w,  Ndim)

      cost = huge(1.d0)

      call system_clock(wall0)

      do sweep = 1, PROBE_SWEEPS
         do i = 1, N_CAND
            k = K_CAND(i)
            if (k > kwide) cycle
            c  = k/2
            tg = probe_time(PROBE_FLUSH, g, xs, ys, v, w, Ndim, k)
            tv = probe_time(PROBE_PANEL, g, xs, ys, v, w, Ndim, c)
            this = tg/real(k, Kind(0.d0)) + 2.d0*tv
            ! Contention only ever adds time, so the least reading is the
            ! least polluted estimator; a mean carries every burst it saw.
            cost(i) = min(cost(i), this)
         enddo
      enddo

      call system_clock(wall1, wall_rate)

      if (wall_rate > 0) probe_seconds = &
      & real(wall1 - wall0, Kind(0.d0))/real(wall_rate, Kind(0.d0))
      probe_cost = cost

      deallocate (g, xs, ys, v, w)

      ! When the probe curve is ~flat, fall back to the formula. Untimed
      ! candidates sit at the largest representable number "huge".
      lo = minval(cost, mask=(cost < huge(1.d0)))
      hi = maxval(cost, mask=(cost < huge(1.d0)))
      if (lo <= 0.d0 .or. hi < PROBE_MARGIN*lo) return

      ! Take the largest candidate within PROBE_MARGIN of the best; overshooting
      ! tends to give better performance on average.
      best = maxval(K_CAND, mask=(cost <= PROBE_MARGIN*lo))

      delay_probe  = best
      delay_source = 'probe'
   end function delay_probe

!-------------------------------------------------------------------------------
!> @brief
!> Fill probe scratch with entries of modulus ~1 and no dominant diagonal.
!>
!> @details
!> Deliberately not an RNG: the probe is a timing, so it would consume a
!> different number of draws on each machine and run, and the Markov chain's
!> generator state must not depend on that. Assumed size, so one routine
!> serves both the matrices and the vectors.
!-------------------------------------------------------------------------------

   subroutine probe_fill(a, n)
      implicit none
      integer, intent(in) :: n
      complex (Kind=Kind(0.d0)), intent(out) :: a(*)
      integer :: i
      do i = 1, n
         a(i) = cmplx(sin(real(i, Kind(0.d0))), cos(real(3*i, Kind(0.d0))), &
         &            Kind(0.d0))
      enddo
   end subroutine probe_fill

!-------------------------------------------------------------------------------
!> @brief
!> Seconds for one call of a probe kernel, repeated until the clock resolves.
!>
!> @details
!> PROBE_FLUSH times one flush ZGEMM('N','T',Ndim,Ndim,n), PROBE_PANEL one
!> panel ZGEMV against n live columns. A clock reporting no rate returns zero,
!> which delay_probe reads as a failed measurement.
!-------------------------------------------------------------------------------

   real (Kind=Kind(0.d0)) function probe_time(kernel, g, xs, ys, v, w, Ndim, n)
      implicit none
      ! kernel is PROBE_FLUSH or PROBE_PANEL; n is the depth k for the flush
      ! and the live column count for the panel. Both arms take every array so
      ! the two calls in delay_probe read alike.
      integer, intent(in) :: kernel, Ndim, n
      complex (Kind=Kind(0.d0)), intent(inout) :: g(Ndim,Ndim), w(Ndim)
      complex (Kind=Kind(0.d0)), intent(in)    :: xs(Ndim,*), ys(Ndim,*), v(*)
      integer :: rep, reps   ! repeats, raised until the clock resolves
      integer (Kind=8) :: c0, c1, rate
      reps = 1
      do
         call system_clock(c0)
         do rep = 1, reps
            if (kernel == PROBE_FLUSH) then
               call ZGEMM('N', 'T', Ndim, Ndim, n, ZONE, xs, Ndim, &
               &          ys, Ndim, ZONE, g, Ndim)
            else
               call ZGEMV('N', Ndim, n, ZONE, xs, Ndim, v, 1, &
               &          ZONE, w, 1)
            endif
         enddo
         call system_clock(c1, rate)
         ! No rate means no measurement, at any repeat count: return the zero
         ! delay_probe reads as a failure rather than escalating to
         ! PROBE_MAX_REPS full-size kernels to be told the same thing.
         if (rate <= 0) then
            probe_time = 0.d0
            return
         endif
         probe_time = real(c1 - c0, Kind(0.d0))/real(rate, Kind(0.d0))
         if (probe_time >= PROBE_MIN_SECONDS .or. reps >= PROBE_MAX_REPS) exit
         reps = reps*4
      enddo
      probe_time = probe_time/real(reps, Kind(0.d0))
   end function probe_time

!-------------------------------------------------------------------------------
!> @brief
!> Live column count for one flavour, for tests and diagnostics.
!-------------------------------------------------------------------------------

   integer function delay_pending(nf)
      implicit none
      integer, intent(in) :: nf
      delay_pending = 0
      if (allocated(ncol)) delay_pending = ncol(nf)
   end function delay_pending

!-------------------------------------------------------------------------------
!> @brief
!> Fail if a factored region is open where the caller cannot cope.
!>
!> @details
!> Wrapgr_PlaceGR and Wrapgr_Random_update read and copy the whole Green's
!> function (GR_st = Gr, Gr = GR_st for the multi-flip restore), which a
!> factored one would silently corrupt. Today neither can run inside the region
!> -- both are reached only outside the sequential vertex loop -- so this
!> asserts an invariant rather than handling a case. A future caller that breaks
!> it stops here instead of quietly reading a stale matrix.
!-------------------------------------------------------------------------------

   subroutine delay_assert_inactive(where)
      implicit none
      character(len=*), intent(in) :: where
      if (delay_active) then
         write(error_unit,*) 'delayed_update: ', trim(where), &
         & ' reached with a factored Green function open'
         Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
   end subroutine delay_assert_inactive

!-------------------------------------------------------------------------------
!> @brief
!> Allocate the panels, kmax + dmax columns wide; no-op when k <= 0.
!>
!> @param[in] dmax Most columns one update appends, maxval(Op_V%N_non_zero).
!> @param[in] k    The depth, from delay_resolve.
!>
!> @details
!> The extra dmax columns hold the append that crosses kmax before its flush.
!-------------------------------------------------------------------------------

   subroutine delay_alloc(Ndim, N_FL, dmax, k)
      implicit none
      integer, intent(in) :: Ndim, N_FL, dmax, k

      ! Recorded even when the delay is off; see the note on the GR dummies.
      ndim_s = Ndim
      nfl_s  = N_FL

      kmax = max(k, 0)
      if (kmax <= 0) return

      panel_w = kmax + dmax
      allocate (xp(Ndim, panel_w, N_FL), yp(Ndim, panel_w, N_FL))
      allocate (ncol(N_FL))
      ncol   = 0
      delay_active = .false.
   end subroutine delay_alloc

   subroutine delay_dealloc()
      implicit none
      if (allocated(xp)) deallocate (xp)
      if (allocated(yp)) deallocate (yp)
      if (allocated(ncol)) deallocate (ncol)
      ! Back to the state delay_alloc found, so that a second allocation cannot
      ! inherit the shape of the first.
      kmax         = 0
      panel_w      = 0
      ndim_s       = 0
      nfl_s        = 0
      delay_active = .false.
   end subroutine delay_dealloc

!-------------------------------------------------------------------------------
!> @brief
!> Open a factored region. No-op when the delay is disabled.
!-------------------------------------------------------------------------------

   subroutine delay_open()
      implicit none
      if (kmax <= 0) return
      ncol   = 0
      delay_active = .true.
   end subroutine delay_open

!-------------------------------------------------------------------------------
!> @brief
!> Flush every flavour into GR and close the region.
!>
!> @details
!> Whatever is still pending has to be paid for: leaving it out would credit the
!> chain with updates it never applied.
!-------------------------------------------------------------------------------

   subroutine delay_close(GR)
      implicit none
      Complex (Kind=Kind(0.d0)), intent(inout) :: GR(ndim_s, ndim_s, nfl_s)
      integer :: nf
      if (.not. delay_active) return
      do nf = 1, nfl_s
         call delay_flush(nf, GR)
      enddo
      delay_active = .false.
   end subroutine delay_close

!-------------------------------------------------------------------------------
!> @brief
!> Apply the pending columns of one flavour to GR: G_stale += X*Y^T.
!-------------------------------------------------------------------------------

   subroutine delay_flush(nf, GR)
      implicit none
      integer, intent(in) :: nf
      Complex (Kind=Kind(0.d0)), intent(inout) :: GR(ndim_s, ndim_s, nfl_s)
      if (kmax <= 0) return
      if (ncol(nf) <= 0) return
      call ZGEMM('N', 'T', ndim_s, ndim_s, ncol(nf), ZONE, xp(1,1,nf), &
      &          ndim_s, yp(1,1,nf), ndim_s, ZONE, GR(1,1,nf), ndim_s)
      ncol(nf) = 0
   end subroutine delay_flush

!-------------------------------------------------------------------------------
!> @brief
!> The d x d block of the *current* Green's function on the operator's support.
!>
!> @details
!> blk(n,m) = G(P(n), P(m)) = G_stale(P(n),P(m)) + sum_c X(P(n),c)*Y(P(m),c).
!>
!> O(d**2 * c), and paid on every proposal including the rejected ones -- this
!> is the cost the delay adds, and the reason a very large k stops paying.
!-------------------------------------------------------------------------------

   subroutine delay_block(nf, GR, P, d, blk, ldb)
      implicit none
      integer, intent(in) :: nf, d, P(d)
      integer, intent(in) :: ldb ! Leading dimension of blk
      Complex (Kind=Kind(0.d0)), intent(in)  :: GR(ndim_s, ndim_s, nfl_s)
      Complex (Kind=Kind(0.d0)), intent(out) :: blk(ldb,*)
      integer :: n, m   ! row and column within the support
      integer :: c
      c = ncol(nf)
      ! A sum over an empty panel is zero, so this covers c = 0 as it stands.
      do m = 1, d
         do n = 1, d
            blk(n,m) = GR(P(n), P(m), nf) &
            &        + sum(xp(P(n),1:c,nf)*yp(P(m),1:c,nf))
         enddo
      enddo
   end subroutine delay_block

!-------------------------------------------------------------------------------
!> @brief
!> d rows of the current Green's function: rows(i,l) = G(P(l), i).
!>
!> @details
!> Stale row plus Y * X(P(l),:)^T, one ZGEMV each. The stale part is a strided
!> read of a column-major matrix -- Ndim cache lines for Ndim*16 bytes -- which
!> is exactly the traffic the delay exists to amortise.
!-------------------------------------------------------------------------------

   subroutine delay_row(nf, GR, P, d, rows)
      implicit none
      integer, intent(in) :: nf, d, P(d)
      Complex (Kind=Kind(0.d0)), intent(in)  :: GR(ndim_s, ndim_s, nfl_s)
      Complex (Kind=Kind(0.d0)), intent(out) :: rows(ndim_s, d)
      ! One row of X gathered off the panel, the ZGEMV's coefficient vector.
      Complex (Kind=Kind(0.d0)) :: tmp(max(ncol(nf),1))
      integer :: i      ! column of G, i.e. position along the row
      integer :: l, c   ! l indexes the support, as in rows(:,l) = G(P(l),:)
      c = ncol(nf)
      do l = 1, d
         do i = 1, ndim_s
            rows(i,l) = GR(P(l), i, nf)
         enddo
         if (c > 0) then
            tmp(1:c) = xp(P(l), 1:c, nf)
            call ZGEMV('N', ndim_s, c, ZONE, yp(1,1,nf), ndim_s, tmp, 1, &
            &          ZONE, rows(1,l), 1)
         endif
      enddo
   end subroutine delay_row

!-------------------------------------------------------------------------------
!> @brief
!> d columns of the current Green's function: cols(i,l) = G(i, P(l)).
!-------------------------------------------------------------------------------

   subroutine delay_col(nf, GR, P, d, cols)
      implicit none
      integer, intent(in) :: nf, d, P(d)
      Complex (Kind=Kind(0.d0)), intent(in)  :: GR(ndim_s, ndim_s, nfl_s)
      Complex (Kind=Kind(0.d0)), intent(out) :: cols(ndim_s, d)
      ! One row of Y gathered off the panel, the ZGEMV's coefficient vector.
      Complex (Kind=Kind(0.d0)) :: tmp(max(ncol(nf),1))
      integer :: l, c   ! l indexes the support, as in cols(:,l) = G(:,P(l))
      c = ncol(nf)
      do l = 1, d
         call ZCOPY(ndim_s, GR(1, P(l), nf), 1, cols(1,l), 1)
         if (c > 0) then
            tmp(1:c) = yp(P(l), 1:c, nf)
            call ZGEMV('N', ndim_s, c, ZONE, xp(1,1,nf), ndim_s, tmp, 1, &
            &          ZONE, cols(1,l), 1)
         endif
      enddo
   end subroutine delay_col

!-------------------------------------------------------------------------------
!> @brief
!> Apply a rank-d update by appending d column pairs; flush when full.
!>
!> @details
!> The caller's update is G += alpha * xcols * ycols^T. The coefficient is
!> folded into the X column here so the flush stays a plain X*Y^T.
!-------------------------------------------------------------------------------

   subroutine delay_append(nf, alpha, xcols, ycols, d, GR)
      implicit none
      integer, intent(in) :: nf, d
      Complex (Kind=Kind(0.d0)), intent(in)    :: alpha   ! update coefficient
      ! The two factors of the caller's rank-d update, one column pair per rank.
      Complex (Kind=Kind(0.d0)), intent(in)    :: xcols(ndim_s, d)
      Complex (Kind=Kind(0.d0)), intent(in)    :: ycols(ndim_s, d)
      Complex (Kind=Kind(0.d0)), intent(inout) :: GR(ndim_s, ndim_s, nfl_s)
      integer :: l      ! which of the d column pairs is being appended
      integer :: c      ! columns already live, so c+l is where l lands
      c = ncol(nf)
      do l = 1, d
         call ZCOPY(ndim_s, xcols(1,l), 1, xp(1, c+l, nf), 1)
         call ZSCAL(ndim_s, alpha, xp(1, c+l, nf), 1)
         call ZCOPY(ndim_s, ycols(1,l), 1, yp(1, c+l, nf), 1)
      enddo
      ncol(nf) = c + d
      if (ncol(nf) >= kmax) call delay_flush(nf, GR)
   end subroutine delay_append

!-------------------------------------------------------------------------------
!> @brief
!> Conjugate the panels with a vertex operator, mirroring Op_Wrapup/Op_Wrapdo.
!>
!> @details
!> A thin pass-through to Op_Wrap_panels, which lives in Operator_mod beside the
!> routines it mirrors so that an edit to one is visible from the other.
!-------------------------------------------------------------------------------

   subroutine delay_wrap(nf, Op, HS_Field, N_Type, nt, updo)
      implicit none
      ! Op, HS_Field, N_Type and nt are the vertex, its field, the wrap variant
      ! and the time slice, exactly as Op_Wrapup and Op_Wrapdo take them; updo
      ! picks which of the two this call mirrors.
      integer, intent(in) :: nf, N_Type, nt
      Type (Operator), intent(in) :: Op
      Complex (Kind=Kind(0.d0)), intent(in) :: HS_Field
      character(len=1), intent(in) :: updo
      if (.not. delay_active) return
      if (ncol(nf) <= 0) return
      call Op_Wrap_panels(xp(1,1,nf), yp(1,1,nf), Op, HS_Field, ndim_s, &
      &                ncol(nf), N_Type, nt, updo)
   end subroutine delay_wrap

end module delayed_update_mod
