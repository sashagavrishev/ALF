! The probe's decision must follow its rules on any cost curve.
!
! delay_probe times kernels and hands the resulting curve to delay_probe_pick.
! The timings cannot be asserted on -- they belong to the machine -- so the
! decision is tested on its own here, against made-up curves over the candidate
! depths [8, 16, 32, 64, 128, 256]:
!
!   a clear minimum         the largest candidate within the 5% margin of it
!   a flat curve            the closed-form fallback, sqrt(2*Ndim) clamped
!   untimed candidates      ignored, the decision taken on the rest
!   nothing timed           the fallback
!   a zero reading          a failed clock, so the fallback

Program DelayedProbePick

   Use delayed_update_mod

   Implicit None

   Real (Kind=Kind(0.D0)), Parameter :: UNTIMED = huge(1.d0)
   Integer :: nfail

   nfail = 0

   Call check("clear minimum", [4.d0, 2.d0, 1.d0, 1.04d0, 1.5d0, 3.d0], &
      &       256, 64, 'probe')
   Call check("flat curve", [1.d0, 1.01d0, 1.02d0, 1.d0, 1.03d0, 1.01d0], &
      &       256, formula(256), 'formula')
   Call check("untimed above Ndim", [2.d0, 1.d0, UNTIMED, UNTIMED, UNTIMED, UNTIMED], &
      &       16, 16, 'probe')
   Call check("nothing timed", [UNTIMED, UNTIMED, UNTIMED, UNTIMED, UNTIMED, UNTIMED], &
      &       256, formula(256), 'formula')
   Call check("zero reading", [0.d0, 2.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
      &       256, formula(256), 'formula')

   If (nfail > 0) Then
      Write (*,*) "FAILURES:", nfail
      Stop 2
   End If

   Write (*,*) "SUCCESS"

Contains

   Subroutine check(name, cost, Ndim, k_expected, source_expected)
      Character (Len=*),      Intent(In) :: name, source_expected
      Real (Kind=Kind(0.D0)), Intent(In) :: cost(:)
      Integer,                Intent(In) :: Ndim, k_expected
      Character (Len=16) :: source
      Integer :: k

      k = delay_probe_pick(cost, Ndim, source)
      If (k /= k_expected .or. trim(source) /= source_expected) Then
         Write (*,*) "ERROR ", name, ": got k =", k, " from ", trim(source), &
            &        "; expected k =", k_expected, " from ", source_expected
         nfail = nfail + 1
      End If
   End Subroutine check

   ! The closed form, restated rather than called: nint(sqrt(2*Ndim)) in [8, 256].
   Integer Function formula(Ndim)
      Integer, Intent(In) :: Ndim
      formula = min(256, max(8, nint(sqrt(2.d0*real(Ndim, Kind(0.D0))))))
   End Function formula

End Program DelayedProbePick