module utils_mod

  use qsort_mod

  implicit none

  interface median
    module procedure median_r4
    module procedure median_r8
  end interface median

  interface variance
    module procedure variance_r4
    module procedure variance_r8
  end interface variance

contains

  real(4) function median_r4(x) result(res)

    real(4), intent(in) :: x(:)

    ! NOTE: We do not deallocate tmp, let OS withdraw the memory.
    real(4), allocatable, save :: tmp(:)

    if (allocated(tmp)) then
      if (size(tmp) < size(x)) then
        deallocate(tmp)
        allocate(tmp(size(x)))
      end if
    else
      allocate(tmp(size(x)))
    end if
    tmp(:size(x)) = x; call qsort(tmp(:size(x)))
    ! NOTE: Plus 1 to avoid 0 index.
    res = tmp(int(size(x) * 0.5) + 1)

  end function median_r4

  real(8) function median_r8(x) result(res)

    real(8), intent(in) :: x(:)

    ! NOTE: We do not deallocate tmp, let OS withdraw the memory.
    real(8), allocatable, save :: tmp(:)

    if (allocated(tmp)) then
      if (size(tmp) < size(x)) then
        deallocate(tmp)
        allocate(tmp(size(x)))
      end if
    else
      allocate(tmp(size(x)))
    end if
    tmp(:size(x)) = x; call qsort(tmp(:size(x)))
    ! NOTE: Plus 1 to avoid 0 index.
    res = tmp(int(size(x) * 0.5) + 1)

  end function median_r8

  real(4) function variance_r4(x) result(res)

    real(4), intent(in) :: x(:)

    real(4) xa

    xa = sum(x) / size(x)
    res = sum((x - xa)**2) / size(x)

  end function variance_r4

  real(8) function variance_r8(x) result(res)

    real(8), intent(in) :: x(:)

    real(8) xa

    xa = sum(x) / size(x)
    res = sum((x - xa)**2) / size(x)

  end function variance_r8

end module utils_mod
