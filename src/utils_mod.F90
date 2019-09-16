module utils_mod

  use qsort_mod

  implicit none

  ! NOTE: We do not deallocate tmp, let OS withdraw the memory.
  real(8), allocatable :: tmp(:)

contains

  real(8) function median(x) result(res)

    real(8), intent(in) :: x(:)

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

  end function median

  real(8) function variance(x) result(res)

    real(8), intent(in) :: x(:)

    real(8) xa

    xa = sum(x) / size(x)
    res = sum((x - xa)**2) / size(x)

  end function variance

end module utils_mod