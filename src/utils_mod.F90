module utils_mod

  use qsort_mod

  implicit none

contains

  real(8) function median(x) result(res)

    real(8), intent(in) :: x(:)

    real(8) xc(size(x))

    xc = x; call qsort(xc)
    res = xc(int(size(x) * 0.5))

  end function median

  real(8) function variance(x) result(res)

    real(8), intent(in) :: x(:)

    real(8) xa

    xa = sum(x) / size(x)
    res = sum((x - xa)**2) / size(x)

  end function variance

end module utils_mod