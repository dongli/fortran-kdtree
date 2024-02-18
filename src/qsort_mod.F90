module qsort_mod

  implicit none

  private

  public qsort

  interface qsort
    module procedure qsort_r4
    module procedure qsort_r8
  end interface qsort

contains

  recursive subroutine qsort_r4(x, left_idx, right_idx)

    real(4), intent(inout) :: x(:)
    integer, intent(in), optional :: left_idx
    integer, intent(in), optional :: right_idx


    integer left_idx_opt, right_idx_opt, i, j
    real(4) tmp, x0

    left_idx_opt  = 1      ; if (present(left_idx )) left_idx_opt  = left_idx
    right_idx_opt = size(x); if (present(right_idx)) right_idx_opt = right_idx

    x0 = x((left_idx_opt + right_idx_opt) / 2)
    i = left_idx_opt; j = right_idx_opt
    do
      do while (x(i) < x0)
        i = i + 1
      end do
      do while (x(j) > x0)
        j = j - 1
      end do
      if (i >= j) exit
      tmp = x(i); x(i) = x(j); x(j) = tmp
      i = i + 1; j = j - 1
    end do

    if (left_idx_opt < i - 1)  call qsort(x, left_idx_opt, i - 1)
    if (j + 1 < right_idx_opt) call qsort(x, j + 1, right_idx_opt)

  end subroutine qsort_r4

  recursive subroutine qsort_r8(x, left_idx, right_idx)

    real(8), intent(inout) :: x(:)
    integer, intent(in), optional :: left_idx
    integer, intent(in), optional :: right_idx


    integer left_idx_opt, right_idx_opt, i, j
    real(8) tmp, x0

    left_idx_opt  = 1      ; if (present(left_idx )) left_idx_opt  = left_idx
    right_idx_opt = size(x); if (present(right_idx)) right_idx_opt = right_idx

    x0 = x((left_idx_opt + right_idx_opt) / 2)
    i = left_idx_opt; j = right_idx_opt
    do
      do while (x(i) < x0)
        i = i + 1
      end do
      do while (x(j) > x0)
        j = j - 1
      end do
      if (i >= j) exit
      tmp = x(i); x(i) = x(j); x(j) = tmp
      i = i + 1; j = j - 1
    end do

    if (left_idx_opt < i - 1)  call qsort(x, left_idx_opt, i - 1)
    if (j + 1 < right_idx_opt) call qsort(x, j + 1, right_idx_opt)

  end subroutine qsort_r8

end module qsort_mod
