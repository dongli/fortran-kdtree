module qsort_mod

  implicit none

  private

  public qsort

contains

  recursive subroutine qsort(x, left_idx_, right_idx_)

    real(8), intent(inout) :: x(:)
    integer, intent(in), optional :: left_idx_
    integer, intent(in), optional :: right_idx_


    integer left_idx, right_idx, i, j
    real(8) tmp, x0

    left_idx  = merge(left_idx_, 1, present(left_idx_))
    right_idx = merge(right_idx_, size(x), present(right_idx_))

    x0 = x((left_idx + right_idx) / 2)
    i = left_idx; j = right_idx
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

    if (left_idx < i - 1)  call qsort(x, left_idx, i - 1)
    if (j + 1 < right_idx) call qsort(x, j + 1, right_idx)

  end subroutine qsort

end module qsort_mod
