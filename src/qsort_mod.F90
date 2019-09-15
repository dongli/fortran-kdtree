module qsort_mod

  implicit none

  private

  public qsort

contains

  recursive subroutine qsort(x, left_idx_, right_idx_)

    real(8), intent(inout) :: x(:)
    integer, intent(in), optional :: left_idx_
    integer, intent(in), optional :: right_idx_

    integer part_idx
    integer left_idx, right_idx

    left_idx  = merge(left_idx_, 1, present(left_idx_))
    right_idx = merge(right_idx_, size(x), present(right_idx_))

    ! Array x should be already sorted.
    if (left_idx >= right_idx) return

    ! Partition the array so that the left elements are smaller than the right ones.
    part_idx = partition(x, left_idx, right_idx)

    call qsort(x, left_idx, part_idx)
    call qsort(x, part_idx + 1, right_idx)

  end subroutine qsort

  integer function partition(x, left_idx, right_idx) result(part_idx)

    real(8), intent(inout) :: x(:)
    integer, intent(inout) :: left_idx
    integer, intent(inout) :: right_idx

    real(8) tmp
    integer pivot_idx, i, j

    pivot_idx = left_idx
    j = -1
    do i = left_idx + 1, right_idx
      if (x(i) < x(pivot_idx) .and. j /= -1) then
        ! Swap the small element with the first large element.
        tmp = x(i); x(i) = x(j); x(j) = tmp
        ! Shift the recorded first large element.
        j = j + 1
      else if (x(i) >= x(pivot_idx) .and. j == -1) then
        ! Record the first element that is larger than or equal to pivot, called first large element.
        j = i
      end if
    end do
    if (j == -1) j = right_idx + 1
    j = j - 1
    if (j /= pivot_idx) then
      ! Swap pivot with the last small element.
      tmp = x(pivot_idx); x(pivot_idx) = x(j); x(j) = tmp
      part_idx = j - 1
    else
      part_idx = pivot_idx
    end if

  end function partition

  subroutine debug_print(x, left_idx, right_idx, part_idx)

    real(8), intent(in) :: x(:)
    integer, intent(in) :: left_idx
    integer, intent(in) :: right_idx
    integer, intent(in), optional :: part_idx

    integer i

    if (size(x) <= 20) then
      if (present(part_idx)) then
        do i = 1, size(x)
          if (i >= left_idx .and. i <= right_idx) then
            write(*, '(F8.2)', advance='no') x(i)
          else
            write(*, '("        ")', advance='no')
          end if
        end do
        write(*, *)
        do i = 1, part_idx - 1
          write(*, '("        ")', advance='no')
        end do
        write(*, '("       ^")')
      else
        do i = 1, size(x)
          write(*, '(A)', advance='no') '--------'
        end do
        write(*, *)
        do i = 1, size(x)
          write(*, '(I8)', advance='no') i
        end do
        write(*, *)
        if (left_idx >= right_idx) then
          do i = 1, size(x)
            write(*, '(F8.2)', advance='no') x(i)
          end do
        else
          do i = 1, size(x)
            if (i >= left_idx .and. i <= right_idx) then
              write(*, '(F8.2)', advance='no') x(i)
            else
              write(*, '("        ")', advance='no')
            end if
          end do
        end if
        write(*, *)
      end if
    end if

  end subroutine debug_print

end module qsort_mod
