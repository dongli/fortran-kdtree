program kdtree_test

  use node_placing_mod
  use kdtree_mod

  implicit none

  integer num_seed
  integer, allocatable :: seed(:)
  type(kdtree_type) kdtree

  call random_seed(num_seed)
  allocate(seed(num_seed))
  seed = 2
  call random_seed(put=seed)

  ! call test_1d()
  call test_2d()

  deallocate(seed)

contains

  subroutine test_1d()

    real(8) x(1,10)
    integer i, ngb_idx(4)

    do i = 1, size(x, 2)
      call random_number(x(1,i))
    end do

    call kdtree%build(x)
    call kdtree%search([0.5d0], ngb_idx)

    write(*, *) ngb_idx

  end subroutine test_1d

  subroutine test_2d()

    real(8), allocatable :: x(:,:)
    real(8) ngb_dist(10)
    integer i, ngb_idx(10)

    call node_placing([0.0d0,1.0d0,0.0d0,1.0d0], radius, x)

    call kdtree%build(x)
    call kdtree%search([0.5d0,0.5d0], ngb_idx, ngb_dist_=ngb_dist)

    do i = 1, size(x, 2)
      if (.not. any(ngb_idx == i)) then
        if (any(norm2(x(:,i) - [0.5d0,0.5d0]) <= ngb_dist)) then
          write(*, *) 'Failed to get neighbor ', i, x(:,i)
        end if
      end if
    end do

    do i = 1, size(ngb_idx)
      write(*, '(I0, ",")', advance='no') ngb_idx(i)
    end do
    write(*, *)

  end subroutine test_2d

  real(8) function radius(x)

    real(8), intent(in) :: x(2)

    if (x(2) < 0.5) then
      radius = 0.1d0
    else
      radius = 0.05d0
    end if

  end function radius

end program kdtree_test