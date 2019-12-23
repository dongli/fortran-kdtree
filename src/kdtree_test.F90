program kdtree_test

  use unit_test
  use node_placing_mod
  use kdtree

  implicit none

  integer num_seed
  integer, allocatable :: seed(:)
  type(kdtree_type) tree

  type(test_suite_type) test_suite

  call test_suite_init('KD-Tree test')

  ! Initialize random seeds.
  call random_seed(num_seed)
  allocate(seed(num_seed))
  seed = 2
  call random_seed(put=seed)

  call test_case_create('Test 1D 1 of 10')
  call test_1d(10, 1)
  call test_case_create('Test 1D 2 of 10')
  call test_1d(10, 2)
  call test_case_create('Test 1D 3 of 10')
  call test_1d(10, 3)
  call test_case_create('Test 1D 4 of 10')
  call test_1d(10, 4)
  call test_case_create('Test 1D 5 of 10')
  call test_1d(10, 5)
  call test_case_create('Test 1D 6 of 10')
  call test_1d(10, 6)
  call test_case_create('Test 1D 7 of 10')
  call test_1d(10, 7)
  
  call test_case_create('Test 2D')

  call test_2d()
  call test_2d_range_search()

  call test_suite_report()

  call test_suite_final()

  deallocate(seed)

contains

  subroutine test_1d(num_point, num_ngb)

    integer, intent(in) :: num_point
    integer, intent(in) :: num_ngb

    real(8), allocatable :: x(:,:)
    real(8), allocatable :: ngb_dist(:)
    integer, allocatable :: ngb_idx(:)
    integer i

    allocate(x(1,num_point))
    allocate(ngb_dist(num_ngb))
    allocate(ngb_idx(num_ngb))

    do i = 1, size(x, 2)
      call random_number(x(1,i))
    end do

    call tree%build(x)
    call tree%search([0.5d0], ngb_idx, ngb_dist_=ngb_dist)

    do i = 1, size(x, 2)
      if (all(ngb_idx /= i)) then
        call assert_true(all(norm2(x(:,i) - [0.5d0]) > ngb_dist), __FILE__, __LINE__)
      end if
    end do

    deallocate(ngb_dist)
    deallocate(ngb_idx)

  end subroutine test_1d

  subroutine test_2d()

    real(8), allocatable :: x(:,:)
    real(8) ngb_dist(10)
    integer i, j, ngb_idx(10), fail_count

    call node_placing([0.0d0,1.0d0,0.0d0,1.0d0], radius, x)

    call tree%build(x)
    call tree%search([0.5d0,0.5d0], ngb_idx, ngb_dist_=ngb_dist)

    do i = 1, size(x, 2)
      if (all(ngb_idx /= i)) then
        call assert_true(all(norm2(x(:,i) - [0.5d0,0.5d0]) > ngb_dist), __FILE__, __LINE__)
      end if
    end do

    fail_count = 0
    do j = 1, size(x, 2)
      call tree%search(x(:,j), ngb_idx, ngb_dist_=ngb_dist)
      do i = 1, size(x, 2)
        if (all(ngb_idx /= i)) then
          if (.not. all(norm2(x(:,i) - x(:,j)) > ngb_dist)) then
            fail_count = fail_count + 1
          end if
        end if
      end do
    end do
    call assert_equal(fail_count, 0)

  end subroutine test_2d

  subroutine test_2d_range_search()

    real(8), allocatable :: x(:,:)
    real(8) :: x0(2) = [0.5d0,0.5d0], r = 0.02d0
    integer i, fail_count
    integer, allocatable :: ngb_idx(:)

    call node_placing([0.0d0,1.0d0,0.0d0,1.0d0], radius, x)

    call tree%build(x)
    call tree%range_search(x0, ngb_idx, r)

    fail_count = 0
    do i = 1, size(x, 2)
      if (norm2(x(:,i) - x0) < r .and. .not. any(ngb_idx == i)) then
        fail_count = fail_count + 1
      end if
    end do
    call assert_equal(fail_count, 0)

    deallocate(ngb_idx)

  end subroutine test_2d_range_search

  real(8) function radius(x)

    real(8), intent(in) :: x(2)

    radius = 0.008d0

  end function radius

end program kdtree_test
