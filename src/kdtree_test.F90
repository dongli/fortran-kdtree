program kdtree_test

  use kdtree_mod

  implicit none

  type(kdtree_type) kdtree
  real(8) x(1,10)
  integer i, ngb(7)

  integer n
  integer, allocatable :: seed(:)

  call random_seed(n)
  allocate(seed(n)); seed = 2
  call random_seed(put=seed)
  do i = 1, size(x, 2)
    call random_number(x(1,i))
    print '(I0, X, F6.4)', i, x(1,i)
  end do

  call kdtree%build(x)

  call kdtree%search([0.5d0], ngb)

end program kdtree_test