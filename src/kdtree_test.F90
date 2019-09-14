program kdtree_test

  use kdtree_mod

  implicit none

  type(kdtree_type) kdtree
  real(8) x(1,10)
  integer i

  do i = 1, size(x, 2)
    call random_number(x(1,i))
  end do

  call kdtree%build(x)

end program kdtree_test