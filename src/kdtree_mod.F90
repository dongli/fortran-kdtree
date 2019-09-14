module kdtree_mod

  use utils_mod
  use node_mod

  implicit none

  type kdtree_type
    type(node_type), pointer :: root_node => null()
  contains
    procedure :: build => kdtree_build
    procedure :: search => kdtree_search
    final :: kdtree_final
  end type kdtree_type

contains

  recursive subroutine kdtree_build(this, x, root_node)

    class(kdtree_type), intent(inout) :: this
    real(8), intent(in) :: x(:,:)
    type(node_type), intent(inout), target, optional :: root_node

    integer num_point, num_dim, part_dim, i
    real(8) xvar(3), max_xvar
    real(8) xmed ! Median coordinate along one dimension
    type(node_type), pointer :: node

    num_dim   = size(x, 1)
    num_point = size(x, 2)

    ! Calculate the variance along each dimension, and choose the dimension with the
    ! largest value as the partition dimension.
    max_xvar = -1
    do i = 1, num_dim
      xvar(i) = variance(x(i,:))
      if (max_xvar < xvar(i)) then
        max_xvar = xvar(i)
        part_dim = i
      end if
    end do

    ! Calculate the median number along the chosen dimension as the partition location.
    xmed = median(x(part_dim,:))

    ! Create tree structures.
    if (present(root_node)) then
      node => root_node
    else
      allocate(this%root_node)
      call this%root_node%init()
      node => this%root_node
      ! Create global index array.
      allocate(node%global_idx_array(num_point))
      do i = 1, num_point
        node%global_idx_array(i) = i
      end do
    end if
    if (num_point == 1) then
      node%x = x(:,1)
      node%global_idx = node%global_idx_array(1)
      call node%discard_arrays()
      ! write(*, '("## leaf node ", I0)') node%id
      ! write(*, '("## point = ", F6.2, X I0)') node%x, node%global_idx
      ! pause
      return
    end if
    call node%create_child_nodes(part_dim, xmed, num_dim, num_point)

    do i = 1, num_point
      if (x(part_dim,i) <= xmed) then
        call node%left %add_point(x(:,i), node%global_idx_array(i))
      else
        call node%right%add_point(x(:,i), node%global_idx_array(i))
      end if
      if (x(part_dim,i) == xmed) then
        node%x = x(:,i)
        node%global_idx = i
      end if
    end do

    ! write(*, '("===== ", I0)') node%id
    ! do i = 1, num_point
    !   write(*, '(F6.2)', advance='no') x(1,i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! if (mod(i, 20) /= 1) write(*, *)
    ! write(*, '("-----")')
    ! write(*, '("xmed = ", F6.2)') node%xmed
    ! write(*, '("cut point = ", F6.2, X, I0)') node%x(1), node%global_idx
    ! write(*, '("left = ", I0)') node%left%id
    ! do i = 1, node%left%num_point
    !   write(*, '(I6)', advance='no') node%left%global_idx_array(i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! if (mod(i, 20) /= 1) write(*, *)
    ! do i = 1, node%left%num_point
    !   write(*, '(F6.2)', advance='no') node%left%x_array(1,i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! if (mod(i, 20) /= 1) write(*, *)
    ! write(*, '("right = ", I0)') node%right%id
    ! do i = 1, node%right%num_point
    !   write(*, '(I6)', advance='no') node%right%global_idx_array(i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! if (mod(i, 20) /= 1) write(*, *)
    ! do i = 1, node%right%num_point
    !   write(*, '(F6.2)', advance='no') node%right%x_array(1,i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! if (mod(i, 20) /= 1) write(*, *)
    ! pause

    ! Clean memory usage.
    call node%discard_arrays()
    call node%left %end_point()
    call node%right%end_point()

    ! Recursively build subtrees.
    call this%build(node%left %x_array, root_node=node%left )
    call this%build(node%right%x_array, root_node=node%right)

  end subroutine kdtree_build

  recursive subroutine kdtree_search(this, x, ngb_idx, root_node)

    ! Firstly, find an initial search path to the leaf node, and calculate the
    ! distance between the node point and the query point.

    class(kdtree_type), intent(in) :: this
    real(8), intent(in) :: x(:)
    integer, intent(in) :: ngb_idx(:)
    type(node_type), intent(in), target, optional :: root_node

    type(node_type), pointer :: node

    if (present(root_node)) then
      node => root_node
    else
      node => this%root_node
    end if
    if (node%num_point == 1) then

    else
      if (x(node%part_dim) <= node%xmed) then
        call this%search(x, ngb_idx, root_node=node%left)
      else
        call this%search(x, ngb_idx, root_node=node%right)
      end if
    end if

  end subroutine kdtree_search

  subroutine kdtree_final(this)

    type(kdtree_type), intent(inout) :: this

    if (associated(this%root_node)) deallocate(this%root_node)

  end subroutine kdtree_final

end module kdtree_mod
