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

  recursive subroutine kdtree_build(this, x, start_node_)

    class(kdtree_type), intent(inout) :: this
    real(8), intent(in) :: x(:,:)
    type(node_type), intent(inout), target, optional :: start_node_

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
    if (present(start_node_)) then
      node => start_node_
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
    call node%create_child_nodes(part_dim, num_dim, num_point)

    do i = 1, num_point
      if (x(part_dim,i) < xmed) then
        call node%left %add_point(x(:,i), node%global_idx_array(i))
      else if (x(part_dim,i) > xmed) then
        call node%right%add_point(x(:,i), node%global_idx_array(i))
      else
        ! Save the cut point.
        node%x = x(:,i)
        node%global_idx = node%global_idx_array(i)
      end if
    end do

    ! write(*, '("===== ", I0)') node%id
    ! do i = 1, num_point
    !   write(*, '(F6.2)', advance='no') x(1,i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! if (mod(i, 20) /= 1) write(*, *)
    ! write(*, '("-----")')
    ! write(*, '("xmed = ", F6.2)') xmed
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
    if (node%left%num_point > 0) then
      call this%build(node%left%x_array, start_node_=node%left)
    else
      deallocate(node%left)
    end if
    if (node%right%num_point > 0) then
      call this%build(node%right%x_array, start_node_=node%right)
    else
      deallocate(node%right)
    end if

  end subroutine kdtree_build

  recursive subroutine kdtree_search(this, x, ngb_idx, start_node_, ngb_dist_, ngb_count_)

    class(kdtree_type), intent(in) :: this
    real(8), intent(in) :: x(:)
    integer, intent(inout) :: ngb_idx(:)
    type(node_type), intent(in), target, optional :: start_node_
    real(8), intent(inout), target, optional :: ngb_dist_(:)
    integer, intent(inout), target, optional :: ngb_count_

    real(8) dist, radius
    integer i, j
    logical replaced, should_enter

    type(node_type), pointer :: node
    real(8), pointer :: ngb_dist(:)
    integer, pointer :: ngb_count

    if (present(start_node_)) then
      node => start_node_
      ngb_dist => ngb_dist_
      ngb_count => ngb_count_
    else
      node => this%root_node
      allocate(ngb_dist(size(ngb_idx)))
      allocate(ngb_count)
      ngb_count = 0
    end if
    dist = norm2(x - node%x)
    ! This acts as a priority queue.
    replaced = .false.
    do i = 1, ngb_count
      if (dist < ngb_dist(i)) then
        ngb_count = min(ngb_count + 1, size(ngb_idx))
        do j = ngb_count, i + 1, -1
          ngb_dist(j) = ngb_dist(j-1)
          ngb_idx (j) = ngb_idx (j-1)
        end do
        ngb_dist(i) = dist
        ngb_idx (i) = node%global_idx
        replaced = .true.
        exit
      end if
    end do
    if (node%part_dim /= 0) then
      ! Leaf node does not have part_dim.
      radius = abs(x(node%part_dim) - node%x(node%part_dim))
    end if
    if (.not. replaced .and. ngb_count < size(ngb_idx)) then
      ngb_count           = ngb_count + 1
      ngb_dist(ngb_count) = dist
      ngb_idx (ngb_count) = node%global_idx
    end if

    ! write(*, '("===== ", I0, X, I0)') node%id, node%global_idx
    ! do i = 1, ngb_count
    !   write(*, '(I8)', advance='no') ngb_idx(i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! write(*, *)
    ! do i = 1, ngb_count
    !   write(*, '(F8.4)', advance='no') ngb_dist(i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! write(*, *)
    ! pause

    if (associated(node%left)) then
      should_enter = .false.
      if (associated(node%right)) then
        if (norm2(x - node%right%x) > radius) then
          should_enter = .true.
        end if
      end if
      if (should_enter .or. x(node%part_dim) < node%x(node%part_dim)) then
          print *, 'Enter left node ', node%left%id
        call this%search(x, ngb_idx, start_node_=node%left, ngb_dist_=ngb_dist, ngb_count_=ngb_count)
      end if
    end if
    if (associated(node%right)) then
      should_enter = .false.
      if (associated(node%left)) then
        if (norm2(x - node%left%x) > radius) then
          should_enter = .true.
        end if
      end if
      if (should_enter .or. x(node%part_dim) > node%x(node%part_dim)) then
        print *, 'Enter right node ', node%right%id
        call this%search(x, ngb_idx, start_node_=node%right, ngb_dist_=ngb_dist, ngb_count_=ngb_count)
      end if
    end if

    if (.not. present(ngb_dist_)) then
      deallocate(ngb_dist)
      deallocate(ngb_count)
    end if

  end subroutine kdtree_search

  subroutine kdtree_final(this)

    type(kdtree_type), intent(inout) :: this

    if (associated(this%root_node)) deallocate(this%root_node)

  end subroutine kdtree_final

end module kdtree_mod
