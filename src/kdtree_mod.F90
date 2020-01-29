module kdtree_mod

  use utils_mod
  use node_mod

  implicit none

  private

  public kdtree_type

  type kdtree_type
    type(node_type), pointer :: root_node => null()
  contains
    procedure :: kdtree_build_1
    procedure :: kdtree_build_2
    generic :: build => kdtree_build_1, kdtree_build_2
    procedure :: search => kdtree_search
    procedure :: range_search => kdtree_range_search
    final :: kdtree_final
  end type kdtree_type

  integer node_access_count

contains

  recursive subroutine kdtree_build_1(this, x, start_node_)

    class(kdtree_type), intent(inout) :: this
    real(8), intent(in) :: x(:,:)
    type(node_type), intent(inout), target, optional :: start_node_

    integer num_point, num_dim, part_dim, i, d
    real(8) xvar, max_xvar
    real(8) xmed ! Median coordinate along one dimension
    type(node_type), pointer :: node

    num_dim   = size(x, 1)
    num_point = size(x, 2)

    ! Calculate the variance along each dimension, and choose the dimension with the
    ! largest value as the partition dimension.
    max_xvar = -1
    do d = 1, num_dim
      xvar = variance(x(d,:))
      if (max_xvar < xvar) then
        max_xvar = xvar
        part_dim = d
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
      node%num_point = num_point
      ! Create global index array.
      allocate(node%global_idx_array(num_point))
      node%global_idx_array = [(i, i = 1, num_point)]
    end if
    if (num_point == 1) then
      ! Reach the leaf node, return back.
      node%x = x(:,1)
      node%global_idx = node%global_idx_array(1)
      call node%discard_arrays()
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

  end subroutine kdtree_build_1

  subroutine kdtree_build_2(this, x, y)

    class(kdtree_type), intent(inout) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(:)

    real(8), allocatable :: xy(:,:)

    if (size(x) /= size(y)) then
      stop '[Error]: kdtree_build: Dimensions of x and y is not the same!'
    end if

    allocate(xy(2,size(x)))

    xy(1,:) = x
    xy(2,:) = y

    call this%build(xy)

    deallocate(xy)

  end subroutine kdtree_build_2

  recursive subroutine kdtree_search(this, x, ngb_idx, mute, start_node_, ngb_dist_, ngb_count_)

    class(kdtree_type), intent(in) :: this
    real(8), intent(in) :: x(:)
    integer, intent(inout) :: ngb_idx(:)
    logical, intent(in), optional :: mute
    type(node_type), intent(in), target, optional :: start_node_
    real(8), intent(inout), target, optional :: ngb_dist_(:)
    integer, intent(inout), target, optional :: ngb_count_

    logical mute_
    real(8) dist
    type(node_type), pointer :: node
    real(8), pointer :: ngb_dist(:)
    integer, pointer :: ngb_count

    if (present(mute)) then
      mute_ = mute
    else
      mute_ = .true.
    end if
    if (present(start_node_)) then
      node => start_node_
    else
      node => this%root_node
    end if
    if (present(ngb_count_)) then
      ngb_count => ngb_count_
    else
      allocate(ngb_count)
      ngb_count = 0
      node_access_count = 0
    end if
    if (present(ngb_dist_)) then
      ngb_dist => ngb_dist_
    else
      allocate(ngb_dist(size(ngb_idx)))
    end if

    node_access_count = node_access_count + 1
    dist = norm2(x - node%x)

    call record_potential_ngb(ngb_count, ngb_idx, ngb_dist, dist, node%global_idx)

    ! Check if the hypersphere with the radius as the distance of the farest neightbor intersects
    ! the splitting hyperplane.
    if (node%part_dim /= 0) then
      dist = abs(x(node%part_dim) - node%x(node%part_dim))
      if (x(node%part_dim) < node%x(node%part_dim)) then
        if (associated(node%left)) then
          call this%search(x, ngb_idx, start_node_=node%left, ngb_dist_=ngb_dist, ngb_count_=ngb_count)
        end if
        if (associated(node%right)) then
          if (dist < ngb_dist(ngb_count)) then
            call this%search(x, ngb_idx, start_node_=node%right, ngb_dist_=ngb_dist, ngb_count_=ngb_count)
          end if
        end if
      else
        if (associated(node%right)) then
          call this%search(x, ngb_idx, start_node_=node%right, ngb_dist_=ngb_dist, ngb_count_=ngb_count)
        end if
        if (associated(node%left)) then
          if (dist < ngb_dist(ngb_count)) then
            call this%search(x, ngb_idx, start_node_=node%left, ngb_dist_=ngb_dist, ngb_count_=ngb_count)
          end if
        end if
      end if
    end if

    if (.not. present(ngb_dist_ )) deallocate(ngb_dist )
    if (.not. present(ngb_count_)) deallocate(ngb_count)
    if (.not. present(start_node_) .and. .not. mute_) then
      write(*, '("Searched ", I0, " of ", I0, " nodes.")') node_access_count, this%root_node%num_point
    end if

  end subroutine kdtree_search

  subroutine kdtree_range_search(this, x, ngb_idx, radius, ngb_dist, ngb_bunch)

    class(kdtree_type), intent(in) :: this
    real(8), intent(in) :: x(:)
    integer, intent(inout), allocatable :: ngb_idx(:)
    real(8), intent(in) :: radius
    real(8), intent(inout), allocatable, optional :: ngb_dist(:)
    integer, intent(in), optional :: ngb_bunch

    integer nb, n, i
    logical finished
    real(8), allocatable :: ngb_dist_(:)
    integer, allocatable :: final_ngb_idx(:)

    if (present(ngb_bunch)) then
      nb = ngb_bunch
    else
      nb = 50
    end if

    if (.not. allocated(ngb_idx)) allocate(ngb_idx(nb))
    allocate(ngb_dist_(size(ngb_idx)))

    ! recursive subroutine kdtree_search(this, x, ngb_idx, mute, start_node_, ngb_dist_, ngb_count_)
    finished = .false.
    do while (.not. finished)
      call this%search(x, ngb_idx, mute=.true., ngb_dist_=ngb_dist_)
      do i = 1, size(ngb_idx)
        if (ngb_dist_(i) > radius) then
          finished = .true.
          exit
        end if
      end do
      if (.not. finished) then
        ! Increase search size by nb.
        n = size(ngb_idx) + nb
        deallocate(ngb_idx  ); allocate(ngb_idx  (n))
        deallocate(ngb_dist_); allocate(ngb_dist_(n))
      end if
    end do

    ! Remove ngb outside range.
    allocate(final_ngb_idx(size(ngb_idx)))
    n = 0
    do i = 1, size(ngb_idx)
      if (ngb_dist_(i) <= radius) then
        n = n + 1
        final_ngb_idx(n) = ngb_idx(i)
      end if
    end do
    deallocate(ngb_idx); allocate(ngb_idx(n))
    do i = 1, n
      ngb_idx(i) = final_ngb_idx(i)
    end do
    deallocate(final_ngb_idx)
    if (present(ngb_dist)) then
      if (allocated(ngb_dist)) deallocate(ngb_dist)
      allocate(ngb_dist(n))
      n = 0
      do i = 1, size(ngb_dist_)
        if (ngb_dist_(i) <= radius) then
          n = n + 1
          ngb_dist(n) = ngb_dist_(i)
        end if
      end do
    end if
    deallocate(ngb_dist_)

  end subroutine kdtree_range_search

  subroutine kdtree_final(this)

    type(kdtree_type), intent(inout) :: this

    if (associated(this%root_node)) deallocate(this%root_node)

  end subroutine kdtree_final

  subroutine record_potential_ngb(ngb_count, ngb_idx, ngb_dist, dist, global_idx)

    integer, intent(inout) :: ngb_count
    integer, intent(inout) :: ngb_idx(:)
    real(8), intent(inout) :: ngb_dist(:)
    real(8), intent(in) :: dist
    integer, intent(in) :: global_idx

    integer i, j
    logical replaced

    ! This acts as a priority queue.
    replaced = .false.
    do i = 1, ngb_count
      if (ngb_idx(i) == global_idx) return
      if (dist < ngb_dist(i)) then
        ngb_count = min(ngb_count + 1, size(ngb_idx))
        do j = ngb_count, i + 1, -1
          ngb_dist(j) = ngb_dist(j-1)
          ngb_idx (j) = ngb_idx (j-1)
        end do
        ngb_dist(i) = dist
        ngb_idx (i) = global_idx
        replaced    = .true.
        exit
      end if
    end do
    if (.not. replaced .and. ngb_count < size(ngb_idx)) then
      ngb_count           = ngb_count + 1
      ngb_dist(ngb_count) = dist
      ngb_idx (ngb_count) = global_idx
    end if

    ! write(*, '("===== ", I0, X, I0)') global_idx
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

  end subroutine record_potential_ngb

end module kdtree_mod
