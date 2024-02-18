module kdtree_mod

  use kdtree_utils_mod
  use node_mod

  implicit none

  private

  public kdtree_type

  type kdtree_type
    type(node_type), pointer :: root_node => null()
  contains
    procedure, private :: kdtree_build_1
    procedure, private :: kdtree_build_2
    generic :: build => kdtree_build_1, kdtree_build_2
    procedure, private :: kdtree_search_r4
    procedure, private :: kdtree_search_r8
    generic :: search => kdtree_search_r4, kdtree_search_r8
    procedure :: range_search => kdtree_range_search
    final :: kdtree_final
  end type kdtree_type

  integer node_access_count

contains

  recursive subroutine kdtree_build_1(this, x, start_node)

    class(kdtree_type), intent(inout) :: this
    real(8), intent(in) :: x(:,:)
    type(node_type), intent(inout), target, optional :: start_node

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
    if (present(start_node)) then
      node => start_node
    else
      allocate(this%root_node)
      call this%root_node%init()
      node => this%root_node
      node%num_point = num_point
      ! Create global index array.
      allocate(node%global_idx_array(num_point))
      do i = 1, num_point
        node%global_idx_array(i) = i
      end do
    end if
    if (num_point == 1) then
      ! Reach the leaf node, return back.
      node%x = x(:,1)
      node%global_idx = node%global_idx_array(1)
      call node%discard_arrays()
      return
    end if
    call node%create_child_nodes(part_dim, num_dim, num_point)

    node%global_idx = -1
    do i = 1, num_point
      if (x(part_dim,i) < xmed) then
        call node%left %add_point(x(:,i), node%global_idx_array(i))
      else if (x(part_dim,i) > xmed) then
        call node%right%add_point(x(:,i), node%global_idx_array(i))
      else
        ! If there is already a point on the cut line, save it to both sides.
        if (node%global_idx /= -1) then
          call node%left %add_point(node%x, node%global_idx)
          call node%right%add_point(node%x, node%global_idx)
        end if
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
      call this%build(node%left%x_array, start_node=node%left)
    else
      deallocate(node%left)
    end if
    if (node%right%num_point > 0) then
      call this%build(node%right%x_array, start_node=node%right)
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

  recursive subroutine kdtree_search_r4(this, x, ngb_idx, mute, start_node, ngb_dist, ngb_count)

    class(kdtree_type), intent(in) :: this
    real(4), intent(in) :: x(:)
    integer, intent(inout) :: ngb_idx(:)
    logical, intent(in), optional :: mute
    type(node_type), intent(in), target, optional :: start_node
    real(4), intent(inout), target, optional :: ngb_dist(:)
    integer, intent(inout), target, optional :: ngb_count

    logical mute_opt
    real(4) dist
    type(node_type), pointer :: node
    real(4), pointer :: ngb_dist_opt(:)
    integer, pointer :: ngb_count_opt

    mute_opt = .true.; if (present(mute)) mute_opt = mute
    if (present(start_node)) then
      node => start_node
    else
      node => this%root_node
    end if
    if (present(ngb_count)) then
      ngb_count_opt => ngb_count
    else
      allocate(ngb_count_opt)
      ngb_count_opt = 0
      node_access_count = 0
    end if
    if (present(ngb_dist)) then
      ngb_dist_opt => ngb_dist
    else
      allocate(ngb_dist_opt(size(ngb_idx)))
    end if

    node_access_count = node_access_count + 1
    dist = norm2(x - node%x)

    call record_potential_ngb_r4(ngb_count_opt, ngb_idx, ngb_dist_opt, dist, node%global_idx)

    ! Check if the hypersphere with the radius as the distance of the farest neightbor intersects
    ! the splitting hyperplane.
    if (node%part_dim /= 0) then
      dist = abs(x(node%part_dim) - node%x(node%part_dim))
      if (x(node%part_dim) < node%x(node%part_dim)) then
        if (associated(node%left)) then
          call this%search(x, ngb_idx, start_node=node%left, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
        end if
        if (associated(node%right)) then
          if (dist < ngb_dist_opt(ngb_count_opt)) then
            call this%search(x, ngb_idx, start_node=node%right, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
          end if
        end if
      else
        if (associated(node%right)) then
          call this%search(x, ngb_idx, start_node=node%right, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
        end if
        if (associated(node%left)) then
          if (dist < ngb_dist_opt(ngb_count_opt)) then
            call this%search(x, ngb_idx, start_node=node%left, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
          end if
        end if
      end if
    end if

    if (.not. present(ngb_dist )) deallocate(ngb_dist_opt)
    if (.not. present(ngb_count)) deallocate(ngb_count_opt)
    if (.not. present(start_node) .and. .not. mute_opt) then
      write(*, '("Searched ", I0, " of ", I0, " nodes.")') node_access_count, this%root_node%num_point
    end if

  end subroutine kdtree_search_r4

  recursive subroutine kdtree_search_r8(this, x, ngb_idx, mute, start_node, ngb_dist, ngb_count)

    class(kdtree_type), intent(in) :: this
    real(8), intent(in) :: x(:)
    integer, intent(inout) :: ngb_idx(:)
    logical, intent(in), optional :: mute
    type(node_type), intent(in), target, optional :: start_node
    real(8), intent(inout), target, optional :: ngb_dist(:)
    integer, intent(inout), target, optional :: ngb_count

    logical mute_opt
    real(8) dist
    type(node_type), pointer :: node
    real(8), pointer :: ngb_dist_opt(:)
    integer, pointer :: ngb_count_opt

    mute_opt = .true.; if (present(mute)) mute_opt = mute
    if (present(start_node)) then
      node => start_node
    else
      node => this%root_node
    end if
    if (present(ngb_count)) then
      ngb_count_opt => ngb_count
    else
      allocate(ngb_count_opt)
      ngb_count_opt = 0
      node_access_count = 0
    end if
    if (present(ngb_dist)) then
      ngb_dist_opt => ngb_dist
    else
      allocate(ngb_dist_opt(size(ngb_idx)))
    end if

    node_access_count = node_access_count + 1
    dist = norm2(x - node%x)

    call record_potential_ngb_r8(ngb_count_opt, ngb_idx, ngb_dist_opt, dist, node%global_idx)

    ! Check if the hypersphere with the radius as the distance of the farest neightbor intersects
    ! the splitting hyperplane.
    if (node%part_dim /= 0) then
      dist = abs(x(node%part_dim) - node%x(node%part_dim))
      if (x(node%part_dim) < node%x(node%part_dim)) then
        if (associated(node%left)) then
          call this%search(x, ngb_idx, start_node=node%left, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
        end if
        if (associated(node%right)) then
          if (dist < ngb_dist_opt(ngb_count_opt)) then
            call this%search(x, ngb_idx, start_node=node%right, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
          end if
        end if
      else
        if (associated(node%right)) then
          call this%search(x, ngb_idx, start_node=node%right, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
        end if
        if (associated(node%left)) then
          if (dist < ngb_dist_opt(ngb_count_opt)) then
            call this%search(x, ngb_idx, start_node=node%left, ngb_dist=ngb_dist_opt, ngb_count=ngb_count_opt)
          end if
        end if
      end if
    end if

    if (.not. present(ngb_dist )) deallocate(ngb_dist_opt)
    if (.not. present(ngb_count)) deallocate(ngb_count_opt)
    if (.not. present(start_node) .and. .not. mute_opt) then
      write(*, '("Searched ", I0, " of ", I0, " nodes.")') node_access_count, this%root_node%num_point
    end if

  end subroutine kdtree_search_r8

  subroutine kdtree_range_search(this, x, ngb_idx, radius, ngb_dist, ngb_bunch)

    class(kdtree_type), intent(in) :: this
    real(8), intent(in) :: x(:)
    integer, intent(inout), allocatable :: ngb_idx(:)
    real(8), intent(in) :: radius
    real(8), intent(inout), allocatable, optional :: ngb_dist(:)
    integer, intent(in), optional :: ngb_bunch

    integer nb, n, i
    logical finished
    real(8), allocatable :: ngb_dist_opt(:)
    integer, allocatable :: final_ngb_idx(:)

    if (present(ngb_bunch)) then
      nb = ngb_bunch
    else
      nb = 50
    end if

    if (.not. allocated(ngb_idx)) allocate(ngb_idx(nb))
    allocate(ngb_dist_opt(size(ngb_idx)))

    finished = .false.
    do while (.not. finished)
      call this%search(x, ngb_idx, mute=.true., ngb_dist=ngb_dist_opt)
      do i = 1, size(ngb_idx)
        if (ngb_dist_opt(i) > radius) then
          finished = .true.
          exit
        end if
      end do
      if (.not. finished) then
        ! Increase search size by nb.
        n = size(ngb_idx) + nb
        deallocate(ngb_idx ); allocate(ngb_idx (n))
        deallocate(ngb_dist_opt); allocate(ngb_dist_opt(n))
      end if
    end do

    ! Remove ngb outside range.
    allocate(final_ngb_idx(size(ngb_idx)))
    n = 0
    do i = 1, size(ngb_idx)
      if (ngb_dist_opt(i) <= radius) then
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
      do i = 1, size(ngb_dist_opt)
        if (ngb_dist_opt(i) <= radius) then
          n = n + 1
          ngb_dist(n) = ngb_dist_opt(i)
        end if
      end do
    end if
    deallocate(ngb_dist_opt)

  end subroutine kdtree_range_search

  subroutine kdtree_final(this)

    type(kdtree_type), intent(inout) :: this

    if (associated(this%root_node)) deallocate(this%root_node)

  end subroutine kdtree_final

  subroutine record_potential_ngb_r4(ngb_count_opt, ngb_idx, ngb_dist, dist, global_idx)

    integer, intent(inout) :: ngb_count_opt
    integer, intent(inout) :: ngb_idx(:)
    real(4), intent(inout) :: ngb_dist(:)
    real(4), intent(in) :: dist
    integer, intent(in) :: global_idx

    integer i, j
    logical replaced

    ! This acts as a priority queue.
    replaced = .false.
    do i = 1, ngb_count_opt
      if (ngb_idx(i) == global_idx) return
      if (dist < ngb_dist(i)) then
        ngb_count_opt = min(ngb_count_opt + 1, size(ngb_idx))
        do j = ngb_count_opt, i + 1, -1
          ngb_dist(j) = ngb_dist(j-1)
          ngb_idx (j) = ngb_idx (j-1)
        end do
        ngb_dist(i) = dist
        ngb_idx (i) = global_idx
        replaced    = .true.
        exit
      end if
    end do
    if (.not. replaced .and. ngb_count_opt < size(ngb_idx)) then
      ngb_count_opt           = ngb_count_opt + 1
      ngb_dist(ngb_count_opt) = dist
      ngb_idx (ngb_count_opt) = global_idx
    end if

    ! write(*, '("===== ", I0, X, I0)') global_idx
    ! do i = 1, ngb_count_opt
    !   write(*, '(I8)', advance='no') ngb_idx(i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! write(*, *)
    ! do i = 1, ngb_count_opt
    !   write(*, '(F8.4)', advance='no') ngb_dist(i)
    !   if (mod(i, 20) == 0) write(*, *)
    ! end do
    ! write(*, *)
    ! pause

  end subroutine record_potential_ngb_r4

  subroutine record_potential_ngb_r8(ngb_count_opt, ngb_idx, ngb_dist, dist, global_idx)

    integer, intent(inout) :: ngb_count_opt
    integer, intent(inout) :: ngb_idx(:)
    real(8), intent(inout) :: ngb_dist(:)
    real(8), intent(in) :: dist
    integer, intent(in) :: global_idx

    integer i, j
    logical replaced

    ! This acts as a priority queue.
    replaced = .false.
    do i = 1, ngb_count_opt
      if (ngb_idx(i) == global_idx) return
      if (dist < ngb_dist(i)) then
        ngb_count_opt = min(ngb_count_opt + 1, size(ngb_idx))
        do j = ngb_count_opt, i + 1, -1
          ngb_dist(j) = ngb_dist(j-1)
          ngb_idx (j) = ngb_idx (j-1)
        end do
        ngb_dist(i) = dist
        ngb_idx (i) = global_idx
        replaced    = .true.
        exit
      end if
    end do
    if (.not. replaced .and. ngb_count_opt < size(ngb_idx)) then
      ngb_count_opt           = ngb_count_opt + 1
      ngb_dist(ngb_count_opt) = dist
      ngb_idx (ngb_count_opt) = global_idx
    end if

  end subroutine record_potential_ngb_r8

end module kdtree_mod
