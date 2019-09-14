module node_mod

  implicit none

  private

  public node_type

  type node_type
    integer :: id = 1
    integer :: part_dim = 0
    real(8) :: xmed
    integer :: num_point = 0
    real(8), allocatable :: x(:,:)
    integer, allocatable :: global_idx(:)
    type(node_type), pointer :: parent  => null()
    type(node_type), pointer :: brother => null()
    type(node_type), pointer :: left    => null()
    type(node_type), pointer :: right   => null()
  contains
    procedure :: init => node_init
    procedure :: create_child_nodes => node_create_child_nodes
    procedure :: add_point => node_add_point
    procedure :: end_point => node_end_point
    final :: node_final
  end type node_type

contains

  subroutine node_init(this, num_dim, max_num_point, parent, brother)

    class(node_type), intent(inout) :: this
    integer, intent(in), optional :: num_dim
    integer, intent(in), optional :: max_num_point
    type(node_type), intent(in), target, optional :: parent
    type(node_type), intent(in), pointer, optional :: brother

    if (allocated(this%x)) then
      deallocate(this%x)
      deallocate(this%global_idx)
    end if
    if (present(num_dim) .and. present(max_num_point)) then
      allocate(this%x(num_dim,max_num_point))
      allocate(this%global_idx(max_num_point))
    end if
    if (present(parent)) this%parent => parent
    if (present(brother)) this%brother => brother

  end subroutine node_init

  subroutine node_create_child_nodes(this, part_dim, xmed, num_dim, max_num_point)

    class(node_type), intent(inout) :: this
    integer, intent(in) :: part_dim
    real(8), intent(in) :: xmed
    integer, intent(in) :: num_dim
    integer, intent(in) :: max_num_point

    this%part_dim = part_dim
    this%xmed     = xmed

    if (associated(this%left )) deallocate(this%left )
    if (associated(this%right)) deallocate(this%right)
    allocate(this%left)
    allocate(this%right)

    this%left %id = this%id * 10
    this%right%id = this%id * 10 + 1

    call this%left %init(num_dim, max_num_point, this, this%right)
    call this%right%init(num_dim, max_num_point, this, this%left )

  end subroutine node_create_child_nodes

  subroutine node_add_point(this, x, global_idx)

    class(node_type), intent(inout) :: this
    real(8), intent(in) :: x(:)
    integer, intent(in) :: global_idx

    this%num_point = this%num_point + 1
    if (this%num_point > size(this%x, 2)) then
      write(*, '("[Error]: ", A, ":", I0, ":", A)') __FILE__, __LINE__, 'Array size is not sufficient!'
      stop 1
    end if
    this%x(:,this%num_point) = x
    this%global_idx(this%num_point) = global_idx

  end subroutine node_add_point

  subroutine node_end_point(this)

    class(node_type), intent(inout) :: this

    real(8), allocatable :: rtmp(:,:)
    integer, allocatable :: itmp(:)
    integer m, n, i, j

    m = size(this%x, 1)
    n = this%num_point

    allocate(rtmp(m,n))
    allocate(itmp(n))

    do j = 1, n
      do i = 1, m
        rtmp(i,j) = this%x(i,j)
      end do
      itmp(j) = this%global_idx(j)
    end do

    deallocate(this%x)
    deallocate(this%global_idx)

    allocate(this%x(m,n))
    allocate(this%global_idx(n))

    this%x = rtmp
    this%global_idx = itmp

    deallocate(rtmp)
    deallocate(itmp)

  end subroutine node_end_point

  subroutine node_final(this)

    type(node_type), intent(inout) :: this

    if (allocated(this%x)) then
      deallocate(this%x)
      deallocate(this%global_idx)
    end if
    if (associated(this%left )) deallocate(this%left )
    if (associated(this%right)) deallocate(this%right)

  end subroutine node_final

end module node_mod
