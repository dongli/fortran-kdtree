module node_mod

  implicit none

  private

  public node_type

  type node_type
    integer :: id = 1
    integer :: part_dim = 0
    real(8), allocatable :: x(:)
    integer global_idx
    integer :: num_point = 0
    ! The following two arrays are used at build time, after that they will be descarded.
    real(8), allocatable :: x_array(:,:)
    integer, allocatable :: global_idx_array(:)
    type(node_type), pointer :: parent  => null()
    type(node_type), pointer :: brother => null()
    type(node_type), pointer :: left    => null()
    type(node_type), pointer :: right   => null()
  contains
    procedure :: init => node_init
    procedure :: create_child_nodes => node_create_child_nodes
    procedure :: add_point => node_add_point
    procedure :: end_point => node_end_point
    procedure :: discard_arrays => node_discard_arrays
    final :: node_final
  end type node_type

contains

  subroutine node_init(this, num_dim, max_num_point, parent, brother)

    class(node_type), intent(inout) :: this
    integer, intent(in), optional :: num_dim
    integer, intent(in), optional :: max_num_point
    type(node_type), intent(in), target, optional :: parent
    type(node_type), intent(in), pointer, optional :: brother

    if (allocated(this%x               )) deallocate(this%x               )
    if (allocated(this%x_array         )) deallocate(this%x_array         )
    if (allocated(this%global_idx_array)) deallocate(this%global_idx_array)
    if (present(num_dim) .and. present(max_num_point)) then
      allocate(this%x               (num_dim              ))
      allocate(this%x_array         (num_dim,max_num_point))
      allocate(this%global_idx_array(        max_num_point))
    end if
    if (present(parent )) this%parent  => parent
    if (present(brother)) this%brother => brother

  end subroutine node_init

  subroutine node_create_child_nodes(this, part_dim, num_dim, max_num_point)

    class(node_type), intent(inout) :: this
    integer, intent(in) :: part_dim
    integer, intent(in) :: num_dim
    integer, intent(in) :: max_num_point

    this%part_dim = part_dim

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
    if (this%num_point > size(this%x_array, 2)) then
      write(*, '("[Error]: ", A, ":", I0, ":", A)') __FILE__, __LINE__, 'Array size is not sufficient!'
      stop 1
    end if
    this%x_array(:,this%num_point)        = x
    this%global_idx_array(this%num_point) = global_idx

  end subroutine node_add_point

  subroutine node_end_point(this)

    class(node_type), intent(inout) :: this

    real(8), allocatable :: rtmp(:,:)
    integer, allocatable :: itmp(:)
    integer m, n, i, j

    if (this%num_point > 0) then
      m = size(this%x_array, 1)
      n = this%num_point

      allocate(rtmp(m,n))
      allocate(itmp(  n))

      do j = 1, n
        do i = 1, m
          rtmp(i,j) = this%x_array(i,j)
        end do
        itmp(j) = this%global_idx_array(j)
      end do

      deallocate(this%x_array         )
      deallocate(this%global_idx_array)

      allocate(this%x_array         (m,n))
      allocate(this%global_idx_array(  n))

      this%x_array          = rtmp
      this%global_idx_array = itmp

      deallocate(rtmp)
      deallocate(itmp)
    end if

  end subroutine node_end_point

  subroutine node_discard_arrays(this)

    ! Discard arrays since they are no longer needed.

    class(node_type), intent(inout) :: this

    if (allocated(this%x_array         )) deallocate(this%x_array         )
    if (allocated(this%global_idx_array)) deallocate(this%global_idx_array)

  end subroutine node_discard_arrays

  subroutine node_final(this)

    type(node_type), intent(inout) :: this

    if (allocated(this%x               )) deallocate(this%x               )
    if (allocated(this%x_array         )) deallocate(this%x_array         )
    if (allocated(this%global_idx_array)) deallocate(this%global_idx_array)
    if (associated(this%left           )) deallocate(this%left            )
    if (associated(this%right          )) deallocate(this%right           )

  end subroutine node_final

end module node_mod
