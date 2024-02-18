module node_placing_mod

  ! This module implement the "node-placing" algorithm to generate 2D grids.
  !
  ! Reference:
  !
  !   - Fornberg, Bengt and Flyer, Natasha, 2015: Fast generation of 2-D node
  !     distributions for mesh-free PDE discretizations. Computer & Mathematics
  !     with Applications, 69, 531-544.

  implicit none

  private

  public node_placing

  interface
    real(8) function node_radius_interface(x)
      real(8), intent(in) :: x(2)
    end function node_radius_interface
  end interface

  real(8), parameter :: pi = atan(1.0d0) * 4.0d0
  integer, parameter :: adj_num_pdp = 5 ! Adjust PDP number 

contains

  subroutine node_placing(box, radius, xy, init_num_node)

    real(8), intent(in) :: box(4)
    procedure(node_radius_interface) radius
    real(8), intent(out), allocatable :: xy(:,:)
    integer, intent(in), optional :: init_num_node

    integer init_num_node_opt ! Initial guessed node number
    integer num_node          ! Node number
    integer num_pdp           ! PDP number
    integer i, im(1), nw
    integer outside_count     ! Counter for PDPs that are outside box
    integer idx_left          ! Leftmost index of PDP within the radius of a grid
    integer idx_right         ! Rightmost index of PDP within the radius of a grid
    real(8) rand, dx
    real(8) r                 ! Radius or resolution of grid
    real(8) ym(1)             ! Minimum y in PDPs
    real(8) d                 ! Distance between placed node and PDPs
    real(8) angle_left
    real(8) angle_right
    real(8) angle
    real(8), allocatable :: pdp(:,:)

    if (allocated(xy)) deallocate(xy)

    ! Set the initial node number.
    init_num_node_opt = 10000; if (present(init_num_node)) init_num_node_opt = init_num_node

    num_node = 0
    num_pdp = init_num_node_opt
    allocate(pdp(2,init_num_node_opt)) ! Allocate more memory to accommodate increasing PDPs.
    allocate(xy (2,init_num_node_opt))

    ! Place initial PDPs along bottom.
    dx = (box(2) - box(1)) / num_pdp
    do i = 1, num_pdp
      pdp(1,i) = box(1) + dx * (i - 0.5)
      call random_number(rand)
      pdp(2,i) = box(3) + 1.0d-4 * rand
    end do

    ! Find the lowest PDP.
    ym = minval(pdp(2,:num_pdp)); im = minloc(pdp(2,:num_pdp))
    do while (ym(1) <= box(4))
      num_node = num_node + 1
      ! Assign the lowest PDP as a node.
      xy(:,num_node) = pdp(:,im(1))
      r = radius(xy(:,num_node))
      idx_left = 0; idx_right = 0
      do i = 1, num_pdp
        d = norm2(pdp(:,i) - xy(:,num_node))
        if (d <= r) then
          if (idx_left == 0) idx_left = i
          idx_right = i
        end if
      end do
      if (idx_left >= im(1)) idx_left = im(1) - 1
      if (idx_left < 1) then
        idx_left = im(1)
        angle_left = pi
      else
        angle_left = atan2(pdp(2,idx_left) - xy(2,num_node), pdp(1,idx_left) - xy(1,num_node))
      end if
      if (idx_right <= im(1)) idx_right = im(1) + 1
      if (idx_right > num_pdp) then
        idx_right = im(1)
        angle_right = 0.0d0
      else
        angle_right = atan2(pdp(2,idx_right) - xy(2,num_node), pdp(1,idx_right) - xy(1,num_node))
      end if
      nw = idx_right - idx_left + 1
      if (nw > adj_num_pdp) then
        ! There are more PDPs contained in the circle.
        do i = idx_left + adj_num_pdp, num_pdp - nw + adj_num_pdp
          pdp(:,i) = pdp(:,i+nw-adj_num_pdp)
        end do
      else if (nw < adj_num_pdp) then
        ! There are less PDPs contained in the circle.
        do i = num_pdp - nw + adj_num_pdp, idx_left + adj_num_pdp, -1
          pdp(:,i) = pdp(:,i+nw-adj_num_pdp)
        end do
      end if
      ! Push PDPs toward top.
      do i = 0, adj_num_pdp - 1
        angle = angle_left - (2 * i + 1) * 0.1 * (angle_left - angle_right)
        pdp(1,idx_left+i) = xy(1,num_node) + r * cos(angle)
        pdp(2,idx_left+i) = xy(2,num_node) + r * sin(angle)
      end do
      num_pdp = num_pdp - nw + adj_num_pdp
      ! Remove PDPs outside box.
      outside_count = count(pdp(1,:) < box(1) .or. pdp(1,:) > box(2))
      if (outside_count > 0) then
        num_pdp = num_pdp - outside_count
        pdp(2,:) = pack(pdp(2,:), pdp(1,:) >= box(1) .and. pdp(1,:) <= box(2))
        pdp(1,:) = pack(pdp(1,:), pdp(1,:) >= box(1) .and. pdp(1,:) <= box(2))
      end if
      ! Find the next lowest PDP.
      ym = minval(pdp(2,:num_pdp)); im = minloc(pdp(2,:num_pdp))
      ! Enlarge array size if necessary.
      if (num_node == init_num_node_opt) then
        init_num_node_opt = init_num_node_opt * 2
        call resize_array(xy, dim=2, new_size=init_num_node_opt)
      end if
    end do

    ! Clean zeros from output array.
    call resize_array(xy, dim=2, new_size=num_node)

    deallocate(pdp)

  end subroutine node_placing

  subroutine resize_array(x, dim, new_size)

    real(8), intent(inout), allocatable :: x(:,:)
    integer, intent(in) :: dim
    integer, intent(in) :: new_size

    real(8), allocatable :: y(:,:)
    integer nx(2), ny(2), i

    nx = shape(x)
    ny = nx
    ny(dim) = new_size

    allocate(y(ny(1),ny(2)))

    if (nx(dim) <= ny(dim)) then
      y(:nx(1),:nx(2)) = x(:nx(1),:nx(2))
    else
      y(:ny(1),:ny(2)) = x(:ny(1),:ny(2))
    end if

    deallocate(x)

    allocate(x(ny(1),ny(2)))

    x = y

    deallocate(y)

  end subroutine resize_array

  subroutine debug_write(num_pdp, pdp, num_node, node)

    integer, intent(in) :: num_pdp
    real(8), intent(in) :: pdp(:,:)
    integer, intent(in) :: num_node
    real(8), intent(in) :: node(:,:)

    character(10) tag
    integer i

    write(tag, '(I0)') num_node

    open(10, file='node_placing.pdp.' // trim(tag) // '.txt')
    do i = 1, num_pdp
      write(10, *) pdp(:,i)
    end do
    close(10)

    open(10, file='node_placing.node.' // trim(tag) // '.txt')
    do i = 1, num_node
      write(10, *) node(:,i)
    end do
    close(10)

  end subroutine debug_write

end module node_placing_mod
