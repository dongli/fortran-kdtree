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

  subroutine node_placing(box, target_num_node, radius, xy)

    real(8), intent(in) :: box(4)
    integer, intent(in) :: target_num_node
    procedure(node_radius_interface) radius
    real(8), intent(out), allocatable :: xy(:,:)

    integer num_node
    integer num_pdp       ! PDP number
    integer i, im(1), nw
    integer outside_count ! Counter for PDPs that are outside box
    integer idx_left      ! Leftmost index of PDP within the radius of a grid
    integer idx_right     ! Rightmost index of PDP within the radius of a grid
    real(8) rand, dx
    real(8) r             ! Radius or resolution of grid
    real(8) ym(1)         ! Minimum y in PDPs
    real(8) d             ! Distance between placed node and PDPs
    real(8) angle_left
    real(8) angle_right
    real(8) angle
    real(8), allocatable :: pdp(:,:)

    if (allocated(xy)) deallocate(xy)

    num_node = 0
    num_pdp = int(sqrt(dble(target_num_node)))
    allocate(pdp(2,int(0.5 * target_num_node))) ! Allocate more memory to accommodate increasing PDPs.
    allocate(xy(2,target_num_node))

    ! Place initial PDPs along bottom.
    dx = (box(2) - box(1)) / num_pdp
    do i = 1, num_pdp
      pdp(1,i) = box(1) + dx * (i - 0.5)
      call random_number(rand)
      pdp(2,i) = box(3) + 1.0d-4 * rand
      ! print *, i, pdp(:,i)
    end do

    ! Find the lowest PDP.
    ym = minval(pdp(2,:num_pdp)); im = minloc(pdp(2,:num_pdp))
    do while (ym(1) <= box(4) .and. num_node < target_num_node)
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
      ! print *, '## idx_left =', idx_left
      ! print *, '## idx_right =', idx_right
      ! print *, '## nw =', nw
      ! print *, '## angle_left =', angle_left
      ! print *, '## angle_right =', angle_right
      if (nw > adj_num_pdp) then
        ! There are more PDPs contained in the circle.
        do i = idx_left + adj_num_pdp, num_pdp - nw + adj_num_pdp
          ! print *, i+nw-adj_num_pdp, '->', i
          pdp(:,i) = pdp(:,i+nw-adj_num_pdp)
        end do
      else if (nw < adj_num_pdp) then
        ! There are less PDPs contained in the circle.
        do i = num_pdp - nw + adj_num_pdp, idx_left + adj_num_pdp, -1
          ! print *, i+nw-adj_num_pdp, '->', i
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
      ! write(*, *) '-----'
      ! write(*, *) '+ node ', num_node, xy(:,num_node)
      ! write(*, *) 'num_pdp =', num_pdp
      ! do i = 1, num_pdp
      !   print *, i, pdp(:,i)
      ! end do
      ! call debug_write(num_pdp, pdp, num_node, xy)
      ! pause
    end do
    call debug_write(num_pdp, pdp, num_node, xy)

    deallocate(pdp)

    ! Clean zeros from output array.
    allocate(pdp(2,num_node))
    pdp = xy(:,:num_node)
    deallocate(xy)
    allocate(xy(2,num_node))
    xy = pdp
    deallocate(pdp)

  end subroutine node_placing

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