# Introduction

KD-Tree is a standard data structure for indexing data, especially in 3D space. It is an extension of Binary-Space Partition (BSP) to more than one dimension. For more information on KD-Tree, please refer Wiki.

This repository is a Fortran implementation of KD-Tree. We need KD-Tree in scientific HPC scenarios, so Fortran is still the right language, and also we need more modern library interfaces.

# Usage

```fortran
  use kdtree_mod

  real(8), allocatable :: x(:,:)       ! num_dim, num_point
  integer, allocatable :: ngb_idx(:,:) ! num_ngb, num_point
  type(kdtree_type) kdtree
  integer i

  ! Allocate x and ngb_idx, and set x accordingly.

  call kdtree%build(x)
  do i = 1, size(x, 2)
    call kdtree%search(x(:,i), ngb_idx(:,i))
  end do

  ! Other works ...
```

# Contributors

- Li Dong <dongli@lasg.iap.ac.cn>