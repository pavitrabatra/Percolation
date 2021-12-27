!!! Go through logic of code once again
!!! Figure out why exactly interface was neededed here but was not needed in previous codes(Read up on modules,contains,interface,subroutine)
!!! Read up on intent
!!!Read up on implementation using modulea

program bernoulli_percolation
implicit none
real*8::p,r
integer::L,i,j,S
integer,dimension(:,:),allocatable:: matrix

interface
        integer function hk(matrix) !result(hk)
                integer,dimension(:,:),allocatable:: matrix
                integer,dimension(:),allocatable :: labels, new_labels !For Labeling used in hk
                integer :: m,n,i,j
                integer :: up, left, label
                end function
end interface

print*, 'Enter the size of the lattice'
read*,L

print*, 'Probability of Occupation'
read*,p

allocate(matrix(L,L))
S=L*L


!initializing lattice
open(10,file='initial_lattice.dat')
do i=1,L
 do j=1,L
  matrix(i,j)=0    
  write(10,*) float(matrix(i,j))
 end do 
end do
close(10)


!Getting a Microstate
open(11, file='Random Microstate.dat')
do i=1,L
 do j=1,L
  call random_number(r)
  if(r>=(1-p))then
          matrix(i,j)=1
  end if
  write(11,*) float(matrix(i,j))
 end do
end do
close(11)

!Identifying Clusters using hk
print*, "The total number of disjoint clusters labelled."
print*, hk(matrix)
end program 

!hk
integer function hk(matrix) !result(hk)
         integer,dimension(:,:),allocatable:: matrix
         integer,dimension(:),allocatable :: labels, new_labels !For Labeling used in hk
         integer :: m,n,i,j
         integer :: up, left, label
 
 INTERFACE
         INTEGER FUNCTION  find(x, labels)
               integer, intent(in) :: x
               integer, dimension(:), intent(inout) :: labels
         END FUNCTION  find
         
         integer function union(x, y, labels)
                 integer, intent(in) :: x, y
                 integer, dimension(:), intent(inout) :: labels
         end function
 END INTERFACE

 m = size(matrix, 1)
 n = size(matrix, 2)
 allocate(labels(m*n/2+1))
 labels = 0
 hk = 0

 do j = 1, n
  do i = 1, m
   if(matrix(i,j) > 0) then
         if(i == 1) then
                 up = 0
         else
                 up = matrix(i-1, j)
         end if
         if(j == 1) then
                 left = 0
         else
                 left = matrix(i, j-1)
         end if
         ! New cluster
         if(up==0 .and. left==0) then
                 hk = hk + 1
                 labels(hk) = hk
                 matrix(i,j) = hk

         ! Site binds clusters
         else if(up > 0 .and. left > 0) then
                 matrix(i,j) = union(up, left, labels)

         ! Only one neighbour 
         else
                 matrix(i,j) = max(up, left)
         end if
   end if
  end do
 end do
 
 allocate(new_labels(m*n/2+1))
 new_labels = 0
 hk = 0

 do j = 1, n
  do i = 1, m
   if(matrix(i,j) > 0) then
          label = find(matrix(i,j), labels)
          if(new_labels(label) == 0) then
                  hk = hk + 1
                  new_labels(label) = hk
          end if
          matrix(i,j) = new_labels(label)
   end if
  end do
 end do
end function hk


!! Union-Find find algorithm:
!! Find the lowest corresponding label.
!! Relabelling is done when necessary.
integer function find(x, labels) !result(y) !!!Figure out if can be done without using interface but using result function
            implicit none
            integer, intent(in) :: x
                !! Label for which to find the lowest corresponding label.
            integer, dimension(:), intent(inout) :: labels
                !! List of labels. **labels(i)** points to the lowest
                !! corresponding label of label **i**.
            integer :: y, z, tmp
            
            find = x
            do while(labels(find) /= find)
                find = labels(find)
            end do

            tmp = x
            do while(labels(tmp) /= tmp)
                z = labels(tmp)
                labels(tmp) = find
                tmp = z
            end do
end function find

!! Union-Find union algorithm:
!! Merge two labels, return the result.
integer function union(x, y, labels) !result(canonical_label) !!!Figure out if can be done without using interface but using result function
    implicit none
    integer, intent(in) :: x, y
    !! Labels to merge.
    integer, dimension(:), intent(inout) :: labels
    !! List of labels.
    INTERFACE
      INTEGER FUNCTION  find(x, labels)
               integer, intent(in) :: x
               integer, dimension(:), intent(inout) :: labels
      END FUNCTION  find
    END INTERFACE
    union = find(y,labels)
    labels(find(x,labels)) = union
end function union
!Identifying Clusters-hk
