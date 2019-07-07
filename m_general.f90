module m_combination
  implicit none
  ! n < 10^12
  integer, parameter :: long = selected_int_kind(12)
  contains

  function factorial(n)
    implicit none
    integer(kind=long) :: factorial
    integer, intent(in) :: n
    integer :: iter

    if (n == 0) factorial = 1
     
    factorial = 1
    do iter = 1, n
      factorial = factorial * iter
    end do

  end function

  function combinations(all,pick)
    implicit none
    integer(kind=long) :: combinations
    integer, intent(in) :: all, pick

    combinations = factorial(all) / &
                    (factorial(pick) * factorial(all - pick))
    
  end function

  subroutine find_all_combination_wrapper(k, n, ncomb, comb)
  ! find all combinations of size n from k numbers
  ! e.g. k = 4, n = 2
  !      all combinations: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
  !      comb = | 1 1 1 2 2 3 |
  !             | 2 3 4 3 4 4 |
    implicit none
    integer, intent(in) :: k, n
    integer(kind=long), intent(in) :: ncomb
    integer, dimension(ncomb,n) :: comb
    !-------------------------------------
    integer, dimension(n) :: tmp
    integer :: counter_comb
    integer :: counter

    counter_comb = 1
    counter = 1
    call find_all_combination(k, n, 1)

    contains
    recursive subroutine find_all_combination(k_in, n_in, left)
      implicit none
      integer :: k_in, n_in, left
      integer :: i
    
      if (n_in==0) then
        comb(counter_comb,:) = tmp
        counter_comb = counter_comb + 1
      else
        do i = left, k_in
          tmp(counter) = i
          counter = counter + 1
          call find_all_combination(k_in,n_in-1,i+1)
          counter = counter - 1
        end do
      end if
    
    end subroutine

  end subroutine

end module

module m_set
  use m_combination, only: long
  implicit none
  real, parameter :: tol = 1.E-4
  contains

    subroutine sort_L(array)
    ! insertion sort
      implicit none
      integer, dimension(:), intent(inout) :: array
      integer :: tmp
      integer(kind=long) :: iter, jter

      do iter = 2, size(array)
        tmp = array(iter)
        jter = iter - 1
        do while(jter >= 1)
          if (array(jter) >= tmp) exit
          array(jter+1) = array(jter)
          jter = jter - 1
        end do
        array(jter+1) = tmp
      end do

    end subroutine

    subroutine sort_S(array)
      implicit none
      real, dimension(:), intent(inout) :: array
      real :: tmp
      integer(kind=long) :: iter, jter

      do iter = 2, size(array)
        tmp = array(iter)
        jter = iter - 1
        do while(jter >= 1)
          if (array(jter) >= (tmp-tol)) exit
          array(jter+1) = array(jter)
          jter = jter - 1
        end do
        array(jter+1) = tmp
      end do
    
    end subroutine

    subroutine set_L(arr, set, size_set)
    ! remove duplicates from sorted array
    ! return a set
      implicit none
      integer, dimension(:), intent(in) :: arr
      integer, dimension(size(arr)), intent(out) :: set
      integer, intent(out) :: size_set
      integer :: iter
      integer(kind=long) :: jter

      if (size(arr) == 1) then
        set = arr
        size_set = 1
        return
      end if

      set(1) = arr(1)
      iter = 1
      do jter = 2, size(arr)
        if (arr(jter) /= set(iter)) then
          iter = iter + 1
          set(iter) = arr(jter)
        end if
      end do
      size_set = iter

    end subroutine

    subroutine set_S(arr, set, size_set)
      implicit none
      real, dimension(:), intent(in) :: arr
      real, dimension(size(arr)), intent(out) :: set
      integer, intent(out) :: size_set
      integer :: iter
      integer(kind=long) :: jter
  
      if (size(arr) == 1) then
        set = arr
        size_set = 1
        return
      end if
  
      set(1) = arr(1)
      iter = 1
      do jter = 2, size(arr)
        if (abs(arr(jter)-set(iter)) > tol) then
          iter = iter + 1
          set(iter) = arr(jter)
        end if
      end do
      size_set = iter
  
    end subroutine

end module