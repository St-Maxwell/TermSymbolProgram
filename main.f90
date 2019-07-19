program main
  use m_shell
  use m_link
  use m_set
  use m_spectra
  use m_combination, only: long
  implicit none
!-----------------------------------------------------------
  integer :: num_shell
  character(len=1), dimension(:), allocatable :: shell_type
  integer, dimension(:), allocatable :: shell_occ
!-----------------------------------------------------------
  type(link), pointer :: head => null()
  type(link), pointer :: tail => null()
  type(link), pointer :: ptr => null()
  class(shell), pointer :: shell_ptr => null()
!-----------------------------------------------------------
  integer(kind=long), dimension(:), allocatable :: shells_num_comb
! the number of combinations should be less than 10^12
  integer(kind=long) :: tot_num_comb
  integer, dimension(:), allocatable :: L_all_comb
  real, dimension(:), allocatable :: S_all_comb
  integer, dimension(:), allocatable :: L_all_comb_sorted
  real, dimension(:), allocatable :: S_all_comb_sorted
  integer, dimension(:), allocatable :: L_shell_comb
  real, dimension(:), allocatable :: S_shell_comb
  integer, dimension(:), allocatable :: L_set
  real, dimension(:), allocatable :: S_set
  integer, dimension(:,:), allocatable :: table
!-----------------------------------------------------------
  integer :: jter
  integer(kind=long) :: kter
  integer(kind=long) :: iter
  integer(kind=long) :: lter
  integer :: iL, iS
!-----------------------------------------------------------
  write(*,*) "------------------------------------------------------"
  write(*,*) "    This program aims for determining term symbols    "
  write(*,*) "               Programmed by St Maxwell               "
  write(*,*) "                    2019-Jul-19th                     "
  write(*,*) "------------------------------------------------------"
  write(*,*)
  write(*,*) "Input the number of shells:"
  read(*,*) num_shell

  allocate(shell_type(num_shell))
  allocate(shell_occ(num_shell))
  allocate(shells_num_comb(num_shell))

  write(*,*) "Input shell symbols (e.g. s p p d):"
  read(*,*) shell_type

  write(*,*) "Input occupation numbers for each shell (e.g. 1 3 1 2):"
  read(*,*) shell_occ

  ! build a linked list
  ! each link receives a shell
  do iter = 1, num_shell
  
    ! creat shell object
    shell_ptr => shell_factory(shell_type(iter), shell_occ(iter))
    shells_num_comb(iter) = shell_ptr%num_comb

    if (.not. associated(head)) then
      allocate(head)
      tail => head
      tail%next => null()
      tail%variable => shell_ptr
    else
      allocate(tail%next)
      tail => tail%next
      tail%next => null()
      tail%variable => shell_ptr
    end if

  end do

  ! obtain the total number of microstates
  ! and allocate working arrays
  tot_num_comb = 1
  do iter = 1, num_shell
    tot_num_comb = tot_num_comb * shells_num_comb(iter)
  end do
  allocate(L_all_comb(tot_num_comb))
  allocate(S_all_comb(tot_num_comb))
  allocate(L_all_comb_sorted(tot_num_comb))
  allocate(S_all_comb_sorted(tot_num_comb))

  L_all_comb = 0
  S_all_comb = 0
  jter = 1
  ptr => head

  do while ( associated(ptr) )
    ! find microstates
    shell_ptr => ptr%variable

    call shell_ptr%get_combination()
    L_shell_comb = shell_ptr%comb(:)%L
    S_shell_comb = shell_ptr%comb(:)%S
    
    ! combinate microstates of each shell
    ! to get all total microstates
    !-------------------------------------------------------
    kter = 1
    do iter = jter+1, num_shell
      kter = kter * shells_num_comb(iter)
    end do
    jter = jter + 1

    lter = 1
    do iter = 1, tot_num_comb, kter

      L_all_comb(iter:iter+kter-1) = &
        L_all_comb(iter:iter+kter-1) + L_shell_comb(lter)
      S_all_comb(iter:iter+kter-1) = &
        S_all_comb(iter:iter+kter-1) + S_shell_comb(lter)

      if (lter == size(L_shell_comb)) then
        lter = 1
        cycle
      end if
      lter = lter + 1

    end do
    !-------------------------------------------------------

    ptr => ptr%next

  end do

  ! obtain all possible values of L & S 
  L_all_comb_sorted = L_all_comb
  S_all_comb_sorted = S_all_comb
  call sort_L(L_all_comb_sorted)
  call sort_S(S_all_comb_sorted)
  L_set = set_of_L(L_all_comb_sorted)
  S_set = set_of_S(S_all_comb_sorted)
  deallocate(L_all_comb_sorted)
  deallocate(S_all_comb_sorted)

  ! build the table
  ! reference: Journal of Chemical Education, 52(2), 1975:87-89
  allocate(table(size(S_set), size(L_set)))
  table = 0
  do iter = 1, tot_num_comb
    do iL = 1, size(L_set)
      if (L_set(iL) == L_all_comb(iter)) exit
    end do
    do iS = 1, size(S_set)
      if (abs(S_set(iS)-S_all_comb(iter)) < tol) exit
    end do
    table(iS, iL) = table(iS, iL) + 1
  end do

  call print_table(table, L_set, S_set)
  call print_term_symbol(table, size(L_set), size(S_set))

  write(*,*) "Finished. Press Enter to exit."
  read(*,*)
  
end program