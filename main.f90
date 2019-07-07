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
  class(*), pointer :: tmp => null()
  type(s_shell), pointer :: s_define_p => null()
  type(p_shell), pointer :: p_define_p => null()
  type(d_shell), pointer :: d_define_p => null()
  type(f_shell), pointer :: f_define_p => null()
  type(g_shell), pointer :: g_define_p => null()
  type(s_shell), pointer :: s_tmp_p => null()
  type(p_shell), pointer :: p_tmp_p => null()
  type(d_shell), pointer :: d_tmp_p => null()
  type(f_shell), pointer :: f_tmp_p => null()
  type(g_shell), pointer :: g_tmp_p => null()
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
  integer, dimension(100) :: L_set
  integer :: unique_L
  real, dimension(100) :: S_set
  integer :: unique_S
  integer, dimension(:,:), allocatable :: table
!-----------------------------------------------------------
  integer :: itrav
  integer(kind=long) :: istep
  integer(kind=long) :: iter
  integer(kind=long) :: icounter
  integer :: iL, iS
!-----------------------------------------------------------
  write(*,*) "------------------------------------------------------"
  write(*,*) " This program is written for determining term symbols "
  write(*,*) "               Programmed by St Maxwell               "
  write(*,*) "                     2019-jul-7th                     "
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

  iter = 1
  do while (.true.)
  ! build a linked list
  ! each link receives a shell
    select case(shell_type(iter))
    case('s')
      allocate(s_define_p)
      call s_define_p%set_nocc(shell_occ(iter))
      call s_define_p%initialize()
      shells_num_comb(iter) = s_define_p%num_comb
      tmp => s_define_p
      s_define_p => null()
    case('p')
      allocate(p_define_p)
      call p_define_p%set_nocc(shell_occ(iter))
      call p_define_p%initialize()
      shells_num_comb(iter) = p_define_p%num_comb
      tmp => p_define_p
      p_define_p => null()
    case('d')
      allocate(d_define_p)
      call d_define_p%set_nocc(shell_occ(iter))
      call d_define_p%initialize()
      shells_num_comb(iter) = d_define_p%num_comb
      tmp => d_define_p
      d_define_p => null()
    case('f')
      allocate(f_define_p)
      call f_define_p%set_nocc(shell_occ(iter))
      call f_define_p%initialize()
      shells_num_comb(iter) = f_define_p%num_comb
      tmp => f_define_p
      f_define_p => null()
    case('g')
      allocate(g_define_p)
      call g_define_p%set_nocc(shell_occ(iter))
      call g_define_p%initialize()
      shells_num_comb(iter) = g_define_p%num_comb
      tmp => g_define_p
      g_define_p => null()
    end select 

    if (.not. associated(head)) then
      allocate(head)
      tail => head
      tail%next => null()
      tail%variable => tmp
    else
      allocate(tail%next)
      tail => tail%next
      tail%next => null()
      tail%variable => tmp
    end if

    if (iter == num_shell) exit
    iter = iter + 1

  end do

  ! obtain the total number of microstates
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
  itrav = 1
  ptr => head

  do while (.true.)
    if (.not. associated(ptr)) exit
    ! initialize shells and find microstates
    tmp => ptr%variable
    select type(tmp)
    type is (s_shell)
      s_tmp_p => tmp
      call s_tmp_p%get_combination()

      allocate(L_shell_comb(s_tmp_p%num_comb))
      allocate(S_shell_comb(s_tmp_p%num_comb))
      L_shell_comb = s_tmp_p%comb(:)%L
      S_shell_comb = s_tmp_p%comb(:)%S

    type is (p_shell)
      p_tmp_p => tmp
      call p_tmp_p%get_combination()

      allocate(L_shell_comb(p_tmp_p%num_comb))
      allocate(S_shell_comb(p_tmp_p%num_comb))
      L_shell_comb = p_tmp_p%comb(:)%L
      S_shell_comb = p_tmp_p%comb(:)%S

    type is (d_shell)
      d_tmp_p => tmp
      call d_tmp_p%get_combination()
      
      allocate(L_shell_comb(d_tmp_p%num_comb))
      allocate(S_shell_comb(d_tmp_p%num_comb))
      L_shell_comb = d_tmp_p%comb(:)%L
      S_shell_comb = d_tmp_p%comb(:)%S

    type is (f_shell)
      f_tmp_p => tmp
      call f_tmp_p%get_combination()
      
      allocate(L_shell_comb(f_tmp_p%num_comb))
      allocate(S_shell_comb(f_tmp_p%num_comb))
      L_shell_comb = f_tmp_p%comb(:)%L
      S_shell_comb = f_tmp_p%comb(:)%S

    type is (g_shell)
      g_tmp_p => tmp
      call g_tmp_p%get_combination()
      
      allocate(L_shell_comb(g_tmp_p%num_comb))
      allocate(S_shell_comb(g_tmp_p%num_comb))
      L_shell_comb = g_tmp_p%comb(:)%L
      S_shell_comb = g_tmp_p%comb(:)%S

    end select
    ! combinate microstates of each shell
    ! to get total microstates
    istep = 1
    do iter = itrav+1, num_shell
      istep = istep * shells_num_comb(iter)
    end do
    itrav = itrav + 1

    icounter = 1
    do iter = 1, tot_num_comb, istep

      L_all_comb(iter:iter+istep-1) = &
        L_all_comb(iter:iter+istep-1) + L_shell_comb(icounter)
      S_all_comb(iter:iter+istep-1) = &
        S_all_comb(iter:iter+istep-1) + S_shell_comb(icounter)

      if (icounter == size(L_shell_comb)) then
        icounter = 1
        cycle
      end if
      icounter = icounter + 1

    end do

    ptr => ptr%next
    deallocate(L_shell_comb)
    deallocate(S_shell_comb)

  end do

  ! obtain all possible values of L & S 
  L_all_comb_sorted = L_all_comb
  S_all_comb_sorted = S_all_comb
  call sort_L(L_all_comb_sorted)
  call sort_S(S_all_comb_sorted)
  call set_L(L_all_comb_sorted, L_set, unique_L)
  call set_S(S_all_comb_sorted, S_set, unique_S)
  deallocate(L_all_comb_sorted)
  deallocate(S_all_comb_sorted)

  ! build the table
  ! reference: Journal of Chemical Education, 52(2), 1975:87-89
  allocate(table(unique_S, unique_L))
  table = 0
  do iter = 1, tot_num_comb
    do iL = 1, unique_L
      if (L_set(iL) == L_all_comb(iter)) exit
    end do
    do iS = 1, unique_S
      if (abs(S_set(iS)-S_all_comb(iter)) < tol) exit
    end do
    table(iS, iL) = table(iS, iL) + 1
  end do

  call print_table(table, L_set, S_set, unique_L, unique_S)
  call print_term_symbol(table, unique_L, unique_S)

  read(*,*)
  
end program