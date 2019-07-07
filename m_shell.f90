module m_shell
  use m_combination
  implicit none
  private
  public s_shell, p_shell, d_shell, f_shell, g_shell

  type :: spin_orb
    integer :: mag_qn
    real :: spin_qn
  end type

  type :: combination
    integer :: L = 0
    real :: S = 0
  end type

  type :: shell
    private
    integer :: num_occ = 0
  end type

  type, extends(shell) :: s_shell
    type(spin_orb), dimension(2) :: orbs
    integer(kind=long) :: num_comb
    type(combination), dimension(:), allocatable :: comb
    contains
      ! set the occupation number of shell
      procedure, public :: set_nocc => set_nocc_s
      ! initialize spin orbitals
      ! calculate the number of microstates
      ! allocate ?_shell%comb
      procedure, public :: initialize => initialize_sub_s
      ! find all possible microstates
      ! calculate L and S of each microstate
      procedure, public :: get_combination => get_combination_s
  end type
    
  type, extends(shell) :: p_shell
    type(spin_orb), dimension(6) :: orbs
    integer(kind=long) :: num_comb
    type(combination), dimension(:), allocatable :: comb
    contains
      procedure, public :: set_nocc => set_nocc_p
      procedure, public :: initialize => initialize_sub_p
      procedure, public :: get_combination => get_combination_p
  end type

  type, extends(shell) :: d_shell
    type(spin_orb), dimension(10) :: orbs
    integer(kind=long) :: num_comb
    type(combination), dimension(:), allocatable :: comb
    contains
      procedure, public :: set_nocc => set_nocc_d
      procedure, public :: initialize => initialize_sub_d
      procedure, public :: get_combination => get_combination_d
  end type

  type, extends(shell) :: f_shell
    type(spin_orb), dimension(14) :: orbs
    integer(kind=long) :: num_comb
    type(combination), dimension(:), allocatable :: comb
    contains
      procedure, public :: set_nocc => set_nocc_f
      procedure, public :: initialize => initialize_sub_f
      procedure, public :: get_combination => get_combination_f
  end type

  type, extends(shell) :: g_shell
    type(spin_orb), dimension(18) :: orbs
    integer(kind=long) :: num_comb
    type(combination), dimension(:), allocatable :: comb
    contains
      procedure, public :: set_nocc => set_nocc_g
      procedure, public :: initialize => initialize_sub_g
      procedure, public :: get_combination => get_combination_g
  end type

  contains

  subroutine set_nocc_s(this, n_occ)
    implicit none
    class(s_shell) :: this
    integer, intent(in) :: n_occ

    if (n_occ > 2) then
      stop "the number of occupation is out of range"
    end if
    this%num_occ = n_occ

  end subroutine

  subroutine set_nocc_p(this, n_occ)
    implicit none
    class(p_shell) :: this
    integer, intent(in) :: n_occ

    if (n_occ > 6) then
      stop "the number of occupation is out of range"
    end if
    this%num_occ = n_occ
    
  end subroutine

  subroutine set_nocc_d(this, n_occ)
    implicit none
    class(d_shell) :: this
    integer, intent(in) :: n_occ

    if (n_occ > 10) then
      stop "the number of occupation is out of range"
    end if
    this%num_occ = n_occ
    
  end subroutine

  subroutine set_nocc_f(this, n_occ)
    implicit none
    class(f_shell) :: this
    integer, intent(in) :: n_occ

    if (n_occ > 14) then
      stop "the number of occupation is out of range"
    end if
    this%num_occ = n_occ
    
  end subroutine

  subroutine set_nocc_g(this, n_occ)
    implicit none
    class(g_shell) :: this
    integer, intent(in) :: n_occ

    if (n_occ > 18) then
      stop "the number of occupation is out of range"
    end if
    this%num_occ = n_occ
    
  end subroutine

  subroutine initialize_sub_s(this)
    implicit none
    class(s_shell) :: this

    if (allocated(this%comb)) deallocate(this%comb)
    this%orbs(:)%mag_qn = [0,0]
    this%orbs(:)%spin_qn = [0.5,-0.5]
    this%num_comb = combinations(2, this%num_occ)
    allocate( this%comb(this%num_comb) )

  end subroutine

  subroutine initialize_sub_p(this)
    implicit none
    class(p_shell) :: this

    if (allocated(this%comb)) deallocate(this%comb)
    this%orbs(:)%mag_qn = [1,1,0,0,-1,-1]
    this%orbs(:)%spin_qn = [0.5,-0.5,0.5,-0.5,0.5,-0.5]
    this%num_comb = combinations(6, this%num_occ)
    allocate( this%comb(this%num_comb) )

  end subroutine

  subroutine initialize_sub_d(this)
    implicit none
    class(d_shell) :: this

    if (allocated(this%comb)) deallocate(this%comb)
    this%orbs(:)%mag_qn = [2,2,1,1,0,0,-1,-1,-2,-2]
    this%orbs(:)%spin_qn = [0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5, &
                           -0.5,0.5,-0.5]
    this%num_comb = combinations(10, this%num_occ)
    allocate( this%comb(this%num_comb) )
                                         
  end subroutine

  subroutine initialize_sub_f(this)
    implicit none
    class(f_shell) :: this

    if (allocated(this%comb)) deallocate(this%comb)
    this%orbs(:)%mag_qn = [3,3,2,2,1,1,0,0,-1,-1,-2,-2,-3,-3]
    this%orbs(:)%spin_qn = [0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5, &
                           -0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5]
    this%num_comb = combinations(14, this%num_occ)
    allocate( this%comb(this%num_comb) )
                             
  end subroutine

  subroutine initialize_sub_g(this)
    implicit none
    class(g_shell) :: this

    if (allocated(this%comb)) deallocate(this%comb)
    this%orbs(:)%mag_qn = [4,4,3,3,2,2,1,1,0,0,-1,-1,-2,-2, &
                          -3,-3,-4,-4]
    this%orbs(:)%spin_qn = [0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5, &
                           -0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5, &
                            0.5,-0.5,0.5,-0.5]
    this%num_comb = combinations(18, this%num_occ)
    allocate( this%comb(this%num_comb) )

  end subroutine

  subroutine get_combination_s(this)
    implicit none
    class(s_shell) :: this
    integer, allocatable, dimension(:,:) :: C
    integer(kind=long) :: i, j

    allocate(C(this%num_comb, this%num_occ))
    call find_all_combination_wrapper(2, this%num_occ, &
                                      this%num_comb, C)

    this%comb(:)%L = 0
    this%comb(:)%S = 0
    do i = 1, this%num_comb
      do j = 1, this%num_occ
        this%comb(i)%L = this%comb(i)%L + this%orbs(C(i,j))%mag_qn
        this%comb(i)%S = this%comb(i)%S + this%orbs(C(i,j))%spin_qn
      end do
    end do

    deallocate(C)

  end subroutine

  subroutine get_combination_p(this)
    implicit none
    class(p_shell) :: this
    integer, allocatable, dimension(:,:) :: C
    integer(kind=long) :: i, j

    allocate(C(this%num_comb, this%num_occ))
    call find_all_combination_wrapper(6, this%num_occ, &
                                      this%num_comb, C)

    this%comb(:)%L = 0
    this%comb(:)%S = 0
    do i = 1, this%num_comb
      do j = 1, this%num_occ
        this%comb(i)%L = this%comb(i)%L + this%orbs(C(i,j))%mag_qn
        this%comb(i)%S = this%comb(i)%S + this%orbs(C(i,j))%spin_qn
      end do
    end do

    deallocate(C)

  end subroutine

  subroutine get_combination_d(this)
    implicit none
    class(d_shell) :: this
    integer, allocatable, dimension(:,:) :: C
    integer(kind=long) :: i, j

    allocate(C(this%num_comb, this%num_occ))
    call find_all_combination_wrapper(10, this%num_occ, &
                                      this%num_comb, C)

    this%comb(:)%L = 0
    this%comb(:)%S = 0
    do i = 1, this%num_comb
      do j = 1, this%num_occ
        this%comb(i)%L = this%comb(i)%L + this%orbs(C(i,j))%mag_qn
        this%comb(i)%S = this%comb(i)%S + this%orbs(C(i,j))%spin_qn
      end do
    end do

    deallocate(C)

  end subroutine

  subroutine get_combination_f(this)
    implicit none
    class(f_shell) :: this
    integer, allocatable, dimension(:,:) :: C
    integer(kind=long) :: i, j

    allocate(C(this%num_comb, this%num_occ))
    call find_all_combination_wrapper(14, this%num_occ, &
                                      this%num_comb, C)

    this%comb(:)%L = 0
    this%comb(:)%S = 0
    do i = 1, this%num_comb
      do j = 1, this%num_occ
        this%comb(i)%L = this%comb(i)%L + this%orbs(C(i,j))%mag_qn
        this%comb(i)%S = this%comb(i)%S + this%orbs(C(i,j))%spin_qn
      end do
    end do

    deallocate(C)

  end subroutine

  subroutine get_combination_g(this)
    implicit none
    class(g_shell) :: this
    integer, allocatable, dimension(:,:) :: C
    integer(kind=long) :: i, j

    allocate(C(this%num_comb, this%num_occ))
    call find_all_combination_wrapper(18, this%num_occ, &
                                      this%num_comb, C)

    this%comb(:)%L = 0
    this%comb(:)%S = 0
    do i = 1, this%num_comb
      do j = 1, this%num_occ
        this%comb(i)%L = this%comb(i)%L + this%orbs(C(i,j))%mag_qn
        this%comb(i)%S = this%comb(i)%S + this%orbs(C(i,j))%spin_qn
      end do
    end do

    deallocate(C)

  end subroutine

end module

module m_spectra
  implicit none
  contains

    subroutine print_table(table, L_set, S_set, uni_L, uni_S)
      implicit none
      integer, dimension(:,:), intent(in) :: table
      integer, dimension(:), intent(in) :: L_set
      real, dimension(:), intent(in) :: S_set
      integer, intent(in) :: uni_L, uni_S
      !-----------------------------------------------------
      character(len=10) :: L_label
      character(len=150) :: table_row
      character(len=:), allocatable :: tab_row_fmt
      character(len=5) :: repeat
      character(len=150) :: bottom_line
      character(len=150) :: S_label
      character(len=:), allocatable :: S_label_fmt      
      !-----------------------------------------------------
      integer :: iter

      ! there are uni_S repeats in a row
      ! format ??I6
      write(repeat,"(I2)") uni_S
      tab_row_fmt = '(' // trim(adjustl(repeat)) // 'I6)'
      ! format ??F6.1
      S_label_fmt = '(' // trim(adjustl(repeat)) // 'F6.1)'

      write(*,*) " "
      write(*,*) "  ML |"
      write(*,*) "     |"

      do iter = 1, uni_L
        write(L_label,"(I4,' |')") L_set(iter)
        write(table_row,fmt=tab_row_fmt) table(:,iter)
        write(*,*) trim(L_label) // trim(table_row)
      end do

      bottom_line = "     |"
      do iter = 1, uni_S+1
        bottom_line = trim(bottom_line) // "______"
      end do
      write(*,*) trim(bottom_line)

      write(S_label,fmt=S_label_fmt) S_set(:uni_S)
      write(*,*) "      " // trim(S_label) // "    MS"

    end subroutine

    subroutine print_term_symbol(table, uni_L, uni_S)
    ! reference: Journal of Chemical Education, 52(2), 1975:87-89
    !       doi: 10.1021/ed052p87
      implicit none
      integer, dimension(:,:), intent(in) :: table
      integer, intent(in) :: uni_L, uni_S
      integer, dimension(:,:), allocatable :: table_out
      integer :: z, w, repeat
      character(len=10) :: symbol
      character(len=10), dimension(50) :: list_symbol
      integer :: num_symbol, remain, nrows
      integer :: iter, jter, kter, lter

      z = (uni_L+1) / 2
      w = ceiling(uni_S/2.0)

      allocate(table_out(z,w))

      table_out = 0
      do iter = 1, z
        do jter = 1, w
          table_out(iter, jter) = table(jter,iter)
        end do
      end do

      num_symbol = 1
      do iter = 1, z
        do while (.not. all(table_out(iter,:) <= 0))
          do jter = 1, w
            if (table_out(iter, jter) <= 0) then
              cycle
            else
              repeat = table_out(iter, jter)
              call output_terms(symbol, (z-iter), &
                                uni_S, jter, repeat) 

              list_symbol(num_symbol) = symbol
              num_symbol = num_symbol + 1
              do kter = iter, z
                do lter = jter, w
                  table_out(kter, lter) = &
                    table_out(kter, lter) - repeat
                end do
              end do
            end if
          end do
        end do
      end do

      write(*,*) 
      write(*,*) "Term Symbols:"
      remain = mod(num_symbol-1, 5)
      nrows = (num_symbol - remain - 1) / 5
      if (nrows /= 0) then
        do iter = 1, nrows
          jter = 5 * iter - 4
          kter = 5 * iter
          write(*,*) list_symbol(jter:kter)
        end do
        write(*,*) list_symbol(kter+1:kter+remain)
      else
        write(*,*) list_symbol(:remain)
      end if

    end subroutine

    subroutine output_terms(symbol, L, S_num, seq, times)
    ! term symbol: (2S+1)L<times>
    ! 2S + 1 = multiplicity, but here we do not obtain it directly via S
      implicit none
      character(len=10) :: symbol
      integer :: L, S_num, seq, times
      integer :: multiplicity
      character :: term

      multiplicity = S_num - 2*seq + 2

      select case(L)
      case(0)
        term = "S"
      case(1)
        term = "P"
      case(2)
        term = "D"
      case(3)
        term = "F"        
      case(4)
        term = "G"
      case(5)
        term = "H"    
      case(6)
        term = "I"
      case(7)
        term = "K"    
      case(8)
        term = "L"
      case(9)
        term = "M"
      case(10)
        term = "N"
      case(11)
        term = "O"
      case(12)
        term = "Q"
      case(13)
        term = "R"        
      case(14)
        term = "T"
      case(15)
        term = "U"    
      case(16)
        term = "V"
      case(17)
        term = "W"    
      case(18)
        term = "X"
      case(19)
        term = "Y" 
      case(20)
        term = "Z"
      case default
        write(term,"(I2)") L
      end select

      select case(times)
      case(1)
        write(symbol,"(I1,A1)") multiplicity, term
      case(2:9)
        write(symbol,"((I1,A1,'<',I1,'>'))") multiplicity, &
                                               term, times
      case(10:)
        write(symbol,"((I1,A1,'<',I2,'>'))") multiplicity, &
                                               term, times
      end select

    end subroutine

end module