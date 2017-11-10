PROGRAM TermSymbolProgram
implicit none
    ! 定义5种角量子数的轨道的轨道角量子数和自旋磁量子数
    integer :: mL_s(2)=(/0,0/)                  
    real :: mS_s(2)=(/0.5,-0.5/)
    integer :: mL_p(6)=(/1,1,0,0,-1,-1/)                  
    real :: mS_p(6)=(/0.5,-0.5,0.5,-0.5,0.5,-0.5/)
    integer :: mL_d(10)=(/2,2,1,1,0,0,-1,-1,-2,-2/)
    real :: mS_d(10)=(/0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5/)
    integer :: mL_f(14)=(/3,3,2,2,1,1,0,0,-1,-1,-2,-2,-3,-3/)                  
    real :: mS_f(14)=(/0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5/)
    integer :: mL_g(18)=(/4,4,3,3,2,2,1,1,0,0,-1,-1,-2,-2,-3,-3,-4,-4/)                  
    real :: mS_g(18)=(/0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5/)
    ! 计算过程中用于储存各角动量的轨道的所有组合形式及后续输出
    integer,allocatable :: ML_s_comb(:),ML_p_comb(:),ML_d_comb(:),ML_f_comb(:),ML_g_comb(:)
    real,allocatable :: MS_s_comb(:),MS_p_comb(:),MS_d_comb(:),MS_f_comb(:),MS_g_comb(:)
    integer,allocatable :: ML(:)
    real,allocatable :: MS(:)
    integer,allocatable :: ML_NET(:)
    real,allocatable :: MS_NET(:)
    integer,allocatable :: TABLE(:,:) 
    integer,allocatable :: TABLE_out(:,:)
    ! 计算中的一些参数，包括一个计算阶乘的函数
    integer :: n1=0
    integer :: n2=0
    integer :: n3=0
    integer :: n4=0
    integer :: n5=0
    integer :: ncomb,ncomb_s,ncomb_p,ncomb_d,ncomb_f,ncomb_g
    real,external :: Factorial
    integer,allocatable :: comb_s(:),comb_p(:),comb_d(:),comb_f(:),comb_g(:)
    ! 大多为临时变量
    integer :: counter_out
    common counter_out
    integer :: i,j,k,l,m,n,counter1,counter2 
    real :: counter3
    integer :: unique_number_L ! used for removing repeated elements
    integer :: unique_number_S ! used for removing repeated elements
    integer :: x,y ! used for TABLE
    integer :: z,w ! used for TABLE_out
    integer :: element
    character(len=60) :: bottom_line 
    
    write(*,*) "This program is written for determining term symbols"
    write(*,*) "Programmed by St Maxwell"
    write(*,*) "Input the number of electrons in each orbital:"
    write(*,*) "(e.g. s1 p1 d1 is noted as 1 1 1 0 0)"
    read(*,*) n1,n2,n3,n4,n5
    
    if(n1/=0) then
        ncomb_s=nint(Factorial(2.0)/(Factorial(float(n1))*Factorial(float(2-n1))))
        allocate(ML_s_comb(ncomb_s))
        allocate(MS_s_comb(ncomb_s))
        allocate(comb_s(0:n1))
        ML_s_comb=0
        MS_s_comb=0
        counter_out=1
        call gen_s(1)
    else
        allocate(ML_s_comb(1))
        allocate(MS_s_comb(1))
        ncomb_s=1
        ML_s_comb=0
        MS_s_comb=0
    end if
    if(n2/=0) then
        ncomb_p=nint(Factorial(6.0)/(Factorial(float(n2))*Factorial(float(6-n2))))
        allocate(ML_p_comb(ncomb_p))
        allocate(MS_p_comb(ncomb_p))
        allocate(comb_p(0:n2))
        ML_p_comb=0
        MS_p_comb=0
        counter_out=1
        call gen_p(1)
    else
        allocate(ML_p_comb(1))
        allocate(MS_p_comb(1))
        ncomb_p=1
        ML_p_comb=0
        MS_p_comb=0
    end if    
    if(n3/=0) then
        ncomb_d=nint(Factorial(10.0)/(Factorial(float(n3))*Factorial(float(10-n3))))
        allocate(ML_d_comb(ncomb_d))
        allocate(MS_d_comb(ncomb_d))
        allocate(comb_d(0:n3))
        ML_d_comb=0
        MS_d_comb=0
        counter_out=1
        call gen_d(1)
    else
        allocate(ML_d_comb(1))
        allocate(MS_d_comb(1))
        ncomb_d=1
        ML_d_comb=0
        MS_d_comb=0
    end if    
    if(n4/=0) then
        ncomb_f=nint(Factorial(14.0)/(Factorial(float(n4))*Factorial(float(14-n4))))
        allocate(ML_f_comb(ncomb_f))
        allocate(MS_f_comb(ncomb_f))
        allocate(comb_f(0:n4))
        ML_f_comb=0
        MS_f_comb=0
        counter_out=1
        call gen_f(1)
    else
        allocate(ML_f_comb(1))
        allocate(MS_f_comb(1))
        ncomb_f=1
        ML_f_comb=0
        MS_f_comb=0
    end if    
    if(n5/=0) then
        ncomb_g=nint(Factorial(18.0)/(Factorial(float(n5))*Factorial(float(18-n5))))
        write(*,*) Factorial(18.0),Factorial(float(n5)),Factorial(float(18-n5)),ncomb_g
        allocate(ML_g_comb(ncomb_g))
        allocate(MS_g_comb(ncomb_g))
        allocate(comb_g(0:n5))
        ML_g_comb=0
        MS_g_comb=0
        counter_out=1
        call gen_g(1)
    else
        allocate(ML_g_comb(1))
        allocate(MS_g_comb(1))
        ncomb_g=1
        ML_g_comb=0
        MS_g_comb=0
    end if    
    
    ncomb=ncomb_s*ncomb_p*ncomb_d*ncomb_f*ncomb_g
    allocate(ML(ncomb))
    allocate(MS(ncomb))
    allocate(ML_NET(ncomb))
    allocate(MS_NET(ncomb))
    counter1=1
    
    do i=1,ncomb_s
        do j=1,ncomb_p
    	      do k=1,ncomb_d
    	  	      do l=1,ncomb_f
    	  	  	      do m=1,ncomb_g
    	  	  	  	      ML(counter1)=ML_s_comb(i)+ML_p_comb(j)+ML_d_comb(k)+ML_f_comb(l)+ML_g_comb(m)
    	  	  	  	      MS(counter1)=MS_s_comb(i)+MS_p_comb(j)+MS_d_comb(k)+MS_f_comb(l)+MS_g_comb(m)
    	  	  	  	      counter1=counter1+1
    	  	  	      end do
    	  	      end do
    	      end do
        end do
    end do
    
    unique_number_L = 1
    ML_NET(1) = ML(1)
    outer_1: do m=2,ncomb
        do n=1,unique_number_L
            if (ML_NET(n) == ML(m)) then
                ! Found a match so start looking again
                cycle outer_1
            end if
        end do
        ! No match found so add it to the output
        unique_number_L = unique_number_L + 1
        ML_NET(unique_number_L) = ML(m)
    end do outer_1
    
    unique_number_S = 1
    MS_NET(1) = MS(1)
    outer_2: do m=2,ncomb
        do n=1,unique_number_S
            if (MS_NET(n) == MS(m)) then
                ! Found a match so start looking again
                cycle outer_2
            end if
        end do
        ! No match found so add it to the output
        unique_number_S = unique_number_S + 1
        MS_NET(unique_number_S) = MS(m)
    end do outer_2   
    
    do i=1, unique_number_L
        counter2 = ML_NET(i)
        do j=i+1, unique_number_L
            if (ML_NET(j) > ML_NET(i)) then
                counter2 = ML_NET(j)
                ML_NET(j) = ML_NET(i)
                ML_NET(i) = counter2
            end if
        end do
    end do
    
    do i=1, unique_number_S
        counter3 = MS_NET(i)
        do j=i+1, unique_number_S
            if (MS_NET(j) > MS_NET(i)) then
                counter3 = MS_NET(j)
                MS_NET(j) = MS_NET(i)
                MS_NET(i) = counter3
            end if
        end do
    end do    
    
    allocate(TABLE(unique_number_S,unique_number_L))
    TABLE=0
    do i=1,ncomb
        x=1
        y=1
        do while(MS_NET(x)/=MS(i))
            x=x+1
        end do
        do while(ML_NET(y)/=ML(i))
            y=y+1
        end do
        TABLE(x,y)=TABLE(x,y)+1
    end do
    
    write(*,*) " "
    write(*,*) "ML |"
    write(*,*) "   |"
    do i=1,unique_number_L
        write(*,20) ML_NET(i)
        write(*,50) "|   "
        if (unique_number_S/=1) then
            write(*,100) TABLE(1:unique_number_S-1,i)
            write(*,"(I2)") TABLE(unique_number_S,i)
        else
            write(*,"(I2)") TABLE(unique_number_S,i)
        end if
    end do
    bottom_line(1:5)="   |_"
    do i=6,(unique_number_S*6+10)
        bottom_line(i:i)="_"
    end do
    write(*,*) bottom_line
    write(*,200) MS_NET(1:unique_number_S)
20  format((I3,' '),\)    
50  format(A2,\)    
100 format((I2,'    '),\)
200 format('     ',(F4.1,'  '),\)
    write(*,*) ' MS'
    write(*,*) " "
    write(*,*) "Term Symbols:"
    
    z=(unique_number_L+1)/2
    w=ceiling(unique_number_S/2.0)
    allocate(TABLE_out(z,w))
    TABLE_out=0
    do i=1,z
        TABLE_out(i,:)=TABLE(:,i)
    end do
    do i=1,z
        do while (.NOT. all(TABLE_out(i,:)<=0))
            do j=1,w
                if (TABLE_out(i,j)<=0) then
                    cycle
                else
                    element=TABLE_out(i,j)
                    call output_terms((z-i),unique_number_S,j,element) 
                    do k=i,z
                        do l=j,w
                            TABLE_out(k,l)=TABLE_out(k,l)-element
                        end do
                    end do
                end if
            end do
        end do
    end do
    
    write(*,*) 
    write(*,*) "Finished, press Enter to exit"
    read(*,*)
    stop
    
    contains
    recursive subroutine gen_s(m)
        implicit none
        integer, intent(in) :: m    
        integer :: l,n
        integer :: counter_in
        common counter_in
        
        if (m > n1) then
            forall(l=1:n1)
                ML_s_comb(counter_in)=ML_s_comb(counter_in)+mL_s(comb_s(l))
                MS_s_comb(counter_in)=MS_s_comb(counter_in)+mS_s(comb_s(l))
            end forall
            counter_in=counter_in+1
        else
            do n = 1, 2
                if ((m == 1) .or. (n > comb_s(m - 1))) then
                    comb_s(m) = n
                    call gen_s(m + 1)
                end if
            end do
        end if
    end subroutine gen_s
    
    recursive subroutine gen_p(m)
        implicit none
        integer, intent(in) :: m    
        integer :: l,n
        integer :: counter_in
        common counter_in
        
        if (m > n2) then
            forall(l=1:n2)
                ML_p_comb(counter_in)=ML_p_comb(counter_in)+mL_p(comb_p(l))
                MS_p_comb(counter_in)=MS_p_comb(counter_in)+mS_p(comb_p(l))
            end forall
            counter_in=counter_in+1
        else
            do n = 1, 6
                if ((m == 1) .or. (n > comb_p(m - 1))) then
                    comb_p(m) = n
                    call gen_p(m + 1)
                end if
            end do
        end if
    end subroutine gen_p
    
    recursive subroutine gen_d(m)
        implicit none
        integer, intent(in) :: m    
        integer :: l,n
        integer :: counter_in
        common counter_in
        
        if (m > n3) then
            forall(l=1:n3)
                ML_d_comb(counter_in)=ML_d_comb(counter_in)+mL_d(comb_d(l))
                MS_d_comb(counter_in)=MS_d_comb(counter_in)+mS_d(comb_d(l))
            end forall
            counter_in=counter_in+1
        else
            do n = 1, 10
                if ((m == 1) .or. (n > comb_d(m - 1))) then
                    comb_d(m) = n
                    call gen_d(m + 1)
                end if
            end do
        end if
    end subroutine gen_d
    
    recursive subroutine gen_f(m)
        implicit none
        integer, intent(in) :: m    
        integer :: l,n
        integer :: counter_in
        common counter_in
        
        if (m > n4) then
            forall(l=1:n4)
                ML_f_comb(counter_in)=ML_f_comb(counter_in)+mL_f(comb_f(l))
                MS_f_comb(counter_in)=MS_f_comb(counter_in)+mS_f(comb_f(l))
            end forall
            counter_in=counter_in+1
        else
            do n = 1, 14
                if ((m == 1) .or. (n > comb_f(m - 1))) then
                    comb_f(m) = n
                    call gen_f(m + 1)
                end if
            end do
        end if
    end subroutine gen_f
    
    recursive subroutine gen_g(m)
        implicit none
        integer, intent(in) :: m    
        integer :: l,n
        integer :: counter_in
        common counter_in
        
        if (m > n5) then
            forall(l=1:n5)
                ML_g_comb(counter_in)=ML_g_comb(counter_in)+mL_g(comb_g(l))
                MS_g_comb(counter_in)=MS_g_comb(counter_in)+mS_g(comb_g(l))
            end forall
            counter_in=counter_in+1
        else
            do n = 1, 18
                if ((m == 1) .or. (n > comb_g(m - 1))) then
                    comb_g(m) = n
                    call gen_g(m + 1)
                end if
            end do
        end if
    end subroutine gen_g

END PROGRAM TermSymbolProgram

real function Factorial(z)
    implicit none
    real :: z,w
    real(kind=8) :: ans
    ans=1
    do w=1,z
  	    ans=ans*w
    end do
    Factorial=ans
end function Factorial
    
subroutine output_terms(L,S_num,seq,times)
    implicit none
    integer :: L,S_num,seq,times
    integer :: multiplicity
    character :: term
    multiplicity=S_num-2*seq+2
    select case(L)
    case(0)
        term="S"
    case(1)
        term="P"
    case(2)
        term="D"
    case(3)
        term="F"        
    case(4)
        term="G"
    case(5)
        term="H"    
    case(6)
        term="I"
    case(7)
        term="K"    
    case(8)
        term="L"
    case(9)
        term="M"
    case(10)
        term="N"
    case(11)
        term="O"
    case(12)
        term="Q"
    case(13)
        term="R"        
    case(14)
        term="T"
    case(15)
        term="U"    
    case(16)
        term="V"
    case(17)
        term="W"    
    case(18)
        term="X"
    case(19)
        term="Y" 
    case(20)
        term="Z"
    case default
        write(term,"(I2)") L
    end select
    if (times==1) then
        write(*,300) multiplicity,term
    else
        write(*,400) multiplicity,term,times
    end if  
300 format((I1,A1," "),\)
400 format((I1,A1,"<",I1,"> "),\)
    return
end subroutine
