program kron_prod_multiplication
    !***********************************************************************
    !SUBMITTED WITH THE MANUSCRIPT:                                         *
    !A TECHNICAL NOTE ON IMPROVING COMPUTATIONAL EFFICIENCY IN LARGE LINEAR * 
    !INVERSE PROBLEMS: AN EXAMPLE FROM CARBON DIOXIDE FLUX ESTIMATION       *
    !VINEET YADAV[1] AND ANNA M. MICHALAK[1]                                *
    ![1]{DEPARTMENT OF GLOBAL ECOLOGY, CARNEGIE INSTITUTION FOR SCIENCE,    *
    !STANFORD, CALIFORNIA, USA 94305}                                       *
    !***********************************************************************
    !WRITTEN BY : VINEET YADAV                                              *  
    !DATE : 2/28/2013                                                       *
    !PURPOSE : COMPUTE HQ FROM THE DIRECT AN INDIRECT METHOD WHERE Q CAN    * 
    !DECOMPOSED AS A KRONECKER PRODUCT OF TWO MATRICES                      *
    !NOTE: THIS IS A COMPLETELY SERIAL FUNCTION                             *
    !***********************************************************************
    !COMPILATION PROCEDURE:                                                 *
    !TESTED WITH COMPILERS IFORT AND GFORTRAN                               *
    !FOR IFORT COMPILE EX:                                                  *
    !ifort -heap-arrays 1024 -o HQ.EXE HQ.f90                               *
    !THIS COMPILATION PROCEDURE IS NECESSARY AS IFORT PUTS TEMPORARY ARRAYS *
    !ON STACK                                                               *
    !FOR GFORTRAN COMPILE EX:                                               *	
    !gfortran -o HQ.EXE HQ.f90                                              *
    !*****************NOTES*************************************************
    !                                                                       *  
    !NOTE 1: THIS CODE DEMONSTRATES MATRIX MULTIPLICATION BETWEEN AN        *
    !ARBITRARY MATRIX H AND AN ARBITRARY COVARIANCE MATRIX Q THAT CAN BE    *
    !EXPRESSED AS A KRONECKER PRODUCT OF TWO MATRICES.                      *
    !                                                                       *  
    !NOTE 2: THIS CODE ASSUMES D AND E MATRIX AS SYMMETRIC MATRIX WHICH IS  * 
    !THE CASE IF Q IS THE COVARIANCE MATRIX                                 *
    !                                                                       *  
    !NOTE 3: THIS IS A COMPLETELY SERIAL PROGRAM WHICH RELIES ON FORTRAN’S  *
    !INTRINSIC FUNCTIONS TO SHOW THE COMPUTATIONAL EFFICIENCY OF THE        *
    !INDIRECT !METHOD IN COMPARISON TO THE DIRECT METHOD AND USES RANDOM    *
    !MATRICES FOR THIS PURPOSE.THE SEED FOR GENERATING RANDOM NUMBERS HAVE  *
    !BEEN KEPT CONSTANT IN THE PROGRAM                                      *
    !                                                                       *  
    !NOTE 4: ALL VARIABLES ARE AS DEFINED IN THE MANUSCRIPT                 * 
    !                                                                       *  
    !NOTE 5: MATLAB PROTOTYPE CODE SHOULD BE USED TO UNDERSTAND THE         * 
    !EQUATIONS IN THE MANUSCRIPT AND THERE IMPLEMENTATION.                  *
    !                                                                       *  
    !NOTE 6: SINCE WE DO NOT HAVE ANY CONTROL OVER HOW MATLAB MANAGES       *
    !MEMORY AND THE ITS EFFICIENCY IN IMPLEMENTING FOR LOOPS IN COMPARISON  * 
    !TO VECTOR OPERATIONS; WE HAVE PROVIDED THE FORTRAN CODE TO PRACTICALLY *
    !SHOW THE COMPUTATIONAL EFFICIENCY OF THE PROPOSED ALGORITHM FOR        * 
    !MATRIX MULTIPLICATION                                                  *
    !                                                                       *  
    !NOTE 7: THE MATRIX MULTIPLICATION IN THIS PROGRAM CAN BE SPED UP BY    *  
    !USING ATLAS/INTEL MKL LIBRARY SUBROUTINE DGEMM. IF YOU WANT TO DO THAT *  
    !JUST UNCOMMENT LINE 298; LINK PROPERLY WITH ATLAS OR MKL LIBRARIES;    *
    !AND ANSWER QUESTION ARE YOU COMPILING WITH ATLAS OR MKL AS YES         *  
    !                                                                       *
    !NOTE 8: INTEL COMPILER GIVES FIRST VALUE GENERATED FROM RANDOM_NUMBER  * 
    !AS ZERO. THIS WOULD HAPPEN IN CASE OF D MATRIX. TO AVOID THIS CALL     *
    !RANDOM_NUMBER FUNCTION TWICE.                                          *
    !***********CODE BEGINS:TYPSET AND DEFINITIONS**************************

    implicit none

    !***********variables***************************************************

    character(len=10) :: col
    character(len=10) :: print_condition
    integer(kind=8), parameter:: rkind=8 ! precision of real numbers can only be 4 or 8
    integer(kind=4), parameter:: ikind=4 ! precision of integers can only be 4 or 8
    integer(kind=ikind):: p ! temporal covariance dimension 1; d matrix in manuscript
    integer(kind=ikind):: r ! spatial covariance dimension 1; e matrix in manuscript
    integer(kind=ikind):: n ! no of observations
    integer(kind=ikind):: row ! a counter for displaying matrices
    integer(kind=ikind):: seed=4 ! fixed random seed
    real(kind=rkind):: start ! variable for calculating time
    real(kind=rkind):: finish ! variable for calculating time
    real(kind=rkind), allocatable, dimension(:,:):: d
    real(kind=rkind), allocatable, dimension(:,:):: e
    real(kind=rkind), allocatable, dimension(:,:):: h
    real(kind=rkind), allocatable, dimension(:,:):: q
    real(kind=rkind), allocatable, dimension(:,:):: hq_indirect
    real(kind=rkind), allocatable, dimension(:,:):: hq_direct

    !**********read input parameters****************************************
    write(*,"(a)"),'********************************************************'
    write(*,"(a)"),'warning: the option of printing matrices is given so'
    write(*,"(a)"),'that user can check results from both methods. do not'
    write(*,"(a)"),'use this functionality to print large arrays which would'
    write(*,"(a)"),'not allow to clearly see the time used in computing '
    write(*,"(a)"),'matrix products from direct and indirect method'
    write(*,"(a)"),'********************************************************'
    write(*,*),' '
    write(*,"(a)"),'now we continue with the program: please enter the responses '
    write(*,*),' '
    write(*,"(a)"),'enter no of rows in h : an integer (n in the manuscript)'
    read (*,*) n
    write(*,"(a)"),'enter no of rows or columns in temporal covariance matrix : an integer (p in the manuscript)'
    read (*,*) p
    write(*,"(a)"),'enter no of rows or columns in spatial covariance matrix :, an integer (r in the manuscript)'
    read (*,*) r
    write(*,"(a)"),'enter yes or no for printing matrices h, q, d and e' 
    read (*,*) print_condition
    print_condition = trim(print_condition)

    !***************print read in information*******************************
    write(*,*),' '
    write(*,"(a)"),'***************************************************'
    write(*,"(a)"),'entered dimensions and choice for printing matrices'
    write(*,"(a)"),'***************************************************'
    write(*,*),' '
    write(*,"(a,i0)"),'no of rows in h: ',n
    write(*,"(a,i0)"),'no of rows in temporal covariance or d matrix: ',p
    write(*,"(a,i0)"),'no of columns in temporal covariance or d matrix: ',p
    write(*,"(a,i0)"),'no of rows in spatial covariance or e matrix: ',r
    write(*,"(a,i0)"),'no of columns in temporal covariance or d matrix: ',r
    write(*,"(a,a)"),'printing matrices: ', print_condition
    write(*,*),' '

    !**********allocation of memory*****************************************

    allocate (d(p,p))
    allocate (e(r,r))
    allocate (h(n,p*r))

    !**********fill arrays with random numbers******************************

    call random_seed(seed)
    call random_number(d) 
    call random_number(d)
    call random_number(e)
    call random_number(h)


    !***********compute hq from the direct method***************************

    allocate (q(p*r,p*r))
    call cpu_time(start)
    call kron(q, d, e)
    call cpu_time(finish)
    write(*,"(a,f12.4)"),'time for kron: ',finish-start
    allocate (hq_direct(n,p*r))

    call cpu_time(start)
    hq_direct = matmul (h,q)
    call cpu_time(finish)

    write(*,"(a,f12.4)"),'time hq from direct method in seconds: ',finish-start

    if (print_condition =='yes') then
        !****************print d ************************
        call print_mat(d,p,p,'temporal covariance')
        !****************print e ************************
        call print_mat(e,r,r,'spatial covariance')
        !****************print h ************************
        call print_mat(h,n,p*r,'h matrix')
        !****************print q**************************
        call print_mat(q,p*r,p*r,'q matrix')
        !****************print hq_direct*******************
        call print_mat(hq_direct,n,p*r,'hq_direct matrix')
    end if

    !***********compute hq from the indirect method**************************

    deallocate(q, hq_direct)
    allocate (hq_indirect(n,p*r))
    call cpu_time(start)
    hq_indirect=hq_ind_ic(p,p,r,r,n,d,e,h) 
    call cpu_time(finish)
    write(*,*),' '
    write(*,"(a,f12.4)"),'time hq from indirect method in seconds: ',finish-start
    write(*,*),' ' 

    if (print_condition =='yes') then
        !****************print hq_indirect*******************
        call print_mat(hq_indirect,n,p*r,'hq_indirect matrix')
    end if


    !*********************functions/subroutines start****************************

    contains
    subroutine kron(k, a, b)
    ! subroutine to compute kronecker product of two matrices
    ! assignment part
    implicit none
    integer(kind=ikind):: i
    integer(kind=ikind):: j
    integer(kind=ikind):: ma
    integer(kind=ikind):: na
    integer(kind=ikind):: mb
    integer(kind=ikind):: nb
    real(kind=rkind), intent(in):: a(:,:)
    real(kind=rkind), intent(in):: b(:,:)
    real(kind=rkind), intent(inout):: k(:,:)
    ! execution part
    ma = ubound(a, 1)
    na = ubound(a, 2)
    mb = ubound(b, 1)
    nb = ubound(b, 2)
    if (size(k,1) /= ma*mb .or. size(k,2) /= na*nb) then
        write(*,*) 'k has invalid size'
    end if
    forall(i=1:ma, j=1:na)
        k(mb*(i-1)+1:mb*i,nb*(j-1)+1:nb*j) = a(i,j)*b
    end forall
    end subroutine kron

    function hq_ind_ic(p,s,r,t,n,d,e,h) result (hq_indirect)
    ! function to compute hq indirectly
    ! note it is assumed that both d and e are square matrix
    ! as is the case in inverse problems
    ! p = number of rows or columns in the temporal covariance
    ! s = number of rows or columns in the temporal covariance
    ! r = number of rows or columns in the spatial covariance
    ! t = number of rows or columns in the spatial covariance
    ! n = total number of observation or the first dimension of h
    ! d = temporal covariance matrix
    ! e = spatial covariance matrix
    ! h = sensitivity matrix
    ! result the multiplication of h*q called as hq_indirect
    implicit none
    integer(kind=ikind),intent(in):: p
    integer(kind=ikind),intent(in):: s
    integer(kind=ikind),intent(in):: r
    integer(kind=ikind),intent(in):: t
    integer(kind=ikind),intent(in):: n
    integer(kind=ikind):: i
    integer(kind=ikind):: j
    integer(kind=ikind):: temp
    real(kind=rkind):: start ! variable for calculating time
    real(kind=rkind):: finish ! variable for calculating time
    real(kind=rkind), intent(in), dimension(:,:):: d(p,s)
    real(kind=rkind), intent(in), dimension(:,:):: e(r,t)
    real(kind=rkind), intent(in), dimension(:,:):: h(n,p*r)
    real(kind=rkind), dimension(:,:):: hqsum(n,r)
    real(kind=rkind), dimension(:,:):: hq_indirect (n,s*t)  
    do i = 1,s
        temp=i*t
        hqsum=0
        do j=1,p
            if (d(j,i) .ne. 0 .and. d(j,i) .ne. 1) then
                hqsum=hqsum+(h(:,1+(r*j)-r:j*r)*d(j,i))
            elseif (d(j,i) .eq. 1) then
                hqsum=hqsum+h(:,1+(r*j)-r:j*r)
            endif
        end do
        hq_indirect(:,temp-(t-1):temp)= matmul(hqsum,e)
    end do
    end function hq_ind_ic

    subroutine print_mat(matrix,rows,columns,name_matrix)
    ! function to print matrices
    ! variable assignment
    character(len=10):: col
    character(len=20):: formt
    character(len=*), intent(in):: name_matrix
    integer(kind=ikind):: n
    integer(kind=ikind), intent(in):: rows
    integer(kind=ikind), intent(in):: columns
    real(kind=rkind), intent(in), dimension(:,:):: matrix
    ! execution part
    write(col,"(i0)") columns
    col = trim(col)
    formt = '('//col//'f7.4))'
    write(*,"(a)"),' '
    write(*,"(a,a)"),'now printing: ', name_matrix
    write(*,"(a)"),' '
    do n = 1, rows
        write(*,formt) matrix(n,:)
    end do
    write(*,*),' '
    end subroutine print_mat

    end program kron_prod_multiplication
