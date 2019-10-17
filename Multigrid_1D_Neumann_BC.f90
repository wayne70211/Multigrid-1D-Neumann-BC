module Multigrid_1D_Neumann_BC
  !use Exact_Solver
  implicit none
contains
  !---------------------------------------------------------------------------!
  ! Multigrid V cycle scheme
  !---------------------------------------------------------------------------!
  subroutine MG_Vcycle(Nx,dx,RHS,U,Level_num,tol,exact_solver,I_cycle)
    implicit none
    integer,intent(in) :: Nx,Level_num,exact_solver
    real*8 ,intent(in) :: dx,tol
    real*8,dimension(0:Nx), intent(in)     :: RHS
    real*8,dimension(0:Nx), intent(inout)  :: U
    integer,dimension(Level_num) :: Level_Nx
    real*8 ,dimension(Level_num) :: Level_dx
    real*8  :: avg,rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,k,Level,iter
    integer, intent(out) :: I_cycle

    ! Defined all data space
    type space
      real*8,allocatable :: data(:)
    end type space

    type (space) UL(Level_num),R(Level_num),F(Level_num),P(Level_num)

    allocate(UL(1)%data(0:Nx),R(1)%data(0:Nx),F(1)%data(0:Nx),P(1)%data(0:Nx))

    Level_Nx(1) = Nx
    Level_dx(1) = dx

    do i=2,Level_num
      k=2**(i-1)
      Level_Nx(i) = Nx / k
      Level_dx(i) = dx * dble(k)
      allocate(UL(i)%data(0:Level_Nx(i)))
      allocate( R(i)%data(0:Level_Nx(i)))
      allocate( F(i)%data(0:Level_Nx(i)))
      allocate( P(i)%data(0:Level_Nx(i)))
    enddo

    Max_Iter   = 1000000     ! Allowed maximum number of outer iteration
    Relax_Iter = 2          ! Number of relaxation for restriction in V-cycle

    !Check the coarsest grid
    if (Level_Nx(Level_num).le.3) then
    write(*,*) Level_num," level is high for ",Nx," grid "
    stop
    end if

    do i=0,Level_Nx(1)
      F(1)%data(i)  = RHS(i)
      UL(1)%data(i) = U(i)
      R(1)%data(i)  = 0.0d0
    enddo

    !Compute initial resitual:
    call Residual(Nx,dx,F(1)%data,UL(1)%data,R(1)%data)
    call L2norm(Nx,R(1)%data,rms0)

    !open(40,file='residual_All_MG6V3.plt')
    !write(40,*) 'variables ="Cycle","Iteration","rms","Residual"'

    iter=1

    do I_cycle=1,Max_Iter

      call Residual(Level_Nx(1),Level_dx(1),F(1)%data,UL(1)%data,R(1)%data)

        ! Check for convergence on finest grid
      call L2norm(Level_Nx(1),R(1)%data,rms)
      !write(40,*) I_cycle,iter,rms,rms/rms0
      !write(*,*) I_cycle,iter,rms,rms/rms0

    !---------------------------------------------------------------------------!
      do Level=1,Level_num-1

        ! Relax
        do i=1,Relax_Iter
          iter=iter+1
          call Relax(Level_Nx(Level),Level_dx(Level),F(Level)%data,UL(Level)%data)
          !call Residual(Level_Nx(Level),Level_dx(Level),F(Level)%data,UL(Level)%data,R(Level)%data)
          !call L2norm(Level_Nx(Level),R(Level)%data,rms)
          !write(*,*) I_cycle,iter,rms,rms/rms0
          !write(40,*) I_cycle,iter,rms,rms/rms0

        end do

        ! Compute residual
        call Residual(Level_Nx(Level),Level_dx(Level),F(Level)%data,UL(Level)%data,R(Level)%data)

        ! Check for convergence on finest grid
        call L2norm(Level_Nx(Level),R(Level)%data,rms)
        !write(40,*) I_cycle,iter,rms,rms/rms0

        if (rms/rms0.le.tol .and. Level .eq. 1) goto 10

        ! Restriction
        call Restriction(Level_Nx(Level),Level_Nx(Level+1),R(Level)%data,F(Level+1)%data)

        ! Gram-Schmidt step
        avg=sum(F(Level+1)%data)/(Level_Nx(Level+1)+1)

        do i=0,Level_Nx(Level+1)
          F(Level+1)%data(i)=F(Level+1)%data(i)-avg
          UL(Level+1)%data(i) = 0.0d0
        end do

      end do
    !---------------------------------------------------------------------------!

      ! Compute residual on coarsest grid
      call Residual(Level_Nx(Level_num),Level_dx(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
      call L2norm(Level_Nx(Level_num),R(Level_num)%data,rmsc)
      !iter=iter+1
      !write(40,*) I_cycle,iter,rmsc,rmsc/rms0

      ! Solve exact solution on coarsest grid
      if (exact_solver .eq. 1) then
        do while (rms/rmsc .gt. tol )!.and. mod(iter,1000) .ne. 0)
          iter=iter+1
          call Relax(Level_Nx(Level_num),Level_dx(Level_num),F(Level_num)%data,UL(Level_num)%data)
          call Residual(Level_Nx(Level_num),Level_dx(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
          ! Check for convergence on smallest grid
          call L2norm(Level_Nx(Level_num),R(Level_num)%data,rms)
          !write(40,*) I_cycle,iter,rms,rms/rms0
          !write(*,*) I_cycle,iter,rms,rms/rms0
        end do

      else if (exact_solver .eq. 2) then
        !call LU_Solver(Level_Nx(Level_num),Level_dx(Level_num),F(Level_num)%data,UL(Level_num)%data)

      else
        !call CG_Solver(Level_Nx(Level_num),Level_dx(Level_num),F(Level_num)%data,UL(Level_num)%data)

      endif

    !---------------------------------------------------------------------------!

      do Level=Level_num-1,1,-1

        ! Prolongation
        call Prolongation(Level_Nx(Level+1),Level_Nx(Level),UL(Level+1)%data,P(Level)%data)

        ! Correct
        do i=0,Level_Nx(Level)
          UL(Level)%data(i) = UL(Level)%data(i) + P(Level)%data(i)
        end do

        ! Relax
        do i=1,Relax_Iter
          iter=iter+1
          call Relax(Level_Nx(Level),Level_dx(Level),F(Level)%data,UL(Level)%data)
        end do

      end do

      ! Gram-Schmidt step
      !avg=sum(F(1)%data)/(Level_Nx(1)+1)
      !do i=0,Level_Nx(1)
      !  F(1)%data(i)=F(1)%data(i)-avg
      !enddo

    end do  ! Outer iteration loop

    10 continue

    do i=0,Nx
      U(i) = UL(1)%data(i)
    end do

    !close(40)

    do i=1,Level_num
      deallocate(UL(i)%data,R(i)%data,F(i)%data,P(i)%data)
    enddo

    return
  end subroutine

  !---------------------------------------------------------------------------!
  !Relaxation formula for Poisson equation
  !Uses GS relaxation
  !Works for Neumann boundary conditions
  !---------------------------------------------------------------------------!
  subroutine Relax(Nx,dx,F,U)
    implicit none
    integer ,intent(in) :: Nx
    real*8  ,intent(in) :: dx
    real*8  ,dimension(0:Nx),intent(in)    :: F
    real*8  ,dimension(0:Nx),intent(inout) :: U
    integer :: i
    real*8  :: avg

    do i=1,Nx-1
      U(i) = 0.5d0*(U(i+1)+U(i-1)-dx*dx*F(i))
    end do

    ! Neumann B.C
    U(0) = U(1) - dx*dx*F(0)
    U(Nx) = U(Nx-1) - dx*dx*F(Nx)

    ! Zero-mean condition
    avg = sum(U)/(Nx+1)
    do i=0,Nx
      U(i) = U(i) - avg
    end do

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Residual formula for Poisson equation
  !---------------------------------------------------------------------------!
  subroutine Residual(Nx,dx,F,U,R)
    implicit none
    integer,intent(in) :: Nx
    real*8 ,intent(in) :: dx
    real*8, dimension(0:Nx),intent(in)  :: U,F
    real*8, dimension(0:Nx),intent(out) :: R

    integer :: i

    do i=1,Nx-1
      R(i) = F(i) - (U(i+1) - 2.0d0*U(i) + U(i-1))/(dx*dx)
    end do

    !Boundary conditions for residuals
    R(0)  = F(0) - (U(1) - U(0))/(dx*dx)
    R(Nx) = F(Nx) - (U(Nx-1) - U(Nx))/(dx*dx)


    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Restriction operators
  !---------------------------------------------------------------------------!
  subroutine Restriction(Nxf,Nxh,R,F)
    implicit none
    integer ,intent(in) :: Nxf,Nxh
    real*8, dimension(0:Nxf),intent(in)  :: R       !on higher grid
    real*8, dimension(0:Nxh),intent(out) :: F       !on lower grid
    integer :: i

    do i=1,Nxh-1
      F(i) = 0.25d0 * (R(2*i-1)+2.0d0*R(2*i)+R(2*i+1))
    end do

    !Boundary conditions
    F(0)   = 0.5d0 * R(0) + 0.25d0 * R(1)
    F(Nxh) = 0.5d0 * R(Nxf) + 0.25d0 * R(Nxf-1)

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Prolongation operator
  !---------------------------------------------------------------------------!
  subroutine Prolongation(Nxh,Nxf,U,P)
    implicit none
    integer ,intent(in) :: Nxf,Nxh
    real*8, dimension(0:Nxh) ,intent(in)  :: U       !on lower grid
    real*8, dimension(0:Nxf) ,intent(out) :: P       !on higher grid
    integer :: i

    do i=0,Nxh-1
      P(2*i)     = U(i)
      P(2*i+1)   = 0.5d0*(U(i)+U(i+1))
    end do

    !Boundary conditions
    P(Nxf)     = U(Nxh)

    return

  end subroutine

  !---------------------------------------------------------------------------!
  ! Compute L2-norm for an array
  !---------------------------------------------------------------------------!
  subroutine L2norm(Nx,R,RMS)
    implicit none
    integer ,intent(in) :: Nx
    real*8, dimension(0:Nx),intent(in) :: R
    real*8  ,intent(out):: RMS
    integer :: i

    RMS=0.0d0

    do i=1,Nx-1
      RMS = RMS + R(i)*R(i)
    end do

    RMS= dsqrt(RMS/dfloat((Nx-1)))

    return
  end subroutine

end module
