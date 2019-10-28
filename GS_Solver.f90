module GS_Solver
  use Multigrid_1D_Neumann_BC
  implicit none
contains
  subroutine GS_Solver(Nx,dx,RHS,U,tol,Iter)
    implicit none
    integer,intent(in) :: Nx
    real*8 ,intent(in) :: dx,tol
    real*8,dimension(0:Nx), intent(in)     :: RHS
    real*8,dimension(0:Nx), intent(inout)  :: U
    real*8,dimension(0:Nx)                 :: R,F
    integer, intent(out) :: Iter
    integer :: i
    real*8  :: rms0,rms

    R = 0.0d0
    F = RHS

    call Residual(Nx,dx,RHS,U,R)
    call L2norm(Nx,R,rms0)
    rms = rms0
    Iter = 0
    open(70,file='residual_All_GS.plt')
    write(70,*) 'variables ="Cycle","Iteration","rms","Residual"'
    
    do while (rms/rms0 .gt. tol )!.and. mod(iter,1000) .ne. 0)
      Iter=Iter+1
      call Relax(Nx,dx,F,U)
      call Residual(Nx,dx,F,U,R)
      ! Check for convergence on smallest grid
      call L2norm(Nx,R,rms)
      write(70,*) iter,iter,rms,rms/rms0
      !write(*,*) I_cycle,iter,rms,rms/rms0
    end do

    close(70)

  end subroutine
end module
