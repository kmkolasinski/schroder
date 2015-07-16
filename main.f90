
program transporter
 use modutils
 use modjed
 use modinip
 use moddevice
 use ifport
 implicit none

! -----------------------------------------------
! Deklaracje zmiennych i tablic
! -----------------------------------------------
 integer :: nx            = 10;
 integer :: ny            = 10;
 integer :: nz            = 10;
 doubleprecision :: dx

! -----------------------------------------------
! Zmienne wczytywane z config.ini
! -----------------------------------------------
call setIniFilename("config.ini")
call getIntValue("Dane","nx",nx)
call getIntValue("Dane","ny",ny)
call getIntValue("Dane","nz",nz)
call getDoubleValue("Dane","dx",atomic_DX)
call getDoubleValue("Dane","m_eff",M_EFF)
call getDoubleValue("Dane","eps",E_MAT)
call getDoubleValue("Dane","g_lan",G_LAN)
call getDoubleValue("Dane","Ef",atomic_Ef)

dx = atomic_DX
call modjed_ustaw_konwersje_jednostek(0.0465D0,12.4D0);

call device_create(nx,ny,nz)

call device_add_tf_dens   (0.0D0,nx*dx,0.0D0,ny*dx,(nz/2-5)*dx,(nz/2-5)*dx,+1.0D-0)
call device_add_const_dens(0.0D0,nx*dx,0.0D0,ny*dx,(nz/2+5)*dx,(nz/2+5)*dx,+1.0D-0)


call device_add_const_pot((nx/2.0-5)*dx,(nx/2.0+5)*dx, &
                          (ny/2-5)*dx,(ny/2+5)*dx,&
                          (nz-1)*dx,nz*dx,-100.0D0)

!call device_add_const_pot(nx/2.0*dx-10,nx/2.0*dx+10,0.0D0,ny/4*dx,(nz-1)*dx,nz*dx,-400.0D0)
!call device_add_const_pot(nx/2.0*dx-10,&
!                         nx/2.0*dx+10,&
!                         ny*dx-ny/4*dx,&
!                         ny*dx,&
!                         (nz-1)*dx,nz*dx,-400.0D0)


DEVICE_FLAGS(:,:,1)  = DFLAG_DIRICHLET
DEVICE_POT(:,:,1)    = 0.0D0
!DEVICE_FLAGS(:,:,nz) = DFLAG_NEUMAN
!DEVICE_FLAGS(:,1,:)  = DFLAG_NEUMAN
!DEVICE_FLAGS(:,ny,:) = DFLAG_NEUMAN
!DEVICE_FLAGS(1,:,:)  = DFLAG_NEUMAN
!DEVICE_FLAGS(nx,:,:) = DFLAG_NEUMAN
call device_solve()
call write_to_file(321,"rho.txt",DEVICE_DENS(:,:,nz/2-4),nx,ny)
call write_to_file(321,"pot.txt",DEVICE_POT(:,:,nz-4),nx,ny)
call modutils_3darray2VTK(DEVICE_POT,dx,"pot3d.vtk")



print*,"finish"
call device_free()

contains



end program transporter
