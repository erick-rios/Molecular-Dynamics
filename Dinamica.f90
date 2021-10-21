module variables

	implicit none
	integer :: nat,npasos,reini
	double precision :: lz,lx,ly
	double precision :: sigma,rho,eps,ulj,rcut, utot, ukin
	double precision :: dt,tempi,temp
	double precision, dimension(:), allocatable :: rx,ry,rz
	double precision, dimension(:), allocatable :: vx,vy,vz
	double precision, dimension(:), allocatable :: fx,fy,fz
	double precision, dimension(:), allocatable :: masa
	double precision, dimension(:), allocatable :: hisrhoz,hisrhoy,hisrhox

end module variables

program MolecularDynamics
	
	use variables
	implicit none

	call LeeDatos
	if (reini == 1) then
		call Guarda(0)
	elseif (reini == 0) then

		call Memoria
		call Posiciones
		call Velocidades
	
	endif

	call Pelicula(0)
	call Perfil(0)
	call Fuerza
	call Repetir
	call Pelicula(2)
	call Perfil(2)
	call Guarda(1)

end program MolecularDynamics


subroutine Memoria
	
	use variables
	implicit none

	allocate(rx(nat), ry(nat), rz(nat))
	allocate(vx(nat), vy(nat), vz(nat))
	allocate(fx(nat), fy(nat), fz(nat))
	allocate(masa(nat))
	masa = 1.0d0

end subroutine Memoria

subroutine LeeDatos
	use variables
	implicit none


	open(1, file = "run.dat", status = "old", action = "read")


		read(1,*) reini
		read(1,*) nat
		read(1,*) rho
		read(1,*) lz
		read(1,*) dt
		read(1,*) sigma
		read(1,*) eps
		read(1,*) npasos
		read(1,*) rcut
		read(1,*) temp

	close(1)

end subroutine LeeDatos

subroutine Posiciones
	use variables
	implicit none 

	integer :: i,j,k,l
	integer :: nx,ny,nz

	lx = sqrt(dble(nat)/(rho*lz))
	ly = lx 

	k = 0

	nx = int(lx/sigma)
	ny = int(ly/sigma)
	nz = int(lz/sigma)

	do j = 1, nz
		do i = 1, ny
			do l = 1, nx
				k = k + 1
				if (k<=nat) then
					write(*,*) k
					rx(k) = dble(l*sigma)
					ry(k) = dble(i*sigma)
					rz(k) = dble(j*sigma)
				endif
			enddo
		enddo
	enddo
	
end subroutine Posiciones



subroutine Velocidades
	
	use variables 
	implicit none

	integer :: i 
	double precision :: rndx, rndy, rndz

	do i = 1, nat

		call random_number(rndx)
		call random_number(rndy)
		call random_number(rndz)

		vx(i) = 2.0d0*rndx - 1.0d0
		vy(i) = 2.0d0*rndy - 1.0d0
		vz(i) = 2.0d0*rndz - 1.0d0

	enddo

end subroutine Velocidades



subroutine Fuerza 

    use variables
    implicit none

    integer :: i,j 
    double precision :: dx, dy, dz
    double precision :: duij, rij

    fx  = 0.0d0
    fy  = 0.0d0
    fz  = 0.0d0
    ulj = 0.0d0

    do i = 1, nat - 1
        do j = i + 1,nat

            dx    = rx(i) - rx(j)
            dy    = ry(i) - ry(j)
            dz    = rz(i) - rz(j)

            !Mínima imágen

            if (dx > lx*0.50d0) then
            	dx = dx - lx
            else 
            	if (dx < -lx*0.50d0) dx = dx + lx
            end if
            if (dy > ly*0.50d0) then
            	dy = dy - ly
            else 
            	if (dy < -ly*0.50d0) dy = dy + ly
            end if 
            if (dz > lz*0.50d0) then
            	dz = dz -lz
            else
            	if (dz < -lz*0.50d0) dz = dz + lz
            endif


            !Termina Minima Imagen

            rij   = sqrt(dx**2 + dy**2 + dz**2)

            if (rij < rcut) then
	            ulj   = ulj + 4.0d0*eps*((sigma/rij)**12-(sigma/rij)**6)
	            duij  = 48.0d0*eps*((sigma/rij)**12 - 0.50d0*(sigma/rij)**6)/rij

	            fx(i) = fx(i) + duij*(dx/rij)
	            fy(i) = fy(i) + duij*(dy/rij)
	            fz(i) = fz(i) + duij*(dz/rij)

	            fx(j) = fx(j) - duij*(dx/rij)
	            fy(j) = fy(j) - duij*(dy/rij)
	            fz(j) = fz(j) - duij*(dz/rij)
	        end if 
        enddo
    enddo

end subroutine Fuerza


subroutine Repetir
	
	use variables
	implicit none

	integer :: paso,i 
	open(9, file = 'TRJ.dat', status = 'replace', action = 'write')

	do paso = 1, npasos
		
		!Velocity-verlet
		
		vx = vx + fx*dt/(2.0d0*masa)
		vy = vy + fy*dt/(2.0d0*masa)
		vz = vz + fz*dt/(2.0d0*masa)

		rx = rx + vx*dt
		ry = ry + vy*dt
		rz = rz + vz*dt

		!Condiciones peridicas de fronteras (PBC)
		
		where (rx < 0.0d0) rx = rx + lx
		where (ry < 0.0d0) ry = ry + ly
		where (rz < 0.0d0) rz = rz + lz

		where (rx > lx) rx = rx - lx
		where (ry > ly) ry = ry - ly
		where (rz > lz) rz = rz - lz

		call Fuerza

		vx = vx + fx*dt/(2.0d0*masa)
		vy = vy + fy*dt/(2.0d0*masa)
		vz = vz + fz*dt/(2.0d0*masa)

		ukin  = sum(0.50d0*masa*vx**2 + 0.50d0*masa*vy**2 + 0.50d0*masa*vz**2)
		tempi = 2.0d0*ukin/dble(3*nat) !Temperatura instantanea
		
		!Termostato: rescalamiento de velocidades
		vx = sqrt(temp/tempi)*vx
		vy = sqrt(temp/tempi)*vy
		vz = sqrt(temp/tempi)*vz


		vx = vx - sum(vx*masa)/dble(nat)
		vy = vy - sum(vy*masa)/dble(nat)
		vz = vz - sum(vz*masa)/dble(nat)

		if (mod(paso,10000) == 0) call Pelicula(1)
		if (mod(paso,100)   == 0) call Perfil(1)
		
		utot = ukin + ulj
		
		if (mod(paso,100)   == 0)	write(*,100)paso, utot/dble(nat), ukin/dble(nat), ulj/dble(nat), tempi
		if (mod(paso,100)   == 0)	write(9,100)paso, utot/dble(nat), ukin/dble(nat), ulj/dble(nat), tempi	

	enddo

	close(9)

	100 format(i10,4f12.5)

end subroutine Repetir 


subroutine Pelicula(flag)
	use variables
	implicit none 

	integer :: i, flag

	if (flag == 0) then
		
		open(2, file = "animacion.dat", status = "replace", action = "write")
		write(2,*)nat
		write(2,200)'Lattice="', lx, 0.0, 0.0, 0.0, ly, 0.0, 0.0, 0.0, lz,'"'
		do i = 1, nat
			write(2,100) 'C', rx(i), ry(i), rz(i)
		enddo

	elseif (flag == 1) then
		
		write(2,*) nat
		write(2,200)'Lattice="', lx, 0.0, 0.0, 0.0, ly, 0.0, 0.0, 0.0, lz,'"'

		do i = 1, nat
			write(2,100) 'C', rx(i),ry(i),rz(i)
		enddo
	else
		close(2)
	endif 

	100 format(a, 3f12.4)
	200 format(a, 9f5.1, a)


end subroutine Pelicula


subroutine Perfil(flag)
    
    use variables
    implicit none

    integer:: i, flag,maxrho,contrho,bin
    double precision::deltar,maxl
    save :: deltar, contrho,maxrho

    if(flag == 0)then

    	deltar = 0.050d0
    	maxl   = max(lx,ly,lz)
    	maxrho = int(lz/deltar)
    	allocate(hisrhoz(0:maxrho))
    	allocate(hisrhoy(0:maxrho))
    	allocate(hisrhox(0:maxrho))

    	hisrhoz = 0.0d0
    	hisrhoy = 0.0d0
    	hisrhox = 0.0d0

	    contrho = 0.0d0

    elseif(flag == 1)then

	    contrho = contrho + 1
    	do i = 1, nat

    		bin          = int(rz(i)/deltar)
    		hisrhoz(bin) = hisrhoz(bin)+1

   			bin          =int(ry(i)/deltar)
    		hisrhoy(bin) = hisrhoy(bin) + 1

    		bin          = int(rx(i)/deltar)
    		hisrhox(bin) = hisrhox(bin) + 1

    	end do

    else
    	open(8, file = 'perfil.dat', status = 'replace', action = 'write')

    	do i=0,maxrho-1
    		
    		write(8,*)dble(i)*deltar,hisrhoz(i)/(lx*ly*deltar*contrho),hisrhoy(i)/(lx*lz*deltar*contrho),hisrhox(i)/(lz*ly*deltar*contrho)

    	end do

    	close(8)

    end if
end subroutine Perfil


subroutine Guarda(flag)
    use variables
    implicit none

    integer:: i,flag

    if(flag==0)then

    	open(9,file='posi.dat',status='old',action='read')
  		read(9,*)nat
  		read(9,*)lx,ly,lz
  		
  		call Memoria
    	
    	do i=1,nat
        	read(9,*)rx(i),ry(i),rz(i)
        	read(9,*)vx(i),vy(i),vz(i)
        end do
        
        close(9)


	elseif(flag==1)then

    	open(9,file='posi.dat',status='replace',action='write')
    	write(9,*)nat
    	write(9,*)lx,ly,lz
    	do i=1,nat
        	write(9,*)rx(i),ry(i),rz(i)
        	write(9,*)vx(i),vy(i),vz(i)
        end do
        close(9)

    endif

end subroutine Guarda