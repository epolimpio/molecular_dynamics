module sub_regions
implicit none
public find_sub_regions

integer :: timeStepCount = 0
integer :: checkBoxPeriod = 100 !timesteps

contains
	subroutine update_sub_regions(r,v,a,maxX,maxY,maxZ)
		integer :: divisions !number of 'slices' through per dimension
		integer :: i,j,k,n,sizeR
		real(8) :: xWidth, yWidth, zWidth
		real (8), dimension(N,3) :: newR
    	real (8), dimension(N,3) :: newV
    	real (8), dimension(N,3) :: newA
		real (8), intent(inout), dimension(:,:) :: r
      	real (8), intent(inout), dimension(:,:) :: v
      	real (8), intent(inout), dimension(:,:) :: a
      	integer, dimension(:,:,:,:), allocatable :: regions

		newR = r
		newV = v
		newA = a

		sizeR = size(r,1)

		!populate the initial subregions with the lists of indices of particles

		xWidth = maxX/real(divisions)
		yWidth = maxY/real(divisions)
		zWidth = maxZ/real(divisions)
		
		do n=1,sizeR+1

			regions[integer(r(n,0)/xWidth),integer(r(n,1)/yWidth),integer(r(n,2)/zWidth)] = & 
			regions[integer(r(n,0)/xWidth),integer(r(n,1)/yWidth),integer(r(n,2)/zWidth)] + n !add n to the list of integers in the region referenced by those indices
		end do


		do i=1,divisions
			do j=1,divisions
				do k=1,divisions
					!for each subregion
						for pt in regions[i,j,k]
							totalForce = [0,0,0]
							for di = -1 to 1
								for dj = -1 to 1
									for dk = -1 to 1
										!for each adjacent subregion
										for pt2 in regions[i+di,j+dj,k+dk]
											totalForce += forceCompare(pt,pt2)!returns a force array ~ [Fx,Fy,Fz] of F on pt1 due to pt2
							!update newR,newV,newA for pt based on totalForce
		

		r = newR !change all pts to new values :: lockstep time update
		v = newV
		a = newA
		timeStepCount = timeStepCount + 1

		if ((timeStepCount % checkBoxPeriod)==0)
			!check each particle to see if it must switch subregion by switching all particles indescriminately to the correct subregion
			for i in 0 to divisions
				for j in 0 to divisions
					for k in 0 to divisions
					!for each subregion
					regions[i,j,k]=[]!empty list of integers

			for n in 0 to r.length()
				regions[int(r[n,0]/xWidth),int(r[n,1]/yWidth),int(r[n,2]/zWidth)] += n !add n to the list of integers in the region referenced by those indices



	end subroutine update_sub_regions
end module