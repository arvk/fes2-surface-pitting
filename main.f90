module commondata
  implicit none
  save

  type :: prol
     integer :: fx,fy,fz,tx,ty,tz
     real*8 :: prob
  end type prol

  type(prol), dimension(:,:), allocatable :: proc1, proc2, proc3
  integer, dimension(:,:,:), allocatable :: posns ! rename to positions
  integer, dimension(:,:,:), allocatable :: initposns ! rename to initial_positions
  integer, dimension(:), allocatable :: noatoms ! rename to number_of_atoms
  integer :: nompa = 7 ! rename to no_of_moves_per_atom
  integer :: cs_x, cs_y, cs_z   ! rename to cell_size_x
  integer :: nospec ! rename to number_of_species
  integer :: noiterations,nomcsteps
  integer :: lx,ly,lz
  integer :: number_of_images = 250
  integer :: temperature
character*3 :: temperature_string
end module commondata


program run_all_temps

include 'mpif.h'

integer :: ierr, num_procs, my_rank

      call MPI_INIT ( ierr )
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)


        if (my_rank .eq. 0) then
           call kmc(298)
        else if (my_rank .eq. 1) then
           call kmc(299)
        else if (my_rank .eq. 2) then
           call kmc(393)
        else if (my_rank .eq. 3) then
           call kmc(394)
        else if (my_rank .eq. 4) then
           call kmc(443)
        else if (my_rank .eq. 5) then
           call kmc(444)
        else if (my_rank .eq. 6) then
           call kmc(483)
        else if (my_rank .eq. 7) then
           call kmc(484)
        else if (my_rank .eq. 8) then
           call kmc(513)
        else if (my_rank .eq. 9) then
           call kmc(514)
        else if (my_rank .eq. 10) then
           call kmc(543)
        else if (my_rank .eq. 11) then
           call kmc(544)
        else if (my_rank .eq. 12) then
           call kmc(573)
        else if (my_rank .eq. 13) then
           call kmc(574)
        else if (my_rank .eq. 14) then
           call kmc(603)
        else if (my_rank .eq. 15) then
           call kmc(604)
        end if

      call MPI_FINALIZE ( ierr )

end program

subroutine kmc(dummy_temperature)
  use commondata
  implicit none

  integer :: atom_id = 0
  integer :: lp, lq, lr, ls ! looping index, rename to loop_p, ...
  real :: rd1, rd2, rd3 ! dummy variables. reduce as much as possible
  integer :: pt ! change to process_type
  real*8 :: p1_prob, p2_prob, p3_prob, p4_prob, p5_prob, p6_prob, p7_prob
  real*8 :: fl_1, fl_2
  integer :: sp1, sp2 ! rename to selected_process_1...
  integer :: rand1, rand2, rand3
  character*12 :: filename
  real*8 :: mc_step_timescale
  integer :: current_mc_step, current_iteration
  real*8 :: get_displacement, total_displacement
  integer :: q1, q2, q3, q4
  integer :: rand_seed_size
  integer, dimension(1) :: rand_seed
  real*8, dimension(:), allocatable :: monomer_s, surface_s
  integer :: xloop, yloop, xl, yl
  real*8 :: neighborhood
  real*8 :: time_of_experiment
  integer :: image_id_mc, image_id_iter
  character*5 :: image_id_mc_string, image_id_iter_string
  integer :: temp_monomer_count
  character(LEN=10) :: date, time
  integer :: dates(8)
  integer :: dummy_temperature

  temperature = dummy_temperature
  write(temperature_string,'(I3)') temperature

  call read_parameters()

  allocate(proc1(sum(noatoms),nompa)); allocate(proc2(sum(noatoms),nompa)); allocate(proc3(sum(noatoms),nompa))
  allocate(monomer_s(nomcsteps))
  allocate(surface_s(nomcsteps))

  call initialize_positions()
  call initialize_plist()

  CALL DATE_AND_TIME(date,time,VALUES=dates)
  rand_seed=dates(1)+dates(2)+dates(3)+dates(5)+dates(6)+dates(7)+dates(8)

  !rand_seed(1) = time()
  rand_seed_size = 1
  call random_seed(size = rand_seed_size)  
  call random_seed(put = rand_seed)

  do current_iteration = 1,noiterations

     time_of_experiment = 0.0d0
     monomer_s = 0.0d0 ; surface_s = 0.0d0

     posns = initposns
     call initialize_plist()

     do current_mc_step = 1,nomcsteps

        p1_prob = sum(proc1%prob) ; p2_prob = sum(proc2%prob) ; p3_prob = sum(proc3%prob)

        call random_number(rd1)

        fl_1 = rd1*(p1_prob+p2_prob+p3_prob)

        if (fl_1 .le. p1_prob) then
           sp1 = 1
!           write(6,*) dummy_temperature, 'CHOSEN PROCESS: 1'
        else if ((fl_1 .gt. p1_prob) .and. (fl_1 .le. p1_prob+p2_prob)) then
           sp1 = 2
!           write(6,*) dummy_temperature, 'CHOSEN PROCESS: 2'
        else
           sp1 = 3
!           write(6,*) dummy_temperature, 'CHOSEN PROCESS: 3'
        end if

        call make_move(sp1)


        call random_number(rd3)

        do while (rd3 .eq. 0 )
           call random_number(rd3)
        end do

        time_of_experiment = time_of_experiment - ((1.0d0/(p1_prob+p2_prob+p3_prob))*log(rd3))



        do rand1 = 1,cs_x
           do rand2 = 1,cs_y
              if (mod(rand1+rand2+2,2).gt.0) then ! Choose S sites
                 if (posns(rand1,rand2,2) .eq. 1) then ! Find monomer vacancies
                    monomer_s(current_mc_step) = monomer_s(current_mc_step) + 1.0d0
                 else ! If the site is not a monomer, i.e. it is a dimer or defect-free
                    surface_s(current_mc_step) = surface_s(current_mc_step) + 1.0d0
                 end if
              end if
           end do
        end do




        if (mod(current_mc_step, ((nomcsteps/number_of_images)-1) ).eq.0) then
           image_id_mc = (current_mc_step/((nomcsteps/number_of_images)-1))
           write(image_id_mc_string,'(I5.5)') image_id_mc
           write(image_id_iter_string,'(I5.5)') current_iteration

           open (unit = 700, file = trim(temperature_string)//'K/forplot.txt')


           write(700,*) "# Experiment Time : ", time_of_experiment, " seconds. "


           temp_monomer_count = 0

           do rand1 = 1,cs_x
              do rand2 = 1,cs_y
                 if (mod(rand1+rand2+2,2).gt.0) then ! Choose S sites
                    if (posns(rand1,rand2,2) .eq. 1) then ! Find monomer vacancies
                       temp_monomer_count = temp_monomer_count + 1
                    end if
                 end if
              end do
           end do


           write(700,*) "# Number of monomer species : ", temp_monomer_count


           do rand1 = 1,cs_x
              do rand2 = 1,cs_y

                 if (posns(rand1,rand2,2) .ne. 0) then

                    write(700,'(I10,I10,F9.3)') rand1, rand2, 2.0d0

                 else

                    neighborhood = 0

                    do xloop = rand1-1,rand1+1
                       do yloop = rand2-1,rand2+1

                          if (xloop<1) then
                             xl = cs_x
                          elseif (xloop>cs_x) then
                             xl = 1
                          else
                             xl = xloop
                          end if

                          if (yloop<1) then
                             yl = cs_y
                          elseif (yloop>cs_y) then
                             yl = 1
                          else
                             yl = yloop
                          end if

                          neighborhood = neighborhood + posns(xl,yl,2)



                       end do
                    end do

                    write(700,'(I10,I10,F9.3)') rand1, rand2, 0.0d0+(neighborhood/9.0d0)

                 end if


              end do
           end do


           close(700)



           call system('cp '//trim(temperature_string)//'K/forplot.txt '//trim(temperature_string)//'K/fp-'//trim(image_id_iter_string)//'-'//trim(image_id_mc_string)//'.txt')

        end if


     end do
  end do



end subroutine kmc



















subroutine make_move(sp1)
  use commondata
  implicit none
  integer :: sp1, rn1, rn2, lp
  integer :: moved_to_x, moved_to_y, moved_to_z
  integer :: moved_from_x, moved_from_y, moved_from_z
  integer :: botval ! PLEASE change variable name
  integer :: lattice_id, mpa
  real*8 :: rd1,cutoff, tripper

  integer :: rand_seed_size
  integer, dimension(1) :: rand_seed

  call random_number(rd1)

  if (sp1 == 1) then

     tripper = 0.0d0
     cutoff = rd1*sum(proc1%prob)

     do lattice_id = 1,sum(noatoms)
        do mpa = 1,nompa

           tripper = tripper + proc1(lattice_id, mpa)%prob

           if (tripper .ge. cutoff) then

              rn1 = lattice_id
              rn2 = mpa

              moved_from_x = proc1(rn1,rn2)%fx ; moved_from_y = proc1(rn1,rn2)%fy ; moved_from_z = proc1(rn1,rn2)%fz
              moved_to_x = proc1(rn1,rn2)%tx ; moved_to_y = proc1(rn1,rn2)%ty ; moved_to_z = proc1(rn1,rn2)%tz

              if (mod(moved_from_x+moved_from_y+moved_from_z,2).gt.0) then
posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)=posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)-1
              else
posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)=posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)-2
              end if

              call rfpls(moved_from_x,moved_from_y,moved_from_z) ; call add_to_plist(moved_from_x,moved_from_y,moved_from_z)
              call rfpls(moved_to_x,moved_to_y,moved_to_z) ; call add_to_plist(moved_to_x,moved_to_y,moved_to_z)
              call upnl(moved_from_x,moved_from_y,moved_from_z) ; call upnl(moved_to_x,moved_to_y,moved_to_z)

              exit
           else

           end if
        end do
        if (tripper .ge. cutoff) exit
     end do


  else if (sp1 == 2) then

     tripper = 0.0d0
     cutoff = rd1*sum(proc2%prob)

     do lattice_id = 1,sum(noatoms)
        do mpa = 1,nompa

           tripper = tripper + proc2(lattice_id, mpa)%prob

           if (tripper .ge. cutoff) then

              rn1 = lattice_id
              rn2 = mpa

              moved_from_x = proc2(rn1,rn2)%fx ; moved_from_y = proc2(rn1,rn2)%fy ; moved_from_z = proc2(rn1,rn2)%fz
              moved_to_x = proc2(rn1,rn2)%tx ; moved_to_y = proc2(rn1,rn2)%ty ; moved_to_z = proc2(rn1,rn2)%tz

              if (mod(moved_from_x+moved_from_y+moved_from_z,2).gt.0) then
posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)=posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)+1
              else
posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)=posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)+2
              end if

              call rfpls(moved_from_x,moved_from_y,moved_from_z) ; call add_to_plist(moved_from_x,moved_from_y,moved_from_z)
              call rfpls(moved_to_x,moved_to_y,moved_to_z) ; call add_to_plist(moved_to_x,moved_to_y,moved_to_z)
              call upnl(moved_from_x,moved_from_y,moved_from_z) ; call upnl(moved_to_x,moved_to_y,moved_to_z)

              exit
           else

           end if
        end do
        if (tripper .ge. cutoff) exit
     end do





  else if (sp1 == 3) then

     tripper = 0.0d0
     cutoff = rd1*sum(proc3%prob)

     do lattice_id = 1,sum(noatoms)
        do mpa = 1,nompa

           tripper = tripper + proc3(lattice_id, mpa)%prob

           if (tripper .ge. cutoff) then

              rn1 = lattice_id
              rn2 = mpa

              moved_from_x = proc3(rn1,rn2)%fx ; moved_from_y = proc3(rn1,rn2)%fy ; moved_from_z = proc3(rn1,rn2)%fz
              moved_to_x = proc3(rn1,rn2)%tx ; moved_to_y = proc3(rn1,rn2)%ty ; moved_to_z = proc3(rn1,rn2)%tz

              if (mod(moved_from_x+moved_from_y+moved_from_z,2).gt.0) then
posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)=posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)+1
posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)=posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)-1
              else
posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)=posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)+2
posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)=posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)-2
              end if

              call rfpls(moved_from_x,moved_from_y,moved_from_z) ; call add_to_plist(moved_from_x,moved_from_y,moved_from_z)
              call rfpls(moved_to_x,moved_to_y,moved_to_z) ; call add_to_plist(moved_to_x,moved_to_y,moved_to_z)
              call upnl(moved_from_x,moved_from_y,moved_from_z) ; call upnl(moved_to_x,moved_to_y,moved_to_z)

              exit
           else


           end if
        end do
        if (tripper .ge. cutoff) exit
     end do

  end if


end subroutine make_move

subroutine initialize_positions()
  use commondata
  implicit none

  integer :: rand_seed_size
  integer, dimension(1) :: rand_seed
  integer :: atom_type, numb_of_atoms, random_integer
  integer :: x_pos, y_pos, z_pos, x_loop, y_loop, z_loop
  real*8 :: random_double

  ! Distribute atoms
  allocate(initposns(cs_x,cs_y,cs_z))
  allocate(posns(cs_x,cs_y,cs_z))

  do x_loop = 1, cs_x
     do y_loop = 1, cs_y
        posns(x_loop,y_loop, cs_z) = 2
        do z_loop = 1, cs_z-1
           posns(x_loop,y_loop,z_loop) = 2
        end do
     end do
  end do


  initposns = posns




end subroutine initialize_positions
subroutine initialize_plist()
  use commondata
  implicit none

  integer :: atom_id, pid
  integer :: x_pos, y_pos, z_pos

  do atom_id = 1, sum(noatoms)
     do pid = 1, nompa

        proc1(atom_id,pid)%fx = 0 ; proc1(atom_id,pid)%fy = 0 ; proc1(atom_id,pid)%fz = 0 
        proc1(atom_id,pid)%tx = 0 ; proc1(atom_id,pid)%ty = 0 ; proc1(atom_id,pid)%tz = 0 
        proc1(atom_id,pid)%prob = 0.0d0

        proc2(atom_id,pid)%fx = 0 ; proc2(atom_id,pid)%fy = 0 ; proc2(atom_id,pid)%fz = 0 
        proc2(atom_id,pid)%tx = 0 ; proc2(atom_id,pid)%ty = 0 ; proc2(atom_id,pid)%tz = 0
        proc2(atom_id,pid)%prob = 0.0d0

        proc3(atom_id,pid)%fx = 0 ; proc3(atom_id,pid)%fy = 0 ; proc3(atom_id,pid)%fz = 0 
        proc3(atom_id,pid)%tx = 0 ; proc3(atom_id,pid)%ty = 0 ; proc3(atom_id,pid)%tz = 0 
        proc3(atom_id,pid)%prob = 0.0d0

     end do
  end do


  atom_id = 0

  do z_pos = 1, cs_z
     do y_pos = 1, cs_y
        do x_pos = 1, cs_x
              call add_to_plist(x_pos,y_pos,z_pos)
        end do
     end do
  end do

end subroutine initialize_plist

subroutine read_parameters()
  use commondata
  implicit none

write(6,*) trim(temperature_string)

  ! Read lattice parameters
  call system("cat "//trim(temperature_string)//"K/param.in | grep ^LATTICE | sed 's/LATTICE//g'| sed 's/=//g' > "//trim(temperature_string)//"K/.toget-lattice")
  open(unit = 7001, file = ""//trim(temperature_string)//"K/.toget-lattice", status = 'old') ; read(7001,*) cs_x, cs_y, cs_z ; close(7001)
write(6,*) 'Hi'
  call system ("rm "//trim(temperature_string)//"K/.toget-lattice")

  ! Read number of species and number of atoms in each species
  call system("cat "//trim(temperature_string)//"K/param.in | grep ^SPECIES | sed 's/SPECIES//g'| sed 's/=//g' | awk '{print NF}' > "//trim(temperature_string)//"K/.toget-numbofspec")
  open(unit = 7002, file = ""//trim(temperature_string)//"K/.toget-numbofspec", status = 'old') ; read(7002,*) nospec ; close(7002)
  call system ("rm "//trim(temperature_string)//"K/.toget-numbofspec")
  allocate(noatoms(nospec))
  call system("cat "//trim(temperature_string)//"K/param.in | grep ^SPECIES | sed 's/SPECIES//g'| sed 's/=//g' > "//trim(temperature_string)//"K/.toget-numbofatoms")
  open(unit = 7003, file = trim(temperature_string)//'K/.toget-numbofatoms', status = 'old') ; read(7003,*) noatoms ; close(7003)
  call system ("rm "//trim(temperature_string)//"K/.toget-numbofatoms")

  ! Read number of steps
  call system("cat "//trim(temperature_string)//"K/param.in | grep ^MCSTEPS | sed 's/MCSTEPS//g'| sed 's/=//g' > "//trim(temperature_string)//"K/.toget-mcsteps")
  open(unit = 7004, file = trim(temperature_string)//'K/.toget-mcsteps', status = 'old') ; read(7004,*) nomcsteps ; close(7004)
  call system ("rm "//trim(temperature_string)//"K/.toget-mcsteps")

  ! Read number of iterations
  call system("cat "//trim(temperature_string)//"K/param.in | grep ^ITERATIONS | sed 's/ITERATIONS//g'| sed 's/=//g' > "//trim(temperature_string)//"K/.toget-iterations")
  open(unit = 7005, file = trim(temperature_string)//'K/.toget-iterations', status = 'old') ; read(7005,*) noiterations ; close(7005)
  call system ("rm "//trim(temperature_string)//"K/.toget-iterations")

end subroutine read_parameters

subroutine add_to_plist(lp,lq,lr)
  use commondata
  implicit none
  integer :: lp,lq,lr, atom_id, pt
  integer :: loop1, loop2 ! PLEASE change variable names
  integer :: xloop, yloop ! PLEASE change variable names
  integer :: xp, yp, zp ! PLEASE change variable names

! Create vacancies
  if (posns(lp,lq,lr) .ne. 0) then 
     pt = 1 ; call make_ptl(pt,lp,lq,lr,lp,lq,lr)
  end if

! Fill vacancies
  if (posns(lp,lq,lr) .ne. 2) then 
     pt = 2 ; call make_ptl(pt,lp,lq,lr,lp,lq,lr)
  end if

! Diffusion
  if (posns(lp,lq,lr) .ne. 2) then 

do xloop = lp-1,lp+1,2
do yloop = lq-1,lq+1,2

if (xloop<1)  then
loop1 = cs_x
elseif (xloop>cs_x)then
loop1 = 1
else
loop1 = xloop
end if

if (yloop<1)  then
loop2 = cs_x
elseif (yloop>cs_x)then
loop2 = 1
else
loop2 = yloop
end if

  if (posns(loop1,loop2,lr) .ne. 0) then 
     pt = 3 ; call make_ptl(pt,lp,lq,lr,loop1,loop2,lr)
  end if

end do
end do

  end if

end subroutine add_to_plist







subroutine make_ptl(pt, fromx, fromy, fromz, tox, toy, toz)
  use commondata
  integer, intent(in) :: pt, fromx, fromy, fromz, tox, toy, toz
  integer :: atom_id
  real*8 :: kbT
  real*8 :: process_barrier
  real*8 :: s_surf, s_bulk, fe_surf, fe_bulk
  real*8 :: neighborhood,neighborhoodfrom,neighborhoodto
  integer :: xloop, yloop

  atom_id = ((fromz-1)*(cs_x*cs_y)) + ((fromy-1)*cs_x) + fromx

  kbT = 0.000086173423d0 * temperature


  if (pt == 1) then 

  s_surf = 0.075d0 ; fe_surf = 0.075d0
  s_bulk = 0.075d0 ; fe_bulk = 0.075d0

     if (mod(fromx+fromy+fromz,2).gt.0) then
        process_barrier = 0.0d0 - (s_surf + ((s_bulk-s_surf)*(cs_z-fromz)))
     else
        process_barrier = 0.0d0 - (fe_surf + ((fe_bulk-fe_surf)*(cs_z-fromz)))
     end if

! Neighborhood sampling
neighborhood = 0.0d0
do xloop = tox-1,tox+1,2
if (xloop<1)  then
neighborhood = neighborhood + posns(cs_x,toy,toz)
elseif (xloop>cs_x)then
neighborhood = neighborhood + posns(1,toy,toz)
else
neighborhood = neighborhood + posns(xloop,toy,toz)
end if
end do

do yloop = toy-1,toy+1,2
if (yloop<1)  then
neighborhood = neighborhood + posns(tox,cs_y,toz)
elseif (yloop>cs_y)then
neighborhood = neighborhood + posns(tox,1,toz)
else
neighborhood = neighborhood + posns(tox,yloop,toz)
end if
end do
! Done neighborhood sampling

 process_barrier = process_barrier * ((neighborhood/8.0d0))

     if (mod(fromx+fromy+fromz,2).gt.0) then
         process_barrier = process_barrier - 1.47d0
     else
         process_barrier = process_barrier - 1.47d0 - 0.08d0
     end if


     do lp = 1, nompa
        if (proc1(atom_id, lp)%fx == 0) then
           proc1(atom_id,lp)%fx = fromx ; proc1(atom_id,lp)%fy = fromy ; proc1(atom_id,lp)%fz = fromz
           proc1(atom_id,lp)%tx = tox ; proc1(atom_id,lp)%ty = toy ; proc1(atom_id,lp)%tz = toz
           proc1(atom_id,lp)%prob = (10**6)*dexp(process_barrier/kbT)*(10**6)
        exit
        end if
     end do










  else if (pt == 2) then 


  s_surf = 0.075d0 ; fe_surf = 0.075d0
  s_bulk = 0.075d0 ; fe_bulk = 0.075d0

     if (mod(fromx+fromy+fromz,2).gt.0) then
        process_barrier = 0.0d0 + (s_surf + ((s_bulk-s_surf)*(cs_z-fromz)))
     else
        process_barrier = 0.0d0 + (fe_surf + ((fe_bulk-fe_surf)*(cs_z-fromz)))
     end if

! Neighborhood sampling
neighborhood = 0.0d0
do xloop = tox-1,tox+1,2
if (xloop<1)  then
neighborhood = neighborhood + posns(cs_x,toy,toz)
elseif (xloop>cs_x)then
neighborhood = neighborhood + posns(1,toy,toz)
else
neighborhood = neighborhood + posns(xloop,toy,toz)
end if
end do

do yloop = toy-1,toy+1,2
if (yloop<1)  then
neighborhood = neighborhood + posns(tox,cs_y,toz)
elseif (yloop>cs_y)then
neighborhood = neighborhood + posns(tox,1,toz)
else
neighborhood = neighborhood + posns(tox,yloop,toz)
end if
end do
! Done neighborhood sampling

 process_barrier = process_barrier * ((neighborhood/8.0d0))

     if (mod(fromx+fromy+fromz,2).gt.0) then
        process_barrier = process_barrier - 1.600
     else
        process_barrier = process_barrier - 1.60d0 + 0.08d0
     end if


     do lp = 1, nompa
        if (proc2(atom_id, lp)%fx == 0) then
           proc2(atom_id,lp)%fx = fromx ; proc2(atom_id,lp)%fy = fromy ; proc2(atom_id,lp)%fz = fromz
           proc2(atom_id,lp)%tx = tox ; proc2(atom_id,lp)%ty = toy ; proc2(atom_id,lp)%tz = toz
           proc2(atom_id,lp)%prob = (10**6)*dexp(process_barrier/kbT)*(10**6)
        exit
        end if
     end do






   else if (pt == 3) then 


! Neighborhood sampling
neighborhood = 0.0d0
do xloop = tox-1,tox+1,2
if (xloop<1)  then
neighborhood = neighborhood + posns(cs_x,toy,toz)
elseif (xloop>cs_x)then
neighborhood = neighborhood + posns(1,toy,toz)
else
neighborhood = neighborhood + posns(xloop,toy,toz)
end if
end do

do yloop = toy-1,toy+1,2
if (yloop<1)  then
neighborhood = neighborhood + posns(tox,cs_y,toz)
elseif (yloop>cs_y)then
neighborhood = neighborhood + posns(tox,1,toz)
else
neighborhood = neighborhood + posns(tox,yloop,toz)
end if
end do

neighborhoodto = neighborhood


neighborhood = 0.0d0
do xloop = fromx-1,fromx+1,2
if (xloop<1)  then
neighborhood = neighborhood + posns(cs_x,fromy,fromz)
elseif (xloop>cs_x)then
neighborhood = neighborhood + posns(1,fromy,fromz)
else
neighborhood = neighborhood + posns(xloop,fromy,fromz)
end if
end do

do yloop = fromy-1,fromy+1,2
if (yloop<1)  then
neighborhood = neighborhood + posns(fromx,cs_y,fromz)
elseif (yloop>cs_y)then
neighborhood = neighborhood + posns(fromx,1,fromz)
else
neighborhood = neighborhood + posns(fromx,yloop,fromz)
end if
end do

neighborhoodfrom = neighborhood

! Done neighborhood sampling

  s_surf = 0.00d0 ; fe_surf = 0.00d0
  s_bulk = 0.00d0 ; fe_bulk = 0.00d0

     if (mod(fromx+fromy+fromz,2).gt.0) then
        process_barrier = 0.0d0 - (s_surf + ((s_bulk-s_surf)*(cs_z-fromz)))
     else
        process_barrier = 0.0d0 - (fe_surf + ((fe_bulk-fe_surf)*(cs_z-fromz)))
     end if

     process_barrier = process_barrier + (0.1125d0 * ((neighborhoodfrom/8.0d0) - (neighborhoodto/8.0d0) )) - 1.35d0

     do lp = 1, nompa
        if (proc3(atom_id, lp)%fx == 0) then
           proc3(atom_id,lp)%fx = fromx ; proc3(atom_id,lp)%fy = fromy ; proc3(atom_id,lp)%fz = fromz
           proc3(atom_id,lp)%tx = tox ; proc3(atom_id,lp)%ty = toy ; proc3(atom_id,lp)%tz = toz 
           proc3(atom_id,lp)%prob = (10**6)*dexp(process_barrier/kbT)*(10**6)
        exit
        end if
     end do

  end if


end subroutine make_ptl








subroutine rfpls(x_pos,y_pos,z_pos)
  use commondata
  implicit none
  integer :: atom_id
  integer :: x_pos, y_pos, z_pos
  integer :: lq

atom_id = ((z_pos-1)*(cs_x*cs_y)) + ((y_pos-1)*(cs_x)) + x_pos

  do lq = 1, nompa
     proc1(atom_id,lq)%fx = 0 ; proc1(atom_id,lq)%fy = 0 ; proc1(atom_id,lq)%fz = 0 
     proc1(atom_id,lq)%tx = 0 ; proc1(atom_id,lq)%ty = 0 ; proc1(atom_id,lq)%tz = 0 ; proc1(atom_id,lq)%prob = 0.0d0
     proc2(atom_id,lq)%fx = 0 ; proc2(atom_id,lq)%fy = 0 ; proc2(atom_id,lq)%fz = 0 
     proc2(atom_id,lq)%tx = 0 ; proc2(atom_id,lq)%ty = 0 ; proc2(atom_id,lq)%tz = 0 ; proc2(atom_id,lq)%prob = 0.0d0
     proc3(atom_id,lq)%fx = 0 ; proc3(atom_id,lq)%fy = 0 ; proc3(atom_id,lq)%fz = 0 
     proc3(atom_id,lq)%tx = 0 ; proc3(atom_id,lq)%ty = 0 ; proc3(atom_id,lq)%tz = 0 ; proc3(atom_id,lq)%prob = 0.0d0
  end do

end subroutine rfpls




subroutine upnl(ma_x,ma_y,ma_z)
  use commondata
  implicit none
  integer :: ma_x,ma_y,ma_z, lp,lq,lr, newlp,newlq,newlr, atom_id, getatomid

  do lp = ma_x-1, ma_x+1
     do lq = ma_y-1, ma_y+1
        do lr = ma_z-1, ma_z+1

           if (lp .lt. 1) then
              newlp = cs_x
           else if (lp .gt. cs_x) then
              newlp = 1
           else
              newlp = lp
           end if

           if (lq .lt. 1) then
              newlq = cs_y
           else if (lq .gt. cs_y) then
              newlq = 1
           else
              newlq = lq
           end if

           if (lr .lt. 1) then
              newlr = cs_z
           else if (lr .gt. cs_z) then
              newlr = 1
           else
              newlr = lr
           end if


           if (posns(newlp, newlq, newlr) .gt. -1) then
              if ((newlp .ne. ma_x) .or. (newlq .ne. ma_y) .or. (newlr .ne. ma_z))  then

                 atom_id = ((newlr-1)*(cs_x*cs_y)) + ((newlq-1)*cs_x) + newlp

                    call rfpls(newlp,newlq,newlr)
                    call add_to_plist(newlp,newlq,newlr)
              
              end if
           end if


        end do
     end do
  end do


end subroutine upnl

  
