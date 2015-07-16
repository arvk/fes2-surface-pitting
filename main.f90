MODULE commondata
  IMPLICIT NONE
  SAVE

  TYPE :: prol
     INTEGER :: fx,fy,fz,tx,ty,tz
     REAL*8 :: prob
  END TYPE prol

  TYPE(prol), DIMENSION(:,:), ALLOCATABLE :: proc1, proc2, proc3
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: posns ! rename to positions
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: initposns ! rename to initial_positions
  INTEGER, DIMENSION(:), ALLOCATABLE :: noatoms ! rename to number_of_atoms
  INTEGER :: nompa = 7 ! rename to no_of_moves_per_atom
  INTEGER :: cs_x, cs_y, cs_z   ! rename to cell_size_x
  INTEGER :: nospec ! rename to number_of_species
  INTEGER :: noiterations,nomcsteps
  INTEGER :: lx,ly,lz
  INTEGER :: number_of_images = 250
  INTEGER :: temperature
  CHARACTER*3 :: temperature_string
END MODULE commondata


PROGRAM run_all_temps

  INCLUDE 'mpif.h'

  INTEGER :: ierr, num_procs, my_rank

  CALL MPI_INIT ( ierr )
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)


  IF (my_rank .EQ. 0) THEN
     CALL kmc(298)
  ELSE IF (my_rank .EQ. 1) THEN
     CALL kmc(299)
  ELSE IF (my_rank .EQ. 2) THEN
     CALL kmc(393)
  ELSE IF (my_rank .EQ. 3) THEN
     CALL kmc(394)
  ELSE IF (my_rank .EQ. 4) THEN
     CALL kmc(443)
  ELSE IF (my_rank .EQ. 5) THEN
     CALL kmc(444)
  ELSE IF (my_rank .EQ. 6) THEN
     CALL kmc(483)
  ELSE IF (my_rank .EQ. 7) THEN
     CALL kmc(484)
  ELSE IF (my_rank .EQ. 8) THEN
     CALL kmc(513)
  ELSE IF (my_rank .EQ. 9) THEN
     CALL kmc(514)
  ELSE IF (my_rank .EQ. 10) THEN
     CALL kmc(543)
  ELSE IF (my_rank .EQ. 11) THEN
     CALL kmc(544)
  ELSE IF (my_rank .EQ. 12) THEN
     CALL kmc(573)
  ELSE IF (my_rank .EQ. 13) THEN
     CALL kmc(574)
  ELSE IF (my_rank .EQ. 14) THEN
     CALL kmc(603)
  ELSE IF (my_rank .EQ. 15) THEN
     CALL kmc(604)
  END IF

  CALL MPI_FINALIZE ( ierr )

END PROGRAM run_all_temps

SUBROUTINE kmc(dummy_temperature)
  USE commondata
  IMPLICIT NONE

  INTEGER :: atom_id = 0
  INTEGER :: lp, lq, lr, ls ! looping index, rename to loop_p, ...
  REAL :: rd1, rd2, rd3 ! dummy variables. reduce as much as possible
  INTEGER :: pt ! change to process_type
  REAL*8 :: p1_prob, p2_prob, p3_prob, p4_prob, p5_prob, p6_prob, p7_prob
  REAL*8 :: fl_1, fl_2
  INTEGER :: sp1, sp2 ! rename to selected_process_1...
  INTEGER :: rand1, rand2, rand3
  CHARACTER*12 :: filename
  REAL*8 :: mc_step_timescale
  INTEGER :: current_mc_step, current_iteration
  REAL*8 :: get_displacement, total_displacement
  INTEGER :: q1, q2, q3, q4
  INTEGER :: rand_seed_size
  INTEGER, DIMENSION(1) :: rand_seed
  REAL*8, DIMENSION(:), ALLOCATABLE :: monomer_s, surface_s
  INTEGER :: xloop, yloop, xl, yl
  REAL*8 :: neighborhood
  REAL*8 :: time_of_experiment
  INTEGER :: image_id_mc, image_id_iter
  CHARACTER*5 :: image_id_mc_string, image_id_iter_string
  INTEGER :: temp_monomer_count
  CHARACTER(LEN=10) :: date, time
  INTEGER :: dates(8)
  INTEGER :: dummy_temperature

  temperature = dummy_temperature
  WRITE(temperature_string,'(I3)') temperature

  CALL read_parameters()

  ALLOCATE(proc1(SUM(noatoms),nompa)); ALLOCATE(proc2(SUM(noatoms),nompa)); ALLOCATE(proc3(SUM(noatoms),nompa))
  ALLOCATE(monomer_s(nomcsteps))
  ALLOCATE(surface_s(nomcsteps))

  CALL initialize_positions()
  CALL initialize_plist()

  CALL DATE_AND_TIME(date,time,VALUES=dates)
  rand_seed=dates(1)+dates(2)+dates(3)+dates(5)+dates(6)+dates(7)+dates(8)

  !rand_seed(1) = time()
  rand_seed_size = 1
  CALL RANDOM_SEED(size = rand_seed_size)
  CALL RANDOM_SEED(put = rand_seed)

  DO current_iteration = 1,noiterations

     time_of_experiment = 0.0d0
     monomer_s = 0.0d0 ; surface_s = 0.0d0

     posns = initposns
     CALL initialize_plist()

     DO current_mc_step = 1,nomcsteps

        p1_prob = SUM(proc1%prob) ; p2_prob = SUM(proc2%prob) ; p3_prob = SUM(proc3%prob)

        CALL RANDOM_NUMBER(rd1)

        fl_1 = rd1*(p1_prob+p2_prob+p3_prob)

        IF (fl_1 .LE. p1_prob) THEN
           sp1 = 1
           !           write(6,*) dummy_temperature, 'CHOSEN PROCESS: 1'
        ELSE IF ((fl_1 .GT. p1_prob) .AND. (fl_1 .LE. p1_prob+p2_prob)) THEN
           sp1 = 2
           !           write(6,*) dummy_temperature, 'CHOSEN PROCESS: 2'
        ELSE
           sp1 = 3
           !           write(6,*) dummy_temperature, 'CHOSEN PROCESS: 3'
        END IF

        CALL make_move(sp1)


        CALL RANDOM_NUMBER(rd3)

        DO WHILE (rd3 .EQ. 0 )
           CALL RANDOM_NUMBER(rd3)
        END DO

        time_of_experiment = time_of_experiment - ((1.0d0/(p1_prob+p2_prob+p3_prob))*LOG(rd3))



        DO rand1 = 1,cs_x
           DO rand2 = 1,cs_y
              IF (MOD(rand1+rand2+2,2).GT.0) THEN ! Choose S sites
                 IF (posns(rand1,rand2,2) .EQ. 1) THEN ! Find monomer vacancies
                    monomer_s(current_mc_step) = monomer_s(current_mc_step) + 1.0d0
                 ELSE ! If the site is not a monomer, i.e. it is a dimer or defect-free
                    surface_s(current_mc_step) = surface_s(current_mc_step) + 1.0d0
                 END IF
              END IF
           END DO
        END DO




        IF (MOD(current_mc_step, ((nomcsteps/number_of_images)-1) ).EQ.0) THEN
           image_id_mc = (current_mc_step/((nomcsteps/number_of_images)-1))
           WRITE(image_id_mc_string,'(I5.5)') image_id_mc
           WRITE(image_id_iter_string,'(I5.5)') current_iteration

           OPEN (unit = 700, file = TRIM(temperature_string)//'K/forplot.txt')


           WRITE(700,*) "# Experiment Time : ", time_of_experiment, " seconds. "


           temp_monomer_count = 0

           DO rand1 = 1,cs_x
              DO rand2 = 1,cs_y
                 IF (MOD(rand1+rand2+2,2).GT.0) THEN ! Choose S sites
                    IF (posns(rand1,rand2,2) .EQ. 1) THEN ! Find monomer vacancies
                       temp_monomer_count = temp_monomer_count + 1
                    END IF
                 END IF
              END DO
           END DO


           WRITE(700,*) "# Number of monomer species : ", temp_monomer_count


           DO rand1 = 1,cs_x
              DO rand2 = 1,cs_y

                 IF (posns(rand1,rand2,2) .NE. 0) THEN

                    WRITE(700,'(I10,I10,F9.3)') rand1, rand2, 2.0d0

                 ELSE

                    neighborhood = 0

                    DO xloop = rand1-1,rand1+1
                       DO yloop = rand2-1,rand2+1

                          IF (xloop<1) THEN
                             xl = cs_x
                          ELSEIF (xloop>cs_x) THEN
                             xl = 1
                          ELSE
                             xl = xloop
                          END IF

                          IF (yloop<1) THEN
                             yl = cs_y
                          ELSEIF (yloop>cs_y) THEN
                             yl = 1
                          ELSE
                             yl = yloop
                          END IF

                          neighborhood = neighborhood + posns(xl,yl,2)



                       END DO
                    END DO

                    WRITE(700,'(I10,I10,F9.3)') rand1, rand2, 0.0d0+(neighborhood/9.0d0)

                 END IF


              END DO
           END DO


           CLOSE(700)



           CALL system('cp '//TRIM(temperature_string)//'K/forplot.txt '//TRIM(temperature_string)//'K/fp-'//TRIM(image_id_iter_string)//'-'//TRIM(image_id_mc_string)//'.txt')

        END IF


     END DO
  END DO



END SUBROUTINE kmc



















SUBROUTINE make_move(sp1)
  USE commondata
  IMPLICIT NONE
  INTEGER :: sp1, rn1, rn2, lp
  INTEGER :: moved_to_x, moved_to_y, moved_to_z
  INTEGER :: moved_from_x, moved_from_y, moved_from_z
  INTEGER :: botval ! PLEASE change variable name
  INTEGER :: lattice_id, mpa
  REAL*8 :: rd1,cutoff, tripper

  INTEGER :: rand_seed_size
  INTEGER, DIMENSION(1) :: rand_seed

  CALL RANDOM_NUMBER(rd1)

  IF (sp1 == 1) THEN

     tripper = 0.0d0
     cutoff = rd1*SUM(proc1%prob)

     DO lattice_id = 1,SUM(noatoms)
        DO mpa = 1,nompa

           tripper = tripper + proc1(lattice_id, mpa)%prob

           IF (tripper .GE. cutoff) THEN

              rn1 = lattice_id
              rn2 = mpa

              moved_from_x = proc1(rn1,rn2)%fx ; moved_from_y = proc1(rn1,rn2)%fy ; moved_from_z = proc1(rn1,rn2)%fz
              moved_to_x = proc1(rn1,rn2)%tx ; moved_to_y = proc1(rn1,rn2)%ty ; moved_to_z = proc1(rn1,rn2)%tz

              IF (MOD(moved_from_x+moved_from_y+moved_from_z,2).GT.0) THEN
                 posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)=posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)-1
              ELSE
                 posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)=posns(proc1(rn1,rn2)%fx,proc1(rn1,rn2)%fy,proc1(rn1,rn2)%fz)-2
              END IF

              CALL rfpls(moved_from_x,moved_from_y,moved_from_z) ; CALL add_to_plist(moved_from_x,moved_from_y,moved_from_z)
              CALL rfpls(moved_to_x,moved_to_y,moved_to_z) ; CALL add_to_plist(moved_to_x,moved_to_y,moved_to_z)
              CALL upnl(moved_from_x,moved_from_y,moved_from_z) ; CALL upnl(moved_to_x,moved_to_y,moved_to_z)

              EXIT
           ELSE

           END IF
        END DO
        IF (tripper .GE. cutoff) EXIT
     END DO


  ELSE IF (sp1 == 2) THEN

     tripper = 0.0d0
     cutoff = rd1*SUM(proc2%prob)

     DO lattice_id = 1,SUM(noatoms)
        DO mpa = 1,nompa

           tripper = tripper + proc2(lattice_id, mpa)%prob

           IF (tripper .GE. cutoff) THEN

              rn1 = lattice_id
              rn2 = mpa

              moved_from_x = proc2(rn1,rn2)%fx ; moved_from_y = proc2(rn1,rn2)%fy ; moved_from_z = proc2(rn1,rn2)%fz
              moved_to_x = proc2(rn1,rn2)%tx ; moved_to_y = proc2(rn1,rn2)%ty ; moved_to_z = proc2(rn1,rn2)%tz

              IF (MOD(moved_from_x+moved_from_y+moved_from_z,2).GT.0) THEN
                 posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)=posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)+1
              ELSE
                 posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)=posns(proc2(rn1,rn2)%fx,proc2(rn1,rn2)%fy,proc2(rn1,rn2)%fz)+2
              END IF

              CALL rfpls(moved_from_x,moved_from_y,moved_from_z) ; CALL add_to_plist(moved_from_x,moved_from_y,moved_from_z)
              CALL rfpls(moved_to_x,moved_to_y,moved_to_z) ; CALL add_to_plist(moved_to_x,moved_to_y,moved_to_z)
              CALL upnl(moved_from_x,moved_from_y,moved_from_z) ; CALL upnl(moved_to_x,moved_to_y,moved_to_z)

              EXIT
           ELSE

           END IF
        END DO
        IF (tripper .GE. cutoff) EXIT
     END DO





  ELSE IF (sp1 == 3) THEN

     tripper = 0.0d0
     cutoff = rd1*SUM(proc3%prob)

     DO lattice_id = 1,SUM(noatoms)
        DO mpa = 1,nompa

           tripper = tripper + proc3(lattice_id, mpa)%prob

           IF (tripper .GE. cutoff) THEN

              rn1 = lattice_id
              rn2 = mpa

              moved_from_x = proc3(rn1,rn2)%fx ; moved_from_y = proc3(rn1,rn2)%fy ; moved_from_z = proc3(rn1,rn2)%fz
              moved_to_x = proc3(rn1,rn2)%tx ; moved_to_y = proc3(rn1,rn2)%ty ; moved_to_z = proc3(rn1,rn2)%tz

              IF (MOD(moved_from_x+moved_from_y+moved_from_z,2).GT.0) THEN
                 posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)=posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)+1
                 posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)=posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)-1
              ELSE
                 posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)=posns(proc3(rn1,rn2)%fx,proc3(rn1,rn2)%fy,proc3(rn1,rn2)%fz)+2
                 posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)=posns(proc3(rn1,rn2)%tx,proc3(rn1,rn2)%ty,proc3(rn1,rn2)%tz)-2
              END IF

              CALL rfpls(moved_from_x,moved_from_y,moved_from_z) ; CALL add_to_plist(moved_from_x,moved_from_y,moved_from_z)
              CALL rfpls(moved_to_x,moved_to_y,moved_to_z) ; CALL add_to_plist(moved_to_x,moved_to_y,moved_to_z)
              CALL upnl(moved_from_x,moved_from_y,moved_from_z) ; CALL upnl(moved_to_x,moved_to_y,moved_to_z)

              EXIT
           ELSE


           END IF
        END DO
        IF (tripper .GE. cutoff) EXIT
     END DO

  END IF


END SUBROUTINE make_move

SUBROUTINE initialize_positions()
  USE commondata
  IMPLICIT NONE

  INTEGER :: rand_seed_size
  INTEGER, DIMENSION(1) :: rand_seed
  INTEGER :: atom_type, numb_of_atoms, random_integer
  INTEGER :: x_pos, y_pos, z_pos, x_loop, y_loop, z_loop
  REAL*8 :: random_double

  ! Distribute atoms
  ALLOCATE(initposns(cs_x,cs_y,cs_z))
  ALLOCATE(posns(cs_x,cs_y,cs_z))

  DO x_loop = 1, cs_x
     DO y_loop = 1, cs_y
        posns(x_loop,y_loop, cs_z) = 2
        DO z_loop = 1, cs_z-1
           posns(x_loop,y_loop,z_loop) = 2
        END DO
     END DO
  END DO


  initposns = posns




END SUBROUTINE initialize_positions
SUBROUTINE initialize_plist()
  USE commondata
  IMPLICIT NONE

  INTEGER :: atom_id, pid
  INTEGER :: x_pos, y_pos, z_pos

  DO atom_id = 1, SUM(noatoms)
     DO pid = 1, nompa

        proc1(atom_id,pid)%fx = 0 ; proc1(atom_id,pid)%fy = 0 ; proc1(atom_id,pid)%fz = 0
        proc1(atom_id,pid)%tx = 0 ; proc1(atom_id,pid)%ty = 0 ; proc1(atom_id,pid)%tz = 0
        proc1(atom_id,pid)%prob = 0.0d0

        proc2(atom_id,pid)%fx = 0 ; proc2(atom_id,pid)%fy = 0 ; proc2(atom_id,pid)%fz = 0
        proc2(atom_id,pid)%tx = 0 ; proc2(atom_id,pid)%ty = 0 ; proc2(atom_id,pid)%tz = 0
        proc2(atom_id,pid)%prob = 0.0d0

        proc3(atom_id,pid)%fx = 0 ; proc3(atom_id,pid)%fy = 0 ; proc3(atom_id,pid)%fz = 0
        proc3(atom_id,pid)%tx = 0 ; proc3(atom_id,pid)%ty = 0 ; proc3(atom_id,pid)%tz = 0
        proc3(atom_id,pid)%prob = 0.0d0

     END DO
  END DO


  atom_id = 0

  DO z_pos = 1, cs_z
     DO y_pos = 1, cs_y
        DO x_pos = 1, cs_x
           CALL add_to_plist(x_pos,y_pos,z_pos)
        END DO
     END DO
  END DO

END SUBROUTINE initialize_plist

SUBROUTINE read_parameters()
  USE commondata
  IMPLICIT NONE

  WRITE(6,*) TRIM(temperature_string)

  ! Read lattice parameters
  CALL system("cat "//TRIM(temperature_string)//"K/param.in | grep ^LATTICE | sed 's/LATTICE//g'| sed 's/=//g' > "//TRIM(temperature_string)//"K/.toget-lattice")
  OPEN(unit = 7001, file = ""//TRIM(temperature_string)//"K/.toget-lattice", status = 'old') ; READ(7001,*) cs_x, cs_y, cs_z ; CLOSE(7001)
  WRITE(6,*) 'Hi'
  CALL system ("rm "//TRIM(temperature_string)//"K/.toget-lattice")

  ! Read number of species and number of atoms in each species
  CALL system("cat "//TRIM(temperature_string)//"K/param.in | grep ^SPECIES | sed 's/SPECIES//g'| sed 's/=//g' | awk '{print NF}' > "//TRIM(temperature_string)//"K/.toget-numbofspec")
  OPEN(unit = 7002, file = ""//TRIM(temperature_string)//"K/.toget-numbofspec", status = 'old') ; READ(7002,*) nospec ; CLOSE(7002)
  CALL system ("rm "//TRIM(temperature_string)//"K/.toget-numbofspec")
  ALLOCATE(noatoms(nospec))
  CALL system("cat "//TRIM(temperature_string)//"K/param.in | grep ^SPECIES | sed 's/SPECIES//g'| sed 's/=//g' > "//TRIM(temperature_string)//"K/.toget-numbofatoms")
  OPEN(unit = 7003, file = TRIM(temperature_string)//'K/.toget-numbofatoms', status = 'old') ; READ(7003,*) noatoms ; CLOSE(7003)
  CALL system ("rm "//TRIM(temperature_string)//"K/.toget-numbofatoms")

  ! Read number of steps
  CALL system("cat "//TRIM(temperature_string)//"K/param.in | grep ^MCSTEPS | sed 's/MCSTEPS//g'| sed 's/=//g' > "//TRIM(temperature_string)//"K/.toget-mcsteps")
  OPEN(unit = 7004, file = TRIM(temperature_string)//'K/.toget-mcsteps', status = 'old') ; READ(7004,*) nomcsteps ; CLOSE(7004)
  CALL system ("rm "//TRIM(temperature_string)//"K/.toget-mcsteps")

  ! Read number of iterations
  CALL system("cat "//TRIM(temperature_string)//"K/param.in | grep ^ITERATIONS | sed 's/ITERATIONS//g'| sed 's/=//g' > "//TRIM(temperature_string)//"K/.toget-iterations")
  OPEN(unit = 7005, file = TRIM(temperature_string)//'K/.toget-iterations', status = 'old') ; READ(7005,*) noiterations ; CLOSE(7005)
  CALL system ("rm "//TRIM(temperature_string)//"K/.toget-iterations")

END SUBROUTINE read_parameters

SUBROUTINE add_to_plist(lp,lq,lr)
  USE commondata
  IMPLICIT NONE
  INTEGER :: lp,lq,lr, atom_id, pt
  INTEGER :: loop1, loop2 ! PLEASE change variable names
  INTEGER :: xloop, yloop ! PLEASE change variable names
  INTEGER :: xp, yp, zp ! PLEASE change variable names

  ! Create vacancies
  IF (posns(lp,lq,lr) .NE. 0) THEN
     pt = 1 ; CALL make_ptl(pt,lp,lq,lr,lp,lq,lr)
  END IF

  ! Fill vacancies
  IF (posns(lp,lq,lr) .NE. 2) THEN
     pt = 2 ; CALL make_ptl(pt,lp,lq,lr,lp,lq,lr)
  END IF

  ! Diffusion
  IF (posns(lp,lq,lr) .NE. 2) THEN

     DO xloop = lp-1,lp+1,2
        DO yloop = lq-1,lq+1,2

           IF (xloop<1)  THEN
              loop1 = cs_x
           ELSEIF (xloop>cs_x)THEN
              loop1 = 1
           ELSE
              loop1 = xloop
           END IF

           IF (yloop<1)  THEN
              loop2 = cs_x
           ELSEIF (yloop>cs_x)THEN
              loop2 = 1
           ELSE
              loop2 = yloop
           END IF

           IF (posns(loop1,loop2,lr) .NE. 0) THEN
              pt = 3 ; CALL make_ptl(pt,lp,lq,lr,loop1,loop2,lr)
           END IF

        END DO
     END DO

  END IF

END SUBROUTINE add_to_plist







SUBROUTINE make_ptl(pt, fromx, fromy, fromz, tox, toy, toz)
  USE commondata
  INTEGER, INTENT(in) :: pt, fromx, fromy, fromz, tox, toy, toz
  INTEGER :: atom_id
  REAL*8 :: kbT
  REAL*8 :: process_barrier
  REAL*8 :: s_surf, s_bulk, fe_surf, fe_bulk
  REAL*8 :: neighborhood,neighborhoodfrom,neighborhoodto
  INTEGER :: xloop, yloop

  atom_id = ((fromz-1)*(cs_x*cs_y)) + ((fromy-1)*cs_x) + fromx

  kbT = 0.000086173423d0 * temperature


  IF (pt == 1) THEN

     s_surf = 0.075d0 ; fe_surf = 0.075d0
     s_bulk = 0.075d0 ; fe_bulk = 0.075d0

     IF (MOD(fromx+fromy+fromz,2).GT.0) THEN
        process_barrier = 0.0d0 - (s_surf + ((s_bulk-s_surf)*(cs_z-fromz)))
     ELSE
        process_barrier = 0.0d0 - (fe_surf + ((fe_bulk-fe_surf)*(cs_z-fromz)))
     END IF

     ! Neighborhood sampling
     neighborhood = 0.0d0
     DO xloop = tox-1,tox+1,2
        IF (xloop<1)  THEN
           neighborhood = neighborhood + posns(cs_x,toy,toz)
        ELSEIF (xloop>cs_x)THEN
           neighborhood = neighborhood + posns(1,toy,toz)
        ELSE
           neighborhood = neighborhood + posns(xloop,toy,toz)
        END IF
     END DO

     DO yloop = toy-1,toy+1,2
        IF (yloop<1)  THEN
           neighborhood = neighborhood + posns(tox,cs_y,toz)
        ELSEIF (yloop>cs_y)THEN
           neighborhood = neighborhood + posns(tox,1,toz)
        ELSE
           neighborhood = neighborhood + posns(tox,yloop,toz)
        END IF
     END DO
     ! Done neighborhood sampling

     process_barrier = process_barrier * ((neighborhood/8.0d0))

     IF (MOD(fromx+fromy+fromz,2).GT.0) THEN
        process_barrier = process_barrier - 1.47d0
     ELSE
        process_barrier = process_barrier - 1.47d0 - 0.08d0
     END IF


     DO lp = 1, nompa
        IF (proc1(atom_id, lp)%fx == 0) THEN
           proc1(atom_id,lp)%fx = fromx ; proc1(atom_id,lp)%fy = fromy ; proc1(atom_id,lp)%fz = fromz
           proc1(atom_id,lp)%tx = tox ; proc1(atom_id,lp)%ty = toy ; proc1(atom_id,lp)%tz = toz
           proc1(atom_id,lp)%prob = (10**6)*dexp(process_barrier/kbT)*(10**6)
           EXIT
        END IF
     END DO










  ELSE IF (pt == 2) THEN


     s_surf = 0.075d0 ; fe_surf = 0.075d0
     s_bulk = 0.075d0 ; fe_bulk = 0.075d0

     IF (MOD(fromx+fromy+fromz,2).GT.0) THEN
        process_barrier = 0.0d0 + (s_surf + ((s_bulk-s_surf)*(cs_z-fromz)))
     ELSE
        process_barrier = 0.0d0 + (fe_surf + ((fe_bulk-fe_surf)*(cs_z-fromz)))
     END IF

     ! Neighborhood sampling
     neighborhood = 0.0d0
     DO xloop = tox-1,tox+1,2
        IF (xloop<1)  THEN
           neighborhood = neighborhood + posns(cs_x,toy,toz)
        ELSEIF (xloop>cs_x)THEN
           neighborhood = neighborhood + posns(1,toy,toz)
        ELSE
           neighborhood = neighborhood + posns(xloop,toy,toz)
        END IF
     END DO

     DO yloop = toy-1,toy+1,2
        IF (yloop<1)  THEN
           neighborhood = neighborhood + posns(tox,cs_y,toz)
        ELSEIF (yloop>cs_y)THEN
           neighborhood = neighborhood + posns(tox,1,toz)
        ELSE
           neighborhood = neighborhood + posns(tox,yloop,toz)
        END IF
     END DO
     ! Done neighborhood sampling

     process_barrier = process_barrier * ((neighborhood/8.0d0))

     IF (MOD(fromx+fromy+fromz,2).GT.0) THEN
        process_barrier = process_barrier - 1.600
     ELSE
        process_barrier = process_barrier - 1.60d0 + 0.08d0
     END IF


     DO lp = 1, nompa
        IF (proc2(atom_id, lp)%fx == 0) THEN
           proc2(atom_id,lp)%fx = fromx ; proc2(atom_id,lp)%fy = fromy ; proc2(atom_id,lp)%fz = fromz
           proc2(atom_id,lp)%tx = tox ; proc2(atom_id,lp)%ty = toy ; proc2(atom_id,lp)%tz = toz
           proc2(atom_id,lp)%prob = (10**6)*dexp(process_barrier/kbT)*(10**6)
           EXIT
        END IF
     END DO






  ELSE IF (pt == 3) THEN


     ! Neighborhood sampling
     neighborhood = 0.0d0
     DO xloop = tox-1,tox+1,2
        IF (xloop<1)  THEN
           neighborhood = neighborhood + posns(cs_x,toy,toz)
        ELSEIF (xloop>cs_x)THEN
           neighborhood = neighborhood + posns(1,toy,toz)
        ELSE
           neighborhood = neighborhood + posns(xloop,toy,toz)
        END IF
     END DO

     DO yloop = toy-1,toy+1,2
        IF (yloop<1)  THEN
           neighborhood = neighborhood + posns(tox,cs_y,toz)
        ELSEIF (yloop>cs_y)THEN
           neighborhood = neighborhood + posns(tox,1,toz)
        ELSE
           neighborhood = neighborhood + posns(tox,yloop,toz)
        END IF
     END DO

     neighborhoodto = neighborhood


     neighborhood = 0.0d0
     DO xloop = fromx-1,fromx+1,2
        IF (xloop<1)  THEN
           neighborhood = neighborhood + posns(cs_x,fromy,fromz)
        ELSEIF (xloop>cs_x)THEN
           neighborhood = neighborhood + posns(1,fromy,fromz)
        ELSE
           neighborhood = neighborhood + posns(xloop,fromy,fromz)
        END IF
     END DO

     DO yloop = fromy-1,fromy+1,2
        IF (yloop<1)  THEN
           neighborhood = neighborhood + posns(fromx,cs_y,fromz)
        ELSEIF (yloop>cs_y)THEN
           neighborhood = neighborhood + posns(fromx,1,fromz)
        ELSE
           neighborhood = neighborhood + posns(fromx,yloop,fromz)
        END IF
     END DO

     neighborhoodfrom = neighborhood

     ! Done neighborhood sampling

     s_surf = 0.00d0 ; fe_surf = 0.00d0
     s_bulk = 0.00d0 ; fe_bulk = 0.00d0

     IF (MOD(fromx+fromy+fromz,2).GT.0) THEN
        process_barrier = 0.0d0 - (s_surf + ((s_bulk-s_surf)*(cs_z-fromz)))
     ELSE
        process_barrier = 0.0d0 - (fe_surf + ((fe_bulk-fe_surf)*(cs_z-fromz)))
     END IF

     process_barrier = process_barrier + (0.1125d0 * ((neighborhoodfrom/8.0d0) - (neighborhoodto/8.0d0) )) - 1.35d0

     DO lp = 1, nompa
        IF (proc3(atom_id, lp)%fx == 0) THEN
           proc3(atom_id,lp)%fx = fromx ; proc3(atom_id,lp)%fy = fromy ; proc3(atom_id,lp)%fz = fromz
           proc3(atom_id,lp)%tx = tox ; proc3(atom_id,lp)%ty = toy ; proc3(atom_id,lp)%tz = toz
           proc3(atom_id,lp)%prob = (10**6)*dexp(process_barrier/kbT)*(10**6)
           EXIT
        END IF
     END DO

  END IF


END SUBROUTINE make_ptl








SUBROUTINE rfpls(x_pos,y_pos,z_pos)
  USE commondata
  IMPLICIT NONE
  INTEGER :: atom_id
  INTEGER :: x_pos, y_pos, z_pos
  INTEGER :: lq

  atom_id = ((z_pos-1)*(cs_x*cs_y)) + ((y_pos-1)*(cs_x)) + x_pos

  DO lq = 1, nompa
     proc1(atom_id,lq)%fx = 0 ; proc1(atom_id,lq)%fy = 0 ; proc1(atom_id,lq)%fz = 0
     proc1(atom_id,lq)%tx = 0 ; proc1(atom_id,lq)%ty = 0 ; proc1(atom_id,lq)%tz = 0 ; proc1(atom_id,lq)%prob = 0.0d0
     proc2(atom_id,lq)%fx = 0 ; proc2(atom_id,lq)%fy = 0 ; proc2(atom_id,lq)%fz = 0
     proc2(atom_id,lq)%tx = 0 ; proc2(atom_id,lq)%ty = 0 ; proc2(atom_id,lq)%tz = 0 ; proc2(atom_id,lq)%prob = 0.0d0
     proc3(atom_id,lq)%fx = 0 ; proc3(atom_id,lq)%fy = 0 ; proc3(atom_id,lq)%fz = 0
     proc3(atom_id,lq)%tx = 0 ; proc3(atom_id,lq)%ty = 0 ; proc3(atom_id,lq)%tz = 0 ; proc3(atom_id,lq)%prob = 0.0d0
  END DO

END SUBROUTINE rfpls




SUBROUTINE upnl(ma_x,ma_y,ma_z)
  USE commondata
  IMPLICIT NONE
  INTEGER :: ma_x,ma_y,ma_z, lp,lq,lr, newlp,newlq,newlr, atom_id, getatomid

  DO lp = ma_x-1, ma_x+1
     DO lq = ma_y-1, ma_y+1
        DO lr = ma_z-1, ma_z+1

           IF (lp .LT. 1) THEN
              newlp = cs_x
           ELSE IF (lp .GT. cs_x) THEN
              newlp = 1
           ELSE
              newlp = lp
           END IF

           IF (lq .LT. 1) THEN
              newlq = cs_y
           ELSE IF (lq .GT. cs_y) THEN
              newlq = 1
           ELSE
              newlq = lq
           END IF

           IF (lr .LT. 1) THEN
              newlr = cs_z
           ELSE IF (lr .GT. cs_z) THEN
              newlr = 1
           ELSE
              newlr = lr
           END IF


           IF (posns(newlp, newlq, newlr) .GT. -1) THEN
              IF ((newlp .NE. ma_x) .OR. (newlq .NE. ma_y) .OR. (newlr .NE. ma_z))  THEN

                 atom_id = ((newlr-1)*(cs_x*cs_y)) + ((newlq-1)*cs_x) + newlp

                 CALL rfpls(newlp,newlq,newlr)
                 CALL add_to_plist(newlp,newlq,newlr)

              END IF
           END IF


        END DO
     END DO
  END DO


END SUBROUTINE upnl
