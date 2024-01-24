PROGRAM HE4_EXCITATION_SPECTRUM


!--- DECLARATIONS OF VARIABLES ----------------------------- 

    IMPLICIT NONE

    ! SYSTEM PARAMETERS kinda particle number, dimension of space, and so forth
    INTEGER, PARAMETER :: Number_of_particles = 25 , spacial_dimension = 2  
    
    !!!! PARTICLE_NUMBER_DENSITY = 0.0625 (A)^(-2) = 25 particles / 400 (A)^(2)

    ! VARIATIONAL PARAMETERS (should be determined through the steepest descent method)
    !! ground state variational parameters
    DOUBLE PRECISION :: c , b , delta , alpha
    !! excited states variational parameters
    DOUBLE PRECISION :: r_0 , w  
    
    ! Positions of Helium atoms
    DOUBLE PRECISION, DIMENSION( spacial_dimension , number_of_particles ) :: RR_present , RR_proposed
    
    ! Random deviations on positions of Helium atoms
    DOUBLE PRECISION :: delta_RR_x , delta_RR_y ! size of box (-0.5 , 0.5) , in Angstrom scale , should be adjusted reffering the acceptance rate
    
    ! Shadow variables
    DOUBLE PRECISION, DIMENSION( spacial_dimension , number_of_particles ) :: SS1_present , SS1_proposed
    DOUBLE PRECISION, DIMENSION( spacial_dimension , number_of_particles ) :: SS2_present , SS2_proposed
    
    ! Random deviations on Shadow variables
    DOUBLE PRECISION :: delta_SS1_x , delta_SS1_y ! size of box (-0.5 , 0.5) , in Angstrom scale , should be adjusted reffering the acceptance rate
    DOUBLE PRECISION :: delta_SS2_x , delta_SS2_y

    ! Probability ratio in the Metropolis algorithm
    DOUBLE PRECISION :: probability_ratio_RR , probability_ratio_SS1 , probability_ratio_SS2
    
    ! Random numbers which generate the initial configurations
    DOUBLE PRECISION :: RND_position , RND_shadows1 , RND_shadows2

    ! Random number for the Metropolis algorithm
    DOUBLE PRECISION :: RND_Metropolis_RR , RND_Metropolis_SS1 , RND_Metropolis_SS2

    ! DUMMY VARIABLES
    INTEGER :: i , j , iteration , axis , particle
    
    ! Monte-Carlo steps
    INTEGER, PARAMETER :: MC_steps = 100000

    ! Acceptance rate
    INTEGER :: N_accepted_RR , N_accepted_SS1 , N_accepted_SS2
    REAL    :: acceptance_rate_RR , acceptance_rate_SS1 , acceptance_rate_SS2

    ! Execution time eastimation
    REAL :: start , end , execution_time

    CALL CPU_TIME(start)




!--- DETERMINATION OF VARIATIONAL PARAMETERS USING THE GROUND STATE SHADOW WAVE FUNCTION


    ! INITIALISATION OF THE LOCATIONS OF HELIUM ATOMS AND THE SHADOW VARIABLES
    RR_present = 0.
    RR_proposed = 0.
    
    SS1_present = 0.
    SS1_proposed = 0.
    
    SS2_present = 0.
    SS2_proposed = 0.


    ! RANDOM INITIAL POSITIONS AND SHADOWS CONFIGURATIONS
    DO particle = 1 , Number_of_particles , +1
        
        DO axis = 1 , spacial_dimension , +1
 
            CALL RANDOM_NUMBER(RND_position)
            CALL RANDOM_NUMBER(RND_shadows1)
            CALL RANDOM_NUMBER(RND_shadows2)

            RR_present(axis , particle)  =  20. * RND_position - 10.    !(-10,10) in Angstrom scale
            SS1_present(axis , particle) =  20. * RND_shadows1 - 10.    !(-10,10) in Angstrom scale
            SS2_present(axis , particle) =  20. * RND_shadows2 - 10.    !(-10,10) in Angstrom scale

        END DO
    
    END DO


    ! INITIALISE VARIATIONAL PARAMETERS
    c = 4.0
    b = 1.0
    delta = 0.05
    alpha = 0.90


    N_accepted_RR = 0
    N_accepted_SS1 = 0
    N_accepted_SS2 = 0


    !=== MAIN LOOP ============================================================
    DO iteration = 1 , MC_steps , +1


        !--- POSITIONAL SPACE METROPOLIS ALGORITHM ----------------------------

        ! Choose a new coordinate with uniform probability in a cube of side delta_RR (-0.5, 0.5) centered at the present coordinate.
        ! The size of box should make the acceptance rate is became less than 0.5.
        particle = 0
        DO particle = 1 , Number_of_particles , +1

            RR_proposed = RR_present
        
            ! (0,1)
            CALL RANDOM_NUMBER(delta_RR_x)
            CALL RANDOM_NUMBER(delta_RR_y)

            ! (-0.5, +0.5)
            delta_RR_x = -0.5 + delta_RR_x
            delta_RR_y = -0.5 + delta_RR_y

            RR_proposed(1, particle) = RR_present(1, particle) + delta_RR_x
            RR_proposed(2, particle) = RR_present(2, particle) + delta_RR_y


            ! Calculation of the probability ratio
            probability_ratio_RR = 1.0

            i = 0
            DO i = 1 , Number_of_particles , +1

                IF ( i == particle ) THEN

                    probability_ratio_RR = probability_ratio_RR &
                        * EXP( -c * (   (separation( RR_proposed( : , particle ) , SS1_present( : , particle ) ))**(2) &
                                      + (separation( RR_proposed( : , particle ) , SS2_present( : , particle ) ))**(2) &
                                      - (separation( RR_present(  : , particle ) , SS1_present( : , particle ) ))**(2) &
                                      - (separation( RR_present(  : , particle ) , SS2_present( : , particle ) ))**(2)   ) )

                ELSE !( i /= particle ) THEN

                    probability_ratio_RR = probability_ratio_RR &
                        * EXP( - (   pseudo_u( separation( RR_proposed( : , particle ) , RR_present( : , i ) ) , b ) &
                                   - pseudo_u( separation( RR_present(  : , particle ) , RR_present( : , i ) ) , b )   ) )

                END IF

            END DO


            ! SHOULD WE NEED TO ACCEPT THE MOVE?

            CALL RANDOM_NUMBER(RND_Metropolis_RR)

            IF ( RND_Metropolis_RR < MIN( 1.0 , probability_ratio_RR ) ) THEN

                ! Accept the move 
                RR_present(1, particle) = RR_proposed(1, particle)
                RR_present(2, particle) = RR_proposed(2, particle)
                
                N_accepted_RR = N_accepted_RR + 1

            ELSE 
                ! Reject the move
                N_accepted_RR = N_accepted_RR + 0
            
            END IF

            !print *, N_accepted_RR
            !print *, probability_ratio_RR


        END DO



        !--- SHADOWS_1 METROPOLIS ALGORITHM ----------------------------

        ! Choose a new coordinate with uniform probability in a cube of side delta_SS1 (-0.5, 0.5) centered at the present coordinate.
        ! The size of box should make the acceptance rate is became less than 0.5.
        particle = 0
        DO particle = 1 , Number_of_particles , +1

            SS1_proposed = SS1_present
        
            ! (0,1)
            CALL RANDOM_NUMBER(delta_SS1_x)
            CALL RANDOM_NUMBER(delta_SS1_y)

            ! (-0.5, +0.5)
            delta_SS1_x = -0.5 + delta_SS1_x
            delta_SS1_y = -0.5 + delta_SS1_y

            SS1_proposed(1, particle) = SS1_present(1, particle) + delta_SS1_x
            SS1_proposed(2, particle) = SS1_present(2, particle) + delta_SS1_y


            ! Calculation of the probability ratio
            probability_ratio_SS1 = 1.0

            i = 0
            DO i = 1 , Number_of_particles , +1

                IF ( i == particle ) THEN

                    probability_ratio_SS1 = probability_ratio_SS1 &
                    * EXP( -c * (   (separation( RR_present( : , particle ) , SS1_proposed( : , particle ) ))**(2) &
                                  - (separation( RR_present( : , particle ) , SS1_present(  : , particle ) ))**(2)   ) )

                ELSE !( i /= particle ) THEN

                    probability_ratio_SS1 = probability_ratio_SS1 &
                * EXP( - (   pseudo_w( separation( SS1_proposed( : , particle ) , SS1_present( : , i ) ) , delta , alpha ) &
                           - pseudo_w( separation( SS1_present(  : , particle ) , SS1_present( : , i ) ) , delta , alpha )   ) )

                END IF

            END DO


            ! SHOULD WE NEED TO ACCEPT THE MOVE?

            CALL RANDOM_NUMBER(RND_Metropolis_SS1)

            IF ( RND_Metropolis_SS1 < MIN( 1. , probability_ratio_SS1 ) ) THEN

                ! Accept the move 
                
                SS1_present(1, particle) = SS1_proposed(1, particle)
                SS1_present(2, particle) = SS1_proposed(2, particle)
                
                N_accepted_SS1 = N_accepted_SS1 + 1

            ELSE 
                ! Reject the move
                N_accepted_SS1 = N_accepted_SS1 + 0
            
            END IF

            !print *, N_accepted_SS1
            !print *, probability_ratio_SS1


        END DO



        !--- SHADOWS_2 METROPOLIS ALGORITHM ----------------------------

        ! Choose a new coordinate with uniform probability in a cube of side delta_SS2 (-0.5, 0.5) centered at the present coordinate.
        ! The size of box should make the acceptance rate is became less than 0.5.
        particle = 0
        DO particle = 1 , Number_of_particles , +1

            SS2_proposed = SS2_present
        
            ! (0,1)
            CALL RANDOM_NUMBER(delta_SS2_x)
            CALL RANDOM_NUMBER(delta_SS2_y)

            ! (-0.1, +0.1)
            delta_SS2_x = -0.5 + delta_SS2_x
            delta_SS2_y = -0.5 + delta_SS2_y

            SS2_proposed(1, particle) = SS2_present(1, particle) + delta_SS2_x
            SS2_proposed(2, particle) = SS2_present(2, particle) + delta_SS2_y


            ! Calculation of the probability ratio
            probability_ratio_SS2 = 1.0

            i = 0
            DO i = 1 , Number_of_particles , +1

                IF ( i == particle ) THEN

                    probability_ratio_SS2 = probability_ratio_SS2 &
                        * EXP( -c * (   (separation( RR_present( : , particle ) , SS2_proposed( : , particle ) ))**(2) &
                                      - (separation( RR_present( : , particle ) , SS2_present(  : , particle ) ))**(2)   ) )

                ELSE !( i /= particle ) THEN

                    probability_ratio_SS2 = probability_ratio_SS2 &
                * EXP( - (   pseudo_w( separation( SS2_proposed( : , particle ) , SS2_present( : , i ) ) , delta , alpha ) &
                           - pseudo_w( separation( SS2_present(  : , particle ) , SS2_present( : , i ) ) , delta , alpha )   ) )

                END IF

            END DO


            ! SHOULD WE NEED TO ACCEPT THE MOVE?

            CALL RANDOM_NUMBER(RND_Metropolis_SS2)

            IF ( RND_Metropolis_SS2 < MIN( 1.0 , probability_ratio_SS2 ) ) THEN

                ! Accept the move 
                
                SS2_present(1, particle) = SS2_proposed(1, particle)
                SS2_present(2, particle) = SS2_proposed(2, particle)
                
                N_accepted_SS2 = N_accepted_SS2 + 1

            ELSE 
                ! Reject the move
                N_accepted_SS2 = N_accepted_SS2 + 0
            
            END IF



        END DO

        !print *, N_accepted_SS2
        !print *, probability_ratio_SS2

        !--------------------------------------------------------------------

    
    END DO



    acceptance_rate_RR  = 1.*N_accepted_RR / ( MC_steps * Number_of_particles)
    acceptance_rate_SS1 = 1.*N_accepted_SS1 / ( MC_steps * Number_of_particles)
    acceptance_rate_SS2 = 1.*N_accepted_SS2 / ( MC_steps * Number_of_particles)


    !PRINT *, acceptance_rate_RR
    PRINT *, "RR  = " , acceptance_rate_RR
    PRINT *, "SS1 = " , acceptance_rate_SS1
    PRINT *, "SS2 = " , acceptance_rate_SS2


    





    ! MONTE-CARLO DO LOOP 바깥쪽에 VARIATIONAL PARAMETER를 조정해서 ENERGY를 MINIMISE하는 
    ! STEEPEST DESCENT METHOD LOOP가 들어와야 할듯.






    CALL CPU_TIME(end)

    execution_time = end - start

    PRINT *, "EXECUTION TIME = ", execution_time , " SEC"





!--- DECLARATIONS OF FUNCTIONS
    
    CONTAINS


    !--- Calculate L2 Norm distance btw two vectors --------
    FUNCTION separation(array1 , array2)

        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(2) :: array1 , array2
        DOUBLE PRECISION               :: separation

        separation = NORM2( array1 - array2 ) !* 10**(-10)  ! Rescale into meter scale from Angstrom scale.

    END FUNCTION
    !-------------------------------------------------------

    
    
    !--- McMillan form -------------------------------------
    FUNCTION pseudo_u(distance , variational) ! input distance is should be in meter scale.

        IMPLICIT NONE 

        DOUBLE PRECISION :: pseudo_u , distance , variational

        pseudo_u = ( variational / distance ) ** ( 5 )

    END FUNCTION pseudo_u
    !-------------------------------------------------------



    !--- Aziz - HFDHE2 potential ---------------------------
    FUNCTION pseudo_w(distance , variational1 , variational2) ! variational1 = delta , variational2 = alpha
                                                    ! input distance is should be in meter scale.
        IMPLICIT NONE
        
        DOUBLE PRECISION            :: pseudo_w , distance , variational1 , variational2
        DOUBLE PRECISION            :: x        
        DOUBLE PRECISION, PARAMETER :: kb = 1.380649D-23
        DOUBLE PRECISION, PARAMETER :: epsilon = kb * 10.8
        DOUBLE PRECISION, PARAMETER :: A = 0.54485046
        DOUBLE PRECISION, PARAMETER :: alphaa = 13.353384
        DOUBLE PRECISION, PARAMETER :: C6 = 1.3732412
        DOUBLE PRECISION, PARAMETER :: C8 = 0.4253785
        DOUBLE PRECISION, PARAMETER :: C10 = 0.178100
        DOUBLE PRECISION, PARAMETER :: D = 1.241314
        DOUBLE PRECISION, PARAMETER :: r_m = 2.9673 !* 10**(-10) ! in meter scale

        x = distance / r_m

        IF ( variational2 * x < D ) THEN
            
            pseudo_w = variational1 * epsilon * ( EXP(- alphaa * (variational2 * x) ) &
                        - ( C6 / (variational2 * x)**6 + C8 / (variational2 * x)**8 + C10 / (variational2 * x)**10 ) &
                        * EXP( -( D/(variational2 * x) - 1 )**2 ) )  

        ELSE 

            pseudo_w = variational1 * epsilon * ( EXP(- alphaa * (variational2 * x) ) &
                        - ( C6 / (variational2 * x)**6 + C8 / (variational2 * x)**8 + C10 / (variational2 * x)**10 ) )

        END IF

    END FUNCTION pseudo_w
    !-------------------------------------------------------


!=========================================================



END PROGRAM HE4_EXCITATION_SPECTRUM