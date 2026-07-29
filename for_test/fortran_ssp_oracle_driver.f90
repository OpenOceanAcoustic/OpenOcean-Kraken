PROGRAM FortranSSPOracle
  USE sspmod
  IMPLICIT NONE
  INTEGER :: i, nOut
  CHARACTER(LEN=1) :: profileType
  CHARACTER(LEN=16) :: argument
  COMPLEX(KIND=8) :: cpOut(5), csOut(5)
  REAL(KIND=8) :: rhoOut(5)

  CALL GET_COMMAND_ARGUMENT(1, argument)
  profileType = argument(1:1)
  IF (INDEX('NCPS', profileType) == 0) &
    ERROR STOP 'type must be N, C, P, or S'

  SSP%NMedia = 1
  SSP%Type = profileType
  SSP%AttenUnit = 'W '
  SSP%Loc(1) = 0
  SSP%NPts(1) = 4
  SSP%Depth(1) = 50.0D0
  SSP%Depth(2) = 100.0D0
  SSP%beta(1) = 0.0D0
  SSP%fT(1) = 100.0D0
  SSP%z(1:4) = [50.0D0, 65.0D0, 80.0D0, 100.0D0]
  SSP%alphaR(1:4) = [1700.0D0, 1750.0D0, 1820.0D0, 1900.0D0]
  SSP%betaR(1:4) = [700.0D0, 760.0D0, 850.0D0, 950.0D0]
  SSP%rho(1:4) = [1.40D0, 1.50D0, 1.70D0, 1.90D0]
  SSP%alphaI(1:4) = [0.20D0, 0.30D0, 0.40D0, 0.50D0]
  SSP%betaI(1:4) = [0.10D0, 0.15D0, 0.20D0, 0.25D0]

  CALL UpdateSSPLoss(100.0D0, 100.0D0)
  nOut = 5
  CALL EvaluateSSP(cpOut, csOut, rhoOut, 1, nOut, &
                   100.0D0, 'TABULATE')
  DO i = 1, nOut
     WRITE(*,'(A1,1X,ES25.17,1X,5(ES25.17,1X))') profileType, &
       50.0D0 + (i-1)*12.5D0, REAL(cpOut(i)), AIMAG(cpOut(i)), &
       REAL(csOut(i)), AIMAG(csOut(i)), rhoOut(i)
  END DO
END PROGRAM FortranSSPOracle
