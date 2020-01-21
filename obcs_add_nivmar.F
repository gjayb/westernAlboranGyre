C $Header: /u/gcmpack/MITgcm/pkg/obcs/obcs_add_nivmar.F,v 1.2 2013/11/09 20:46:52 jc Exp $
C $Name:  $

#include "OBCS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: OBCS_ADD_TIDES

C     !INTERFACE:
      SUBROUTINE OBCS_ADD_NIVMAR( myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE OBCS_ADD_NIVMAR
C     | o Modify OB normal flow to add storm surge velocity
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "OBCS_PARAMS.h"
#include "OBCS_GRID.h"
#include "OBCS_FIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP
cjc(Defining RAMPA
      _RL RAMPA
cjc)
#ifdef ALLOW_OBCS

C     !FUNCTIONS:

C     !LOCAL VARIABLES:
C     bi, bj       :: tile indices
C     i,j,k        :: loop indices
C     iB, jB       :: local index of open boundary
C     msgBuf       :: Informational/error message buffer
      INTEGER bi, bj
      INTEGER i, j, k, iB, jB
      INTEGER td

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('OBCS_ADD_NIVMAR',myThid)
#endif

cjc(Calculating RAMPA
      RAMPA=1.
      IF ( Select_Rampa ) THEN
        IF ( myTime .GE. DtRampa ) THEN
             RAMPA = 1.
        ELSE
             RAMPA = myTime/DtRampa
        ENDIF
      ENDIF
cjc)
C--   Add nivmar currents:
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_OBCS_EAST
        IF ( tileHasOBE(bi,bj) ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           iB = OB_Ie(j,bi,bj)
           IF ( iB.NE.OB_indexNone ) THEN
             OBEu(j,k,bi,bj) = OBEu(j,k,bi,bj) +
     &              RAMPA*maskW(iB,j,k,bi,bj)*
     &              OBEuNiv(j,bi,bj) 
           ENDIF
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_OBCS_EAST */

#ifdef ALLOW_OBCS_WEST
        IF ( tileHasOBW(bi,bj) ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           iB = OB_Iw(j,bi,bj)
           IF ( iB.NE.OB_indexNone ) THEN
             OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj) +
     &              RAMPA*maskW(1+iB,j,k,bi,bj)* 
     &              OBWuNiv(j,bi,bj)
           ENDIF
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_OBCS_WEST */

       ENDDO
      ENDDO

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('OBCS_ADD_NIVMAR',myThid)
#endif

#endif /* ALLOW_OBCS */

      RETURN
      END
