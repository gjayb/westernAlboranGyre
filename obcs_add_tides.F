C $Header: /u/gcmpack/MITgcm/pkg/obcs/obcs_add_tides.F,v 1.2 2012/11/15 20:46:52 dimitri Exp $
C $Name:  $

#include "OBCS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: OBCS_ADD_TIDES

C     !INTERFACE:
      SUBROUTINE OBCS_ADD_TIDES( myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE OBCS_ADD_TIDES
C     | o Modify OB normal flow to add tidal forcing
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
cjc(Defining RAMPA, sigmma, arg, Amp
      _RL RAMPA, sigmma, arg, Amp
cjc)
#ifdef ALLOW_OBCS
#ifdef ALLOW_OBCS_TIDES

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
      IF (debugMode) CALL DEBUG_ENTER('OBCS_ADD_TIDES',myThid)
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
C--   Add tidal currents:
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_OBCS_EAST
        IF ( tileHasOBE(bi,bj) ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           iB = OB_Ie(j,bi,bj)
           IF ( iB.NE.OB_indexNone ) THEN
            DO td=1,tidalComponents
cjc(Including tides East
             sigmma = 2.D0*PI/tidalPeriod(td)
             arg = sigmma*myTime - OBEph(j,td,bi,bj) + 
     &             (V0cn(td)+ucn(td))*2.D0*PI 
             Amp = fcn(td)*OBEam(j,td,bi,bj)
            
             OBEu(j,k,bi,bj) = OBEu(j,k,bi,bj) +
     &                         RAMPA*maskW(iB,j,k,bi,bj)*
     &                         Amp*COS(arg)
cjc)
c-orig              OBEu(j,k,bi,bj) = OBEu(j,k,bi,bj) +
c-orig     &        maskW(iB,j,k,bi,bj) * OBEam(j,td,bi,bj) *
c-orig     &        COS( 2.D0 * PI * (myTime-OBEph(j,td,bi,bj)) /
c-orig     &             tidalPeriod(td) )
cjc)
            ENDDO
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
            DO td=1,tidalComponents
cjc(Including tides West
             sigmma = 2.D0*PI/tidalPeriod(td)
             arg = sigmma*myTime - OBWph(j,td,bi,bj) +
     &             (V0cn(td)+ucn(td))*2.D0*PI
             Amp = fcn(td)*OBWam(j,td,bi,bj)

             OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj) +
     &                         RAMPA*maskW(1+iB,j,k,bi,bj)*
     &                         Amp*COS(arg)
cjc) 
c-orig                    OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj) +
c-orig     &              maskW(1+iB,j,k,bi,bj) * OBWam(j,td,bi,bj) *
c-orig     &              COS( 2.D0 * PI * (myTime-OBWph(j,td,bi,bj)) /
c-orig     &                   tidalPeriod(td) )
            ENDDO
           ENDIF
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_OBCS_WEST */

#ifdef ALLOW_OBCS_NORTH
        IF ( tileHasOBN(bi,bj) ) THEN
         DO k=1,Nr
          DO i=1-OLx,sNx+OLx
           jB = OB_Jn(i,bi,bj)
           IF ( jB.NE.OB_indexNone ) THEN
            DO td=1,tidalComponents
cjc(Including tides North
             sigmma = 2.D0*PI/tidalPeriod(td)
             arg = sigmma*myTime - OBNph(i,td,bi,bj) +
     &             (V0cn(td)+ucn(td))*2.D0*PI
             Amp = fcn(td)*OBNam(i,td,bi,bj)

             OBNv(i,k,bi,bj) = OBNv(i,k,bi,bj) +
     &                         RAMPA*maskS(i,jB,k,bi,bj)*
     &                         Amp*COS(arg)
cjc)
c-orig             OBNv(i,k,bi,bj) = OBNv(i,k,bi,bj) +
c-orig     &              maskS(i,jB,k,bi,bj) * OBNam(i,td,bi,bj) *
c-orig     &              COS( 2.D0 * PI * (myTime-OBNph(i,td,bi,bj)) /
c-orig     &                   tidalPeriod(td) )
            ENDDO
           ENDIF
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_OBCS_NORTH */

#ifdef ALLOW_OBCS_SOUTH
        IF ( tileHasOBS(bi,bj) ) THEN
         DO k=1,Nr
          DO i=1-OLx,sNx+OLx
           jB = OB_Js(i,bi,bj)
           IF ( jB.NE.OB_indexNone ) THEN
            DO td=1,tidalComponents
cjc(Inserting tides South
             sigmma = 2.D0*PI/tidalPeriod(td)
             arg = sigmma*myTime - OBSph(i,td,bi,bj) +
     &             (V0cn(td)+ucn(td))*2.D0*PI
             Amp = fcn(td)*OBSam(i,td,bi,bj)

             OBSv(i,k,bi,bj) = OBSv(i,k,bi,bj) +
     &                         RAMPA*maskS(i,1+jB,k,bi,bj)*
     &                         Amp*COS(arg)
cjc)
c-orig                OBSv(i,k,bi,bj) = OBSv(i,k,bi,bj) +
c-orig     &              maskS(i,1+jB,k,bi,bj)* OBSam(i,td,bi,bj) *
c-orig     &              COS( 2.D0 * PI * (myTime-OBSph(i,td,bi,bj)) /
c-orig     &                   tidalPeriod(td) )
            ENDDO
           ENDIF
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_OBCS_SOUTH */

       ENDDO
      ENDDO

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('OBCS_ADD_TIDES',myThid)
#endif

#endif /* ALLOW_OBCS_TIDES */
#endif /* ALLOW_OBCS */

      RETURN
      END
