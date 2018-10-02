

  SUBROUTINE meridional_boundary_tz(field,locations,nt,nz,ny,nx,obc)

    IMPLICIT NONE

    REAL(4),DIMENSION(nt,nz,ny,nx),INTENT(in) :: field
    INTEGER,DIMENSION(ny),INTENT(in) :: locations
    INTEGER,INTENT(in) :: nt,nz,ny,nx
    REAL(4),DIMENSION(nt,nz,ny),INTENT(out) :: obc
    INTEGER :: jj,jk,jt


    DO jj=1,ny
       DO jk=1,nz
          DO jt=1,nt
             IF (locations(jj) == 0) THEN
                ! land
                obc(jt,jk,jj) = 0
             ELSE
                obc(jt,jk,jj) = field(jt,jk,jj,locations(jj))
             ENDIF
          ENDDO
       ENDDO
    ENDDO



  END SUBROUTINE
