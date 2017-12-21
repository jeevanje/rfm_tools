      PROGRAM HITBIN
C
C VERSION (update VERSID)
C     14SEP17 AD 2.31 Read ISO in .par file as Z1 instead of I1
C     15MAR16 AD 2.30 Change Class10 (CH4) Vib Level assignments
C                     Correction to Strength calc for old HITRAN data
C                     Conform to gfortran compilation with -Wall
C     21AUG15 AD 2.22 Fix bug detecting old format files
C                     Check for read errors and print offending record.
C     03FEB14 AD 2.21 Fix bug with copying CO2 rotational ID
C     10JAN14 AD 2.20 Allow for additional non-HITRAN molecules >#47
C     24OCT13 AD 2.10 Correction: Change ISO=0 from input to 10 for output
C                     (currently only applies to CO2 isotopologue 838)
C     22AUG13 AD 2.01 Limit number of warning messages within IDXVIB
C     26JUL13 AD 2.00 Simplified code - remove some options from old version
C                     Directly convert new or old HITRAN format files
C     30JUN09 AD 1.20 HITRAN2008: allow for SBROAD in F5.4 or F5.3 (2008)
C     09MAR02 AD 1.10 DUPCHK: if two lines have equal status change from 
C                       termination condition to warning.
C                       Remove redundant MAXNEW parameter.
C     09SEP01 AD 1.00 Original. Based on HITLIN.
C
C DESCRIPTION
C     Convert HITRAN sequential ASCII file to direct access binary file
C     Note: the record length of the unformatted output file depends
C           on the no of bytes per word of the machine being used. 
C           System specific RECL is used in the output binary file 
C           open statement
C
C     Variables for each transition record (100 characters in input records)
C              Input  Output Description
C             Old/New
C      LSTAT     -      I 4   Priority of transition information.
C      MOL      I2      I 8   Molecule number.
C      ISO      I1      I 12  Isotope number (=1 most abundant, 2=second etc)
C      WNO     F12.6    D 20  Line frequency [cm-1].
C      STR     D10.3    R 24  Line strength  I=[cm-1./(molec.cm-2)] @ 296K.
C                                          O=[cm-1./(kg.moles.cm-2)] @296K
C                           (Scale input by Avogadro No.to avoid underflows.)
C      TPR     E10.3    R 28  Transition probab. [Debyes2] (old) or Einstein A
C      ABR      F5.4    R 32  Air-broadened halfwidth  (HWHM) [cm-1/atm] @ 296K.
C      SBR    F5.4/F5.3 R 36  Self-broadened halfwidth (HWHM) [cm-1/atm] @ 296K.
C      ELS     F10.4    R 40   Lower-state energy [cm-1].
C      ABC      F4.2    R 44   Coefficient of temperature dependance of ABROAD
C      TSP      F8.6    R 48   Transition shift due to pressure 
C      IUV    I3/A15    I 52    Upper state global (=Vib) quanta index / ID
C      ILV    I3/A15    I 56    Lower state global quanta index / ID
C      IUR    A9/A15   A9 65   Upper state local (=Rot) quanta. 
C      ILR    A9/A15   A9 74   Lower state local quanta. 
C      SPARE9          A9 83   Spare (used to be quality and reference info)
C      IFP       -      I 87    Forward pointer on data line.
C      SPARE1          A1 88   Spare (used to be empty byte at end of record)
C
C     LSTAT values for different types of binary file record
C        -7 forward pointer block
C        -2 file termination record
C         0 file header record
C         4 other header records
C        10 line transition
C
      IMPLICIT NONE
C
      EXTERNAL
     &  IDXVIB ! Translate new HITRAN C*15 Vib. ID to old integer index
     &, IRECL  ! Return no.bytes assumed for RECL parameter in OPEN statement
      INTEGER IDXVIB
      INTEGER IRECL
C
C LOCAL CONSTANTS
      INTEGER IFORM               ! format version# for HITRAN output
        PARAMETER ( IFORM = 1 )
      INTEGER LUNHIT              ! LUN for ASCII (input) file
        PARAMETER ( LUNHIT = 1 )
      INTEGER LUNBIN              ! LUN for new binary (output) file
        PARAMETER ( LUNBIN = 2 )
      INTEGER MAXNFP              ! Max no.of forward pointers per record
        PARAMETER ( MAXNFP = 14 ) ! Don't change this
      INTEGER NFPREC              ! No. forward pointer records in a f.p.block
        PARAMETER ( NFPREC = 4 )  
      INTEGER MAXMOL              ! Max HITRAN molecule index allowed
        PARAMETER ( MAXMOL = NFPREC*MAXNFP ) ! Change NFPREC to increase MAXMOL
      INTEGER NRECFP              ! No of records between f.p. blocks
        PARAMETER ( NRECFP = 200 )
      REAL AVOG                   ! Avogadro's number (*1e-26)
        PARAMETER ( AVOG = 6.02214199 )
      DOUBLE PRECISION WNOBIG     ! Larger Wno than likely to occur in HITRAN
        PARAMETER ( WNOBIG = 1.0D6 )  
C
C LOCAL VARIABLES
      INTEGER IFP             ! Forward pointer read from binary file
      INTEGER IFPREC          ! Record# within forward pointer block
      INTEGER ILV             ! Lower level vibrational ID (old format HITRAN)
      INTEGER IMOL            ! Molecule counter
      INTEGER IOFF            ! Offset for molecule# in f.p. records
      INTEGER IOSVAL          ! Saved value of IOSTAT
      INTEGER IRCMOL(MAXMOL)  ! Next record# containing transition for each mol
      INTEGER IREC            ! Record# in new binary file
      INTEGER IREC1           ! First record# in new file after headers
      INTEGER IREC2           ! Last record# in new file (termination rec)
      INTEGER IRECFP          ! First record# in new file after last f.p.block
      INTEGER ISO             ! Isotope ID 
      INTEGER IUV             ! Upper level vibrational ID (old format HITRAN)
      INTEGER JREC            ! Secondary record counter
      INTEGER LSTAT           ! Type of binary file record
      INTEGER MOL             ! Molecule ID 
      INTEGER RECLEN          ! Record length of binary files
      LOGICAL OLDFMT          ! T=old HITRAN format (<2004), F=new format
      REAL    ABC             ! Air Broad. Temp. Coeff 
      REAL    ABR             ! Air Broadened halfwidth
      REAL    ELS             ! Energy of the lower state
      REAL    SBR             ! Self broadened halfwidth
      REAL    STR             ! Line strength
      REAL    TPR             ! Transition Probability
      REAL    TSP             ! Pressure shift
      DOUBLE PRECISION DSTR   ! Line strength allowing for < 1.0E-38
      DOUBLE PRECISION WNO    ! Wavenumber 
      DOUBLE PRECISION WNOUPP ! Last wavenumber written/read
      DOUBLE PRECISION WNORQ1 ! Lowest wavenumber to select from ASCII file
      DOUBLE PRECISION WNORQ2 ! Highest wavenumber to select from ASCII file
      CHARACTER*80 FILNAM     ! User-input name of input and output files
      CHARACTER*48 HEAD48     ! User header for binary file
      CHARACTER*84 HEADR2     ! HITBIN header record for binary file
      CHARACTER*9  ILR        ! Lower level rotational ID (old format HITRAN)
      CHARACTER*9  IUR        ! Upper level rotational ID (old format HITRAN)
      CHARACTER*15 JLR        ! Lower level rotational ID (new format HITRAN)
      CHARACTER*15 JLV        ! Lower level vibrational ID (new format HITRAN)
      CHARACTER*15 JUR        ! Upper level rotational ID (new format HITRAN)
      CHARACTER*15 JUV        ! Upper level vibrational ID (new format HITRAN)
      CHARACTER*80  RECORD    ! User input for wavenumber range
      CHARACTER*160 REC160    ! 160-character record for old/new HITRAN format
      CHARACTER*1  SPARE1     ! Spare byte in output record
      CHARACTER*9  SPARE9     ! Spare bytes in output record (was C*9 acc.index
C                               and data references)
      CHARACTER*4  VERSID     ! Program Version identifier
C
      DATA SPARE9 / '         ' /    ! 9 spaces
      DATA SPARE1 / ' ' /            ! 1 space
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      VERSID = '2.31'
      WRITE ( *, '(A,A,A,I2)' ) 'R-HITBIN: Running HITBIN v', VERSID, 
     &                   ' generating output Format#', IFORM
C
C Get HITRAN (ASCII) input file and open
      WRITE ( *, '(A$)' ) 'Input HITRAN file: '
      READ ( *, '(A)' ) FILNAM
      OPEN ( UNIT=LUNHIT, FILE=FILNAM, STATUS='OLD' )
C Check format
      OLDFMT = .FALSE.
      REC160 = ' '
      READ ( LUNHIT, '(A)', ERR=900, IOSTAT=IOSVAL ) REC160
      OLDFMT = REC160(102:) .EQ. ' '
      IF ( OLDFMT ) WRITE (*,*) 
     &  'File identifed as old (pre 2004) format'
      REWIND ( LUNHIT )
C
C Get wavenumber range for ASCII file (default = use all)
      WRITE ( *, '(A$)' ) 'Wavenumber range (cm-1) [<CR>=all]: '
      READ ( *, '(A)' ) RECORD
      IF ( RECORD .EQ. ' ' ) THEN
        WNORQ1 = -1.0D0
        WNORQ2 = WNOBIG
      ELSE
        READ ( RECORD, * ) WNORQ1, WNORQ2
        WRITE ( *, '(A)' ) 
     &    'I-HITBIN: Finding first record of ASCII file...'
        WNO = -1.0D0
        DO WHILE ( WNO .LT. WNORQ1 )
          READ ( LUNHIT, '(3X,F12.6)', END=900 ) WNO
        END DO
        IF ( WNO .GT. WNORQ2 ) GOTO 900           ! no overlap with reqd range
        BACKSPACE ( LUNHIT )       
      END IF
C
C Get name of the new binary file to be created
      WRITE ( *, '(A$)' ) 'New binary file: '
      READ ( *, '(A)' ) FILNAM
C Get record length for new binary file
      RECLEN = 22 * IRECL ( LUNBIN )
      WRITE ( *, * ) 'Assuming appropriate RECL value is:', RECLEN
      OPEN ( UNIT=LUNBIN, FILE=FILNAM, STATUS='NEW', ACCESS='DIRECT',
     &       RECL=RECLEN )
C
C Get header for new file 
      WRITE ( *, '(A$)' ) 'Header for new file (up to 48 chars): '
      READ ( *, '(A)' ) HEAD48
C
C Write HITBIN header (record length and version ID) to rec#2
      WRITE ( HEADR2, '(A,I3,A,A)' ) 'RECL=', RECLEN, 
     &  ' Converted to binary format by HITBIN v.', VERSID
      WRITE ( LUNBIN, REC=2 ) 4, HEADR2
C
      IREC = 3
      IREC1 = IREC
C
      WRITE ( *, '(A)' ) 'I-HITBIN: Writing new binary file...'
C Begin output data with an empty forward pointer block
      DO JREC = 1, NFPREC 
        WRITE ( LUNBIN, REC=IREC ) -7
        IREC = IREC + 1
      END DO
      IRECFP = IREC 
C
C Repeat for each record, taking lowest wavenumber 
      DO WHILE ( WNO .LE. WNORQ2 ) 
C Check if forward pointer block required
        IF ( IREC - IRECFP .EQ. NRECFP ) THEN          ! insert f.p. block
          DO JREC = 1, NFPREC  
            WRITE ( LUNBIN, REC=IREC ) -7
            IREC = IREC + 1
          END DO
          IRECFP = IREC
        END IF

C Copy record from ASCII file to new file
C To avoid underflow problems, STR read as double precision and scaled by
C Avogadro's number before converting to single precision 
C Note that format descriptors Fw.d (eg F5.2) do not generally represent the
C actual HITRAN format (eg F5.3) but simply set to avoid warning messages from
C some compilers (eg ifort) which recommend w ge d+3 
        READ ( LUNHIT, '(A)', END=800, ERR=901, IOSTAT=IOSVAL ) REC160
        IF ( OLDFMT ) THEN 
          READ ( REC160, 1000, ERR=902, IOSTAT=IOSVAL ) 
     &      MOL, ISO, WNO, DSTR, TPR,
     &      ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR
 1000     FORMAT( I2, I1, F12.6, F10.3, E10.3, 
     &            f5.2, F5.2, F10.4, F4.1, F8.5, 2I3, 2A9 )
        ELSE
          READ ( REC160, 1001, ERR=902, IOSTAT=IOSVAL ) 
     &      MOL, ISO, WNO, DSTR, TPR, 
     &      ABR, SBR, ELS, ABC, TSP, JUV, JLV, JUR, JLR
 1001     FORMAT( I2, Z1, F12.6, F10.3, E10.3, 
     &            F5.2, F5.2, F10.4, F4.1, F8.5, 4A15 )
       
C Convert new format C*15 codes for upper,lower vibrational level into old
C format integer indices required by RFM if particular Vib levels have to be
C identified (user-selected or non-LTE).
          IUV = IDXVIB ( MOL, JUV )
          ILV = IDXVIB ( MOL, JLV )
C Reformat local quantum numbers for CO2 lines (required by RFM for line-mixing)
C The first 9 characters of the HITRAN04 C*15 strings seems to match the 
C old C*9 fields except that in the old format the characters are shifted one
C space left and the new format has an extra character 'e' or 'f' which needs
C to be set blank (this was fixed in v2.21 of this code)
C No idea what this will do for local quantum numbers for other molecules (!)
C but hopefully this information won't be required
          IUR = JUR(2:9)//' '
          ILR = JLR(2:9)//' '
        ENDIF
        STR = SNGL ( DSTR * AVOG * 1.0E26 ) 
C
C HITRAN 2012 has 10 CO2 isotopologues with '0' used for #10. Use 10 in bin file
        IF ( ISO .EQ. 0 ) ISO = 10
C
        IF ( WNO .LE. WNORQ2 ) THEN 
          WRITE ( LUNBIN, REC=IREC ) 10, MOL, ISO, WNO, STR, TPR, ! 10=LSTAT
     &       ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR, 
     &       SPARE9, 0, SPARE1   ! 0 = forward ptr
          IREC = IREC + 1
          WNOUPP = WNO
C
          IF ( MOD ( IREC, 100000 ) .EQ. 0 ) 
     &      WRITE ( *, '(A,I8,A,F12.6)' ) 
     &      'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
        END IF 
      END DO
  800 CONTINUE
C
C Save record after last data record
      IREC2 = IREC
C
C Write header record
      WRITE ( LUNBIN, REC=1 ) 0, IFORM, IREC1, IREC2, HEAD48
C
C Write termination record
      WRITE ( LUNBIN, REC=IREC2 ) -2, 0, 0, WNOUPP, (0,IMOL=1,14)
      WRITE ( *, '(A,I8,A,F12.6)' )
     &  'I-HITBIN: Last Record#', IREC2, ' Wavenumber=', WNOUPP
C
C Calculate forward pointers
      WRITE ( *, '(A)' ) 'I-HITBIN: Writing forward pointers...'
      DO IMOL = 1, MAXMOL
        IRCMOL(IMOL) = IREC2
      END DO
      IFPREC = NFPREC
      DO IREC = IREC2-1, IREC1, -1 
        READ ( LUNBIN, REC=IREC ) LSTAT
        IF ( LSTAT .EQ. -7 ) THEN
          IOFF = ( IFPREC - 1 ) * MAXNFP
          WRITE ( LUNBIN, REC=IREC ) -7, IOFF+1, 0, WNOUPP, 
     &      ( IRCMOL(IMOL)-IREC, IMOL = IOFF+1, IOFF+MAXNFP ) 
          IF ( IFPREC .EQ. 1 ) THEN
            IFPREC = NFPREC
          ELSE
            IFPREC = IFPREC - 1
          END IF
        ELSE
          READ ( LUNBIN, REC=IREC ) LSTAT, MOL, ISO, WNO, STR, TPR, 
     &      ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR, SPARE9, IFP,
     &      SPARE1
          IFP = IRCMOL(MOL) - IREC 
          IRCMOL(MOL) = IREC
          WNOUPP = WNO
          WRITE ( LUNBIN, REC=IREC ) LSTAT, MOL, ISO, WNO, STR, TPR, 
     &      ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR, SPARE9, IFP,
     &      SPARE1
        END IF
        IF ( MOD ( IREC, 100000 ) .EQ. 0 ) 
     &    WRITE ( *, '(A,I8,A,F12.6)' ) 
     &    'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
      END DO
C
      STOP  'R-HITBIN: Successful completion'
C
 900  CONTINUE
      STOP 'F-HITBIN: HITRAN file does not overlap required wno.range'
C
 901  CONTINUE
      WRITE ( *, '(A,I11)' ) 'F-HITBIN: Error reading record '//
     &  'from HITRAN file. IOSTAT=', IOSVAL
      WRITE ( *, '(A)' ) 'Last valid record was:'
      IF ( OLDFMT ) THEN
        WRITE ( *, '(A)' ) REC160(1:100)
      ELSE
        WRITE ( *, '(A)' ) REC160
      END IF
      STOP
C
 902  CONTINUE
      WRITE ( *, '(A,I11)' ) 'F-HITBIN: Error decoding following '//
     &  'record from HITRAN file. IOSTAT=', IOSVAL
      IF ( OLDFMT ) THEN
        WRITE ( *, '(A)' ) REC160(1:100)
      ELSE
        WRITE ( *, '(A)' ) REC160
      END IF
      STOP
C
      END
C
      INTEGER FUNCTION IRECL ( LUN )
C
C VERSION
C     25JUL13  AD  Original.
C
C DESCRIPTION
C     Return no.bytes assumed for RECL parameter in OPEN statement
C     RECL specifies record lengths either in bytes or units of 4 bytes,
C     depending on the compiler and compilation options.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER  LUN  !  I  Spare LUN used for testing (freed after use)
C
C LOCAL VARIABLES
      INTEGER  I, J ! dummy variables written to file
      INTEGER  IOS  ! saved value of IOSTAT
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      I = 0
      J = 0
      IRECL = 1
      OPEN ( UNIT=LUN, RECL=2*IRECL, ACCESS='DIRECT', STATUS='SCRATCH' )
      WRITE ( LUN, REC=1, IOSTAT=IOS ) I, J
      CLOSE ( LUN )
      IF ( IOS .EQ. 0 ) RETURN

      IRECL = 4
      OPEN ( UNIT=LUN, RECL=2*IRECL, ACCESS='DIRECT', STATUS='SCRATCH' )
      WRITE ( LUN, REC=1, IOSTAT=IOS ) I, J
      CLOSE ( LUN )
      IF ( IOS .EQ. 0 ) RETURN
C
      STOP 'F-IRECL: Unable to establish RECL definition'
C
      END
      INTEGER FUNCTION IDXVIB ( IDXMOL, STRVIB )
C
C VERSION
C     15MAR16  AD Change Class10 assignments to match CH4 Vib Temps
C                 Remove local IVIB - set IDXVIB directly
C                 Remove <TAB> characters - replace with spaces
C     22AUG13  AD Limit number of warning messages
C     26JUL13  AD Add new molec 43-47 & vib states found in HITRAN2012
C     30JUN09  AD Add new HITRAN molec#40-42
C     19DEC07  AD Set IVIB locally, and set unidentified levels to 999
C     24MAR06  AD Original.
C
C DESCRIPTION
C     Translate HITRAN C*15 Global Quantum ID to old integer index
C     Based on a program supplied by Javier Martin-Torres (NASA LaRC).
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IDXMOL  !  HITRAN code for molecule
      CHARACTER*15 STRVIB  !  HITRAN04 field defining vibration level
C
C LOCAL CONSTANTS
      INTEGER MAXMOL       !  .GE. Max number of different molecules
        PARAMETER ( MAXMOL = 47 )
      INTEGER MAXWRN       !  Max number of warning messages per molecules
        PARAMETER ( MAXWRN = 10 )
      INTEGER MAXNEW       !  Max number of new molecules beyond MAXMOL
        PARAMETER ( MAXNEW = 10 )
C
C LOCAL VARIABLES
      INTEGER INEW           ! Counter for new molecules
      INTEGER NNEW           ! No. of new molecules found
      INTEGER NWRN(MAXMOL) ! No. warnings per molecule.
      INTEGER NEWLST(MAXNEW) ! List of IDXMOL values for new molecules
C
C DATA STATEMENTS
      DATA NWRN / MAXMOL * 0 /      
      DATA NNEW / 0 /
      DATA NEWLST / MAXNEW * 0 /
      SAVE NWRN, NNEW, NEWLST
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IDXVIB = 0
C 
      IF ( IDXMOL .GT. MAXMOL ) THEN
        DO INEW = 1, NNEW
          IF ( IDXMOL .EQ. NEWLST(INEW) ) RETURN
        END DO
        IF ( NNEW .LT. MAXNEW ) THEN
          NNEW = NNEW + 1
          NEWLST(NNEW) = IDXMOL
          WRITE (*,*) 
     &      'W-IDXVIB: Setting IDXVIB=0 for unidentified IDXMOL=',IDXMOL
        END IF
        RETURN
      END IF
C
C Class 1: Diatomic Molecules
      IF ( IDXMOL .EQ.  5 .OR.                                    ! CO
     &     IDXMOL .EQ. 14 .OR.                                    ! HF
     &     IDXMOL .EQ. 15 .OR.                                    ! HCl
     &     IDXMOL .EQ. 16 .OR.                                    ! HBr
     &     IDXMOL .EQ. 17 .OR.                                    ! HI
     &     IDXMOL .EQ. 22 .OR.                                    ! N2
     &     IDXMOL .EQ. 36 .OR.                                    ! NO+
     &     IDXMOL .EQ. 45 .OR.                                    ! H2
     &     IDXMOL .EQ. 46      ) THEN                             ! CS
        IF ( STRVIB .EQ. '              0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '              1' ) IDXVIB = 2
        IF ( STRVIB .EQ. '              2' ) IDXVIB = 3
        IF ( STRVIB .EQ. '              3' ) IDXVIB = 4
        IF ( STRVIB .EQ. '              4' ) IDXVIB = 5
        IF ( STRVIB .EQ. '              5' ) IDXVIB = 6
        IF ( STRVIB .EQ. '              6' ) IDXVIB = 7
        IF ( STRVIB .EQ. '              7' ) IDXVIB = 8
        IF ( STRVIB .EQ. '              8' ) IDXVIB = 9
        IF ( STRVIB .EQ. '              9' ) IDXVIB = 10
        IF ( STRVIB .EQ. '             10' ) IDXVIB = 11
        IF ( STRVIB .EQ. '             11' ) IDXVIB = 12
        IF ( STRVIB .EQ. '             12' ) IDXVIB = 13
        IF ( STRVIB .EQ. '             13' ) IDXVIB = 14
        IF ( STRVIB .EQ. '             14' ) IDXVIB = 15
        IF ( STRVIB .EQ. '             15' ) IDXVIB = 16
        IF ( STRVIB .EQ. '             16' ) IDXVIB = 17
        IF ( STRVIB .EQ. '             17' ) IDXVIB = 18
        IF ( STRVIB .EQ. '             18' ) IDXVIB = 19
        IF ( STRVIB .EQ. '             19' ) IDXVIB = 20
        IF ( STRVIB .EQ. '             20' ) IDXVIB = 21
C new levels for HF, HCl found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '             21' ) IDXVIB = 22
        IF ( STRVIB .EQ. '             22' ) IDXVIB = 23
        IF ( STRVIB .EQ. '             23' ) IDXVIB = 24
        IF ( STRVIB .EQ. '             24' ) IDXVIB = 25
        IF ( STRVIB .EQ. '             25' ) IDXVIB = 26
        IF ( STRVIB .EQ. '             26' ) IDXVIB = 27
C
C Class 2: Diatomic molecules with different electonic levels
C HITRAN 2012 uses different spacing to previous versions so allow for both
      ELSE IF ( IDXMOL .EQ.  7 ) THEN                              ! O2
        IF ( STRVIB .EQ. '            X 0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '       X      0' ) IDXVIB = 1    ! HITRAN 2012
        IF ( STRVIB .EQ. '            X 1' ) IDXVIB = 2
        IF ( STRVIB .EQ. '       X      1' ) IDXVIB = 2    ! HITRAN 2012
        IF ( STRVIB .EQ. '            a 0' ) IDXVIB = 3
        IF ( STRVIB .EQ. '       a      0' ) IDXVIB = 3    ! HITRAN 2012
        IF ( STRVIB .EQ. '            a 1' ) IDXVIB = 4
        IF ( STRVIB .EQ. '       a      1' ) IDXVIB = 4    ! HITRAN 2012
        IF ( STRVIB .EQ. '            b 0' ) IDXVIB = 5
        IF ( STRVIB .EQ. '       b      0' ) IDXVIB = 5    ! HITRAN 2012
        IF ( STRVIB .EQ. '            b 1' ) IDXVIB = 6
        IF ( STRVIB .EQ. '       b      1' ) IDXVIB = 6    ! HITRAN 2012
        IF ( STRVIB .EQ. '            b 2' ) IDXVIB = 7
        IF ( STRVIB .EQ. '       b      2' ) IDXVIB = 7    ! HITRAN 2012
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 8
        IF ( STRVIB .EQ. '            X 2' ) IDXVIB = 9
        IF ( STRVIB .EQ. '            B 0' ) IDXVIB = 10
        IF ( STRVIB .EQ. '            B 1' ) IDXVIB = 11
        IF ( STRVIB .EQ. '            B 2' ) IDXVIB = 12
        IF ( STRVIB .EQ. '            B 3' ) IDXVIB = 13
        IF ( STRVIB .EQ. '            B 4' ) IDXVIB = 14
        IF ( STRVIB .EQ. '            B 5' ) IDXVIB = 15
        IF ( STRVIB .EQ. '            B 6' ) IDXVIB = 16
        IF ( STRVIB .EQ. '            B 7' ) IDXVIB = 17
        IF ( STRVIB .EQ. '            B 8' ) IDXVIB = 18
        IF ( STRVIB .EQ. '            B 9' ) IDXVIB = 19
        IF ( STRVIB .EQ. '           B 10' ) IDXVIB = 20
        IF ( STRVIB .EQ. '           B 11' ) IDXVIB = 21
        IF ( STRVIB .EQ. '           B 12' ) IDXVIB = 22
        IF ( STRVIB .EQ. '           B 13' ) IDXVIB = 23
        IF ( STRVIB .EQ. '           B 14' ) IDXVIB = 24
        IF ( STRVIB .EQ. '           B 15' ) IDXVIB = 25
        IF ( STRVIB .EQ. '           B 16' ) IDXVIB = 26
        IF ( STRVIB .EQ. '           B 17' ) IDXVIB = 27
        IF ( STRVIB .EQ. '           B 18' ) IDXVIB = 28
        IF ( STRVIB .EQ. '           B 19' ) IDXVIB = 29
C
C Class 3: Diatomic molecules with Pi-doublet electronic state
      ELSE IF ( IDXMOL .EQ.  8 .OR.                                ! NO
     &          IDXMOL .EQ. 13 .OR.                                ! OH
     &          IDXMOL .EQ. 18      ) THEN                         ! ClO
        IF ( STRVIB .EQ. '       X3/2   0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '       X3/2   1' ) IDXVIB = 2
        IF ( STRVIB .EQ. '       X3/2   2' ) IDXVIB = 3
        IF ( STRVIB .EQ. '       X3/2   3' ) IDXVIB = 4
        IF ( STRVIB .EQ. '       X3/2   4' ) IDXVIB = 5
        IF ( STRVIB .EQ. '       X3/2   5' ) IDXVIB = 6
        IF ( STRVIB .EQ. '       X3/2   6' ) IDXVIB = 7
        IF ( STRVIB .EQ. '       X3/2   7' ) IDXVIB = 8
        IF ( STRVIB .EQ. '       X3/2   8' ) IDXVIB = 9
        IF ( STRVIB .EQ. '       X3/2   9' ) IDXVIB = 10
        IF ( STRVIB .EQ. '       X1/2   0' ) IDXVIB = 11
        IF ( STRVIB .EQ. '       X1/2   1' ) IDXVIB = 12
        IF ( STRVIB .EQ. '       X1/2   2' ) IDXVIB = 13
        IF ( STRVIB .EQ. '       X1/2   3' ) IDXVIB = 14
        IF ( STRVIB .EQ. '       X1/2   4' ) IDXVIB = 15
        IF ( STRVIB .EQ. '       X1/2   5' ) IDXVIB = 16
        IF ( STRVIB .EQ. '       X1/2   6' ) IDXVIB = 17
        IF ( STRVIB .EQ. '       X1/2   7' ) IDXVIB = 18
        IF ( STRVIB .EQ. '       X1/2   8' ) IDXVIB = 19
        IF ( STRVIB .EQ. '       X1/2   9' ) IDXVIB = 20
        IF ( STRVIB .EQ. '       X1/2  10' ) IDXVIB = 21
        IF ( STRVIB .EQ. '       X1/2  11' ) IDXVIB = 22
        IF ( STRVIB .EQ. '       X1/2  12' ) IDXVIB = 23
        IF ( STRVIB .EQ. '       X3/2  10' ) IDXVIB = 24
        IF ( STRVIB .EQ. '       X3/2  11' ) IDXVIB = 25
        IF ( STRVIB .EQ. '       X3/2  12' ) IDXVIB = 26
        IF ( STRVIB .EQ. '       A1     0' ) IDXVIB = 27
        IF ( STRVIB .EQ. '       A1     1' ) IDXVIB = 28
        IF ( STRVIB .EQ. '       A1     2' ) IDXVIB = 29
        IF ( STRVIB .EQ. '       A1     3' ) IDXVIB = 30
        IF ( STRVIB .EQ. '       A2     0' ) IDXVIB = 31
        IF ( STRVIB .EQ. '       A2     1' ) IDXVIB = 32
        IF ( STRVIB .EQ. '       A2     2' ) IDXVIB = 33
        IF ( STRVIB .EQ. '       A2     3' ) IDXVIB = 34
        IF ( STRVIB .EQ. '       X3/2  13' ) IDXVIB = 35
        IF ( STRVIB .EQ. '       X3/2  14' ) IDXVIB = 36
        IF ( STRVIB .EQ. '       X1/2  13' ) IDXVIB = 37
        IF ( STRVIB .EQ. '       X1/2  14' ) IDXVIB = 38
C
C Class 4: Linear triatomic molecules
      ELSE IF ( IDXMOL .EQ.  4 .OR.                                 ! N2O
     &          IDXMOL .EQ. 19 .OR.                                 ! OCS
     &          IDXMOL .EQ. 23      ) THEN                          ! HCN
        IF ( STRVIB .EQ. '        0 0 0 0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '        0 1 1 0' ) IDXVIB = 2
        IF ( STRVIB .EQ. '        0 2 0 0' ) IDXVIB = 3
        IF ( STRVIB .EQ. '        0 2 2 0' ) IDXVIB = 4
        IF ( STRVIB .EQ. '        1 0 0 0' ) IDXVIB = 5
        IF ( STRVIB .EQ. '        0 3 1 0' ) IDXVIB = 6
        IF ( STRVIB .EQ. '        0 3 3 0' ) IDXVIB = 7
        IF ( STRVIB .EQ. '        1 1 1 0' ) IDXVIB = 8
        IF ( STRVIB .EQ. '        0 4 0 0' ) IDXVIB = 9
        IF ( STRVIB .EQ. '        0 4 2 0' ) IDXVIB = 10
        IF ( STRVIB .EQ. '        1 2 0 0' ) IDXVIB = 11
        IF ( STRVIB .EQ. '        1 2 2 0' ) IDXVIB = 12
        IF ( STRVIB .EQ. '        2 0 0 0' ) IDXVIB = 13
        IF ( STRVIB .EQ. '        0 0 0 1' ) IDXVIB = 14
        IF ( STRVIB .EQ. '        0 5 1 0' ) IDXVIB = 15
        IF ( STRVIB .EQ. '        1 3 1 0' ) IDXVIB = 16
        IF ( STRVIB .EQ. '        1 3 3 0' ) IDXVIB = 17
        IF ( STRVIB .EQ. '        2 1 1 0' ) IDXVIB = 18
        IF ( STRVIB .EQ. '        0 1 1 1' ) IDXVIB = 19
        IF ( STRVIB .EQ. '        1 4 0 0' ) IDXVIB = 20
        IF ( STRVIB .EQ. '        1 4 2 0' ) IDXVIB = 21
        IF ( STRVIB .EQ. '        2 2 0 0' ) IDXVIB = 22
        IF ( STRVIB .EQ. '        2 2 2 0' ) IDXVIB = 23
        IF ( STRVIB .EQ. '        3 0 0 0' ) IDXVIB = 24
        IF ( STRVIB .EQ. '        0 2 0 1' ) IDXVIB = 25
        IF ( STRVIB .EQ. '        0 2 2 1' ) IDXVIB = 26
        IF ( STRVIB .EQ. '        1 0 0 1' ) IDXVIB = 27
        IF ( STRVIB .EQ. '        2 3 1 0' ) IDXVIB = 28
        IF ( STRVIB .EQ. '        3 1 1 0' ) IDXVIB = 29
        IF ( STRVIB .EQ. '        0 3 1 1' ) IDXVIB = 30
        IF ( STRVIB .EQ. '        0 3 3 1' ) IDXVIB = 31
        IF ( STRVIB .EQ. '        1 1 1 1' ) IDXVIB = 32
        IF ( STRVIB .EQ. '        4 0 0 0' ) IDXVIB = 33
        IF ( STRVIB .EQ. '        3 2 0 0' ) IDXVIB = 34
        IF ( STRVIB .EQ. '        2 0 0 1' ) IDXVIB = 35
        IF ( STRVIB .EQ. '        1 2 0 1' ) IDXVIB = 36
        IF ( STRVIB .EQ. '        1 2 2 1' ) IDXVIB = 37
        IF ( STRVIB .EQ. '        0 0 0 2' ) IDXVIB = 38
        IF ( STRVIB .EQ. '        2 1 1 1' ) IDXVIB = 39
        IF ( STRVIB .EQ. '        0 1 1 2' ) IDXVIB = 40
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 41
        IF ( STRVIB .EQ. '        0 6 0 0' ) IDXVIB = 42
        IF ( STRVIB .EQ. '        0 6 2 0' ) IDXVIB = 43
        IF ( STRVIB .EQ. '        0 4 4 0' ) IDXVIB = 44
        IF ( STRVIB .EQ. '        0 5 3 0' ) IDXVIB = 45
        IF ( STRVIB .EQ. '        0 4 4 1' ) IDXVIB = 46
        IF ( STRVIB .EQ. '        0 4 2 1' ) IDXVIB = 47
        IF ( STRVIB .EQ. '        0 4 0 1' ) IDXVIB = 48
        IF ( STRVIB .EQ. '        1 5 1 0' ) IDXVIB = 49
        IF ( STRVIB .EQ. '        1 5 3 0' ) IDXVIB = 50
        IF ( STRVIB .EQ. '        2 3 3 0' ) IDXVIB = 51
        IF ( STRVIB .EQ. '        0 5 3 1' ) IDXVIB = 52
        IF ( STRVIB .EQ. '        0 5 1 1' ) IDXVIB = 53
        IF ( STRVIB .EQ. '        1 0 0 2' ) IDXVIB = 54
        IF ( STRVIB .EQ. '        1 3 1 1' ) IDXVIB = 55
        IF ( STRVIB .EQ. '        1 3 3 1' ) IDXVIB = 56
        IF ( STRVIB .EQ. '        0 7 3 0' ) IDXVIB = 57
        IF ( STRVIB .EQ. '        0 7 1 0' ) IDXVIB = 58
        IF ( STRVIB .EQ. '        1 6 0 0' ) IDXVIB = 59
        IF ( STRVIB .EQ. '        2 4 0 0' ) IDXVIB = 60
        IF ( STRVIB .EQ. '        2 4 2 0' ) IDXVIB = 61
        IF ( STRVIB .EQ. '        4 1 1 0' ) IDXVIB = 62
        IF ( STRVIB .EQ. '        3 2 2 0' ) IDXVIB = 63
        IF ( STRVIB .EQ. '        0 2 2 2' ) IDXVIB = 64
        IF ( STRVIB .EQ. '        0 2 0 2' ) IDXVIB = 65
        IF ( STRVIB .EQ. '        1 4 0 1' ) IDXVIB = 66
        IF ( STRVIB .EQ. '        1 4 2 1' ) IDXVIB = 67
        IF ( STRVIB .EQ. '        2 2 0 1' ) IDXVIB = 68
        IF ( STRVIB .EQ. '        3 0 0 1' ) IDXVIB = 69
        IF ( STRVIB .EQ. '        2 5 1 0' ) IDXVIB = 70
        IF ( STRVIB .EQ. '        4 2 0 0' ) IDXVIB = 71
        IF ( STRVIB .EQ. '        3 3 1 0' ) IDXVIB = 72
        IF ( STRVIB .EQ. '        0 6 2 1' ) IDXVIB = 73
        IF ( STRVIB .EQ. '        1 1 1 2' ) IDXVIB = 74
        IF ( STRVIB .EQ. '        2 3 1 1' ) IDXVIB = 75
        IF ( STRVIB .EQ. '        3 1 1 1' ) IDXVIB = 76
        IF ( STRVIB .EQ. '        3 4 0 0' ) IDXVIB = 77
        IF ( STRVIB .EQ. '        5 0 0 0' ) IDXVIB = 78
        IF ( STRVIB .EQ. '        0 1 1 3' ) IDXVIB = 79
        IF ( STRVIB .EQ. '        0 0 0 3' ) IDXVIB = 80
        IF ( STRVIB .EQ. '        2 0 0 2' ) IDXVIB = 81
        IF ( STRVIB .EQ. '        3 2 0 1' ) IDXVIB = 82
        IF ( STRVIB .EQ. '        4 0 0 1' ) IDXVIB = 83
        IF ( STRVIB .EQ. '        1 0 0 3' ) IDXVIB = 84
        IF ( STRVIB .EQ. '        2 2 2 1' ) IDXVIB = 85
        IF ( STRVIB .EQ. '        0 9 1 0' ) IDXVIB = 86
        IF ( STRVIB .EQ. '        0 7 1 1' ) IDXVIB = 87
        IF ( STRVIB .EQ. '        0 2 0 2' ) IDXVIB = 88
        IF ( STRVIB .EQ. '        0 5 1 1' ) IDXVIB = 89
        IF ( STRVIB .EQ. '        1 1 1 2' ) IDXVIB = 90
        IF ( STRVIB .EQ. '        3 4 2 0' ) IDXVIB = 91
C new levels for OCS found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '        0 3 1 2' ) IDXVIB = 92
        IF ( STRVIB .EQ. '        0 4 0 2' ) IDXVIB = 93
        IF ( STRVIB .EQ. '        0 4 2 2' ) IDXVIB = 94
        IF ( STRVIB .EQ. '        1 2 0 2' ) IDXVIB = 95
        IF ( STRVIB .EQ. '        1 2 2 2' ) IDXVIB = 96
C
C Class 5: Linear triatomic molecules with large Fermi resonance
      ELSE IF ( IDXMOL .EQ.  2 ) THEN                                ! CO2
        IF ( STRVIB .EQ. '       0 0 0 01' ) IDXVIB = 1
        IF ( STRVIB .EQ. '       0 1 1 01' ) IDXVIB = 2
        IF ( STRVIB .EQ. '       1 0 0 02' ) IDXVIB = 3
        IF ( STRVIB .EQ. '       0 2 2 01' ) IDXVIB = 4
        IF ( STRVIB .EQ. '       1 0 0 01' ) IDXVIB = 5
        IF ( STRVIB .EQ. '       1 1 1 02' ) IDXVIB = 6
        IF ( STRVIB .EQ. '       0 3 3 01' ) IDXVIB = 7
        IF ( STRVIB .EQ. '       1 1 1 01' ) IDXVIB = 8
        IF ( STRVIB .EQ. '       0 0 0 11' ) IDXVIB = 9
        IF ( STRVIB .EQ. '       2 0 0 03' ) IDXVIB = 10
        IF ( STRVIB .EQ. '       1 2 2 02' ) IDXVIB = 11
        IF ( STRVIB .EQ. '       2 0 0 02' ) IDXVIB = 12
        IF ( STRVIB .EQ. '       0 4 4 01' ) IDXVIB = 13
        IF ( STRVIB .EQ. '       1 2 2 01' ) IDXVIB = 14
        IF ( STRVIB .EQ. '       2 0 0 01' ) IDXVIB = 15
        IF ( STRVIB .EQ. '       0 1 1 11' ) IDXVIB = 16
        IF ( STRVIB .EQ. '       2 1 1 03' ) IDXVIB = 17
        IF ( STRVIB .EQ. '       1 3 3 02' ) IDXVIB = 18
        IF ( STRVIB .EQ. '       2 1 1 02' ) IDXVIB = 19
        IF ( STRVIB .EQ. '       0 5 5 01' ) IDXVIB = 20
        IF ( STRVIB .EQ. '       1 3 3 01' ) IDXVIB = 21
        IF ( STRVIB .EQ. '       2 1 1 01' ) IDXVIB = 22
        IF ( STRVIB .EQ. '       1 0 0 12' ) IDXVIB = 23
        IF ( STRVIB .EQ. '       0 2 2 11' ) IDXVIB = 24
        IF ( STRVIB .EQ. '       1 0 0 11' ) IDXVIB = 25
        IF ( STRVIB .EQ. '       3 0 0 04' ) IDXVIB = 26
        IF ( STRVIB .EQ. '       2 2 2 03' ) IDXVIB = 27
        IF ( STRVIB .EQ. '       1 4 4 02' ) IDXVIB = 28
        IF ( STRVIB .EQ. '       3 0 0 03' ) IDXVIB = 29
        IF ( STRVIB .EQ. '       2 2 2 02' ) IDXVIB = 30
        IF ( STRVIB .EQ. '       0 6 6 01' ) IDXVIB = 31
        IF ( STRVIB .EQ. '       3 0 0 02' ) IDXVIB = 32
        IF ( STRVIB .EQ. '       1 4 4 01' ) IDXVIB = 33
        IF ( STRVIB .EQ. '       2 2 2 01' ) IDXVIB = 34
        IF ( STRVIB .EQ. '       3 0 0 01' ) IDXVIB = 35
        IF ( STRVIB .EQ. '       1 1 1 12' ) IDXVIB = 36
        IF ( STRVIB .EQ. '       0 3 3 11' ) IDXVIB = 37
        IF ( STRVIB .EQ. '       1 1 1 11' ) IDXVIB = 38
        IF ( STRVIB .EQ. '       0 0 0 21' ) IDXVIB = 39
        IF ( STRVIB .EQ. '       3 1 1 04' ) IDXVIB = 40
        IF ( STRVIB .EQ. '       3 1 1 03' ) IDXVIB = 41
        IF ( STRVIB .EQ. '       3 1 1 02' ) IDXVIB = 42
        IF ( STRVIB .EQ. '       2 0 0 13' ) IDXVIB = 43
        IF ( STRVIB .EQ. '       1 2 2 12' ) IDXVIB = 44
        IF ( STRVIB .EQ. '       2 3 3 01' ) IDXVIB = 45
        IF ( STRVIB .EQ. '       3 1 1 01' ) IDXVIB = 46
        IF ( STRVIB .EQ. '       0 4 4 11' ) IDXVIB = 47
        IF ( STRVIB .EQ. '       2 0 0 12' ) IDXVIB = 48
        IF ( STRVIB .EQ. '       1 2 2 11' ) IDXVIB = 49
        IF ( STRVIB .EQ. '       2 0 0 11' ) IDXVIB = 50
        IF ( STRVIB .EQ. '       0 1 1 21' ) IDXVIB = 51
        IF ( STRVIB .EQ. '       4 0 0 04' ) IDXVIB = 52
        IF ( STRVIB .EQ. '       3 2 2 03' ) IDXVIB = 53
        IF ( STRVIB .EQ. '       2 1 1 13' ) IDXVIB = 54
        IF ( STRVIB .EQ. '       4 0 0 02' ) IDXVIB = 55
        IF ( STRVIB .EQ. '       1 3 3 12' ) IDXVIB = 56
        IF ( STRVIB .EQ. '       0 5 5 11' ) IDXVIB = 57
        IF ( STRVIB .EQ. '       2 1 1 12' ) IDXVIB = 58
        IF ( STRVIB .EQ. '       1 3 3 11' ) IDXVIB = 59
        IF ( STRVIB .EQ. '       2 1 1 11' ) IDXVIB = 60
        IF ( STRVIB .EQ. '       1 0 0 22' ) IDXVIB = 61
        IF ( STRVIB .EQ. '       0 2 2 21' ) IDXVIB = 62
        IF ( STRVIB .EQ. '       1 0 0 21' ) IDXVIB = 63
        IF ( STRVIB .EQ. '       3 0 0 14' ) IDXVIB = 64
        IF ( STRVIB .EQ. '       2 2 2 13' ) IDXVIB = 65
        IF ( STRVIB .EQ. '       1 4 4 12' ) IDXVIB = 66
        IF ( STRVIB .EQ. '       4 1 1 02' ) IDXVIB = 67
        IF ( STRVIB .EQ. '       3 0 0 13' ) IDXVIB = 68
        IF ( STRVIB .EQ. '       0 6 6 11' ) IDXVIB = 69
        IF ( STRVIB .EQ. '       2 2 2 12' ) IDXVIB = 70
        IF ( STRVIB .EQ. '       3 0 0 12' ) IDXVIB = 71
        IF ( STRVIB .EQ. '       4 1 1 01' ) IDXVIB = 72
        IF ( STRVIB .EQ. '       1 4 4 11' ) IDXVIB = 73
        IF ( STRVIB .EQ. '       2 2 2 11' ) IDXVIB = 74
        IF ( STRVIB .EQ. '       3 0 0 11' ) IDXVIB = 75
        IF ( STRVIB .EQ. '       1 1 1 22' ) IDXVIB = 76
        IF ( STRVIB .EQ. '       0 3 3 21' ) IDXVIB = 77
        IF ( STRVIB .EQ. '       1 1 1 21' ) IDXVIB = 78
        IF ( STRVIB .EQ. '       0 0 0 31' ) IDXVIB = 79
        IF ( STRVIB .EQ. '       3 1 1 14' ) IDXVIB = 80
        IF ( STRVIB .EQ. '       2 3 3 13' ) IDXVIB = 81
        IF ( STRVIB .EQ. '       3 1 1 13' ) IDXVIB = 82
        IF ( STRVIB .EQ. '       2 3 3 12' ) IDXVIB = 83
        IF ( STRVIB .EQ. '       3 1 1 12' ) IDXVIB = 84
        IF ( STRVIB .EQ. '       1 5 5 11' ) IDXVIB = 85
        IF ( STRVIB .EQ. '       2 0 0 23' ) IDXVIB = 86
        IF ( STRVIB .EQ. '       2 3 3 11' ) IDXVIB = 87
        IF ( STRVIB .EQ. '       1 2 2 22' ) IDXVIB = 88
        IF ( STRVIB .EQ. '       3 1 1 11' ) IDXVIB = 89
        IF ( STRVIB .EQ. '       2 0 0 22' ) IDXVIB = 90
        IF ( STRVIB .EQ. '       1 2 2 21' ) IDXVIB = 91
        IF ( STRVIB .EQ. '       2 0 0 21' ) IDXVIB = 92
        IF ( STRVIB .EQ. '       0 1 1 31' ) IDXVIB = 93
        IF ( STRVIB .EQ. '       4 0 0 15' ) IDXVIB = 94
        IF ( STRVIB .EQ. '       3 2 2 14' ) IDXVIB = 95
        IF ( STRVIB .EQ. '       4 0 0 14' ) IDXVIB = 96
        IF ( STRVIB .EQ. '       3 2 2 13' ) IDXVIB = 97
        IF ( STRVIB .EQ. '       4 0 0 13' ) IDXVIB = 98
        IF ( STRVIB .EQ. '       5 1 1 02' ) IDXVIB = 99
        IF ( STRVIB .EQ. '       3 2 2 12' ) IDXVIB = 100
        IF ( STRVIB .EQ. '       4 0 0 12' ) IDXVIB = 101
        IF ( STRVIB .EQ. '       2 1 1 23' ) IDXVIB = 102
        IF ( STRVIB .EQ. '       3 2 2 11' ) IDXVIB = 103
        IF ( STRVIB .EQ. '       2 1 1 22' ) IDXVIB = 104
        IF ( STRVIB .EQ. '       4 0 0 11' ) IDXVIB = 105
        IF ( STRVIB .EQ. '       2 1 1 21' ) IDXVIB = 106
        IF ( STRVIB .EQ. '       1 0 0 32' ) IDXVIB = 107
        IF ( STRVIB .EQ. '       0 2 2 31' ) IDXVIB = 108
        IF ( STRVIB .EQ. '       1 0 0 31' ) IDXVIB = 109
        IF ( STRVIB .EQ. '       4 1 1 14' ) IDXVIB = 110
        IF ( STRVIB .EQ. '       4 1 1 13' ) IDXVIB = 111
        IF ( STRVIB .EQ. '       4 1 1 12' ) IDXVIB = 112
        IF ( STRVIB .EQ. '       1 1 1 32' ) IDXVIB = 113
        IF ( STRVIB .EQ. '       0 3 3 31' ) IDXVIB = 114
        IF ( STRVIB .EQ. '       1 1 1 31' ) IDXVIB = 115
        IF ( STRVIB .EQ. '       2 0 0 33' ) IDXVIB = 116
        IF ( STRVIB .EQ. '       1 2 2 32' ) IDXVIB = 117
        IF ( STRVIB .EQ. '       2 0 0 32' ) IDXVIB = 118
        IF ( STRVIB .EQ. '       1 2 2 31' ) IDXVIB = 119
        IF ( STRVIB .EQ. '       2 0 0 31' ) IDXVIB = 120
        IF ( STRVIB .EQ. '       2 1 1 33' ) IDXVIB = 121
        IF ( STRVIB .EQ. '       2 1 1 32' ) IDXVIB = 122
        IF ( STRVIB .EQ. '       2 1 1 31' ) IDXVIB = 123
        IF ( STRVIB .EQ. '       2 3 3 03' ) IDXVIB = 124
        IF ( STRVIB .EQ. '       1 5 5 02' ) IDXVIB = 125
        IF ( STRVIB .EQ. '       2 3 3 02' ) IDXVIB = 126
        IF ( STRVIB .EQ. '       0 7 7 01' ) IDXVIB = 127
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 128
        IF ( STRVIB .EQ. '       1 0 0 41' ) IDXVIB = 129
        IF ( STRVIB .EQ. '       1 0 0 51' ) IDXVIB = 130
        IF ( STRVIB .EQ. '       1 0 0 52' ) IDXVIB = 131
        IF ( STRVIB .EQ. '       0 0 0 51' ) IDXVIB = 132
C new levels found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '       0 0 0 00' ) IDXVIB = 133
        IF ( STRVIB .EQ. '       0 0 0 41' ) IDXVIB = 134
        IF ( STRVIB .EQ. '       0 1 1 41' ) IDXVIB = 135
        IF ( STRVIB .EQ. '       0 1 1 51' ) IDXVIB = 136
        IF ( STRVIB .EQ. '       0 2 2 41' ) IDXVIB = 137
        IF ( STRVIB .EQ. '       0 2 2 51' ) IDXVIB = 138
        IF ( STRVIB .EQ. '       0 4 4 21' ) IDXVIB = 139
        IF ( STRVIB .EQ. '       0 4 4 31' ) IDXVIB = 140
        IF ( STRVIB .EQ. '       0 5 5 21' ) IDXVIB = 141
        IF ( STRVIB .EQ. '       0 5 5 31' ) IDXVIB = 142
        IF ( STRVIB .EQ. '       0 7 7 11' ) IDXVIB = 143
        IF ( STRVIB .EQ. '       0 8 8 01' ) IDXVIB = 144
        IF ( STRVIB .EQ. '       0 8 8 11' ) IDXVIB = 145
        IF ( STRVIB .EQ. '       0 9 9 01' ) IDXVIB = 146
        IF ( STRVIB .EQ. '       1 0 0 42' ) IDXVIB = 147
        IF ( STRVIB .EQ. '       1 1 1 41' ) IDXVIB = 148
        IF ( STRVIB .EQ. '       1 1 1 42' ) IDXVIB = 149
        IF ( STRVIB .EQ. '       1 1 1 51' ) IDXVIB = 150
        IF ( STRVIB .EQ. '       1 1 1 52' ) IDXVIB = 151
        IF ( STRVIB .EQ. '       1 2 2 51' ) IDXVIB = 152
        IF ( STRVIB .EQ. '       1 2 2 52' ) IDXVIB = 153
        IF ( STRVIB .EQ. '       1 3 3 21' ) IDXVIB = 154
        IF ( STRVIB .EQ. '       1 3 3 22' ) IDXVIB = 155
        IF ( STRVIB .EQ. '       1 3 3 31' ) IDXVIB = 156
        IF ( STRVIB .EQ. '       1 3 3 32' ) IDXVIB = 157
        IF ( STRVIB .EQ. '       1 4 4 21' ) IDXVIB = 158
        IF ( STRVIB .EQ. '       1 4 4 22' ) IDXVIB = 159
        IF ( STRVIB .EQ. '       1 4 4 31' ) IDXVIB = 160
        IF ( STRVIB .EQ. '       1 4 4 32' ) IDXVIB = 161
        IF ( STRVIB .EQ. '       1 5 5 01' ) IDXVIB = 162
        IF ( STRVIB .EQ. '       1 5 5 12' ) IDXVIB = 163
        IF ( STRVIB .EQ. '       1 6 6 01' ) IDXVIB = 164
        IF ( STRVIB .EQ. '       1 6 6 02' ) IDXVIB = 165
        IF ( STRVIB .EQ. '       1 6 6 11' ) IDXVIB = 166
        IF ( STRVIB .EQ. '       1 6 6 12' ) IDXVIB = 167
        IF ( STRVIB .EQ. '       1 7 7 01' ) IDXVIB = 168
        IF ( STRVIB .EQ. '       1 7 7 02' ) IDXVIB = 169
        IF ( STRVIB .EQ. '       1 7 7 11' ) IDXVIB = 170
        IF ( STRVIB .EQ. '       1 7 7 12' ) IDXVIB = 171
        IF ( STRVIB .EQ. '       2 0 0 41' ) IDXVIB = 172
        IF ( STRVIB .EQ. '       2 0 0 42' ) IDXVIB = 173
        IF ( STRVIB .EQ. '       2 0 0 51' ) IDXVIB = 174
        IF ( STRVIB .EQ. '       2 0 0 52' ) IDXVIB = 175
        IF ( STRVIB .EQ. '       2 0 0 53' ) IDXVIB = 176
        IF ( STRVIB .EQ. '       2 2 2 21' ) IDXVIB = 177
        IF ( STRVIB .EQ. '       2 2 2 22' ) IDXVIB = 178
        IF ( STRVIB .EQ. '       2 2 2 23' ) IDXVIB = 179
        IF ( STRVIB .EQ. '       2 2 2 31' ) IDXVIB = 180
        IF ( STRVIB .EQ. '       2 2 2 32' ) IDXVIB = 181
        IF ( STRVIB .EQ. '       2 2 2 33' ) IDXVIB = 182
        IF ( STRVIB .EQ. '       2 3 3 21' ) IDXVIB = 183
        IF ( STRVIB .EQ. '       2 3 3 22' ) IDXVIB = 184
        IF ( STRVIB .EQ. '       2 3 3 23' ) IDXVIB = 185
        IF ( STRVIB .EQ. '       2 3 3 31' ) IDXVIB = 186
        IF ( STRVIB .EQ. '       2 3 3 32' ) IDXVIB = 187
        IF ( STRVIB .EQ. '       2 3 3 33' ) IDXVIB = 188
        IF ( STRVIB .EQ. '       2 4 4 01' ) IDXVIB = 189
        IF ( STRVIB .EQ. '       2 4 4 02' ) IDXVIB = 190
        IF ( STRVIB .EQ. '       2 4 4 03' ) IDXVIB = 191
        IF ( STRVIB .EQ. '       2 4 4 11' ) IDXVIB = 192
        IF ( STRVIB .EQ. '       2 4 4 12' ) IDXVIB = 193
        IF ( STRVIB .EQ. '       2 4 4 13' ) IDXVIB = 194
        IF ( STRVIB .EQ. '       2 5 5 01' ) IDXVIB = 195
        IF ( STRVIB .EQ. '       2 5 5 02' ) IDXVIB = 196
        IF ( STRVIB .EQ. '       2 5 5 03' ) IDXVIB = 197
        IF ( STRVIB .EQ. '       2 5 5 11' ) IDXVIB = 198
        IF ( STRVIB .EQ. '       2 5 5 12' ) IDXVIB = 199
        IF ( STRVIB .EQ. '       2 5 5 13' ) IDXVIB = 200
        IF ( STRVIB .EQ. '       2 6 6 01' ) IDXVIB = 201
        IF ( STRVIB .EQ. '       2 6 6 02' ) IDXVIB = 202
        IF ( STRVIB .EQ. '       2 6 6 11' ) IDXVIB = 203
        IF ( STRVIB .EQ. '       2 6 6 12' ) IDXVIB = 204
        IF ( STRVIB .EQ. '       2 6 6 13' ) IDXVIB = 205
        IF ( STRVIB .EQ. '       3 0 0 21' ) IDXVIB = 206
        IF ( STRVIB .EQ. '       3 0 0 22' ) IDXVIB = 207
        IF ( STRVIB .EQ. '       3 0 0 23' ) IDXVIB = 208
        IF ( STRVIB .EQ. '       3 0 0 24' ) IDXVIB = 209
        IF ( STRVIB .EQ. '       3 0 0 31' ) IDXVIB = 210
        IF ( STRVIB .EQ. '       3 0 0 32' ) IDXVIB = 211
        IF ( STRVIB .EQ. '       3 0 0 33' ) IDXVIB = 212
        IF ( STRVIB .EQ. '       3 0 0 34' ) IDXVIB = 213
        IF ( STRVIB .EQ. '       3 1 1 21' ) IDXVIB = 214
        IF ( STRVIB .EQ. '       3 1 1 22' ) IDXVIB = 215
        IF ( STRVIB .EQ. '       3 1 1 23' ) IDXVIB = 216
        IF ( STRVIB .EQ. '       3 1 1 24' ) IDXVIB = 217
        IF ( STRVIB .EQ. '       3 1 1 31' ) IDXVIB = 218
        IF ( STRVIB .EQ. '       3 1 1 32' ) IDXVIB = 219
        IF ( STRVIB .EQ. '       3 1 1 33' ) IDXVIB = 220
        IF ( STRVIB .EQ. '       3 1 1 34' ) IDXVIB = 221
        IF ( STRVIB .EQ. '       3 2 2 01' ) IDXVIB = 222
        IF ( STRVIB .EQ. '       3 2 2 02' ) IDXVIB = 223
        IF ( STRVIB .EQ. '       3 2 2 04' ) IDXVIB = 224
        IF ( STRVIB .EQ. '       3 3 3 01' ) IDXVIB = 225
        IF ( STRVIB .EQ. '       3 3 3 02' ) IDXVIB = 226
        IF ( STRVIB .EQ. '       3 3 3 03' ) IDXVIB = 227
        IF ( STRVIB .EQ. '       3 3 3 04' ) IDXVIB = 228
        IF ( STRVIB .EQ. '       3 3 3 11' ) IDXVIB = 229
        IF ( STRVIB .EQ. '       3 3 3 12' ) IDXVIB = 230
        IF ( STRVIB .EQ. '       3 3 3 13' ) IDXVIB = 231
        IF ( STRVIB .EQ. '       3 3 3 14' ) IDXVIB = 232
        IF ( STRVIB .EQ. '       3 4 4 01' ) IDXVIB = 233
        IF ( STRVIB .EQ. '       3 4 4 02' ) IDXVIB = 234
        IF ( STRVIB .EQ. '       3 4 4 03' ) IDXVIB = 235
        IF ( STRVIB .EQ. '       3 4 4 04' ) IDXVIB = 236
        IF ( STRVIB .EQ. '       3 4 4 11' ) IDXVIB = 237
        IF ( STRVIB .EQ. '       3 4 4 12' ) IDXVIB = 238
        IF ( STRVIB .EQ. '       3 4 4 13' ) IDXVIB = 239
        IF ( STRVIB .EQ. '       3 4 4 14' ) IDXVIB = 240
        IF ( STRVIB .EQ. '       3 5 5 01' ) IDXVIB = 241
        IF ( STRVIB .EQ. '       3 5 5 02' ) IDXVIB = 242
        IF ( STRVIB .EQ. '       3 5 5 12' ) IDXVIB = 243
        IF ( STRVIB .EQ. '       3 5 5 13' ) IDXVIB = 244
        IF ( STRVIB .EQ. '       4 0 0 01' ) IDXVIB = 245
        IF ( STRVIB .EQ. '       4 0 0 03' ) IDXVIB = 246
        IF ( STRVIB .EQ. '       4 0 0 05' ) IDXVIB = 247
        IF ( STRVIB .EQ. '       4 0 0 21' ) IDXVIB = 248
        IF ( STRVIB .EQ. '       4 0 0 22' ) IDXVIB = 249
        IF ( STRVIB .EQ. '       4 0 0 23' ) IDXVIB = 250
        IF ( STRVIB .EQ. '       4 0 0 24' ) IDXVIB = 251
        IF ( STRVIB .EQ. '       4 0 0 33' ) IDXVIB = 252
        IF ( STRVIB .EQ. '       4 0 0 34' ) IDXVIB = 253
        IF ( STRVIB .EQ. '       4 1 1 03' ) IDXVIB = 254
        IF ( STRVIB .EQ. '       4 1 1 04' ) IDXVIB = 255
        IF ( STRVIB .EQ. '       4 1 1 05' ) IDXVIB = 256
        IF ( STRVIB .EQ. '       4 1 1 11' ) IDXVIB = 257
        IF ( STRVIB .EQ. '       4 1 1 15' ) IDXVIB = 258
        IF ( STRVIB .EQ. '       4 1 1 25' ) IDXVIB = 259
        IF ( STRVIB .EQ. '       4 2 2 01' ) IDXVIB = 260
        IF ( STRVIB .EQ. '       4 2 2 02' ) IDXVIB = 261
        IF ( STRVIB .EQ. '       4 2 2 03' ) IDXVIB = 262
        IF ( STRVIB .EQ. '       4 2 2 04' ) IDXVIB = 263
        IF ( STRVIB .EQ. '       4 2 2 05' ) IDXVIB = 264
        IF ( STRVIB .EQ. '       4 2 2 11' ) IDXVIB = 265
        IF ( STRVIB .EQ. '       4 2 2 12' ) IDXVIB = 266
        IF ( STRVIB .EQ. '       4 2 2 13' ) IDXVIB = 267
        IF ( STRVIB .EQ. '       4 2 2 14' ) IDXVIB = 268
        IF ( STRVIB .EQ. '       4 2 2 15' ) IDXVIB = 269
        IF ( STRVIB .EQ. '       4 3 3 01' ) IDXVIB = 270
        IF ( STRVIB .EQ. '       4 3 3 02' ) IDXVIB = 271
        IF ( STRVIB .EQ. '       4 3 3 03' ) IDXVIB = 272
        IF ( STRVIB .EQ. '       4 3 3 04' ) IDXVIB = 273
        IF ( STRVIB .EQ. '       4 3 3 11' ) IDXVIB = 274
        IF ( STRVIB .EQ. '       4 3 3 12' ) IDXVIB = 275
        IF ( STRVIB .EQ. '       4 3 3 13' ) IDXVIB = 276
        IF ( STRVIB .EQ. '       4 3 3 14' ) IDXVIB = 277
        IF ( STRVIB .EQ. '       4 3 3 15' ) IDXVIB = 278
        IF ( STRVIB .EQ. '       4 4 4 01' ) IDXVIB = 279
        IF ( STRVIB .EQ. '       4 4 4 02' ) IDXVIB = 280
        IF ( STRVIB .EQ. '       4 4 4 03' ) IDXVIB = 281
        IF ( STRVIB .EQ. '       5 0 0 01' ) IDXVIB = 282
        IF ( STRVIB .EQ. '       5 0 0 02' ) IDXVIB = 283
        IF ( STRVIB .EQ. '       5 0 0 03' ) IDXVIB = 284
        IF ( STRVIB .EQ. '       5 0 0 04' ) IDXVIB = 285
        IF ( STRVIB .EQ. '       5 0 0 05' ) IDXVIB = 286
        IF ( STRVIB .EQ. '       5 0 0 06' ) IDXVIB = 287
        IF ( STRVIB .EQ. '       5 0 0 11' ) IDXVIB = 288
        IF ( STRVIB .EQ. '       5 0 0 12' ) IDXVIB = 289
        IF ( STRVIB .EQ. '       5 0 0 13' ) IDXVIB = 290
        IF ( STRVIB .EQ. '       5 0 0 14' ) IDXVIB = 291
        IF ( STRVIB .EQ. '       5 0 0 15' ) IDXVIB = 292
        IF ( STRVIB .EQ. '       5 0 0 16' ) IDXVIB = 293
        IF ( STRVIB .EQ. '       5 1 1 01' ) IDXVIB = 294
        IF ( STRVIB .EQ. '       5 1 1 03' ) IDXVIB = 295
        IF ( STRVIB .EQ. '       5 1 1 04' ) IDXVIB = 296
        IF ( STRVIB .EQ. '       5 1 1 05' ) IDXVIB = 297
        IF ( STRVIB .EQ. '       5 1 1 06' ) IDXVIB = 298
        IF ( STRVIB .EQ. '       5 1 1 11' ) IDXVIB = 299
        IF ( STRVIB .EQ. '       5 1 1 12' ) IDXVIB = 300
        IF ( STRVIB .EQ. '       5 1 1 13' ) IDXVIB = 301
        IF ( STRVIB .EQ. '       5 1 1 14' ) IDXVIB = 302
        IF ( STRVIB .EQ. '       5 1 1 15' ) IDXVIB = 303
        IF ( STRVIB .EQ. '       5 1 1 16' ) IDXVIB = 304
        IF ( STRVIB .EQ. '       5 2 2 01' ) IDXVIB = 305
        IF ( STRVIB .EQ. '       5 2 2 02' ) IDXVIB = 306
        IF ( STRVIB .EQ. '       5 2 2 03' ) IDXVIB = 307
        IF ( STRVIB .EQ. '       5 2 2 04' ) IDXVIB = 308
        IF ( STRVIB .EQ. '       5 2 2 05' ) IDXVIB = 309
        IF ( STRVIB .EQ. '       5 2 2 11' ) IDXVIB = 310
        IF ( STRVIB .EQ. '       5 2 2 12' ) IDXVIB = 311
        IF ( STRVIB .EQ. '       5 2 2 13' ) IDXVIB = 312
        IF ( STRVIB .EQ. '       5 2 2 14' ) IDXVIB = 313
        IF ( STRVIB .EQ. '       5 2 2 15' ) IDXVIB = 314
        IF ( STRVIB .EQ. '       5 3 3 02' ) IDXVIB = 315
        IF ( STRVIB .EQ. '       5 4 4 03' ) IDXVIB = 316
        IF ( STRVIB .EQ. '       6 0 0 01' ) IDXVIB = 317
        IF ( STRVIB .EQ. '       6 0 0 02' ) IDXVIB = 318
        IF ( STRVIB .EQ. '       6 0 0 03' ) IDXVIB = 319
        IF ( STRVIB .EQ. '       6 0 0 04' ) IDXVIB = 320
        IF ( STRVIB .EQ. '       6 0 0 06' ) IDXVIB = 321
        IF ( STRVIB .EQ. '       6 0 0 12' ) IDXVIB = 322
        IF ( STRVIB .EQ. '       6 0 0 13' ) IDXVIB = 323
        IF ( STRVIB .EQ. '       6 0 0 14' ) IDXVIB = 324
        IF ( STRVIB .EQ. '       6 0 0 15' ) IDXVIB = 325
        IF ( STRVIB .EQ. '       6 0 0 16' ) IDXVIB = 326
        IF ( STRVIB .EQ. '       6 0 0 17' ) IDXVIB = 327
        IF ( STRVIB .EQ. '       6 1 1 01' ) IDXVIB = 328
        IF ( STRVIB .EQ. '       6 1 1 02' ) IDXVIB = 329
        IF ( STRVIB .EQ. '       6 1 1 03' ) IDXVIB = 330
        IF ( STRVIB .EQ. '       6 1 1 04' ) IDXVIB = 331
        IF ( STRVIB .EQ. '       6 1 1 14' ) IDXVIB = 332
        IF ( STRVIB .EQ. '       6 1 1 15' ) IDXVIB = 333
C
C Class 6: Non-linear triatomic molecules
      ELSE IF ( IDXMOL .EQ.  1 .OR.                                ! H2O
     &          IDXMOL .EQ.  3 .OR.                                  ! O3
     &          IDXMOL .EQ.  9 .OR.                                  ! SO2
     &          IDXMOL .EQ. 10 .OR.                                  ! NO2
     &          IDXMOL .EQ. 21 .OR.                                  ! HOCl
     &          IDXMOL .EQ. 31 .OR.                                  ! H2S
     &          IDXMOL .EQ. 33 .OR.                                  ! HO2
     &          IDXMOL .EQ. 37      ) THEN                         ! HOBr
        IF ( STRVIB .EQ. '          0 0 0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '          0 1 0' ) IDXVIB = 2
        IF ( STRVIB .EQ. '          0 2 0' ) IDXVIB = 3
        IF ( STRVIB .EQ. '          1 0 0' ) IDXVIB = 4
        IF ( STRVIB .EQ. '          0 0 1' ) IDXVIB = 5
        IF ( STRVIB .EQ. '          0 3 0' ) IDXVIB = 6
        IF ( STRVIB .EQ. '          1 1 0' ) IDXVIB = 7
        IF ( STRVIB .EQ. '          0 1 1' ) IDXVIB = 8
        IF ( STRVIB .EQ. '          0 4 0' ) IDXVIB = 9
        IF ( STRVIB .EQ. '          1 2 0' ) IDXVIB = 10
        IF ( STRVIB .EQ. '          0 2 1' ) IDXVIB = 11
        IF ( STRVIB .EQ. '          2 0 0' ) IDXVIB = 12
        IF ( STRVIB .EQ. '          1 0 1' ) IDXVIB = 13
        IF ( STRVIB .EQ. '          0 0 2' ) IDXVIB = 14
        IF ( STRVIB .EQ. '          1 3 0' ) IDXVIB = 15
        IF ( STRVIB .EQ. '          0 3 1' ) IDXVIB = 16
        IF ( STRVIB .EQ. '          2 1 0' ) IDXVIB = 17
        IF ( STRVIB .EQ. '          1 1 1' ) IDXVIB = 18
        IF ( STRVIB .EQ. '          0 1 2' ) IDXVIB = 19
        IF ( STRVIB .EQ. '          0 4 1' ) IDXVIB = 20
        IF ( STRVIB .EQ. '          2 2 0' ) IDXVIB = 21
        IF ( STRVIB .EQ. '          1 2 1' ) IDXVIB = 22
        IF ( STRVIB .EQ. '          0 2 2' ) IDXVIB = 23
        IF ( STRVIB .EQ. '          3 0 0' ) IDXVIB = 24
        IF ( STRVIB .EQ. '          2 0 1' ) IDXVIB = 25
        IF ( STRVIB .EQ. '          1 0 2' ) IDXVIB = 26
        IF ( STRVIB .EQ. '          0 0 3' ) IDXVIB = 27
        IF ( STRVIB .EQ. '          1 3 1' ) IDXVIB = 28
        IF ( STRVIB .EQ. '          3 1 0' ) IDXVIB = 29
        IF ( STRVIB .EQ. '          2 1 1' ) IDXVIB = 30
        IF ( STRVIB .EQ. '          1 1 2' ) IDXVIB = 31
        IF ( STRVIB .EQ. '          0 1 3' ) IDXVIB = 32
        IF ( STRVIB .EQ. '          1 4 1' ) IDXVIB = 33
        IF ( STRVIB .EQ. '          0 4 2' ) IDXVIB = 34
        IF ( STRVIB .EQ. '          3 2 0' ) IDXVIB = 35
        IF ( STRVIB .EQ. '          2 2 1' ) IDXVIB = 36
        IF ( STRVIB .EQ. '          3 0 1' ) IDXVIB = 37
        IF ( STRVIB .EQ. '          2 0 2' ) IDXVIB = 38
        IF ( STRVIB .EQ. '          1 2 2' ) IDXVIB = 39
        IF ( STRVIB .EQ. '          0 2 3' ) IDXVIB = 40
        IF ( STRVIB .EQ. '          4 0 0' ) IDXVIB = 41
        IF ( STRVIB .EQ. '          1 0 3' ) IDXVIB = 42
        IF ( STRVIB .EQ. '          0 0 4' ) IDXVIB = 43
        IF ( STRVIB .EQ. '          1 5 1' ) IDXVIB = 44
        IF ( STRVIB .EQ. '          3 3 0' ) IDXVIB = 45
        IF ( STRVIB .EQ. '          2 3 1' ) IDXVIB = 46
        IF ( STRVIB .EQ. '          2 1 2' ) IDXVIB = 47
        IF ( STRVIB .EQ. '          3 1 1' ) IDXVIB = 48
        IF ( STRVIB .EQ. '          4 1 0' ) IDXVIB = 49
        IF ( STRVIB .EQ. '          1 1 3' ) IDXVIB = 50
        IF ( STRVIB .EQ. '          3 2 1' ) IDXVIB = 51
        IF ( STRVIB .EQ. '          2 2 2' ) IDXVIB = 52
        IF ( STRVIB .EQ. '          3 0 2' ) IDXVIB = 53
        IF ( STRVIB .EQ. '          4 0 1' ) IDXVIB = 54
        IF ( STRVIB .EQ. '          4 2 0' ) IDXVIB = 55
        IF ( STRVIB .EQ. '          1 2 3' ) IDXVIB = 56
        IF ( STRVIB .EQ. '          5 0 0' ) IDXVIB = 57
        IF ( STRVIB .EQ. '          2 0 3' ) IDXVIB = 58
        IF ( STRVIB .EQ. '          1 0 4' ) IDXVIB = 59
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 60
        IF ( STRVIB .EQ. '          3 3 1' ) IDXVIB = 61
        IF ( STRVIB .EQ. '          2 1 3' ) IDXVIB = 62
        IF ( STRVIB .EQ. '          3 1 2' ) IDXVIB = 63
        IF ( STRVIB .EQ. '          4 1 1' ) IDXVIB = 64
        IF ( STRVIB .EQ. '          3 0 3' ) IDXVIB = 65
        IF ( STRVIB .EQ. '          4 0 2' ) IDXVIB = 66
        IF ( STRVIB .EQ. '          4 0 3' ) IDXVIB = 67
        IF ( STRVIB .EQ. '          4 2 1' ) IDXVIB = 68
        IF ( STRVIB .EQ. '          5 0 1' ) IDXVIB = 69
        IF ( STRVIB .EQ. '          3 1 3' ) IDXVIB = 70
        IF ( STRVIB .EQ. '          4 1 2' ) IDXVIB = 71
        IF ( STRVIB .EQ. '          2 3 2' ) IDXVIB = 72
        IF ( STRVIB .EQ. '          0 5 0' ) IDXVIB = 73
        IF ( STRVIB .EQ. '          0 6 0' ) IDXVIB = 74
        IF ( STRVIB .EQ. '          0 7 0' ) IDXVIB = 75
        IF ( STRVIB .EQ. '          0 3 2' ) IDXVIB = 76
        IF ( STRVIB .EQ. '          0 5 1' ) IDXVIB = 77
        IF ( STRVIB .EQ. '          0 6 1' ) IDXVIB = 78
        IF ( STRVIB .EQ. '          0 8 0' ) IDXVIB = 79
        IF ( STRVIB .EQ. '          1 4 0' ) IDXVIB = 80
        IF ( STRVIB .EQ. '          1 5 0' ) IDXVIB = 81
        IF ( STRVIB .EQ. '          0 3 3' ) IDXVIB = 82
        IF ( STRVIB .EQ. '          0 3 4' ) IDXVIB = 83
        IF ( STRVIB .EQ. '          0 4 3' ) IDXVIB = 84
        IF ( STRVIB .EQ. '          0 5 3' ) IDXVIB = 85
        IF ( STRVIB .EQ. '          0 6 1' ) IDXVIB = 86
        IF ( STRVIB .EQ. '          0 6 3' ) IDXVIB = 87
        IF ( STRVIB .EQ. '          0 7 1' ) IDXVIB = 88
        IF ( STRVIB .EQ. '          1 1 5' ) IDXVIB = 89
        IF ( STRVIB .EQ. '          1 3 2' ) IDXVIB = 90
        IF ( STRVIB .EQ. '          1 3 3' ) IDXVIB = 91
        IF ( STRVIB .EQ. '          1 4 2' ) IDXVIB = 92
        IF ( STRVIB .EQ. '          1 6 0' ) IDXVIB = 93
        IF ( STRVIB .EQ. '          1 7 0' ) IDXVIB = 94
        IF ( STRVIB .EQ. '          2 2 3' ) IDXVIB = 95
        IF ( STRVIB .EQ. '          2 4 0' ) IDXVIB = 96
        IF ( STRVIB .EQ. '          2 4 1' ) IDXVIB = 97
        IF ( STRVIB .EQ. '          3 2 2' ) IDXVIB = 98
        IF ( STRVIB .EQ. '          3 4 0' ) IDXVIB = 99
        IF ( STRVIB .EQ. '          3 4 1' ) IDXVIB = 100
        IF ( STRVIB .EQ. '          4 3 0' ) IDXVIB = 101
        IF ( STRVIB .EQ. '          4 3 1' ) IDXVIB = 102
        IF ( STRVIB .EQ. '          5 1 0' ) IDXVIB = 103
        IF ( STRVIB .EQ. '          5 1 1' ) IDXVIB = 104
        IF ( STRVIB .EQ. '          5 2 0' ) IDXVIB = 105
        IF ( STRVIB .EQ. '          6 0 0' ) IDXVIB = 106
        IF ( STRVIB .EQ. '          6 0 1' ) IDXVIB = 107
        IF ( STRVIB .EQ. '          6 1 0' ) IDXVIB = 108
        IF ( STRVIB .EQ. '          6 1 1' ) IDXVIB = 109
        IF ( STRVIB .EQ. '          6 2 0' ) IDXVIB = 110
        IF ( STRVIB .EQ. '          7 0 0' ) IDXVIB = 111
        IF ( STRVIB .EQ. '          7 0 1' ) IDXVIB = 112
        IF ( STRVIB .EQ. '          8 0 0' ) IDXVIB = 113
        IF ( STRVIB .EQ. '          0 5 0' ) IDXVIB = 114
        IF ( STRVIB .EQ. '          0 6 0' ) IDXVIB = 115
C new levels for H2O found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '          0 0 5' ) IDXVIB = 116
        IF ( STRVIB .EQ. '          0 0 6' ) IDXVIB = 117
        IF ( STRVIB .EQ. '          0 0 7' ) IDXVIB = 118
        IF ( STRVIB .EQ. '          0 1 4' ) IDXVIB = 119
        IF ( STRVIB .EQ. '          0 1 5' ) IDXVIB = 120
        IF ( STRVIB .EQ. '          0 2 4' ) IDXVIB = 121
        IF ( STRVIB .EQ. '          0 2 5' ) IDXVIB = 122
        IF ( STRVIB .EQ. '          0 3 5' ) IDXVIB = 123
        IF ( STRVIB .EQ. '          0 4 4' ) IDXVIB = 124
        IF ( STRVIB .EQ. '          0 4 5' ) IDXVIB = 125
        IF ( STRVIB .EQ. '          0 5 2' ) IDXVIB = 126
        IF ( STRVIB .EQ. '          0 6 2' ) IDXVIB = 127
        IF ( STRVIB .EQ. '          0 7 2' ) IDXVIB = 128
        IF ( STRVIB .EQ. '          0 7 3' ) IDXVIB = 129
        IF ( STRVIB .EQ. '          0 8 1' ) IDXVIB = 130
        IF ( STRVIB .EQ. '          0 8 2' ) IDXVIB = 131
        IF ( STRVIB .EQ. '          0 9 0' ) IDXVIB = 132
        IF ( STRVIB .EQ. '          0 9 1' ) IDXVIB = 133
        IF ( STRVIB .EQ. '          0 9 2' ) IDXVIB = 134
        IF ( STRVIB .EQ. '          0 9 3' ) IDXVIB = 135
        IF ( STRVIB .EQ. '          010 0' ) IDXVIB = 136
        IF ( STRVIB .EQ. '          010 1' ) IDXVIB = 137
        IF ( STRVIB .EQ. '          011 0' ) IDXVIB = 138
        IF ( STRVIB .EQ. '          011 1' ) IDXVIB = 139
        IF ( STRVIB .EQ. '          012 0' ) IDXVIB = 140
        IF ( STRVIB .EQ. '          012 1' ) IDXVIB = 141
        IF ( STRVIB .EQ. '          013 0' ) IDXVIB = 142
        IF ( STRVIB .EQ. '          014 0' ) IDXVIB = 143
        IF ( STRVIB .EQ. '          015 0' ) IDXVIB = 144
        IF ( STRVIB .EQ. '          1 0 5' ) IDXVIB = 145
        IF ( STRVIB .EQ. '          1 0 6' ) IDXVIB = 146
        IF ( STRVIB .EQ. '          1 1 4' ) IDXVIB = 147
        IF ( STRVIB .EQ. '          1 2 4' ) IDXVIB = 148
        IF ( STRVIB .EQ. '          1 2 5' ) IDXVIB = 149
        IF ( STRVIB .EQ. '          1 3 4' ) IDXVIB = 150
        IF ( STRVIB .EQ. '          1 4 3' ) IDXVIB = 151
        IF ( STRVIB .EQ. '          1 5 2' ) IDXVIB = 152
        IF ( STRVIB .EQ. '          1 5 3' ) IDXVIB = 153
        IF ( STRVIB .EQ. '          1 5 4' ) IDXVIB = 154
        IF ( STRVIB .EQ. '          1 6 1' ) IDXVIB = 155
        IF ( STRVIB .EQ. '          1 6 2' ) IDXVIB = 156
        IF ( STRVIB .EQ. '          1 6 3' ) IDXVIB = 157
        IF ( STRVIB .EQ. '          1 7 1' ) IDXVIB = 158
        IF ( STRVIB .EQ. '          1 8 0' ) IDXVIB = 159
        IF ( STRVIB .EQ. '          1 8 1' ) IDXVIB = 160
        IF ( STRVIB .EQ. '          1 8 2' ) IDXVIB = 161
        IF ( STRVIB .EQ. '          1 9 0' ) IDXVIB = 162
        IF ( STRVIB .EQ. '          1 9 1' ) IDXVIB = 163
        IF ( STRVIB .EQ. '          110 0' ) IDXVIB = 164
        IF ( STRVIB .EQ. '          110 1' ) IDXVIB = 165
        IF ( STRVIB .EQ. '          111 0' ) IDXVIB = 166
        IF ( STRVIB .EQ. '          111 1' ) IDXVIB = 167
        IF ( STRVIB .EQ. '          114 0' ) IDXVIB = 168
        IF ( STRVIB .EQ. '          2 0 4' ) IDXVIB = 169
        IF ( STRVIB .EQ. '          2 1 4' ) IDXVIB = 170
        IF ( STRVIB .EQ. '          2 2 4' ) IDXVIB = 171
        IF ( STRVIB .EQ. '          2 3 0' ) IDXVIB = 172
        IF ( STRVIB .EQ. '          2 3 3' ) IDXVIB = 173
        IF ( STRVIB .EQ. '          2 4 2' ) IDXVIB = 174
        IF ( STRVIB .EQ. '          2 4 3' ) IDXVIB = 175
        IF ( STRVIB .EQ. '          2 5 0' ) IDXVIB = 176
        IF ( STRVIB .EQ. '          2 5 1' ) IDXVIB = 177
        IF ( STRVIB .EQ. '          2 5 3' ) IDXVIB = 178
        IF ( STRVIB .EQ. '          2 6 0' ) IDXVIB = 179
        IF ( STRVIB .EQ. '          2 6 1' ) IDXVIB = 180
        IF ( STRVIB .EQ. '          2 7 0' ) IDXVIB = 181
        IF ( STRVIB .EQ. '          2 7 1' ) IDXVIB = 182
        IF ( STRVIB .EQ. '          2 7 2' ) IDXVIB = 183
        IF ( STRVIB .EQ. '          2 8 0' ) IDXVIB = 184
        IF ( STRVIB .EQ. '          2 8 1' ) IDXVIB = 185
        IF ( STRVIB .EQ. '          2 9 0' ) IDXVIB = 186
        IF ( STRVIB .EQ. '          210 0' ) IDXVIB = 187
        IF ( STRVIB .EQ. '          210 1' ) IDXVIB = 188
        IF ( STRVIB .EQ. '          211 0' ) IDXVIB = 189
        IF ( STRVIB .EQ. '          3 2 3' ) IDXVIB = 190
        IF ( STRVIB .EQ. '          3 3 2' ) IDXVIB = 191
        IF ( STRVIB .EQ. '          3 3 3' ) IDXVIB = 192
        IF ( STRVIB .EQ. '          3 4 2' ) IDXVIB = 193
        IF ( STRVIB .EQ. '          3 5 0' ) IDXVIB = 194
        IF ( STRVIB .EQ. '          3 5 1' ) IDXVIB = 195
        IF ( STRVIB .EQ. '          3 5 2' ) IDXVIB = 196
        IF ( STRVIB .EQ. '          3 6 0' ) IDXVIB = 197
        IF ( STRVIB .EQ. '          3 6 1' ) IDXVIB = 198
        IF ( STRVIB .EQ. '          3 7 0' ) IDXVIB = 199
        IF ( STRVIB .EQ. '          3 7 1' ) IDXVIB = 200
        IF ( STRVIB .EQ. '          3 8 0' ) IDXVIB = 201
        IF ( STRVIB .EQ. '          3 8 1' ) IDXVIB = 202
        IF ( STRVIB .EQ. '          4 2 2' ) IDXVIB = 203
        IF ( STRVIB .EQ. '          4 3 2' ) IDXVIB = 204
        IF ( STRVIB .EQ. '          4 4 0' ) IDXVIB = 205
        IF ( STRVIB .EQ. '          4 4 1' ) IDXVIB = 206
        IF ( STRVIB .EQ. '          4 5 0' ) IDXVIB = 207
        IF ( STRVIB .EQ. '          4 5 1' ) IDXVIB = 208
        IF ( STRVIB .EQ. '          4 6 0' ) IDXVIB = 209
        IF ( STRVIB .EQ. '          4 7 0' ) IDXVIB = 210
        IF ( STRVIB .EQ. '          5 0 2' ) IDXVIB = 211
        IF ( STRVIB .EQ. '          5 1 2' ) IDXVIB = 212
        IF ( STRVIB .EQ. '          5 2 1' ) IDXVIB = 213
        IF ( STRVIB .EQ. '          5 3 0' ) IDXVIB = 214
        IF ( STRVIB .EQ. '          5 3 1' ) IDXVIB = 215
        IF ( STRVIB .EQ. '          5 5 0' ) IDXVIB = 216
        IF ( STRVIB .EQ. '          6 2 1' ) IDXVIB = 217
        IF ( STRVIB .EQ. '          6 3 0' ) IDXVIB = 218
        IF ( STRVIB .EQ. '          6 4 0' ) IDXVIB = 219
        IF ( STRVIB .EQ. '          7 1 0' ) IDXVIB = 220
        IF ( STRVIB .EQ. '         -2-2-2' ) IDXVIB = 221
C new levels for O3 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '       I  1 0 5' ) IDXVIB = 222
        IF ( STRVIB .EQ. '       I  1 2 4' ) IDXVIB = 223
        IF ( STRVIB .EQ. '       I  2 0 5' ) IDXVIB = 224
        IF ( STRVIB .EQ. '       I  2 2 3' ) IDXVIB = 225
        IF ( STRVIB .EQ. '       I  2 3 3' ) IDXVIB = 226
        IF ( STRVIB .EQ. '      II  1 0 5' ) IDXVIB = 227
        IF ( STRVIB .EQ. '      II  1 2 4' ) IDXVIB = 228
        IF ( STRVIB .EQ. '      II  2 2 3' ) IDXVIB = 229
        IF ( STRVIB .EQ. '      II  2 3 3' ) IDXVIB = 230
C
C Class 7: Linear tetratomic molecules
      ELSE IF ( IDXMOL .EQ. 26 ) THEN                                ! C2H2
        IF ( STRVIB .EQ. ' 0 0 0 0 1 1   ' ) IDXVIB = 1
        IF ( STRVIB .EQ. ' 0 0 0 0 0 0+  ' ) IDXVIB = 2
        IF ( STRVIB .EQ. ' 0 0 1 0 0 0+  ' ) IDXVIB = 3
        IF ( STRVIB .EQ. ' 1 0 1 0 0 0+  ' ) IDXVIB = 4
        IF ( STRVIB .EQ. ' 0 0 0 0 1 1  u' ) IDXVIB = 5
        IF ( STRVIB .EQ. ' 0 0 0 0 0 0+ g' ) IDXVIB = 6
        IF ( STRVIB .EQ. ' 0 0 0 0 3 1  u' ) IDXVIB = 7
        IF ( STRVIB .EQ. ' 0 0 0 1 1 0+ u' ) IDXVIB = 8
        IF ( STRVIB .EQ. ' 0 0 0 2 1 1 1u' ) IDXVIB = 9
        IF ( STRVIB .EQ. ' 0 0 0 2 1 1 2u' ) IDXVIB = 10
        IF ( STRVIB .EQ. ' 0 0 1 0 0 0+ u' ) IDXVIB = 11
        IF ( STRVIB .EQ. ' 0 1 0 1 1 0+ u' ) IDXVIB = 12
        IF ( STRVIB .EQ. ' 1 0 1 0 0 0+ u' ) IDXVIB = 13
        IF ( STRVIB .EQ. ' 1 1 0 1 1 0+ u' ) IDXVIB = 14
        IF ( STRVIB .EQ. ' 0 0 0 0 2 0+ g' ) IDXVIB = 15
        IF ( STRVIB .EQ. ' 0 0 0 0 2 2  g' ) IDXVIB = 16
        IF ( STRVIB .EQ. ' 0 0 0 0 4 0+ g' ) IDXVIB = 17
        IF ( STRVIB .EQ. ' 0 0 0 0 4 2  g' ) IDXVIB = 18
        IF ( STRVIB .EQ. ' 0 0 0 2 2 0+2g' ) IDXVIB = 19
        IF ( STRVIB .EQ. ' 0 0 0 2 2 0- g' ) IDXVIB = 20
        IF ( STRVIB .EQ. ' 0 0 0 2 2 2 2g' ) IDXVIB = 21
        IF ( STRVIB .EQ. ' 0 1 0 1 0 1  g' ) IDXVIB = 22
        IF ( STRVIB .EQ. ' 1 0 1 0 1 1  g' ) IDXVIB = 23
        IF ( STRVIB .EQ. ' 0 0 0 1 0 1  g' ) IDXVIB = 24
        IF ( STRVIB .EQ. ' 0 0 0 1 1 0- u' ) IDXVIB = 25
        IF ( STRVIB .EQ. ' 0 0 0 1 1 2  u' ) IDXVIB = 26
        IF ( STRVIB .EQ. ' 0 0 0 1 3 0+ u' ) IDXVIB = 27
        IF ( STRVIB .EQ. ' 0 0 0 1 3 0- u' ) IDXVIB = 28
        IF ( STRVIB .EQ. ' 0 0 0 1 3 2 1u' ) IDXVIB = 29
        IF ( STRVIB .EQ. ' 0 0 0 1 3 2 2u' ) IDXVIB = 30
        IF ( STRVIB .EQ. ' 0 0 0 3 1 0+ u' ) IDXVIB = 31
        IF ( STRVIB .EQ. ' 0 0 0 3 1 0- u' ) IDXVIB = 32
        IF ( STRVIB .EQ. ' 0 0 0 3 1 2 1u' ) IDXVIB = 33
        IF ( STRVIB .EQ. ' 0 0 0 3 1 2 2u' ) IDXVIB = 34
        IF ( STRVIB .EQ. ' 0 1 0 0 1 1  u' ) IDXVIB = 35
        IF ( STRVIB .EQ. ' 1 0 1 1 0 1  u' ) IDXVIB = 36
C new levels for C2H2 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. ' 0 0 0 0 2 0+  ' ) IDXVIB = 37
        IF ( STRVIB .EQ. ' 0 0 0 0 2 2   ' ) IDXVIB = 38
        IF ( STRVIB .EQ. ' 0 0 0 0 3 1   ' ) IDXVIB = 39
        IF ( STRVIB .EQ. ' 0 0 0 0 3 3   ' ) IDXVIB = 40
        IF ( STRVIB .EQ. ' 0 0 0 1 0 1   ' ) IDXVIB = 41
        IF ( STRVIB .EQ. ' 0 0 0 1 1 0+  ' ) IDXVIB = 42
        IF ( STRVIB .EQ. ' 0 0 0 1 1 0-  ' ) IDXVIB = 43
        IF ( STRVIB .EQ. ' 0 0 0 1 1 2   ' ) IDXVIB = 44
        IF ( STRVIB .EQ. ' 0 0 0 1 2 1   ' ) IDXVIB = 45
        IF ( STRVIB .EQ. ' 0 0 0 1 2 1 1g' ) IDXVIB = 46
        IF ( STRVIB .EQ. ' 0 0 0 1 2 1 2g' ) IDXVIB = 47
        IF ( STRVIB .EQ. ' 0 0 0 1 2 3   ' ) IDXVIB = 48
        IF ( STRVIB .EQ. ' 0 0 0 1 3 2  u' ) IDXVIB = 49
        IF ( STRVIB .EQ. ' 0 0 0 2 0 0+  ' ) IDXVIB = 50
        IF ( STRVIB .EQ. ' 0 0 0 2 0 0+ g' ) IDXVIB = 51
        IF ( STRVIB .EQ. ' 0 0 0 2 0 2   ' ) IDXVIB = 52
        IF ( STRVIB .EQ. ' 0 0 0 2 0 2  g' ) IDXVIB = 53
        IF ( STRVIB .EQ. ' 0 0 0 2 1 1   ' ) IDXVIB = 54
        IF ( STRVIB .EQ. ' 0 0 0 2 1 3   ' ) IDXVIB = 55
        IF ( STRVIB .EQ. ' 0 0 0 2 1 3  u' ) IDXVIB = 56
        IF ( STRVIB .EQ. ' 0 0 0 2 2 0+ g' ) IDXVIB = 57
        IF ( STRVIB .EQ. ' 0 0 0 3 0 1   ' ) IDXVIB = 58
        IF ( STRVIB .EQ. ' 0 0 0 3 0 3   ' ) IDXVIB = 59
        IF ( STRVIB .EQ. ' 0 0 0 3 1 2  u' ) IDXVIB = 60
        IF ( STRVIB .EQ. ' 0 0 1 0 1 1  g' ) IDXVIB = 61
        IF ( STRVIB .EQ. ' 0 0 1 0 2 0+ u' ) IDXVIB = 62
        IF ( STRVIB .EQ. ' 0 0 1 0 2 2  u' ) IDXVIB = 63
        IF ( STRVIB .EQ. ' 0 0 1 1 0 1  u' ) IDXVIB = 64
        IF ( STRVIB .EQ. ' 0 0 1 1 1 0+ g' ) IDXVIB = 65
        IF ( STRVIB .EQ. ' 0 0 1 1 1 0- g' ) IDXVIB = 66
        IF ( STRVIB .EQ. ' 0 0 1 1 1 2  g' ) IDXVIB = 67
        IF ( STRVIB .EQ. ' 0 0 1 2 0 0+ u' ) IDXVIB = 68
        IF ( STRVIB .EQ. ' 0 0 1 2 0 2  u' ) IDXVIB = 69
        IF ( STRVIB .EQ. ' 0 0 1 3 0 1  u' ) IDXVIB = 70
        IF ( STRVIB .EQ. ' 0 0 2 0 0 0+ g' ) IDXVIB = 71
        IF ( STRVIB .EQ. ' 0 0 2 0 1 1  u' ) IDXVIB = 72
        IF ( STRVIB .EQ. ' 0 0 2 1 0 1  g' ) IDXVIB = 73
        IF ( STRVIB .EQ. ' 0 0 2 1 1 0+ u' ) IDXVIB = 74
        IF ( STRVIB .EQ. ' 0 0 2 2 0 2  g' ) IDXVIB = 75
        IF ( STRVIB .EQ. ' 0 0 3 0 0 0+ u' ) IDXVIB = 76
        IF ( STRVIB .EQ. ' 0 0 3 1 0 1  u' ) IDXVIB = 77
        IF ( STRVIB .EQ. ' 0 1 0 0 0 0+ g' ) IDXVIB = 78
        IF ( STRVIB .EQ. ' 0 1 0 0 3 1  u' ) IDXVIB = 79
        IF ( STRVIB .EQ. ' 0 1 0 1 2 1 2g' ) IDXVIB = 80
        IF ( STRVIB .EQ. ' 0 1 0 1 3 0+2u' ) IDXVIB = 81
        IF ( STRVIB .EQ. ' 0 1 0 1 3 2 2u' ) IDXVIB = 82
        IF ( STRVIB .EQ. ' 0 1 0 2 1 1 2u' ) IDXVIB = 83
        IF ( STRVIB .EQ. ' 0 1 0 2 2 0- g' ) IDXVIB = 84
        IF ( STRVIB .EQ. ' 0 1 0 2 2 2 2g' ) IDXVIB = 85
        IF ( STRVIB .EQ. ' 0 1 0 2 3 1 3u' ) IDXVIB = 86
        IF ( STRVIB .EQ. ' 0 1 0 3 1 0+ u' ) IDXVIB = 87
        IF ( STRVIB .EQ. ' 0 1 0 3 1 0+2u' ) IDXVIB = 88
        IF ( STRVIB .EQ. ' 0 1 0 3 1 2 2u' ) IDXVIB = 89
        IF ( STRVIB .EQ. ' 0 1 0 4 1 1  u' ) IDXVIB = 90
        IF ( STRVIB .EQ. ' 0 1 0 4 1 1 2u' ) IDXVIB = 91
        IF ( STRVIB .EQ. ' 0 1 1 0 0 0+ u' ) IDXVIB = 92
        IF ( STRVIB .EQ. ' 0 1 1 0 2 0+ u' ) IDXVIB = 93
        IF ( STRVIB .EQ. ' 0 1 1 1 0 1  u' ) IDXVIB = 94
        IF ( STRVIB .EQ. ' 0 1 1 2 0 0+ u' ) IDXVIB = 95
        IF ( STRVIB .EQ. ' 0 2 0 1 1 0+ u' ) IDXVIB = 96
        IF ( STRVIB .EQ. ' 0 2 0 1 3 0+ u' ) IDXVIB = 97
        IF ( STRVIB .EQ. ' 0 2 0 2 1 1 2u' ) IDXVIB = 98
        IF ( STRVIB .EQ. ' 0 2 0 3 1 0+ u' ) IDXVIB = 99
        IF ( STRVIB .EQ. ' 1 0 0 0 0 0+ g' ) IDXVIB = 100
        IF ( STRVIB .EQ. ' 1 0 0 0 1 1  u' ) IDXVIB = 101
        IF ( STRVIB .EQ. ' 1 0 0 0 2 0+ g' ) IDXVIB = 102
        IF ( STRVIB .EQ. ' 1 0 0 0 2 2  g' ) IDXVIB = 103
        IF ( STRVIB .EQ. ' 1 0 0 0 3 1  u' ) IDXVIB = 104
        IF ( STRVIB .EQ. ' 1 0 0 1 0 1  g' ) IDXVIB = 105
        IF ( STRVIB .EQ. ' 1 0 0 1 1 0+ u' ) IDXVIB = 106
        IF ( STRVIB .EQ. ' 1 0 0 1 1 0- u' ) IDXVIB = 107
        IF ( STRVIB .EQ. ' 1 0 0 1 1 2  u' ) IDXVIB = 108
        IF ( STRVIB .EQ. ' 1 0 0 1 2 1  g' ) IDXVIB = 109
        IF ( STRVIB .EQ. ' 1 0 0 2 1 1  u' ) IDXVIB = 110
        IF ( STRVIB .EQ. ' 1 0 0 2 1 1 1u' ) IDXVIB = 111
        IF ( STRVIB .EQ. ' 1 0 1 0 2 0+ u' ) IDXVIB = 112
        IF ( STRVIB .EQ. ' 1 0 1 0 2 2  u' ) IDXVIB = 113
        IF ( STRVIB .EQ. ' 1 0 1 1 1 0+ g' ) IDXVIB = 114
        IF ( STRVIB .EQ. ' 1 0 1 1 1 0- g' ) IDXVIB = 115
        IF ( STRVIB .EQ. ' 1 0 1 1 1 2  g' ) IDXVIB = 116
        IF ( STRVIB .EQ. ' 1 0 1 2 0 0+ u' ) IDXVIB = 117
        IF ( STRVIB .EQ. ' 1 0 1 2 0 2  u' ) IDXVIB = 118
        IF ( STRVIB .EQ. ' 1 1 0 1 2 1 2g' ) IDXVIB = 119
        IF ( STRVIB .EQ. ' 1 1 0 2 0 0+ g' ) IDXVIB = 120
        IF ( STRVIB .EQ. ' 1 1 0 2 1 1  u' ) IDXVIB = 121
        IF ( STRVIB .EQ. ' 1 1 0 2 1 1 1u' ) IDXVIB = 122
        IF ( STRVIB .EQ. ' 1 1 0 2 1 1 2u' ) IDXVIB = 123
        IF ( STRVIB .EQ. ' 1 1 1 0 0 0+ u' ) IDXVIB = 124
        IF ( STRVIB .EQ. ' 1 1 1 2 0 0+ u' ) IDXVIB = 125
        IF ( STRVIB .EQ. ' 1 2 0 1 1 0+ u' ) IDXVIB = 126
        IF ( STRVIB .EQ. ' 2 0 0 0 0 0+ g' ) IDXVIB = 127
        IF ( STRVIB .EQ. ' 2 0 0 0 1 1  u' ) IDXVIB = 128
        IF ( STRVIB .EQ. ' 2 0 0 1 0 1  g' ) IDXVIB = 129
        IF ( STRVIB .EQ. ' 2 0 1 0 0 0+ u' ) IDXVIB = 130

C
C Class 8: Pyramidal tetratomic molecules
      ELSE IF ( IDXMOL .EQ. 11 .OR.                                ! NH3
     &          IDXMOL .EQ. 28      ) THEN                         ! PH3
        IF ( STRVIB .EQ. '      0 0 0 0  ' ) IDXVIB = 1
        IF ( STRVIB .EQ. '      0 1 0 0  ' ) IDXVIB = 2
        IF ( STRVIB .EQ. '      0 2 0 0  ' ) IDXVIB = 3
        IF ( STRVIB .EQ. '      0 0 0 1  ' ) IDXVIB = 4
        IF ( STRVIB .EQ. '      0 1 0 1  ' ) IDXVIB = 5
        IF ( STRVIB .EQ. '      0 0 0 2  ' ) IDXVIB = 6
        IF ( STRVIB .EQ. '      0 0 1 0  ' ) IDXVIB = 7
        IF ( STRVIB .EQ. '      1 0 0 0  ' ) IDXVIB = 8
        IF ( STRVIB .EQ. '      0 0 0 0 a' ) IDXVIB = 9
        IF ( STRVIB .EQ. '      0 1 0 0 a' ) IDXVIB = 10
        IF ( STRVIB .EQ. '      0 2 0 0 a' ) IDXVIB = 11
        IF ( STRVIB .EQ. '      0 0 0 1 a' ) IDXVIB = 12
        IF ( STRVIB .EQ. '      0 0 0 0 s' ) IDXVIB = 13
        IF ( STRVIB .EQ. '      0 1 0 0 s' ) IDXVIB = 14
        IF ( STRVIB .EQ. '      0 2 0 0 s' ) IDXVIB = 15
        IF ( STRVIB .EQ. '      0 0 0 1 s' ) IDXVIB = 16
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 17
        IF ( STRVIB .EQ. '      0 3 0 0 s' ) IDXVIB = 18
        IF ( STRVIB .EQ. '      0 1 0 1 s' ) IDXVIB = 19
        IF ( STRVIB .EQ. '      0 1 0 1 a' ) IDXVIB = 20
        IF ( STRVIB .EQ. '      0 3 0 0 a' ) IDXVIB = 21
        IF ( STRVIB .EQ. '      0 2 0 1 s' ) IDXVIB = 22
        IF ( STRVIB .EQ. '      0 0 0 2As' ) IDXVIB = 23
        IF ( STRVIB .EQ. '      0 0 0 2Aa' ) IDXVIB = 24
        IF ( STRVIB .EQ. '      0 0 0 2Es' ) IDXVIB = 25
        IF ( STRVIB .EQ. '      0 0 0 2Ea' ) IDXVIB = 26
        IF ( STRVIB .EQ. '      1 0 0 0 s' ) IDXVIB = 27
        IF ( STRVIB .EQ. '      1 0 0 0 a' ) IDXVIB = 28
        IF ( STRVIB .EQ. '      0 0 1 0 s' ) IDXVIB = 29
        IF ( STRVIB .EQ. '      0 0 1 0 a' ) IDXVIB = 30
        IF ( STRVIB .EQ. '      0 4 0 0 s' ) IDXVIB = 31
        IF ( STRVIB .EQ. '      0 2 0 1 a' ) IDXVIB = 32
        IF ( STRVIB .EQ. '      0 3 0 1 s' ) IDXVIB = 33
        IF ( STRVIB .EQ. '      0 4 0 0 a' ) IDXVIB = 34
        IF ( STRVIB .EQ. '      0 1 0 2As' ) IDXVIB = 35
        IF ( STRVIB .EQ. '      0 1 0 2Es' ) IDXVIB = 36
        IF ( STRVIB .EQ. '      0 1 0 2Aa' ) IDXVIB = 37
        IF ( STRVIB .EQ. '      0 1 0 2Ea' ) IDXVIB = 38
        IF ( STRVIB .EQ. '      1 1 0 0 s' ) IDXVIB = 39
        IF ( STRVIB .EQ. '      1 1 0 0 a' ) IDXVIB = 40
        IF ( STRVIB .EQ. '      0 1 1 0 s' ) IDXVIB = 41
        IF ( STRVIB .EQ. '      0 1 1 0 a' ) IDXVIB = 42
        IF ( STRVIB .EQ. '      0 3 0 1 a' ) IDXVIB = 43
        IF ( STRVIB .EQ. '      1 0 0 1 s' ) IDXVIB = 44
        IF ( STRVIB .EQ. '      1 0 0 1 a' ) IDXVIB = 45
        IF ( STRVIB .EQ. '      0 0 1 1 s' ) IDXVIB = 46
        IF ( STRVIB .EQ. '      0 0 1 1 a' ) IDXVIB = 47
C new levels for NH3 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '    0000 00 0  ' ) IDXVIB = 48
        IF ( STRVIB .EQ. ' 0000 00 0     ' ) IDXVIB = 49
        IF ( STRVIB .EQ. ' 0000 00 0 A1'' ' ) IDXVIB = 50
        IF ( STRVIB .EQ. ' 0000 00 0 A2" ' ) IDXVIB = 51
        IF ( STRVIB .EQ. ' 0000 00 0 A2'' ' ) IDXVIB = 52
        IF ( STRVIB .EQ. ' 0001 01 1     ' ) IDXVIB = 53
        IF ( STRVIB .EQ. ' 0001 01 1 E"  ' ) IDXVIB = 54
        IF ( STRVIB .EQ. ' 0001 01 1 E''  ' ) IDXVIB = 55
        IF ( STRVIB .EQ. ' 0002 00 0     ' ) IDXVIB = 56
        IF ( STRVIB .EQ. ' 0002 00 0 A1'' ' ) IDXVIB = 57
        IF ( STRVIB .EQ. ' 0002 00 0 A2" ' ) IDXVIB = 58
        IF ( STRVIB .EQ. ' 0002 02 2     ' ) IDXVIB = 59
        IF ( STRVIB .EQ. ' 0002 02 2 E"  ' ) IDXVIB = 60
        IF ( STRVIB .EQ. ' 0002 02 2 E''  ' ) IDXVIB = 61
        IF ( STRVIB .EQ. ' 0003 01 1     ' ) IDXVIB = 62
        IF ( STRVIB .EQ. ' 0003 01 1 E"  ' ) IDXVIB = 63
        IF ( STRVIB .EQ. ' 0003 01 1 E''  ' ) IDXVIB = 64
        IF ( STRVIB .EQ. ' 0003 03 3     ' ) IDXVIB = 65
        IF ( STRVIB .EQ. ' 0003 03 3 A2'' ' ) IDXVIB = 66
        IF ( STRVIB .EQ. ' 0010 10 1     ' ) IDXVIB = 67
        IF ( STRVIB .EQ. ' 0010 10 1 E"  ' ) IDXVIB = 68
        IF ( STRVIB .EQ. ' 0010 10 1 E''  ' ) IDXVIB = 69
        IF ( STRVIB .EQ. ' 0011 11 0 A1" ' ) IDXVIB = 70
        IF ( STRVIB .EQ. ' 0011 11 2     ' ) IDXVIB = 71
        IF ( STRVIB .EQ. ' 0011 11 2 A1" ' ) IDXVIB = 72
        IF ( STRVIB .EQ. ' 0011 11 2 A1'' ' ) IDXVIB = 73
        IF ( STRVIB .EQ. ' 0011 11 2 A2" ' ) IDXVIB = 74
        IF ( STRVIB .EQ. ' 0011 11 2 A2'' ' ) IDXVIB = 75
        IF ( STRVIB .EQ. ' 0011 11 2 E"  ' ) IDXVIB = 76
        IF ( STRVIB .EQ. ' 0011 11 2 E''  ' ) IDXVIB = 77
        IF ( STRVIB .EQ. ' 0012 10 1 E"  ' ) IDXVIB = 78
        IF ( STRVIB .EQ. ' 0012 10 1 E''  ' ) IDXVIB = 79
        IF ( STRVIB .EQ. ' 0012 11 1 E"  ' ) IDXVIB = 80
        IF ( STRVIB .EQ. ' 0012 12 1 E"  ' ) IDXVIB = 81
        IF ( STRVIB .EQ. ' 0012 12 1 E''  ' ) IDXVIB = 82
        IF ( STRVIB .EQ. ' 0012 12 3     ' ) IDXVIB = 83
        IF ( STRVIB .EQ. ' 0020 00 0 A1'' ' ) IDXVIB = 84
        IF ( STRVIB .EQ. ' 0020 00 0 A2" ' ) IDXVIB = 85
        IF ( STRVIB .EQ. ' 0020 20 2 E"  ' ) IDXVIB = 86
        IF ( STRVIB .EQ. ' 0020 20 2 E''  ' ) IDXVIB = 87
        IF ( STRVIB .EQ. ' 0100 00 0     ' ) IDXVIB = 88
        IF ( STRVIB .EQ. ' 0100 00 0 A1'' ' ) IDXVIB = 89
        IF ( STRVIB .EQ. ' 0100 00 0 A2" ' ) IDXVIB = 90
        IF ( STRVIB .EQ. ' 0101 01 1     ' ) IDXVIB = 91
        IF ( STRVIB .EQ. ' 0101 01 1 E"  ' ) IDXVIB = 92
        IF ( STRVIB .EQ. ' 0101 01 1 E''  ' ) IDXVIB = 93
        IF ( STRVIB .EQ. ' 0102 00 0 A1'' ' ) IDXVIB = 94
        IF ( STRVIB .EQ. ' 0102 02 2 E''  ' ) IDXVIB = 95
        IF ( STRVIB .EQ. ' 0110 10 1     ' ) IDXVIB = 96
        IF ( STRVIB .EQ. ' 0110 10 1 E"  ' ) IDXVIB = 97
        IF ( STRVIB .EQ. ' 0110 10 1 E''  ' ) IDXVIB = 98
        IF ( STRVIB .EQ. ' 0200 00 0     ' ) IDXVIB = 99
        IF ( STRVIB .EQ. ' 0200 00 0 A1'' ' ) IDXVIB = 100
        IF ( STRVIB .EQ. ' 0200 00 0 A2" ' ) IDXVIB = 101
        IF ( STRVIB .EQ. ' 0201 01 1 E"  ' ) IDXVIB = 102
        IF ( STRVIB .EQ. ' 0201 01 1 E''  ' ) IDXVIB = 103
        IF ( STRVIB .EQ. ' 0202 00 0     ' ) IDXVIB = 104
        IF ( STRVIB .EQ. ' 0202 00 0 A1'' ' ) IDXVIB = 105
        IF ( STRVIB .EQ. ' 0202 00 0 A2" ' ) IDXVIB = 106
        IF ( STRVIB .EQ. ' 0202 02 2     ' ) IDXVIB = 107
        IF ( STRVIB .EQ. ' 0202 02 2 E"  ' ) IDXVIB = 108
        IF ( STRVIB .EQ. ' 0203 01 1 E"  ' ) IDXVIB = 109
        IF ( STRVIB .EQ. ' 0203 01 1 E''  ' ) IDXVIB = 110
        IF ( STRVIB .EQ. ' 0210 10 1     ' ) IDXVIB = 111
        IF ( STRVIB .EQ. ' 0210 10 1 E''  ' ) IDXVIB = 112
        IF ( STRVIB .EQ. ' 0300 00 0     ' ) IDXVIB = 113
        IF ( STRVIB .EQ. ' 0300 00 0 A1'' ' ) IDXVIB = 114
        IF ( STRVIB .EQ. ' 0300 00 0 A2" ' ) IDXVIB = 115
        IF ( STRVIB .EQ. ' 0301 01 1 E"  ' ) IDXVIB = 116
        IF ( STRVIB .EQ. ' 0400 00 0 A1'' ' ) IDXVIB = 117
        IF ( STRVIB .EQ. ' 0401 01 1     ' ) IDXVIB = 118
        IF ( STRVIB .EQ. ' 0401 01 1 E''  ' ) IDXVIB = 119
        IF ( STRVIB .EQ. ' 1000 00 0     ' ) IDXVIB = 120
        IF ( STRVIB .EQ. ' 1000 00 0 A1'' ' ) IDXVIB = 121
        IF ( STRVIB .EQ. ' 1000 00 0 A2" ' ) IDXVIB = 122
        IF ( STRVIB .EQ. ' 1001 01 1     ' ) IDXVIB = 123
        IF ( STRVIB .EQ. ' 1001 01 1 E"  ' ) IDXVIB = 124
        IF ( STRVIB .EQ. ' 1001 01 1 E''  ' ) IDXVIB = 125
        IF ( STRVIB .EQ. ' 1002 02 2 E"  ' ) IDXVIB = 126
        IF ( STRVIB .EQ. ' 1002 02 2 E''  ' ) IDXVIB = 127
        IF ( STRVIB .EQ. ' 1010 10 1 E"  ' ) IDXVIB = 128
        IF ( STRVIB .EQ. ' 1010 10 1 E''  ' ) IDXVIB = 129
        IF ( STRVIB .EQ. ' 1100 00 0     ' ) IDXVIB = 130
        IF ( STRVIB .EQ. ' 1100 00 0 A1'' ' ) IDXVIB = 131
        IF ( STRVIB .EQ. ' 1100 00 0 A2" ' ) IDXVIB = 132
        IF ( STRVIB .EQ. ' 1101 01 1 E''  ' ) IDXVIB = 133
        IF ( STRVIB .EQ. ' 1200 00 0     ' ) IDXVIB = 134
        IF ( STRVIB .EQ. ' 1200 00 0 A1'' ' ) IDXVIB = 135
C new levels for PH3 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '      0 0 1 1  ' ) IDXVIB = 136
        IF ( STRVIB .EQ. '      0 1 0 2  ' ) IDXVIB = 137
        IF ( STRVIB .EQ. '      0 1 1 0  ' ) IDXVIB = 138
        IF ( STRVIB .EQ. '      0 2 0 1  ' ) IDXVIB = 139
        IF ( STRVIB .EQ. '      0 3 0 0  ' ) IDXVIB = 140
        IF ( STRVIB .EQ. '      0 4 0 0  ' ) IDXVIB = 141
        IF ( STRVIB .EQ. '      1 0 0 1  ' ) IDXVIB = 142
        IF ( STRVIB .EQ. '      1 1 0 0  ' ) IDXVIB = 143
C
C Class 9: Non-linear tetratomic molecules 
      ELSE IF ( IDXMOL .EQ. 20 .OR.                                  ! H2CO
     &          IDXMOL .EQ. 25 .OR.                                  ! H2O2
     &          IDXMOL .EQ. 29      ) THEN                           ! COF2
        IF ( STRVIB .EQ. '    0 0 0 0 0 0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '    0 0 0 0 0 2' ) IDXVIB = 2
        IF ( STRVIB .EQ. '    0 0 1 1 0 0' ) IDXVIB = 3
        IF ( STRVIB .EQ. '    0 0 1 0 0 1' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    1 0 0 0 0 0' ) IDXVIB = 5
        IF ( STRVIB .EQ. '    0 0 0 0 1 0' ) IDXVIB = 6
        IF ( STRVIB .EQ. '    0 1 0 1 0 0' ) IDXVIB = 7
        IF ( STRVIB .EQ. '    0 1 0 0 0 1' ) IDXVIB = 8
        IF ( STRVIB .EQ. '    0 0 0 0 0 1' ) IDXVIB = 9
        IF ( STRVIB .EQ. '    0 1 0 0 0 0' ) IDXVIB = 10
        IF ( STRVIB .EQ. '    0 0 0 1 0 0' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 2 0 0 0 0' ) IDXVIB = 12
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 13
        IF ( STRVIB .EQ. '    0 0 2 0 0 1' ) IDXVIB = 14
        IF ( STRVIB .EQ. '    0 0 0 0 2 0' ) IDXVIB = 15
        IF ( STRVIB .EQ. '    0 0 000 0 0' ) IDXVIB = 16
        IF ( STRVIB .EQ. '    0 0 000 0 1' ) IDXVIB = 17
        IF ( STRVIB .EQ. '    0 0 001 0 0' ) IDXVIB = 18
        IF ( STRVIB .EQ. '    0 0 002 0 0' ) IDXVIB = 19
        IF ( STRVIB .EQ. '    0 0 011 0 0' ) IDXVIB = 20
        IF ( STRVIB .EQ. '    0 0 012 0 0' ) IDXVIB = 21
        IF ( STRVIB .EQ. '    0 0 021 0 0' ) IDXVIB = 22
        IF ( STRVIB .EQ. '    0 0 022 0 0' ) IDXVIB = 23
        IF ( STRVIB .EQ. '    0 0 031 0 0' ) IDXVIB = 24
        IF ( STRVIB .EQ. '    0 0 032 0 0' ) IDXVIB = 25
        IF ( STRVIB .EQ. '    0 0 101 0 0' ) IDXVIB = 26
        IF ( STRVIB .EQ. '    0 0 102 0 0' ) IDXVIB = 27
        IF ( STRVIB .EQ. '    0 0 111 0 0' ) IDXVIB = 28
        IF ( STRVIB .EQ. '    0 0 112 0 0' ) IDXVIB = 29
        IF ( STRVIB .EQ. '    0 0 003 0 0' ) IDXVIB = 30
        IF ( STRVIB .EQ. '    0 0 004 0 0' ) IDXVIB = 31
        IF ( STRVIB .EQ. '    0 0 013 0 0' ) IDXVIB = 32
        IF ( STRVIB .EQ. '    0 0 014 0 0' ) IDXVIB = 33
        IF ( STRVIB .EQ. '    0 0 023 0 0' ) IDXVIB = 34
        IF ( STRVIB .EQ. '    0 0 024 0 0' ) IDXVIB = 35
        IF ( STRVIB .EQ. '    0 0 033 0 0' ) IDXVIB = 36
        IF ( STRVIB .EQ. '    0 0 034 0 0' ) IDXVIB = 37
        IF ( STRVIB .EQ. '    0 0 103 0 0' ) IDXVIB = 38
        IF ( STRVIB .EQ. '    0 0 104 0 0' ) IDXVIB = 39
C new levels for H2O2 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '    0 0 001 0 1' ) IDXVIB = 40
        IF ( STRVIB .EQ. '    0 0 002 0 1' ) IDXVIB = 41
        IF ( STRVIB .EQ. '    0 0 003 0 1' ) IDXVIB = 42
        IF ( STRVIB .EQ. '    0 0 004 0 1' ) IDXVIB = 43
        IF ( STRVIB .EQ. '    0 0 011 0 1' ) IDXVIB = 44
        IF ( STRVIB .EQ. '    0 0 012 0 1' ) IDXVIB = 45
        IF ( STRVIB .EQ. '    0 0 013 0 1' ) IDXVIB = 46
        IF ( STRVIB .EQ. '    0 0 014 0 1' ) IDXVIB = 47
        IF ( STRVIB .EQ. '    0 0 021 0 1' ) IDXVIB = 48
        IF ( STRVIB .EQ. '    0 0 022 0 1' ) IDXVIB = 49
        IF ( STRVIB .EQ. '    0 0 023 0 1' ) IDXVIB = 50
        IF ( STRVIB .EQ. '    0 0 024 0 1' ) IDXVIB = 51
        IF ( STRVIB .EQ. '    0 0 113 0 0' ) IDXVIB = 52
        IF ( STRVIB .EQ. '    0 0 114 0 0' ) IDXVIB = 53
        IF ( STRVIB .EQ. '    0 0 131 0 0' ) IDXVIB = 54
        IF ( STRVIB .EQ. '    0 0 132 0 0' ) IDXVIB = 55
        IF ( STRVIB .EQ. '    0 0 133 0 0' ) IDXVIB = 56
        IF ( STRVIB .EQ. '    0 0 134 0 0' ) IDXVIB = 57
        IF ( STRVIB .EQ. '    0 1 001 0 0' ) IDXVIB = 58
        IF ( STRVIB .EQ. '    0 1 002 0 0' ) IDXVIB = 59
        IF ( STRVIB .EQ. '    0 1 003 0 0' ) IDXVIB = 60
        IF ( STRVIB .EQ. '    0 1 004 0 0' ) IDXVIB = 61
        IF ( STRVIB .EQ. '    0 1 011 0 0' ) IDXVIB = 62
        IF ( STRVIB .EQ. '    0 1 012 0 0' ) IDXVIB = 63
        IF ( STRVIB .EQ. '    0 1 013 0 0' ) IDXVIB = 64
        IF ( STRVIB .EQ. '    0 1 014 0 0' ) IDXVIB = 65
C new levels for COF2 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '    0 0 0 0 1 1' ) IDXVIB = 66
        IF ( STRVIB .EQ. '    0 0 1 0 0 0' ) IDXVIB = 67
C new levels for H2CO found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '    0 0 0 1 0 1' ) IDXVIB = 68
        IF ( STRVIB .EQ. '    0 0 0 2 0 0' ) IDXVIB = 69
        IF ( STRVIB .EQ. '    0 0 2 0 0 0' ) IDXVIB = 70
        IF ( STRVIB .EQ. '    0 1 1 0 0 0' ) IDXVIB = 71
C
C Class 10: Pentatomic or greater polyatomic molecules
      ELSE IF ( IDXMOL .EQ.  6 .OR.                                ! CH4
     &          IDXMOL .EQ. 12 .OR.                                ! HNO3
     &          IDXMOL .EQ. 24 .OR.                                ! CH3Cl
     &          IDXMOL .EQ. 27 .OR.                                ! C2H6
     &          IDXMOL .EQ. 30 .OR.                                ! SF6
     &          IDXMOL .EQ. 32 .OR.                                ! HCOOH
     &          IDXMOL .EQ. 35 .OR.                                ! ClONO2
     &          IDXMOL .EQ. 38 .OR.                                ! C2H4
     &          IDXMOL .EQ. 39 .OR.                                ! CH3OH
     &          IDXMOL .EQ. 40 .OR.                                ! CH3Br
     &          IDXMOL .EQ. 41 .OR.                                ! CH3CN
     &          IDXMOL .EQ. 42      ) THEN                         ! CF4
        IF ( STRVIB .EQ. '         GROUND' ) IDXVIB = 1
        IF ( STRVIB .EQ. '             V1' ) IDXVIB = 2
        IF ( STRVIB .EQ. '             V2' ) IDXVIB = 3
        IF ( STRVIB .EQ. '             V4' ) IDXVIB = 4
        IF ( STRVIB .EQ. '             V5' ) IDXVIB = 5
        IF ( STRVIB .EQ. '             V9' ) IDXVIB = 6
        IF ( STRVIB .EQ. '            2V5' ) IDXVIB = 7
        IF ( STRVIB .EQ. '            2V9' ) IDXVIB = 8
        IF ( STRVIB .EQ. '            3V6' ) IDXVIB = 9
        IF ( STRVIB .EQ. '            3V9' ) IDXVIB = 10
        IF ( STRVIB .EQ. '          V5+V9' ) IDXVIB = 11
        IF ( STRVIB .EQ. '               ' ) IDXVIB = 12
        IF ( STRVIB .EQ. '             V6' ) IDXVIB = 13
        IF ( STRVIB .EQ. '             V3' ) IDXVIB = 14
        IF ( STRVIB .EQ. '            2V6' ) IDXVIB = 15
        IF ( STRVIB .EQ. '             V7' ) IDXVIB = 16
        IF ( STRVIB .EQ. '             V8' ) IDXVIB = 17
        IF ( STRVIB .EQ. '          V8+V9' ) IDXVIB = 18
        IF ( STRVIB .EQ. '          V3+V6' ) IDXVIB = 19
        IF ( STRVIB .EQ. '            2V3' ) IDXVIB = 20
        IF ( STRVIB .EQ. '          V5+V6' ) IDXVIB = 21
        IF ( STRVIB .EQ. '          V3+V5' ) IDXVIB = 22
        IF ( STRVIB .EQ. '          V4+V9' ) IDXVIB = 23
        IF ( STRVIB .EQ. '            V10' ) IDXVIB = 24
        IF ( STRVIB .EQ. '            V11' ) IDXVIB = 25
        IF ( STRVIB .EQ. '         V2+V12' ) IDXVIB = 26
        IF ( STRVIB .EQ. '       2V10+V12' ) IDXVIB = 27
        IF ( STRVIB .EQ. '         V9+V10' ) IDXVIB = 28
        IF ( STRVIB .EQ. '            V12' ) IDXVIB = 29
        IF ( STRVIB .EQ. '           2V12' ) IDXVIB = 30
        IF ( STRVIB .EQ. '         V6+V12' ) IDXVIB = 31
        IF ( STRVIB .EQ. '         V7+V12' ) IDXVIB = 32
        IF ( STRVIB .EQ. '         V8+V12' ) IDXVIB = 33
        IF ( STRVIB .EQ. '        V8+2V12' ) IDXVIB = 34
        IF ( STRVIB .EQ. '           3V12' ) IDXVIB = 35
        IF ( STRVIB .EQ. '           4V12' ) IDXVIB = 36
C 15Mar16 modified to use CH4 assignments from IAA for vib temps which only
C use first four indices
C 0 0 0 0 assume ground state
        IF ( STRVIB .EQ. '    0 0 0 0    ' ) IDXVIB = 1
        IF ( STRVIB .EQ. '    0 0 0 0 1A1' ) IDXVIB = 1
C 0 1 0 0 CH4 Energy 1310.7606 cm-1
        IF ( STRVIB .EQ. '    0 0 0 1 1F2' ) IDXVIB = 2
C 0 1 0 0 CH4 Energy 1533.3320 cm-1
        IF ( STRVIB .EQ. '    0 1 0 0 1 E' ) IDXVIB = 3
        IF ( STRVIB .EQ. '    0 1 0 0 1E ' ) IDXVIB = 3
C 0 0 0 2 CH4 Energy 2608.6899 cm-1
        IF ( STRVIB .EQ. '    0 0 0 2 1 E' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    0 0 0 2 1A1' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    0 0 0 2 1F2' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    0 0 0 2 1 E' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    0 0 0 2 1A1' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    0 0 0 2 1F2' ) IDXVIB = 4
        IF ( STRVIB .EQ. '    0 0 0 2 1E ' ) IDXVIB = 4
C 0 1 0 1 CH4 Energy 2838.2600 cm-1
        IF ( STRVIB .EQ. '    0 1 0 1 1F1' ) IDXVIB = 6
        IF ( STRVIB .EQ. '    0 1 0 1 1F2' ) IDXVIB = 6
C 1 0 0 0 CH4 Energy 2916.4829 cm-1
        IF ( STRVIB .EQ. '    1 0 0 0 1A1' ) IDXVIB = 7
C 0 0 1 0 CH4 Energy 3019.4971 cm-1
        IF ( STRVIB .EQ. '    0 0 1 0 1F2' ) IDXVIB = 8
C 0 2 0 0 CH4 Energy 3064.4800 cm-1
        IF ( STRVIB .EQ. '    0 2 0 0 1 E' ) IDXVIB = 9
        IF ( STRVIB .EQ. '    0 2 0 0 1A1' ) IDXVIB = 9
        IF ( STRVIB .EQ. '    0 2 0 0 1E ' ) IDXVIB = 9
C 1 0 0 1 CH4 Energy 4223.4600 cm-1
        IF ( STRVIB .EQ. '    1 0 0 1 1F2' ) IDXVIB = 10
C 0 0 1 1 CH4 Energy 4321.6699 cm-1
        IF ( STRVIB .EQ. '    0 0 1 1 1 E' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1A1' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1F1' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1F2' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1 E' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1A1' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1F1' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1F2' ) IDXVIB = 11
        IF ( STRVIB .EQ. '    0 0 1 1 1E ' ) IDXVIB = 11
C 0 1 1 0 CH4 Energy 4540.6499 cm-1
        IF ( STRVIB .EQ. '    0 1 1 0 1F1' ) IDXVIB = 12
        IF ( STRVIB .EQ. '    0 1 1 0 1F2' ) IDXVIB = 12
C 0 0 2 0 CH4 Energy 6024.2798 cm-1
        IF ( STRVIB .EQ. '    0 0 2 0 1F2' ) IDXVIB = 13
        IF ( STRVIB .EQ. '    0 0 2 0  A1' ) IDXVIB = 13
        IF ( STRVIB .EQ. '    0 0 2 0  E ' ) IDXVIB = 13
        IF ( STRVIB .EQ. '    0 0 2 0  F2' ) IDXVIB = 13
C
        IF ( STRVIB .EQ. '    0 0 0 3 1A1' ) IDXVIB = 43
        IF ( STRVIB .EQ. '    0 0 0 3 1F1' ) IDXVIB = 44
        IF ( STRVIB .EQ. '    0 0 0 3 1F2' ) IDXVIB = 45
        IF ( STRVIB .EQ. '    0 0 0 3 2F2' ) IDXVIB = 46
        IF ( STRVIB .EQ. '    0 0 0 4    ' ) IDXVIB = 47
        IF ( STRVIB .EQ. '    0 0 1 2    ' ) IDXVIB = 53
        IF ( STRVIB .EQ. '    0 0 1 2 1F2' ) IDXVIB = 54
        IF ( STRVIB .EQ. '    0 0 3 0    ' ) IDXVIB = 56
        IF ( STRVIB .EQ. '    0 1 0 2 1 E' ) IDXVIB = 60
        IF ( STRVIB .EQ. '    0 1 0 2 1A1' ) IDXVIB = 61
        IF ( STRVIB .EQ. '    0 1 0 2 1A2' ) IDXVIB = 62
        IF ( STRVIB .EQ. '    0 1 0 2 1F1' ) IDXVIB = 63
        IF ( STRVIB .EQ. '    0 1 0 2 1F2' ) IDXVIB = 64
        IF ( STRVIB .EQ. '    0 1 0 2 2 E' ) IDXVIB = 65
        IF ( STRVIB .EQ. '    0 1 2 0    ' ) IDXVIB = 68
        IF ( STRVIB .EQ. '    0 2 0 1 1F1' ) IDXVIB = 71
        IF ( STRVIB .EQ. '    0 2 0 1 1F2' ) IDXVIB = 72
        IF ( STRVIB .EQ. '    0 2 0 1 2F2' ) IDXVIB = 73
        IF ( STRVIB .EQ. '    0 3 0 0    ' ) IDXVIB = 74
        IF ( STRVIB .EQ. '    0 3 0 0 1 E' ) IDXVIB = 75
        IF ( STRVIB .EQ. '    0 3 0 0 1A1' ) IDXVIB = 76
        IF ( STRVIB .EQ. '    0 3 0 0 1A2' ) IDXVIB = 77
        IF ( STRVIB .EQ. '    1 1 0 0 1 E' ) IDXVIB = 80
        IF ( STRVIB .EQ. '    0 0 0 3 1A1' ) IDXVIB = 84
        IF ( STRVIB .EQ. '    0 0 0 3 1F1' ) IDXVIB = 85
        IF ( STRVIB .EQ. '    0 0 0 3 1F2' ) IDXVIB = 86
        IF ( STRVIB .EQ. '    0 0 0 3 2F2' ) IDXVIB = 87
C new levels for CH4 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '            3V2' ) IDXVIB = 92
        IF ( STRVIB .EQ. '            3V3' ) IDXVIB = 93
        IF ( STRVIB .EQ. '          V2+V3' ) IDXVIB = 94
        IF ( STRVIB .EQ. '          V2+V5' ) IDXVIB = 95
        IF ( STRVIB .EQ. '          V2+V6' ) IDXVIB = 96
        IF ( STRVIB .EQ. '         2V2 E ' ) IDXVIB = 97
        IF ( STRVIB .EQ. '         2V3+V6' ) IDXVIB = 98
        IF ( STRVIB .EQ. '         3V5 A1' ) IDXVIB = 99
        IF ( STRVIB .EQ. '         3V5 E ' ) IDXVIB = 100
        IF ( STRVIB .EQ. '         V3+2V6' ) IDXVIB = 101
        IF ( STRVIB .EQ. '         V5+2V6' ) IDXVIB = 102
        IF ( STRVIB .EQ. '       V1+V2+V6' ) IDXVIB = 103
        IF ( STRVIB .EQ. '       V1+V5 E ' ) IDXVIB = 104
        IF ( STRVIB .EQ. '       V1+V6 E ' ) IDXVIB = 105
        IF ( STRVIB .EQ. '       V2+V4+V6' ) IDXVIB = 106
        IF ( STRVIB .EQ. '       V3+V4 E ' ) IDXVIB = 107
        IF ( STRVIB .EQ. '       V3+V5+V6' ) IDXVIB = 108
        IF ( STRVIB .EQ. '       V4+V5 A1' ) IDXVIB = 109
        IF ( STRVIB .EQ. '       V4+V5 A2' ) IDXVIB = 110
        IF ( STRVIB .EQ. '       V4+V5 E ' ) IDXVIB = 111
        IF ( STRVIB .EQ. '       V4+V6 A1' ) IDXVIB = 112
        IF ( STRVIB .EQ. '       V4+V6 A2' ) IDXVIB = 113
        IF ( STRVIB .EQ. '       V4+V6 E ' ) IDXVIB = 114
        IF ( STRVIB .EQ. '      2V3+V5 E ' ) IDXVIB = 115
        IF ( STRVIB .EQ. '      2V5+V6 A1' ) IDXVIB = 116
        IF ( STRVIB .EQ. '      2V5+V6 A2' ) IDXVIB = 117
        IF ( STRVIB .EQ. '      2V5+V6 E ' ) IDXVIB = 118
        IF ( STRVIB .EQ. '      V2+2V5+V6' ) IDXVIB = 119
        IF ( STRVIB .EQ. '      V2+2V6 A1' ) IDXVIB = 120
        IF ( STRVIB .EQ. '      V2+2V6 E ' ) IDXVIB = 121
        IF ( STRVIB .EQ. '      V3+2V5 A1' ) IDXVIB = 122
        IF ( STRVIB .EQ. '      V3+2V5 E ' ) IDXVIB = 123
        IF ( STRVIB .EQ. '    0 0 0 5  A1' ) IDXVIB = 125
        IF ( STRVIB .EQ. '    0 0 0 5  E ' ) IDXVIB = 126
        IF ( STRVIB .EQ. '    0 0 0 5  F1' ) IDXVIB = 127
        IF ( STRVIB .EQ. '    0 0 0 5  F2' ) IDXVIB = 128
        IF ( STRVIB .EQ. '    0 0 1 2 1A1' ) IDXVIB = 130
        IF ( STRVIB .EQ. '    0 0 1 2 1E ' ) IDXVIB = 131
        IF ( STRVIB .EQ. '    0 0 1 2 1F1' ) IDXVIB = 132
        IF ( STRVIB .EQ. '    0 1 0 2 1E ' ) IDXVIB = 137
        IF ( STRVIB .EQ. '    0 1 0 2 2E ' ) IDXVIB = 138
        IF ( STRVIB .EQ. '    0 1 0 3 1F1' ) IDXVIB = 139
        IF ( STRVIB .EQ. '    0 1 0 4  A1' ) IDXVIB = 140
        IF ( STRVIB .EQ. '    0 1 0 4  A2' ) IDXVIB = 141
        IF ( STRVIB .EQ. '    0 1 0 4  E ' ) IDXVIB = 142
        IF ( STRVIB .EQ. '    0 1 0 4  F1' ) IDXVIB = 143
        IF ( STRVIB .EQ. '    0 1 0 4  F2' ) IDXVIB = 144
        IF ( STRVIB .EQ. '    0 1 1 1    ' ) IDXVIB = 145
        IF ( STRVIB .EQ. '    0 1 1 1  A1' ) IDXVIB = 146
        IF ( STRVIB .EQ. '    0 1 1 1  A2' ) IDXVIB = 147
        IF ( STRVIB .EQ. '    0 1 1 1  E ' ) IDXVIB = 148
        IF ( STRVIB .EQ. '    0 1 1 1  F1' ) IDXVIB = 149
        IF ( STRVIB .EQ. '    0 1 1 1  F2' ) IDXVIB = 150
        IF ( STRVIB .EQ. '    0 1 1 1 1A1' ) IDXVIB = 151
        IF ( STRVIB .EQ. '    0 1 1 1 1A2' ) IDXVIB = 152
        IF ( STRVIB .EQ. '    0 1 1 1 1E ' ) IDXVIB = 153
        IF ( STRVIB .EQ. '    0 1 1 1 1F1' ) IDXVIB = 154
        IF ( STRVIB .EQ. '    0 1 1 1 1F2' ) IDXVIB = 155
        IF ( STRVIB .EQ. '    0 2 0 2 1A1' ) IDXVIB = 157
        IF ( STRVIB .EQ. '    0 2 0 2 1E ' ) IDXVIB = 158
        IF ( STRVIB .EQ. '    0 2 0 2 1F1' ) IDXVIB = 159
        IF ( STRVIB .EQ. '    0 2 0 2 1F2' ) IDXVIB = 160
        IF ( STRVIB .EQ. '    0 2 1 0    ' ) IDXVIB = 161
        IF ( STRVIB .EQ. '    0 2 1 0  F1' ) IDXVIB = 162
        IF ( STRVIB .EQ. '    0 2 1 0  F2' ) IDXVIB = 163
        IF ( STRVIB .EQ. '    0 3 0 0 1E ' ) IDXVIB = 164
        IF ( STRVIB .EQ. '    0 3 0 1  F1' ) IDXVIB = 165
        IF ( STRVIB .EQ. '    0 3 0 1  F2' ) IDXVIB = 166
        IF ( STRVIB .EQ. '    0 3 0 1 1F1' ) IDXVIB = 167
        IF ( STRVIB .EQ. '    0 3 0 1 1F2' ) IDXVIB = 168
        IF ( STRVIB .EQ. '    1 0 0 2 1F2' ) IDXVIB = 169
        IF ( STRVIB .EQ. '    1 0 0 3  A1' ) IDXVIB = 170
        IF ( STRVIB .EQ. '    1 0 0 3  F2' ) IDXVIB = 171
        IF ( STRVIB .EQ. '    1 0 1 0  F2' ) IDXVIB = 172
        IF ( STRVIB .EQ. '    1 0 1 0 1F2' ) IDXVIB = 173
        IF ( STRVIB .EQ. '    1 1 0 0 1E ' ) IDXVIB = 174
        IF ( STRVIB .EQ. '    1 1 0 1 1F2' ) IDXVIB = 175
        IF ( STRVIB .EQ. '    1 2 0 0  E ' ) IDXVIB = 176
        IF ( STRVIB .EQ. '    1 2 0 0  F1' ) IDXVIB = 177
        IF ( STRVIB .EQ. '    2 0 0 0  A1' ) IDXVIB = 178
        IF ( STRVIB .EQ. '    2 0 0 0 1A1' ) IDXVIB = 179
C new levels for HNO3 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '          V5+V7' ) IDXVIB = 180
        IF ( STRVIB .EQ. '          V6+V7' ) IDXVIB = 181
C new levels for C2H6 found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '            2V4' ) IDXVIB = 182
        IF ( STRVIB .EQ. '            3V4' ) IDXVIB = 183
        IF ( STRVIB .EQ. '          V4+V6' ) IDXVIB = 184
        IF ( STRVIB .EQ. '          V4+V8' ) IDXVIB = 185
        IF ( STRVIB .EQ. '         2V4+V9' ) IDXVIB = 186
        IF ( STRVIB .EQ. '         V4+V12' ) IDXVIB = 187
        IF ( STRVIB .EQ. '        2V4+V12' ) IDXVIB = 188
C new level for HCOOH found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '          V6+V9' ) IDXVIB = 189
C new level for CH3Cl found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '         2V3+V5' ) IDXVIB = 190
      ELSE IF ( IDXMOL .EQ. 43 ) THEN                               ! C4H2
C new levels found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. ' 000000000 e g ' ) IDXVIB = 1
        IF ( STRVIB .EQ. ' 000000000 e u ' ) IDXVIB = 2
        IF ( STRVIB .EQ. ' 000000001 e g ' ) IDXVIB = 3
        IF ( STRVIB .EQ. ' 000000001 e u ' ) IDXVIB = 4
        IF ( STRVIB .EQ. ' 000000001 f g ' ) IDXVIB = 5
        IF ( STRVIB .EQ. ' 000000001 f u ' ) IDXVIB = 6
        IF ( STRVIB .EQ. ' 000000002 e g ' ) IDXVIB = 7
        IF ( STRVIB .EQ. ' 000000002 e u ' ) IDXVIB = 8
        IF ( STRVIB .EQ. ' 000000002 f g ' ) IDXVIB = 9
        IF ( STRVIB .EQ. ' 000000002 f u ' ) IDXVIB = 10
        IF ( STRVIB .EQ. ' 000000003 e g ' ) IDXVIB = 11
        IF ( STRVIB .EQ. ' 000000003 e u ' ) IDXVIB = 12
        IF ( STRVIB .EQ. ' 000000003 f g ' ) IDXVIB = 13
        IF ( STRVIB .EQ. ' 000000003 f u ' ) IDXVIB = 14
        IF ( STRVIB .EQ. ' 000000004 e g ' ) IDXVIB = 15
        IF ( STRVIB .EQ. ' 000000004 e u ' ) IDXVIB = 16
        IF ( STRVIB .EQ. ' 000000004 f g ' ) IDXVIB = 17
        IF ( STRVIB .EQ. ' 000000004 f u ' ) IDXVIB = 18
        IF ( STRVIB .EQ. ' 000000005 e g ' ) IDXVIB = 19
        IF ( STRVIB .EQ. ' 000000005 e u ' ) IDXVIB = 20
        IF ( STRVIB .EQ. ' 000000005 f g ' ) IDXVIB = 21
        IF ( STRVIB .EQ. ' 000000005 f u ' ) IDXVIB = 22
        IF ( STRVIB .EQ. ' 000000006 e g ' ) IDXVIB = 23
        IF ( STRVIB .EQ. ' 000000006 e u ' ) IDXVIB = 24
        IF ( STRVIB .EQ. ' 000000006 f g ' ) IDXVIB = 25
        IF ( STRVIB .EQ. ' 000000006 f u ' ) IDXVIB = 26
        IF ( STRVIB .EQ. ' 000000007 e g ' ) IDXVIB = 27
        IF ( STRVIB .EQ. ' 000000007 e u ' ) IDXVIB = 28
        IF ( STRVIB .EQ. ' 000000007 f g ' ) IDXVIB = 29
        IF ( STRVIB .EQ. ' 000000007 f u ' ) IDXVIB = 30
        IF ( STRVIB .EQ. ' 000000008 e g ' ) IDXVIB = 31
        IF ( STRVIB .EQ. ' 000000008 e u ' ) IDXVIB = 32
        IF ( STRVIB .EQ. ' 000000008 f u ' ) IDXVIB = 33
        IF ( STRVIB .EQ. ' 000000010 e g ' ) IDXVIB = 34
        IF ( STRVIB .EQ. ' 000000010 e u ' ) IDXVIB = 35
        IF ( STRVIB .EQ. ' 000000010 f g ' ) IDXVIB = 36
        IF ( STRVIB .EQ. ' 000000010 f u ' ) IDXVIB = 37
        IF ( STRVIB .EQ. ' 000000011 e g ' ) IDXVIB = 38
        IF ( STRVIB .EQ. ' 000000011 e u ' ) IDXVIB = 39
        IF ( STRVIB .EQ. ' 000000011 f g ' ) IDXVIB = 40
        IF ( STRVIB .EQ. ' 000000011 f u ' ) IDXVIB = 41
        IF ( STRVIB .EQ. ' 000000012 e g ' ) IDXVIB = 42
        IF ( STRVIB .EQ. ' 000000012 e u ' ) IDXVIB = 43
        IF ( STRVIB .EQ. ' 000000012 f g ' ) IDXVIB = 44
        IF ( STRVIB .EQ. ' 000000012 f u ' ) IDXVIB = 45
        IF ( STRVIB .EQ. ' 000000013 e g ' ) IDXVIB = 46
        IF ( STRVIB .EQ. ' 000000013 e u ' ) IDXVIB = 47
        IF ( STRVIB .EQ. ' 000000013 f g ' ) IDXVIB = 48
        IF ( STRVIB .EQ. ' 000000013 f u ' ) IDXVIB = 49
        IF ( STRVIB .EQ. ' 000000014 e g ' ) IDXVIB = 50
        IF ( STRVIB .EQ. ' 000000014 e u ' ) IDXVIB = 51
        IF ( STRVIB .EQ. ' 000000014 f g ' ) IDXVIB = 52
        IF ( STRVIB .EQ. ' 000000014 f u ' ) IDXVIB = 53
        IF ( STRVIB .EQ. ' 000000015 e g ' ) IDXVIB = 54
        IF ( STRVIB .EQ. ' 000000015 e u ' ) IDXVIB = 55
        IF ( STRVIB .EQ. ' 000000015 f g ' ) IDXVIB = 56
        IF ( STRVIB .EQ. ' 000000015 f u ' ) IDXVIB = 57
        IF ( STRVIB .EQ. ' 000000016 e g ' ) IDXVIB = 58
        IF ( STRVIB .EQ. ' 000000016 e u ' ) IDXVIB = 59
        IF ( STRVIB .EQ. ' 000000016 f g ' ) IDXVIB = 60
        IF ( STRVIB .EQ. ' 000000016 f u ' ) IDXVIB = 61
        IF ( STRVIB .EQ. ' 000000020 e g ' ) IDXVIB = 62
        IF ( STRVIB .EQ. ' 000000020 e u ' ) IDXVIB = 63
        IF ( STRVIB .EQ. ' 000000020 f g ' ) IDXVIB = 64
        IF ( STRVIB .EQ. ' 000000020 f u ' ) IDXVIB = 65
        IF ( STRVIB .EQ. ' 000000021 e g ' ) IDXVIB = 66
        IF ( STRVIB .EQ. ' 000000021 e u ' ) IDXVIB = 67
        IF ( STRVIB .EQ. ' 000000021 f g ' ) IDXVIB = 68
        IF ( STRVIB .EQ. ' 000000021 f u ' ) IDXVIB = 69
        IF ( STRVIB .EQ. ' 000000022 e g ' ) IDXVIB = 70
        IF ( STRVIB .EQ. ' 000000022 e u ' ) IDXVIB = 71
        IF ( STRVIB .EQ. ' 000000022 f g ' ) IDXVIB = 72
        IF ( STRVIB .EQ. ' 000000022 f u ' ) IDXVIB = 73
        IF ( STRVIB .EQ. ' 000000023 e g ' ) IDXVIB = 74
        IF ( STRVIB .EQ. ' 000000023 e u ' ) IDXVIB = 75
        IF ( STRVIB .EQ. ' 000000023 f g ' ) IDXVIB = 76
        IF ( STRVIB .EQ. ' 000000023 f u ' ) IDXVIB = 77
        IF ( STRVIB .EQ. ' 000000030 e g ' ) IDXVIB = 78
        IF ( STRVIB .EQ. ' 000000030 e u ' ) IDXVIB = 79
        IF ( STRVIB .EQ. ' 000000030 f g ' ) IDXVIB = 80
        IF ( STRVIB .EQ. ' 000000030 f u ' ) IDXVIB = 81
        IF ( STRVIB .EQ. ' 000000100 e g ' ) IDXVIB = 82
        IF ( STRVIB .EQ. ' 000000100 e u ' ) IDXVIB = 83
        IF ( STRVIB .EQ. ' 000000100 f g ' ) IDXVIB = 84
        IF ( STRVIB .EQ. ' 000000100 f u ' ) IDXVIB = 85
        IF ( STRVIB .EQ. ' 000000101 e g ' ) IDXVIB = 86
        IF ( STRVIB .EQ. ' 000000101 e u ' ) IDXVIB = 87
        IF ( STRVIB .EQ. ' 000000101 f g ' ) IDXVIB = 88
        IF ( STRVIB .EQ. ' 000000101 f u ' ) IDXVIB = 89
        IF ( STRVIB .EQ. ' 000000102 e g ' ) IDXVIB = 90
        IF ( STRVIB .EQ. ' 000000102 e u ' ) IDXVIB = 91
        IF ( STRVIB .EQ. ' 000000102 f g ' ) IDXVIB = 92
        IF ( STRVIB .EQ. ' 000000102 f u ' ) IDXVIB = 93
        IF ( STRVIB .EQ. ' 000000103 e g ' ) IDXVIB = 94
        IF ( STRVIB .EQ. ' 000000103 e u ' ) IDXVIB = 95
        IF ( STRVIB .EQ. ' 000000103 f g ' ) IDXVIB = 96
        IF ( STRVIB .EQ. ' 000000103 f u ' ) IDXVIB = 97
        IF ( STRVIB .EQ. ' 000000105 e g ' ) IDXVIB = 98
        IF ( STRVIB .EQ. ' 000000105 e u ' ) IDXVIB = 99
        IF ( STRVIB .EQ. ' 000000105 f g ' ) IDXVIB = 100
        IF ( STRVIB .EQ. ' 000000105 f u ' ) IDXVIB = 101
        IF ( STRVIB .EQ. ' 000000106 e g ' ) IDXVIB = 102
        IF ( STRVIB .EQ. ' 000000106 e u ' ) IDXVIB = 103
        IF ( STRVIB .EQ. ' 000000106 f g ' ) IDXVIB = 104
        IF ( STRVIB .EQ. ' 000000106 f u ' ) IDXVIB = 105
        IF ( STRVIB .EQ. ' 000000110 e g ' ) IDXVIB = 106
        IF ( STRVIB .EQ. ' 000000110 e u ' ) IDXVIB = 107
        IF ( STRVIB .EQ. ' 000000110 f g ' ) IDXVIB = 108
        IF ( STRVIB .EQ. ' 000000110 f u ' ) IDXVIB = 109
        IF ( STRVIB .EQ. ' 000000111 e g ' ) IDXVIB = 110
        IF ( STRVIB .EQ. ' 000000111 e u ' ) IDXVIB = 111
        IF ( STRVIB .EQ. ' 000000111 f g ' ) IDXVIB = 112
        IF ( STRVIB .EQ. ' 000000111 f u ' ) IDXVIB = 113
        IF ( STRVIB .EQ. ' 000000112 e g ' ) IDXVIB = 114
        IF ( STRVIB .EQ. ' 000000112 e u ' ) IDXVIB = 115
        IF ( STRVIB .EQ. ' 000000112 f g ' ) IDXVIB = 116
        IF ( STRVIB .EQ. ' 000000112 f u ' ) IDXVIB = 117
        IF ( STRVIB .EQ. ' 000000113 e g ' ) IDXVIB = 118
        IF ( STRVIB .EQ. ' 000000113 e u ' ) IDXVIB = 119
        IF ( STRVIB .EQ. ' 000000113 f g ' ) IDXVIB = 120
        IF ( STRVIB .EQ. ' 000000113 f u ' ) IDXVIB = 121
        IF ( STRVIB .EQ. ' 000000120 e g ' ) IDXVIB = 122
        IF ( STRVIB .EQ. ' 000000120 e u ' ) IDXVIB = 123
        IF ( STRVIB .EQ. ' 000000120 f g ' ) IDXVIB = 124
        IF ( STRVIB .EQ. ' 000000120 f u ' ) IDXVIB = 125
        IF ( STRVIB .EQ. ' 000000200 e g ' ) IDXVIB = 126
        IF ( STRVIB .EQ. ' 000000200 e u ' ) IDXVIB = 127
        IF ( STRVIB .EQ. ' 000000200 f g ' ) IDXVIB = 128
        IF ( STRVIB .EQ. ' 000000200 f u ' ) IDXVIB = 129
        IF ( STRVIB .EQ. ' 000000201 e g ' ) IDXVIB = 130
        IF ( STRVIB .EQ. ' 000000201 e u ' ) IDXVIB = 131
        IF ( STRVIB .EQ. ' 000000201 f g ' ) IDXVIB = 132
        IF ( STRVIB .EQ. ' 000000201 f u ' ) IDXVIB = 133
        IF ( STRVIB .EQ. ' 000000202 e g ' ) IDXVIB = 134
        IF ( STRVIB .EQ. ' 000000202 e u ' ) IDXVIB = 135
        IF ( STRVIB .EQ. ' 000000202 f g ' ) IDXVIB = 136
        IF ( STRVIB .EQ. ' 000000202 f u ' ) IDXVIB = 137
        IF ( STRVIB .EQ. ' 000000210 e g ' ) IDXVIB = 138
        IF ( STRVIB .EQ. ' 000000210 e u ' ) IDXVIB = 139
        IF ( STRVIB .EQ. ' 000000210 f g ' ) IDXVIB = 140
        IF ( STRVIB .EQ. ' 000000210 f u ' ) IDXVIB = 141
        IF ( STRVIB .EQ. ' 000000211 e g ' ) IDXVIB = 142
        IF ( STRVIB .EQ. ' 000000211 e u ' ) IDXVIB = 143
        IF ( STRVIB .EQ. ' 000000211 f g ' ) IDXVIB = 144
        IF ( STRVIB .EQ. ' 000000211 f u ' ) IDXVIB = 145
        IF ( STRVIB .EQ. ' 000000212 e g ' ) IDXVIB = 146
        IF ( STRVIB .EQ. ' 000000212 e u ' ) IDXVIB = 147
        IF ( STRVIB .EQ. ' 000000212 f g ' ) IDXVIB = 148
        IF ( STRVIB .EQ. ' 000000212 f u ' ) IDXVIB = 149
        IF ( STRVIB .EQ. ' 000001000 e g ' ) IDXVIB = 150
        IF ( STRVIB .EQ. ' 000001000 e u ' ) IDXVIB = 151
        IF ( STRVIB .EQ. ' 000001000 f g ' ) IDXVIB = 152
        IF ( STRVIB .EQ. ' 000001000 f u ' ) IDXVIB = 153
        IF ( STRVIB .EQ. ' 000001001 e g ' ) IDXVIB = 154
        IF ( STRVIB .EQ. ' 000001001 e u ' ) IDXVIB = 155
        IF ( STRVIB .EQ. ' 000001001 f g ' ) IDXVIB = 156
        IF ( STRVIB .EQ. ' 000001001 f u ' ) IDXVIB = 157
        IF ( STRVIB .EQ. ' 000001002 e g ' ) IDXVIB = 158
        IF ( STRVIB .EQ. ' 000001002 e u ' ) IDXVIB = 159
        IF ( STRVIB .EQ. ' 000001002 f g ' ) IDXVIB = 160
        IF ( STRVIB .EQ. ' 000001002 f u ' ) IDXVIB = 161
        IF ( STRVIB .EQ. ' 000001003 e g ' ) IDXVIB = 162
        IF ( STRVIB .EQ. ' 000001003 e u ' ) IDXVIB = 163
        IF ( STRVIB .EQ. ' 000001003 f g ' ) IDXVIB = 164
        IF ( STRVIB .EQ. ' 000001003 f u ' ) IDXVIB = 165
        IF ( STRVIB .EQ. ' 000001005 f u ' ) IDXVIB = 166
        IF ( STRVIB .EQ. ' 000001006 e g ' ) IDXVIB = 167
        IF ( STRVIB .EQ. ' 000001006 e u ' ) IDXVIB = 168
        IF ( STRVIB .EQ. ' 000001006 f g ' ) IDXVIB = 169
        IF ( STRVIB .EQ. ' 000001010 e g ' ) IDXVIB = 170
        IF ( STRVIB .EQ. ' 000001010 e u ' ) IDXVIB = 171
        IF ( STRVIB .EQ. ' 000001010 f g ' ) IDXVIB = 172
        IF ( STRVIB .EQ. ' 000001010 f u ' ) IDXVIB = 173
        IF ( STRVIB .EQ. ' 000001011 e g ' ) IDXVIB = 174
        IF ( STRVIB .EQ. ' 000001011 e u ' ) IDXVIB = 175
        IF ( STRVIB .EQ. ' 000001011 f g ' ) IDXVIB = 176
        IF ( STRVIB .EQ. ' 000001011 f u ' ) IDXVIB = 177
        IF ( STRVIB .EQ. ' 000001012 e g ' ) IDXVIB = 178
        IF ( STRVIB .EQ. ' 000001012 e u ' ) IDXVIB = 179
        IF ( STRVIB .EQ. ' 000001012 f g ' ) IDXVIB = 180
        IF ( STRVIB .EQ. ' 000001012 f u ' ) IDXVIB = 181
        IF ( STRVIB .EQ. ' 000001013 e g ' ) IDXVIB = 182
        IF ( STRVIB .EQ. ' 000001013 e u ' ) IDXVIB = 183
        IF ( STRVIB .EQ. ' 000001013 f g ' ) IDXVIB = 184
        IF ( STRVIB .EQ. ' 000001013 f u ' ) IDXVIB = 185
        IF ( STRVIB .EQ. ' 000001020 e g ' ) IDXVIB = 186
        IF ( STRVIB .EQ. ' 000001020 e u ' ) IDXVIB = 187
        IF ( STRVIB .EQ. ' 000001020 f g ' ) IDXVIB = 188
        IF ( STRVIB .EQ. ' 000001020 f u ' ) IDXVIB = 189
        IF ( STRVIB .EQ. ' 000001100 e g ' ) IDXVIB = 190
        IF ( STRVIB .EQ. ' 000001100 e u ' ) IDXVIB = 191
        IF ( STRVIB .EQ. ' 000001100 f g ' ) IDXVIB = 192
        IF ( STRVIB .EQ. ' 000001100 f u ' ) IDXVIB = 193
        IF ( STRVIB .EQ. ' 000001110 e g ' ) IDXVIB = 194
        IF ( STRVIB .EQ. ' 000001110 e u ' ) IDXVIB = 195
        IF ( STRVIB .EQ. ' 000001110 f g ' ) IDXVIB = 196
        IF ( STRVIB .EQ. ' 000001110 f u ' ) IDXVIB = 197
        IF ( STRVIB .EQ. ' 000002000 e g ' ) IDXVIB = 198
        IF ( STRVIB .EQ. ' 000002000 e u ' ) IDXVIB = 199
        IF ( STRVIB .EQ. ' 000002000 f g ' ) IDXVIB = 200
        IF ( STRVIB .EQ. ' 000002000 f u ' ) IDXVIB = 201
        IF ( STRVIB .EQ. ' 000002001 e g ' ) IDXVIB = 202
        IF ( STRVIB .EQ. ' 000002001 e u ' ) IDXVIB = 203
        IF ( STRVIB .EQ. ' 000002001 f g ' ) IDXVIB = 204
        IF ( STRVIB .EQ. ' 000002001 f u ' ) IDXVIB = 205
        IF ( STRVIB .EQ. ' 000002002 e g ' ) IDXVIB = 206
        IF ( STRVIB .EQ. ' 000002002 e u ' ) IDXVIB = 207
        IF ( STRVIB .EQ. ' 000002002 f g ' ) IDXVIB = 208
        IF ( STRVIB .EQ. ' 000002002 f u ' ) IDXVIB = 209
        IF ( STRVIB .EQ. ' 000002003 e g ' ) IDXVIB = 210
        IF ( STRVIB .EQ. ' 000002003 e u ' ) IDXVIB = 211
        IF ( STRVIB .EQ. ' 000002003 f g ' ) IDXVIB = 212
        IF ( STRVIB .EQ. ' 000002003 f u ' ) IDXVIB = 213
        IF ( STRVIB .EQ. ' 000002010 e g ' ) IDXVIB = 214
        IF ( STRVIB .EQ. ' 000002010 e u ' ) IDXVIB = 215
        IF ( STRVIB .EQ. ' 000002010 f g ' ) IDXVIB = 216
        IF ( STRVIB .EQ. ' 000002010 f u ' ) IDXVIB = 217
        IF ( STRVIB .EQ. ' 000002100 e g ' ) IDXVIB = 218
        IF ( STRVIB .EQ. ' 000002100 e u ' ) IDXVIB = 219
        IF ( STRVIB .EQ. ' 000002100 f g ' ) IDXVIB = 220
        IF ( STRVIB .EQ. ' 000002100 f u ' ) IDXVIB = 221
        IF ( STRVIB .EQ. ' 000003000 e g ' ) IDXVIB = 222
        IF ( STRVIB .EQ. ' 000003000 e u ' ) IDXVIB = 223
        IF ( STRVIB .EQ. ' 000003000 f g ' ) IDXVIB = 224
        IF ( STRVIB .EQ. ' 000003000 f u ' ) IDXVIB = 225
        IF ( STRVIB .EQ. ' 001000000 e g ' ) IDXVIB = 226
        IF ( STRVIB .EQ. ' 001000000 e u ' ) IDXVIB = 227
        IF ( STRVIB .EQ. ' 001000001 e g ' ) IDXVIB = 228
        IF ( STRVIB .EQ. ' 001000001 e u ' ) IDXVIB = 229
        IF ( STRVIB .EQ. ' 001000001 f g ' ) IDXVIB = 230
        IF ( STRVIB .EQ. ' 001000001 f u ' ) IDXVIB = 231
        IF ( STRVIB .EQ. ' 001000002 e g ' ) IDXVIB = 232
        IF ( STRVIB .EQ. ' 001000002 e u ' ) IDXVIB = 233
        IF ( STRVIB .EQ. ' 001000002 f g ' ) IDXVIB = 234
        IF ( STRVIB .EQ. ' 001000002 f u ' ) IDXVIB = 235
        IF ( STRVIB .EQ. ' 001000003 e g ' ) IDXVIB = 236
        IF ( STRVIB .EQ. ' 001000003 e u ' ) IDXVIB = 237
        IF ( STRVIB .EQ. ' 001000003 f g ' ) IDXVIB = 238
        IF ( STRVIB .EQ. ' 001000003 f u ' ) IDXVIB = 239
        IF ( STRVIB .EQ. ' 001000004 e g ' ) IDXVIB = 240
        IF ( STRVIB .EQ. ' 001000004 e u ' ) IDXVIB = 241
        IF ( STRVIB .EQ. ' 001000004 f g ' ) IDXVIB = 242
        IF ( STRVIB .EQ. ' 001000004 f u ' ) IDXVIB = 243
        IF ( STRVIB .EQ. ' 001000010 e g ' ) IDXVIB = 244
        IF ( STRVIB .EQ. ' 001000010 e u ' ) IDXVIB = 245
        IF ( STRVIB .EQ. ' 001000010 f g ' ) IDXVIB = 246
        IF ( STRVIB .EQ. ' 001000010 f u ' ) IDXVIB = 247
        IF ( STRVIB .EQ. ' 001000011 e g ' ) IDXVIB = 248
        IF ( STRVIB .EQ. ' 001000011 e u ' ) IDXVIB = 249
        IF ( STRVIB .EQ. ' 001000011 f g ' ) IDXVIB = 250
        IF ( STRVIB .EQ. ' 001000011 f u ' ) IDXVIB = 251
        IF ( STRVIB .EQ. ' 001000012 e g ' ) IDXVIB = 252
        IF ( STRVIB .EQ. ' 001000012 e u ' ) IDXVIB = 253
        IF ( STRVIB .EQ. ' 001000012 f g ' ) IDXVIB = 254
        IF ( STRVIB .EQ. ' 001000012 f u ' ) IDXVIB = 255
        IF ( STRVIB .EQ. ' 001000100 e g ' ) IDXVIB = 256
        IF ( STRVIB .EQ. ' 001000100 e u ' ) IDXVIB = 257
        IF ( STRVIB .EQ. ' 001000100 f g ' ) IDXVIB = 258
        IF ( STRVIB .EQ. ' 001000100 f u ' ) IDXVIB = 259
        IF ( STRVIB .EQ. ' 001000101 e g ' ) IDXVIB = 260
        IF ( STRVIB .EQ. ' 001000101 e u ' ) IDXVIB = 261
        IF ( STRVIB .EQ. ' 001000101 f g ' ) IDXVIB = 262
        IF ( STRVIB .EQ. ' 001000101 f u ' ) IDXVIB = 263
        IF ( STRVIB .EQ. ' 001000102 e g ' ) IDXVIB = 264
        IF ( STRVIB .EQ. ' 001000102 e u ' ) IDXVIB = 265
        IF ( STRVIB .EQ. ' 001000102 f g ' ) IDXVIB = 266
        IF ( STRVIB .EQ. ' 001000102 f u ' ) IDXVIB = 267
        IF ( STRVIB .EQ. ' 001000200 e g ' ) IDXVIB = 268
        IF ( STRVIB .EQ. ' 001000200 e u ' ) IDXVIB = 269
        IF ( STRVIB .EQ. ' 001000200 f g ' ) IDXVIB = 270
        IF ( STRVIB .EQ. ' 001000200 f u ' ) IDXVIB = 271
        IF ( STRVIB .EQ. ' 001001000 e g ' ) IDXVIB = 272
        IF ( STRVIB .EQ. ' 001001000 e u ' ) IDXVIB = 273
        IF ( STRVIB .EQ. ' 001001000 f g ' ) IDXVIB = 274
        IF ( STRVIB .EQ. ' 001001000 f u ' ) IDXVIB = 275
        IF ( STRVIB .EQ. ' 001001001 e g ' ) IDXVIB = 276
        IF ( STRVIB .EQ. ' 001001001 e u ' ) IDXVIB = 277
        IF ( STRVIB .EQ. ' 001001001 f g ' ) IDXVIB = 278
        IF ( STRVIB .EQ. ' 001001001 f u ' ) IDXVIB = 279
        IF ( STRVIB .EQ. ' 001001002 e g ' ) IDXVIB = 280
        IF ( STRVIB .EQ. ' 001001002 e u ' ) IDXVIB = 281
        IF ( STRVIB .EQ. ' 001001002 f g ' ) IDXVIB = 282
        IF ( STRVIB .EQ. ' 001001002 f u ' ) IDXVIB = 283
        IF ( STRVIB .EQ. ' 002000000 e g ' ) IDXVIB = 284
        IF ( STRVIB .EQ. ' 002000000 e u ' ) IDXVIB = 285
        IF ( STRVIB .EQ. ' 002000001 e g ' ) IDXVIB = 286
        IF ( STRVIB .EQ. ' 002000001 e u ' ) IDXVIB = 287
        IF ( STRVIB .EQ. ' 002000001 f g ' ) IDXVIB = 288
        IF ( STRVIB .EQ. ' 002000001 f u ' ) IDXVIB = 289
      ELSE IF ( IDXMOL .EQ. 44 ) THEN                               ! HC3N
C new levels found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. '  0000000 0 0 0' ) IDXVIB = 1
        IF ( STRVIB .EQ. '  0000001 0 0 1' ) IDXVIB = 2
        IF ( STRVIB .EQ. '  0000002 0 0 0' ) IDXVIB = 3
        IF ( STRVIB .EQ. '  0000002 0 0 2' ) IDXVIB = 4
        IF ( STRVIB .EQ. '  0000003 0 0 1' ) IDXVIB = 5
        IF ( STRVIB .EQ. '  0000003 0 0 3' ) IDXVIB = 6
        IF ( STRVIB .EQ. '  0000004 0 0 0' ) IDXVIB = 7
        IF ( STRVIB .EQ. '  0000004 0 0 2' ) IDXVIB = 8
        IF ( STRVIB .EQ. '  0000004 0 0 4' ) IDXVIB = 9
        IF ( STRVIB .EQ. '  0000005 0 0 1' ) IDXVIB = 10
        IF ( STRVIB .EQ. '  0000005 0 0 3' ) IDXVIB = 11
        IF ( STRVIB .EQ. '  0000005 0 0 5' ) IDXVIB = 12
        IF ( STRVIB .EQ. '  0000006 0 0 0' ) IDXVIB = 13
        IF ( STRVIB .EQ. '  0000006 0 0 2' ) IDXVIB = 14
        IF ( STRVIB .EQ. '  0000006 0 0 4' ) IDXVIB = 15
        IF ( STRVIB .EQ. '  0000006 0 0 6' ) IDXVIB = 16
        IF ( STRVIB .EQ. '  0000007 0 0 1' ) IDXVIB = 17
        IF ( STRVIB .EQ. '  0000007 0 0 3' ) IDXVIB = 18
        IF ( STRVIB .EQ. '  0000007 0 0 5' ) IDXVIB = 19
        IF ( STRVIB .EQ. '  0000007 0 0 7' ) IDXVIB = 20
        IF ( STRVIB .EQ. '  0000008 0 0 0' ) IDXVIB = 21
        IF ( STRVIB .EQ. '  0000008 0 0 2' ) IDXVIB = 22
        IF ( STRVIB .EQ. '  0000008 0 0 4' ) IDXVIB = 23
        IF ( STRVIB .EQ. '  0000008 0 0 6' ) IDXVIB = 24
        IF ( STRVIB .EQ. '  0000009 0 0 1' ) IDXVIB = 25
        IF ( STRVIB .EQ. '  0000009 0 0 3' ) IDXVIB = 26
        IF ( STRVIB .EQ. '  0000009 0 0 5' ) IDXVIB = 27
        IF ( STRVIB .EQ. '  0000009 0 0 7' ) IDXVIB = 28
        IF ( STRVIB .EQ. '  0000010 0 1 0' ) IDXVIB = 29
        IF ( STRVIB .EQ. '  0000011 0 1 1' ) IDXVIB = 30
        IF ( STRVIB .EQ. '  0000011 0 1-1' ) IDXVIB = 31
        IF ( STRVIB .EQ. '  0000012 0 1 0' ) IDXVIB = 32
        IF ( STRVIB .EQ. '  0000012 0 1 2' ) IDXVIB = 33
        IF ( STRVIB .EQ. '  0000012 0-1 2' ) IDXVIB = 34
        IF ( STRVIB .EQ. '  0000013 0 1 1' ) IDXVIB = 35
        IF ( STRVIB .EQ. '  0000013 0 1 3' ) IDXVIB = 36
        IF ( STRVIB .EQ. '  0000013 0 1-1' ) IDXVIB = 37
        IF ( STRVIB .EQ. '  0000013 0-1 3' ) IDXVIB = 38
        IF ( STRVIB .EQ. '  0000014 0 1 0' ) IDXVIB = 39
        IF ( STRVIB .EQ. '  0000014 0 1 2' ) IDXVIB = 40
        IF ( STRVIB .EQ. '  0000014 0 1 4' ) IDXVIB = 41
        IF ( STRVIB .EQ. '  0000014 0-1 2' ) IDXVIB = 42
        IF ( STRVIB .EQ. '  0000014 0-1 4' ) IDXVIB = 43
        IF ( STRVIB .EQ. '  0000015 0 1 1' ) IDXVIB = 44
        IF ( STRVIB .EQ. '  0000015 0 1 3' ) IDXVIB = 45
        IF ( STRVIB .EQ. '  0000015 0 1 5' ) IDXVIB = 46
        IF ( STRVIB .EQ. '  0000015 0 1-1' ) IDXVIB = 47
        IF ( STRVIB .EQ. '  0000015 0-1 3' ) IDXVIB = 48
        IF ( STRVIB .EQ. '  0000015 0-1 5' ) IDXVIB = 49
        IF ( STRVIB .EQ. '  0000016 0 1 0' ) IDXVIB = 50
        IF ( STRVIB .EQ. '  0000016 0 1 2' ) IDXVIB = 51
        IF ( STRVIB .EQ. '  0000016 0 1 4' ) IDXVIB = 52
        IF ( STRVIB .EQ. '  0000016 0 1 6' ) IDXVIB = 53
        IF ( STRVIB .EQ. '  0000017 0 1 1' ) IDXVIB = 54
        IF ( STRVIB .EQ. '  0000017 0 1 3' ) IDXVIB = 55
        IF ( STRVIB .EQ. '  0000017 0 1 5' ) IDXVIB = 56
        IF ( STRVIB .EQ. '  0000017 0 1 7' ) IDXVIB = 57
        IF ( STRVIB .EQ. '  0000017 0 1-1' ) IDXVIB = 58
        IF ( STRVIB .EQ. '  0000020 0 0 0' ) IDXVIB = 59
        IF ( STRVIB .EQ. '  0000020 0 2 0' ) IDXVIB = 60
        IF ( STRVIB .EQ. '  0000021 0 0 1' ) IDXVIB = 61
        IF ( STRVIB .EQ. '  0000021 0 2 1' ) IDXVIB = 62
        IF ( STRVIB .EQ. '  0000021 0 2-1' ) IDXVIB = 63
        IF ( STRVIB .EQ. '  0000022 0 0 0' ) IDXVIB = 64
        IF ( STRVIB .EQ. '  0000022 0 0 2' ) IDXVIB = 65
        IF ( STRVIB .EQ. '  0000022 0 2 0' ) IDXVIB = 66
        IF ( STRVIB .EQ. '  0000022 0 2 2' ) IDXVIB = 67
        IF ( STRVIB .EQ. '  0000022 0 2-2' ) IDXVIB = 68
        IF ( STRVIB .EQ. '  0000023 0 0 1' ) IDXVIB = 69
        IF ( STRVIB .EQ. '  0000023 0 0 3' ) IDXVIB = 70
        IF ( STRVIB .EQ. '  0000023 0 2 1' ) IDXVIB = 71
        IF ( STRVIB .EQ. '  0000023 0 2 3' ) IDXVIB = 72
        IF ( STRVIB .EQ. '  0000023 0 2-1' ) IDXVIB = 73
        IF ( STRVIB .EQ. '  0000023 0-2 3' ) IDXVIB = 74
        IF ( STRVIB .EQ. '  0000024 0 0 0' ) IDXVIB = 75
        IF ( STRVIB .EQ. '  0000024 0 0 2' ) IDXVIB = 76
        IF ( STRVIB .EQ. '  0000024 0 0 4' ) IDXVIB = 77
        IF ( STRVIB .EQ. '  0000024 0 2 0' ) IDXVIB = 78
        IF ( STRVIB .EQ. '  0000024 0 2 2' ) IDXVIB = 79
        IF ( STRVIB .EQ. '  0000024 0 2 4' ) IDXVIB = 80
        IF ( STRVIB .EQ. '  0000024 0 2-2' ) IDXVIB = 81
        IF ( STRVIB .EQ. '  0000025 0 0 1' ) IDXVIB = 82
        IF ( STRVIB .EQ. '  0000025 0 0 3' ) IDXVIB = 83
        IF ( STRVIB .EQ. '  0000025 0 0 5' ) IDXVIB = 84
        IF ( STRVIB .EQ. '  0000025 0 2 1' ) IDXVIB = 85
        IF ( STRVIB .EQ. '  0000025 0 2 3' ) IDXVIB = 86
        IF ( STRVIB .EQ. '  0000025 0 2 5' ) IDXVIB = 87
        IF ( STRVIB .EQ. '  0000025 0 2-1' ) IDXVIB = 88
        IF ( STRVIB .EQ. '  0000030 0 1 0' ) IDXVIB = 89
        IF ( STRVIB .EQ. '  0000030 0 3 0' ) IDXVIB = 90
        IF ( STRVIB .EQ. '  0000031 0 1 1' ) IDXVIB = 91
        IF ( STRVIB .EQ. '  0000031 0 1-1' ) IDXVIB = 92
        IF ( STRVIB .EQ. '  0000031 0 3 1' ) IDXVIB = 93
        IF ( STRVIB .EQ. '  0000031 0 3-1' ) IDXVIB = 94
        IF ( STRVIB .EQ. '  0000032 0 1 0' ) IDXVIB = 95
        IF ( STRVIB .EQ. '  0000032 0 1 2' ) IDXVIB = 96
        IF ( STRVIB .EQ. '  0000032 0 3 0' ) IDXVIB = 97
        IF ( STRVIB .EQ. '  0000032 0 3 2' ) IDXVIB = 98
        IF ( STRVIB .EQ. '  0000032 0 3-2' ) IDXVIB = 99
        IF ( STRVIB .EQ. '  0000033 0 1 1' ) IDXVIB = 100
        IF ( STRVIB .EQ. '  0000033 0 1 3' ) IDXVIB = 101
        IF ( STRVIB .EQ. '  0000033 0 1-1' ) IDXVIB = 102
        IF ( STRVIB .EQ. '  0000033 0 3 1' ) IDXVIB = 103
        IF ( STRVIB .EQ. '  0000033 0 3 3' ) IDXVIB = 104
        IF ( STRVIB .EQ. '  0000033 0 3-1' ) IDXVIB = 105
        IF ( STRVIB .EQ. '  0000033 0 3-3' ) IDXVIB = 106
        IF ( STRVIB .EQ. '  0000040 0 0 0' ) IDXVIB = 107
        IF ( STRVIB .EQ. '  0000040 0 2 0' ) IDXVIB = 108
        IF ( STRVIB .EQ. '  0000040 0 4 0' ) IDXVIB = 109
        IF ( STRVIB .EQ. '  0000041 0 0 1' ) IDXVIB = 110
        IF ( STRVIB .EQ. '  0000041 0 2 1' ) IDXVIB = 111
        IF ( STRVIB .EQ. '  0000041 0 2-1' ) IDXVIB = 112
        IF ( STRVIB .EQ. '  0000041 0 4 1' ) IDXVIB = 113
        IF ( STRVIB .EQ. '  0000041 0 4-1' ) IDXVIB = 114
        IF ( STRVIB .EQ. '  0000100 1 0 0' ) IDXVIB = 115
        IF ( STRVIB .EQ. '  0000101 1 0 1' ) IDXVIB = 116
        IF ( STRVIB .EQ. '  0000101 1 0-1' ) IDXVIB = 117
        IF ( STRVIB .EQ. '  0000102 1 0 0' ) IDXVIB = 118
        IF ( STRVIB .EQ. '  0000102 1 0 2' ) IDXVIB = 119
        IF ( STRVIB .EQ. '  0000102-1 0 2' ) IDXVIB = 120
        IF ( STRVIB .EQ. '  0000103 1 0 1' ) IDXVIB = 121
        IF ( STRVIB .EQ. '  0000103 1 0 3' ) IDXVIB = 122
        IF ( STRVIB .EQ. '  0000103 1 0-1' ) IDXVIB = 123
        IF ( STRVIB .EQ. '  0000103-1 0 3' ) IDXVIB = 124
        IF ( STRVIB .EQ. '  0000104 1 0 0' ) IDXVIB = 125
        IF ( STRVIB .EQ. '  0000104 1 0 2' ) IDXVIB = 126
        IF ( STRVIB .EQ. '  0000104 1 0 4' ) IDXVIB = 127
        IF ( STRVIB .EQ. '  0000104-1 0 2' ) IDXVIB = 128
        IF ( STRVIB .EQ. '  0000104-1 0 4' ) IDXVIB = 129
        IF ( STRVIB .EQ. '  0000105 1 0 1' ) IDXVIB = 130
        IF ( STRVIB .EQ. '  0000105 1 0 3' ) IDXVIB = 131
        IF ( STRVIB .EQ. '  0000105 1 0 5' ) IDXVIB = 132
        IF ( STRVIB .EQ. '  0000105 1 0-1' ) IDXVIB = 133
        IF ( STRVIB .EQ. '  0000105-1 0 3' ) IDXVIB = 134
        IF ( STRVIB .EQ. '  0000105-1 0 5' ) IDXVIB = 135
        IF ( STRVIB .EQ. '  0000106 1 0 0' ) IDXVIB = 136
        IF ( STRVIB .EQ. '  0000106 1 0 2' ) IDXVIB = 137
        IF ( STRVIB .EQ. '  0000106 1 0 4' ) IDXVIB = 138
        IF ( STRVIB .EQ. '  0000106 1 0 6' ) IDXVIB = 139
        IF ( STRVIB .EQ. '  0000106-1 0 2' ) IDXVIB = 140
        IF ( STRVIB .EQ. '  0000106-1 0 4' ) IDXVIB = 141
        IF ( STRVIB .EQ. '  0000106-1 0 6' ) IDXVIB = 142
        IF ( STRVIB .EQ. '  0000110 1 1 0' ) IDXVIB = 143
        IF ( STRVIB .EQ. '  0000110 1-1 0' ) IDXVIB = 144
        IF ( STRVIB .EQ. '  0000111 1 1 1' ) IDXVIB = 145
        IF ( STRVIB .EQ. '  0000111 1 1-1' ) IDXVIB = 146
        IF ( STRVIB .EQ. '  0000111 1-1 1' ) IDXVIB = 147
        IF ( STRVIB .EQ. '  0000111-1 1 1' ) IDXVIB = 148
        IF ( STRVIB .EQ. '  0000112 1 1 0' ) IDXVIB = 149
        IF ( STRVIB .EQ. '  0000112 1 1 2' ) IDXVIB = 150
        IF ( STRVIB .EQ. '  0000112 1 1-2' ) IDXVIB = 151
        IF ( STRVIB .EQ. '  0000112 1-1 0' ) IDXVIB = 152
        IF ( STRVIB .EQ. '  0000112 1-1 2' ) IDXVIB = 153
        IF ( STRVIB .EQ. '  0000112-1 1 2' ) IDXVIB = 154
        IF ( STRVIB .EQ. '  0000113 1 1 1' ) IDXVIB = 155
        IF ( STRVIB .EQ. '  0000113 1 1 3' ) IDXVIB = 156
        IF ( STRVIB .EQ. '  0000113 1 1-1' ) IDXVIB = 157
        IF ( STRVIB .EQ. '  0000113-1 1 1' ) IDXVIB = 158
        IF ( STRVIB .EQ. '  0000113-1 1 3' ) IDXVIB = 159
        IF ( STRVIB .EQ. '  0000114 1 1 0' ) IDXVIB = 160
        IF ( STRVIB .EQ. '  0000114 1 1 2' ) IDXVIB = 161
        IF ( STRVIB .EQ. '  0000114 1 1 4' ) IDXVIB = 162
        IF ( STRVIB .EQ. '  0000114 1 1-2' ) IDXVIB = 163
        IF ( STRVIB .EQ. '  0000114-1 1 2' ) IDXVIB = 164
        IF ( STRVIB .EQ. '  0000114-1 1 4' ) IDXVIB = 165
        IF ( STRVIB .EQ. '  0000120 1 0 0' ) IDXVIB = 166
        IF ( STRVIB .EQ. '  0000120 1 2 0' ) IDXVIB = 167
        IF ( STRVIB .EQ. '  0000120-1 2 0' ) IDXVIB = 168
        IF ( STRVIB .EQ. '  0000121 1 0 1' ) IDXVIB = 169
        IF ( STRVIB .EQ. '  0000121 1 0-1' ) IDXVIB = 170
        IF ( STRVIB .EQ. '  0000121 1 2 1' ) IDXVIB = 171
        IF ( STRVIB .EQ. '  0000121 1 2-1' ) IDXVIB = 172
        IF ( STRVIB .EQ. '  0000121-1 2 1' ) IDXVIB = 173
        IF ( STRVIB .EQ. '  0000122 1 0 0' ) IDXVIB = 174
        IF ( STRVIB .EQ. '  0000122 1 0 2' ) IDXVIB = 175
        IF ( STRVIB .EQ. '  0000122 1 2 0' ) IDXVIB = 176
        IF ( STRVIB .EQ. '  0000122 1 2 2' ) IDXVIB = 177
        IF ( STRVIB .EQ. '  0000122 1 2-2' ) IDXVIB = 178
        IF ( STRVIB .EQ. '  0000122-1 0 2' ) IDXVIB = 179
        IF ( STRVIB .EQ. '  0000122-1 2 0' ) IDXVIB = 180
        IF ( STRVIB .EQ. '  0000122-1 2 2' ) IDXVIB = 181
        IF ( STRVIB .EQ. '  0000130 1 1 0' ) IDXVIB = 182
        IF ( STRVIB .EQ. '  0000130 1 3 0' ) IDXVIB = 183
        IF ( STRVIB .EQ. '  0000130-1 3 0' ) IDXVIB = 184
        IF ( STRVIB .EQ. '  0000200 0 0 0' ) IDXVIB = 185
        IF ( STRVIB .EQ. '  0000200 2 0 0' ) IDXVIB = 186
        IF ( STRVIB .EQ. '  0000201 0 0 1' ) IDXVIB = 187
        IF ( STRVIB .EQ. '  0000201 2 0 1' ) IDXVIB = 188
        IF ( STRVIB .EQ. '  0000201 2 0-1' ) IDXVIB = 189
        IF ( STRVIB .EQ. '  0000202 0 0 0' ) IDXVIB = 190
        IF ( STRVIB .EQ. '  0000202 0 0 2' ) IDXVIB = 191
        IF ( STRVIB .EQ. '  0000202 2 0 0' ) IDXVIB = 192
        IF ( STRVIB .EQ. '  0000202 2 0 2' ) IDXVIB = 193
        IF ( STRVIB .EQ. '  0000202 2 0-2' ) IDXVIB = 194
        IF ( STRVIB .EQ. '  0000203 0 0 1' ) IDXVIB = 195
        IF ( STRVIB .EQ. '  0000203 0 0 3' ) IDXVIB = 196
        IF ( STRVIB .EQ. '  0000203 2 0 1' ) IDXVIB = 197
        IF ( STRVIB .EQ. '  0000203 2 0 3' ) IDXVIB = 198
        IF ( STRVIB .EQ. '  0000203 2 0-1' ) IDXVIB = 199
        IF ( STRVIB .EQ. '  0000203-2 0 3' ) IDXVIB = 200
        IF ( STRVIB .EQ. '  0000210 0 1 0' ) IDXVIB = 201
        IF ( STRVIB .EQ. '  0000210 2 1 0' ) IDXVIB = 202
        IF ( STRVIB .EQ. '  0000211 0 1 1' ) IDXVIB = 203
        IF ( STRVIB .EQ. '  0000211 0 1-1' ) IDXVIB = 204
        IF ( STRVIB .EQ. '  0000211 2 1 1' ) IDXVIB = 205
        IF ( STRVIB .EQ. '  0000211 2 1-1' ) IDXVIB = 206
        IF ( STRVIB .EQ. '  0000300 1 0 0' ) IDXVIB = 207
        IF ( STRVIB .EQ. '  0000300 3 0 0' ) IDXVIB = 208
        IF ( STRVIB .EQ. '  0001000 0 0 0' ) IDXVIB = 209
        IF ( STRVIB .EQ. '  0001001 0 0 1' ) IDXVIB = 210
        IF ( STRVIB .EQ. '  0001002 0 0 0' ) IDXVIB = 211
        IF ( STRVIB .EQ. '  0001002 0 0 2' ) IDXVIB = 212
        IF ( STRVIB .EQ. '  0001003 0 0 1' ) IDXVIB = 213
        IF ( STRVIB .EQ. '  0001003 0 0 3' ) IDXVIB = 214
        IF ( STRVIB .EQ. '  0001004 0 0 0' ) IDXVIB = 215
        IF ( STRVIB .EQ. '  0001004 0 0 2' ) IDXVIB = 216
        IF ( STRVIB .EQ. '  0001004 0 0 4' ) IDXVIB = 217
        IF ( STRVIB .EQ. '  0001005 0 0 1' ) IDXVIB = 218
        IF ( STRVIB .EQ. '  0001005 0 0 3' ) IDXVIB = 219
        IF ( STRVIB .EQ. '  0001005 0 0 5' ) IDXVIB = 220
        IF ( STRVIB .EQ. '  0001010 0 1 0' ) IDXVIB = 221
        IF ( STRVIB .EQ. '  0001011 0 1 1' ) IDXVIB = 222
        IF ( STRVIB .EQ. '  0001011 0 1-1' ) IDXVIB = 223
        IF ( STRVIB .EQ. '  0001012 0 1 0' ) IDXVIB = 224
        IF ( STRVIB .EQ. '  0001012 0 1 2' ) IDXVIB = 225
        IF ( STRVIB .EQ. '  0001013 0 1 1' ) IDXVIB = 226
        IF ( STRVIB .EQ. '  0001013 0 1 3' ) IDXVIB = 227
        IF ( STRVIB .EQ. '  0001013 0 1-1' ) IDXVIB = 228
        IF ( STRVIB .EQ. '  0001020 0 0 0' ) IDXVIB = 229
        IF ( STRVIB .EQ. '  0001020 0 2 0' ) IDXVIB = 230
        IF ( STRVIB .EQ. '  0001021 0 0 1' ) IDXVIB = 231
        IF ( STRVIB .EQ. '  0001021 0 2 1' ) IDXVIB = 232
        IF ( STRVIB .EQ. '  0001021 0 2-1' ) IDXVIB = 233
        IF ( STRVIB .EQ. '  0001100 1 0 0' ) IDXVIB = 234
        IF ( STRVIB .EQ. '  0001101 1 0 1' ) IDXVIB = 235
        IF ( STRVIB .EQ. '  0001101 1 0-1' ) IDXVIB = 236
        IF ( STRVIB .EQ. '  0001102 1 0 0' ) IDXVIB = 237
        IF ( STRVIB .EQ. '  0001102 1 0 2' ) IDXVIB = 238
        IF ( STRVIB .EQ. '  0001102-1 0 2' ) IDXVIB = 239
        IF ( STRVIB .EQ. '  0001110 1 1 0' ) IDXVIB = 240
        IF ( STRVIB .EQ. '  0002000 0 0 0' ) IDXVIB = 241
        IF ( STRVIB .EQ. '  0002001 0 0 1' ) IDXVIB = 242
      ELSE IF ( IDXMOL .EQ. 47 ) THEN                               ! SO3
C new levels found in HITRAN2012 - arbitrarily assigned IDXVIB values
        IF ( STRVIB .EQ. ' 0 0 0 0 0 0A1''' ) IDXVIB = 1
        IF ( STRVIB .EQ. ' 0 0 0 0 1 1E'' ' ) IDXVIB = 2
        IF ( STRVIB .EQ. ' 0 0 0 0 2 0A1''' ) IDXVIB = 3
        IF ( STRVIB .EQ. ' 0 0 0 0 2 2E'' ' ) IDXVIB = 4
        IF ( STRVIB .EQ. ' 0 0 1 1 0 0E'' ' ) IDXVIB = 5
        IF ( STRVIB .EQ. ' 0 0 2 2 0 0E'' ' ) IDXVIB = 6
        IF ( STRVIB .EQ. ' 0 1 0 0 0 0A2"' ) IDXVIB = 7
        IF ( STRVIB .EQ. ' 0 1 0 0 1 1E" ' ) IDXVIB = 8
        IF ( STRVIB .EQ. ' 0 2 0 0 0 0A1''' ) IDXVIB = 9
        IF ( STRVIB .EQ. ' 1 0 0 0 0 0A1''' ) IDXVIB = 10
C 
C Monatomic molecules
      ELSE IF ( IDXMOL .EQ. 34 ) THEN                               ! O
        IF ( STRVIB .EQ. '              0' ) IDXVIB = 1
      ELSE
        WRITE ( *, * ) 'F-IDXVIB: Unrecognised value IDXMOL=', IDXMOL
        STOP
      END IF
C
C Check that level has been identified
      IF ( IDXVIB .EQ. 0 ) THEN
        IF ( NWRN(IDXMOL) .LT. MAXWRN ) THEN
          WRITE ( *, * ) 'W-IDXVIB: Unrecognised Vibrational Quanta=''', 
     &      STRVIB, ''' for IDXMOL=', IDXMOL
          NWRN(IDXMOL) = NWRN(IDXMOL) + 1
          IF ( NWRN(IDXMOL) .EQ. MAXWRN ) WRITE ( *, * )
     &      'Further such warnings for this molecule suppressed'
        END IF
        IDXVIB = 999
      END IF
C
      END
