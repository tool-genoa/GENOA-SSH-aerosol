      subroutine  ssh_creac(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization for gas-phase reactions.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     NR: number of reactions.
C     S: stoichiometric matrix
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'
      character *500 mot(nbmot)
      CHARACTER *12 NAM(5)
      dimension imot(nbmot),inam(5)



c
      nr=nr+1
c
      if (nr.gt.nrmax) then
         write(*,*)'ERROR: bad dimension for nr>nrmax.'
         stop 1
      endif

c     Check reactants
c
      molec(nr)=1
      nam(1)(1:imot(1)) = mot(1)(1:imot(1))
      inam(1)=imot(1)
      icurseur=1
 100  icurseur = icurseur+1
      if (mot(icurseur)(1:imot(icurseur)).eq.'+') then
         icurseur = icurseur+1
         molec(nr)=molec(nr)+1
         nam(molec(nr))(1:imot(icurseur))=
     &        mot(icurseur)(1:imot(icurseur))
         inam(molec(nr))=imot(icurseur)
         goto 100
      endif
c
      if (molec(nr).gt.3) then
         write(*,*)'ERROR: more than 3 reactants.'
         stop 1
      endif
c
      do JE = 1 ,molec(nr)
         indic=0
         do IE = 1 ,nesp(1)
            if (inam(je).eq.inom(ie))  then
               if (nAM(JE)(1:inam(je)).eq.NOM(IE)(1:inom(ie))) then
                  jer(je,nr)=ie
                  s(ie,nr) = s(ie,nr) -1.D0
                  indic=1
               endif
            endif
         enddo
c
         if (indic.eq.0) then
            write(6,*)'ERROR: the following reactant is unknown ',
     &           nam(je)(1:inam(je))
            stop 1
         endif
      enddo
C
      call ssh_WW90(nr,molec(nr),jer(1,nr),jer(2,nr),jer(3,nr))
      call ssh_DW90(nr,molec(nr),jer(1,nr),jer(2,nr),jer(3,nr))
c
c     Check products
c
      if ((mot(icurseur)(1:imot(icurseur)).ne.'>').and.
     &     (mot(icurseur)(1:imot(icurseur)).ne.'->')) then
         write(*,*)'ERROR: > or -> expected '
         stop 1
      endif
      icurseur = icurseur+1
c
 200  stoieir=1.D0
      if (icurseur.lt.nmot) then
         if ((mot(icurseur+1)(1:imot(icurseur+1)).ne.'+').and.
     &        (mot(icurseur+1)(1:imot(icurseur+1)).ne.'-')) then
            call ssh_reel(stoieir,mot(icurseur),imot(icurseur))
            icurseur = icurseur+1
         endif
      endif
      nam(1)(1:imot(icurseur))=
     &     mot(icurseur)(1:imot(icurseur))
c
      indic=0
      do IE = 1 ,nesp(1)
         if ((imot(icurseur).eq.inom(ie)).and.
     &        (NAM(1)(1:imot(icurseur)).eq.NOM(IE)(1:inom(ie)))) then
            s(ie,nr)=s(ie,nr)+stoIEIR
            indic=1
         endif
      enddo
!      if (indic.eq.0) then
!zhizhao        write(6,*)'WARNING: product unknown or reaction without product'
!zhizhao     &        ,nam(1)(1:imot(icurseur))!zhizhao
!      endif
c
      icurseur = icurseur+1
      if (icurseur.le.nmot) then
         if ((mot(icurseur)(1:imot(icurseur)).eq.'+').or.
     &        (mot(icurseur)(1:imot(icurseur)).eq.'-')) then
            icurseur = icurseur+1
            goto 200
         else
            write(*,*)'ERROR: + or - expected.'
            stop 1
         endif
      endif
c
c     Check gas/liquid
c
      do i=1,nesp(1)
         if ((s(i,nr).ne.0.).and.(indaq(i).ne.indaqr(nr)))
     &        then
            write(*,*)'ERROR: the phases are not coherent'
            write(*,*)'reaction ',nr,' species ',nom(i)
            stop 1
         endif
      enddo
      return
      end

C------------------------------------------------------------------------
      subroutine  ssh_kinreac(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization for kinetics of gas-phase reactions.
C     The different routines associated to the kinetic laws are
C     called.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     BP(.,NR): coefficient for kinetic rates.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'

      character *500 mot(nbmot)
      dimension imot(nbmot)
      dimension b(ntabphotmax)

c
c     Arrhenius' law
c     Gas-phase only
      i=2

 100  continue
      if (mot(i)(1:4).eq.'ARR1') then
         nb(nr)=1
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_WK190(nr,bp(1,nr))
      elseif (mot(i)(1:4).eq.'ARR2') then
         nb(nr)=2
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_WK290(nr,bp(1,nr),bp(2,nr))
      elseif (mot(i)(1:4).eq.'ARR3') then
         nb(nr)=3
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_WK390(nr,bp(1,nr),bp(2,nr),bp(3,nr))
      elseif (mot(i)(1:5).eq.'ARRC2') then
         nb(nr)=2
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         iprecalc(nr)=2
      elseif (mot(i)(1:5).eq.'ARRC3') then
         nb(nr)=3
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         iprecalc(nr)=3

c     PHOT: Photolysis.
      elseif (mot(i)(1:4).eq.'PHOT') then
         nb(nr)=10
         nrphot = nrphot + 1
         if (ntabphot.eq.0) then
            write(*,*)
     $           'ERROR: tabulated angles not defined for photolysis'
            write(*,*)'Needs to be defined by SET TABULATION ...'
            stop 1
         endif
c     Read according increasing or decreasing sequence
         if (ireversetab.eq.0) then
            do j=1,ntabphot
               call ssh_reel(bp(j,nr),mot(i+j),imot(i+j))
               b(j)=bp(j,nr)
            enddo
         else
            do j=1,ntabphot
               jj=ntabphot+1-j
               call ssh_reel(bp(jj,nr),mot(i+j),imot(i+j))
               b(jj)=bp(jj,nr)
            enddo
         endif
c
         call ssh_WPHOT90(nr,ntabphot,b,tabphot,0)

c     TROE: TROE/Fall-off for the general case.
      elseif (mot(i)(1:5).eq.'TROE4') then
         nb(nr)=4
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_WTROE90 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),0.6d0)

      elseif (mot(i)(1:5).eq.'TROE5') then
         nb(nr)=4
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_WTROE90 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),
     &        bp(5,nr))

      elseif (mot(i)(1:5).eq.'TROE7') then
         nb(nr)=4
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_reel(bp(6,nr),mot(i+6),imot(i+6))
         call ssh_reel(bp(7,nr),mot(i+7),imot(i+7))
         call ssh_WTROE790 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),
     &        bp(5,nr),bp(6,nr),bp(7,nr))

c     TROE12: TROE/Fall-off for MOCA.
      elseif (mot(i)(1:6).eq.'TROE10') then
         nb(nr)=4
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_reel(bp(6,nr),mot(i+6),imot(i+6))
         call ssh_reel(bp(7,nr),mot(i+7),imot(i+7))
         call ssh_reel(bp(8,nr),mot(i+8),imot(i+8))
         call ssh_reel(bp(9,nr),mot(i+9),imot(i+9))
         call ssh_reel(bp(10,nr),mot(i+10),imot(i+10))
         call ssh_WTROE1090 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),
     &        bp(5,nr),bp(6,nr),bp(7,nr),bp(8,nr),bp(9,nr),bp(10,nr))

c     CVAR/MOCA: temperature-dependent stoichiometry.
      elseif (mot(i)(1:5).eq.'CVAR') then
         nb(nr)=5
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_reel(bp(6,nr),mot(i+6),imot(i+6))
         call ssh_reel(bp(7,nr),mot(i+7),imot(i+7))
         call ssh_WCV90(nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr),
     &        bp(6,nr),bp(7,nr))

c     RCFE: Reactions Calculated From Equilibria.
      elseif (mot(i)(1:4).eq.'RCFE') then
         nb(nr)=8
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_reel(bp(6,nr),mot(i+6),imot(i+6))
         call ssh_WRCFE90 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),
     &        bp(5,nr),bp(6,nr))

C     MCM add for Complex rate coefficients in MCM v3.3.1
      elseif (mot(i)(1:4).eq.'MCM1') then
         nb(nr)=4
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_reel(bp(6,nr),mot(i+6),imot(i+6))
         call ssh_reel(bp(7,nr),mot(i+7),imot(i+7))
         call ssh_reel(bp(8,nr),mot(i+8),imot(i+8))
         call ssh_WKMCM1 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),
     &        bp(5,nr),bp(6,nr),bp(7,nr),bp(8,nr))

      elseif (mot(i)(1:4).eq.'MCM2') then
         nb(nr)=4
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_reel(bp(4,nr),mot(i+4),imot(i+4))
         call ssh_reel(bp(5,nr),mot(i+5),imot(i+5))
         call ssh_reel(bp(6,nr),mot(i+6),imot(i+6))
         call ssh_reel(bp(7,nr),mot(i+7),imot(i+7))
         call ssh_reel(bp(8,nr),mot(i+8),imot(i+8))
         call ssh_WKMCM2 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),
     &        bp(5,nr),bp(6,nr),bp(7,nr),bp(8,nr))

c     GENERAL: IUPAC from MCM photolysis genoa
      elseif (mot(i)(1:4).eq.'MCM3') then
         nb(nr)=10
         nrphot = nrphot + 1
         call ssh_reel(bp(1,nr),mot(i+1),imot(i+1))
         call ssh_reel(bp(2,nr),mot(i+2),imot(i+2))
         call ssh_reel(bp(3,nr),mot(i+3),imot(i+3))
         call ssh_WKMCM3 (nr,bp(1,nr),bp(2,nr),bp(3,nr))

c     SPEC: specific reactions.
      elseif (mot(i)(1:4).eq.'SPEC') then
         nb(nr)=5
         call ssh_entier(ispebp(nr),mot(i+1),imot(i+1))

c modif: YK(2010/02/15)
c!zhizhao         write(*,*) 'For some specific reactionconstants',mechanism_name
c         write(*,*) mechanism_name
c         if (chem_mechanism.eq.1) then
         if (mechanism_name .eq. "racm  ") then
            call ssh_WSPEC_RACM90 (nr,ispebp(nr))
C!zhizhao            write(*,*) '   expressions in RACM mechanism are used'
c         elseif (chem_mechanism.eq.2) then
         elseif (mechanism_name .eq. "cb05  " .or.
     &           mechanism_name .eq. "cb05en" .or.
     &           mechanism_name .eq. "cb05v0") then
            call ssh_WSPEC_CB0590 (nr,ispebp(nr))
C!zhizhao            write(*,*) '   expressions in CB05 mechanism are used'
c         elseif (chem_mechanism.eq.3) then
         elseif (mechanism_name .eq. "racm2 " .or.
     &           mechanism_name .eq. "racm2-2020" .or.
     &           mechanism_name .eq. "racm2-2020rad") then
            call ssh_WSPEC_RACM290 (nr,ispebp(nr))
C!zhizhao            write(*,*) '   expressions in RACM2 mechanism are used'
	 elseif (mechanism_name .eq. "MELCHIOR2") then
C!zhizhao            write(*,*) nr, ispebp(nr)
	    call ssh_WSPEC_MELCHIOR2 (nr, ispebp(nr))
C!zhizhao	    write(*,*) '   expressions in MELCHIOR2 mechanism are used'
         elseif (mechanism_name .eq. "MCM") then
            call ssh_WSPEC_MCM (nr,ispebp(nr))
C!zhizhao         else
C!zhizhao            write(*,*) nr, ispebp(nr)
C!zhizhao            write(*,*) 'Warning: specific reaction expression
C!zhizhao     &   is not given for ', mechanism_name
         endif

c end of modif

!!!!!!!!!! ADDED BY VICTOR2020 :
c     HETERO: heterogeneous reactions.
      elseif (mot(i)(1:6).eq.'HETERO') then
         nb(nr)=5
         call ssh_entier(ihetero(nr),mot(i+1),imot(i+1))
         call ssh_WHETERO90 (nr,ihetero(nr))
!!!!!!!!!! 

c     EXTRA: specific reaction with corrected factors
C     O3 -> 2. OH with corrected photolysis
      elseif (mot(i)(1:5).eq.'EXTRA') then
         nb(nr)=10
         if (ntabphot.eq.0) then
            write(*,*)'ERROR: tabulation not given for photolysis'
            write(*,*)'Needs to be defined by SET TABULATION ...'
            stop 1
         endif
c     Read according increasing or decreasing sequence
         if (ireversetab.eq.0) then
            do j=1,ntabphot
               call ssh_reel(bp(j,nr),mot(i+j),imot(i+j))
               b(j)=bp(j,nr)
            enddo
         else
            do j=1,ntabphot
               jj=ntabphot+1-j
               call ssh_reel(bp(jj,nr),mot(i+j),imot(i+j))
               b(jj)=bp(jj,nr)
            enddo
         endif
         call ssh_WPHOT90(nr,ntabphot,b,tabphot,1)

c     Third body.
      elseif (mot(i)(1:2).eq.'TB') then
         if (mot(i+1)(1:1).eq.'M') then
            ittb(nr)=1
         elseif (mot(i+1)(1:2).eq.'O2') then
            ittb(nr)=2
         elseif (mot(i+1)(1:2).eq.'N2') then
            ittb(nr)=3
         elseif (mot(i+1)(1:3).eq.'H2O') then
            ittb(nr)=4
         elseif ((mot(i+1)(1:2).eq.'H2').and.(imot(i+1).eq.2)) then
            ittb(nr)=5
         elseif (mot(i+1)(1:3).eq.'RO2') then !zhizhao add RO2 in TB
            ittb(nr)=6
         ! zhizhao add the possibility to add branch ratio (6 digits)
         else
            call ssh_reel(ratio,mot(i+1),imot(i+1))
            !read(mot(i+1)(1:6),*,IOSTAT=e) ratio
            !if (e.ne.0) print*,'can not read TB ratio',nr
            ittb(nr)=7
         !zhizhao
         !else
         !   write(*,*)'ERROR: syntax for Third Body'
         !   write(*,*)'M, O2, N2, H20 or H2 expected.'
         !   stop 1
         endif
         i=i+2
      else
         write(*,*)
     $        'ERROR: unknown syntax for KINETIC definition ',mot(2)
         stop 1
      endif

c     Modification BS/KS 21/05/2002
c     Case of a third body reaction: need for kinetics.

      if ((ittb(nr).ne.0).and.(nb(nr).eq.0)) goto 100
      if (ittb(nr).ne.0) call ssh_WTB90(nr,ittb(nr),ratio)!zhizhao

c     Update the chemical production term and the Jacobian matrix.
      call ssh_WFJ90(s,nr,jer)

c     Update the production and loss terms (P-Lc formulation)
      call ssh_WPL90(s,nr,jer)

c     Conversion mol/l -> molec/cm3.
      if (indaqr(nr).eq.1) then
         if(molec(nr).eq.2) bp(1,nr)=bp(1,nr)*1.0d3/av
         if(molec(nr).eq.3) bp(1,nr)=bp(1,nr)*1.0d6/av**2
      endif
c
      return
      end
C------------------------------------------------------------------------
      subroutine  ssh_cdis(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization for ionic dissociations.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     NR: number of reactions.
C     S: stoichiometric matrix
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'

      character *500 mot(nbmot)
      CHARACTER *12 NAM(5)
      dimension imot(nbmot),inam(5)
c
      nequil=nequil+1
      if (nequil.gt.nionx) then
         write(*,*)'ERROR: bad dimension for nionx=',nionx
         stop 1
      endif

c     Reactant.
      nam(1)(1:imot(1))=mot(1)(1:imot(1))
      inam(1)=imot(1)
c
      indic=0
      do IE = 1 ,nesp(1)
         if ((inam(1).eq.inom(ie)).and.
     &        (NAM(1)(1:inam(1)).eq.NOM(IE)(1:inom(ie)))) then
            indic=ie
            iesp(nequil,1) = ie
         endif
      enddo
      if (indic.eq.0) then
         write(6,*)'ERROR: unknown species ',nam(1)
         stop 1
      endif
      if (indaq(indic).eq.0) then
         write(6,*)'ERROR: this is a gaseous species ',nam(1)
         stop 1
      endif
c
      icurseur=1
      icurseur = icurseur+1
      if (mot(icurseur)(1:imot(icurseur)).ne.'=') then
         write(*,*)'ERROR: syntax = expected.'
         stop 1
      endif
c
      icurseur = icurseur+1

c     Products.
      nneq(nequil)=2
      inam(nneq(nequil))=imot(icurseur)
      nam(nneq(nequil))(1:imot(icurseur))=
     &     mot(icurseur)(1:imot(icurseur))
c
 200  stoeq=1.D0

c     BS 2003: to be checked !
      if (icurseur.lt.nmot) then
         if ((mot(icurseur+1)(1:imot(icurseur+1)).ne.'+').and.
     &        (mot(icurseur+1)(1:imot(icurseur+1)).ne.'-')) then
            call ssh_reel(stoeq,mot(icurseur),imot(icurseur))
            icurseur = icurseur+1
         endif
      endif
c
      nam(nneq(nequil))(1:imot(icurseur))=
     &     mot(icurseur)(1:imot(icurseur))
      inam(nneq(nequil))=imot(icurseur)

c
      write(*,*)'Ionic species: ',
     &     nam(nneq(nequil))(1:inam(nneq(nequil)))
c
      indic=0
      do IE = 1 ,nesp(1)
         if ((inam(nneq(nequil)).eq.inom(ie)).and.
     &        (NAM(nneq(nequil))(1:imot(icurseur)).eq.
     $        NOM(IE)(1:inom(ie))))
     &        then
            seqion(IE,nequil) = seqion(IE,nequil) + stoeq
            iesp(nequil,nneq(nequil)) = ie
            indic=1
            if (indaq(ie).eq.0) then
               write(*,*)'ERROR: gaseous species ',nom(ie),
     &              ' in dissociation ',nequil
               stop 1
            endif
         endif
      enddo
      if (indic.eq.0) then
         write(6,*)'ERROR: unknown species in dissociation ',
     &        nam(nneq(nequil))
         stop 1
      endif
c
      icurseur = icurseur+1
c
      if (icurseur.le.nmot) then
         if ((mot(icurseur)(1:imot(icurseur)).eq.'+').or.
     &        (mot(icurseur)(1:imot(icurseur)).eq.'-'))  then
            icurseur = icurseur + 1
            nneq(nequil) = nneq(nequil) + 1
            goto 200
         else
            write(*,*)'ERROR: syntax + expected.'
            stop 1
         endif
      endif

c     Check lumping for more than 4 ions.
      if (nneq(nequil).gt.3) then
         write(*,*)'ERROR: more than 2 products in dissociation.
     &        See routine "lump"'
         stop 1
      endif
c
      return
      end
C------------------------------------------------------------------------
      subroutine  ssh_kindis(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization of kinetics for ionic dissociations.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'

      character *500 mot(nbmot)
      dimension imot(nbmot)
c
      if (mot(2)(1:5).eq.'ARRC2') then
         call ssh_reel(xk1(nequil),mot(3),imot(3))
         call ssh_reel(xk2(nequil),mot(4),imot(4))
      else
         write(*,*)nom(2)(1:4),' ','ERROR: syntax.'
         stop 1
      endif
      return
      end
C------------------------------------------------------------------------
      subroutine  ssh_chenry(mot,imot,nmot,ieq)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization for Henry's equilibrium.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'
      character *12 nam(5)
      character *500 mot(nbmot)
      dimension imot(nbmot),inam(5)
c
c     Henry gas --> liq
c     Creation of Henry liq --> gas.

      nr=nr+1
c
      if (nr.gt.nrmax) then
         write(*,*)'ERROR: Bad dimension for nr>nrmax'
         stop 1
      endif

c     Reactant.

      molec(nr)=1
c
      nam(1)(1:imot(1)) = mot(1)(1:imot(1))
      inam(1)=imot(1)
      indic=0
      do IE = 1 ,nesp(1)
         if (inam(1).eq.inom(ie))  then
            if (nAM(1)(1:inam(1)).eq.nom(IE)(1:inom(ie))) then
               indr=ie
               jer(1,nr) = IE
               s(ie,nr) = -1.D0
               indic=1
            endif
         endif
      enddo
c
      if (indic.eq.0) then
         write(6,*)'ERROR: unknown species ',nam(1)
         stop 1
      endif

c     Product.

      nam(1)(1:imot(3)) = mot(3)(1:imot(3))
      inam(1)=imot(3)
      indic=0
      do IE = 1 ,nesp(1)
         if (inam(1).eq.inom(ie))  then
            if (nAM(1)(1:inam(1)).eq.nom(IE)(1:inom(ie))) then
               indp=ie
               s(ie,nr) = 1.D0
               indic=1
            endif
         endif
      enddo
c
      if (indic.eq.0) then
         write(6,*)'ERROR: unknown species ',nam(1)
         stop 1
      endif
      indaqr(nr)=indaq(jer(1,nr))
c
c     Creation reaction L->G
c
      nr=nr+1
      if (nr.gt.nrmax) then
         write(*,*)'ERROR: Bad dimension for nr>nrmax.'
         stop 1
      endif
c
      molec(nr)=1
c
      s(indp,nr)=-1.
      s(indr,nr)=1.D0
      jep(nr)=indr
      jer(1,nr)=indp
c
c     Lumping.
c
      iheq(nr-1)=ieq
      iheq(nr)=ieq
      j1=jer(1,nr)
      j2=jer(1,nr-1)
      ihreac(j1)=nr
      ihreac(j2)=nr-1
c
      if (indaq(jer(1,nr)).ne.1) then
         write(*,*)'ERROR: phases are not coherent for'
         write(*,*)'Henry ',nr
         write(*,*)'Species ',nom(i),' has to be aqueous.'
         stop 1
      endif
      if (indaq(jep(nr)).ne.0) then
         write(*,*)'ERROR phases are not coherent for'
         write(*,*)'Henry ',nr
         write(*,*)'Species ',nom(i),' has to be gaseous.'
         stop 1
      endif
      indaqr(nr)=1
c
      return
      end
C------------------------------------------------------------------------
      subroutine  ssh_kinhenry(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization for kinetics of Henry's equilibrium.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'

      character *500 mot(nbmot)
      dimension imot(nbmot)
c
c     Reaction G ->L
c
      nr=nr-1
c
      nb(nr)=6
      i=jer(1,nr)
c
      if (rmol(i).eq.0.) then
         ilettre=1
c
         do while (ilettre.le.inom(i))
            ncomp=0
            if (nom(i)(ilettre:ilettre).eq.'H') then
               ncomp=1
            elseif (nom(i)(ilettre:ilettre).eq.'C') then
               ncomp=12
            elseif (nom(i)(ilettre:ilettre).eq.'O') then
               ncomp=16
            elseif (nom(i)(ilettre:ilettre).eq.'N') then
               ncomp=14
            elseif (nom(i)(ilettre:ilettre).eq.'S') then
               ncomp=32
            else
               write(*,*)'ERROR: species ',nom(i),' unknown molar mass.'
               stop 1
            endif
            nn=1
            ilet=ilettre+1
            if (inom(i).ge.ilet) then
               if (nom(i)(ilet:ilet).eq.'2') nn=2
               if (nom(i)(ilet:ilet).eq.'3') nn=3
               if (nom(i)(ilet:ilet).eq.'4') nn=4
               if (nom(i)(ilet:ilet).eq.'5') nn=5
               if (nom(i)(ilet:ilet).eq.'6') nn=6
            endif
            rmol(i)=rmol(i)+ncomp*nn
            if (nn.gt.1) ilettre=ilettre+1
            ilettre=ilettre+1
         enddo
         write(*,*)'Computed molar mass for ',nom(i),'=',rmol(i)
      else
         write(*,*)'Read molar mass for ',nom(i),'=',rmol(i)
      endif
c
      if (rmol(i).eq.0.) then
         write(*,*)'ERROR: Species ',nom(i),' unknown molar mass.'
         stop 1
      endif
c
c     Reaction L->G
c
      nr=nr+1
      nb(nr)=7
      if (mot(2)(1:5).eq.'ARRC2') then
         call ssh_reel(bp(1,nr),mot(3),imot(3))
         call ssh_reel(bp(2,nr),mot(4),imot(4))
         call ssh_reel(bp(3,nr),mot(5),imot(5))
         iprecalc(nr)=2
      elseif (mot(2)(1:4).eq.'ARR1') then
         call ssh_reel(bp(1,nr),mot(3),imot(3))
         bp(2,nr)=0.D0
      elseif (mot(2)(1:4).eq.'ARR2') then
         call ssh_reel(bp(1,nr),mot(3),imot(3))
         call ssh_reel(bp(2,nr),mot(4),imot(4))
      else
         write(*,*)'ERROR: syntax kinetics Henry ',nr-1
         stop 1
      endif
c
      return
      end
C------------------------------------------------------------------------
      subroutine ssh_initphase
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization for multiphase description.
C     IPHASE = 1 : multiphase description.
C     IPHASE = 2 : gas-phase.
C     IPHASE = 3 : aqueous-phase.
C
C     IRMONODI: index reaction in one phase --> multiphase.
C
C     Notice that Henry's reactions are only to be taken into account for
C     the multiphase case.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      include 'ficcom'
c
      do ir=1,nr
c
         ip=indaqr(ir)+2
         nrp(1)=nrp(1)+1
         if ((nb(ir).ne.6).and.(nb(ir).ne.7)) then
            nrp(ip)=nrp(ip)+1
            irmonodi(nrp(ip),ip)=ir
         endif
c
         if (molec(ir).eq.1) then
            nrmol1(1)=nrmol1(1)+1
            imolec1(nrmol1(1),1)=nrp(1)
            if ((nb(ir).ne.6).and.(nb(ir).ne.7)) then
               nrmol1(ip)=nrmol1(ip)+1
               imolec1(nrmol1(ip),ip)=ir
            endif
c
         elseif (molec(ir).eq.2) then
            nrmol2(1)=nrmol2(1)+1
            imolec2(nrmol2(1),1)=nrp(1)
            nrmol2(ip)=nrmol2(ip)+1
            imolec2(nrmol2(ip),ip)=ir
c
         elseif (molec(ir).eq.3) then
            nrmol3(1)=nrmol3(1)+1
            imolec3(nrmol3(1),1)=nrp(1)
            nrmol3(ip)=nrmol3(ip)+1
            imolec3(nrmol3(ip),ip)=ir
         endif
c
         if (nb(ir).eq.1) then
            narr1(1)=narr1(1)+1
            iarr1(narr1(1),1)=nrp(1)
            narr1(ip)=narr1(ip)+1
            iarr1(narr1(ip),ip)=ir
         elseif (nb(ir).eq.2) then
            narr2(1)=narr2(1)+1
            iarr2(narr2(1),1)=nrp(1)
            narr2(ip)=narr2(ip)+1
            iarr2(narr2(ip),ip)=ir
         elseif (nb(ir).eq.3) then
            narr3(1)=narr3(1)+1
            iarr3(narr3(1),1)=nrp(1)
            narr3(ip)=narr3(ip)+1
            iarr3(narr3(ip),ip)=ir
         elseif (nb(ir).eq.4) then
            narr4(1)=narr4(1)+1
            iarr4(narr4(1),1)=nrp(1)
            narr4(ip)=narr4(ip)+1
            iarr4(narr4(ip),ip)=ir
         elseif (nb(ir).eq.5) then
            narr5(1)=narr5(1)+1
            iarr5(narr5(1),1)=nrp(1)
            narr5(ip)=narr5(ip)+1
            iarr5(narr5(ip),ip)=ir
         elseif (nb(ir).eq.6) then
            narr6(1)=narr6(1)+1
            iarr6(narr6(1),1)=nrp(1)
         elseif (nb(ir).eq.7) then
            narr7(1)=narr7(1)+1
            iarr7(narr7(1),1)=nrp(1)
         elseif (nb(ir).eq.8) then
            narr8(1)=narr8(1)+1
            iarr8(narr8(1),1)=nrp(1)
            narr8(ip)=narr8(ip)+1
            iarr8(narr8(ip),ip)=ir
         endif
c
         if (indaqr(ir).eq.1) then
            if (molec(ir).eq.2) then
               naq2(1)=naq2(1)+1
               inaq2(naq2(1),1)=ir
               naq2(3)=naq2(3)+1
               inaq2(naq2(3),1)=ir
            endif
         endif
c
         if (ittb(ir).ne.0) then
            nthird(1)=nthird(1)+1
            indthird(nthird(1),1)=nrp(1)
            nthird(ip)=nthird(ip)+1
            indthird(nthird(ip),ip)=ir
         endif
c
      enddo
c
c     write(*,*)' Total number of reactions         =',nrp(1)
C!zhizhao      write(*,*)' Number of gas-phase reactions     =',nrp(2)
C!zhizhao      write(*,*)' Number of aqueous-phase reactions =',nrp(3)
C!zhizhao      write(*,*)
C!zhizhao     &     ' Nbr of Henry reversible reactions =',nrp(1)-nrp(2)-nrp(3)
C!zhizhao      write(*,*)' Third body reactions              =',nthird(1)
C!zhizhao      write(*,*)'#############################################'
      return
      end

