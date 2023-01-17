C*************************************************************
C     Molecular descriptors
C       (1) for protein sequences
C       (2) protein-protein sequence pairs
C       (3) ligands
C       (4) protein sequence-ligand
C****************************************************************
C     input.dat: one or more protein sequences in Fasta format
C*************************************************************
C     input-ligand.sdf: one or more ligands in SDF format
      ! call main('input-protein.fasta')
      ! end
      subroutine run(protein_fasta_file, 
     & output_dir, output_dir_len,
     & config_dir, config_dir_len)
      implicit double precision(a-h,o-z)
      character, intent(in)::protein_fasta_file*100
      character, intent(in)::output_dir*100
      character, intent(in)::config_dir*300
      integer, intent(in)::output_dir_len
      integer, intent(in)::config_dir_len
      character seq*10000, name*100,line*101
      character molname*40
      common/aaindex/aadb(1000,20)
      dimension iad(1000),iappn(1000)
      dimension xtop(1000)
      dimension nacp1(1000), nacp2(1000)
      dimension idacp1(1000,30),idacp2(10000,30)
C
      open(unit=1,
     & file=config_dir(1:config_dir_len) // 'input-param.dat',
     & status='old')
      open(unit=3,
     & file=config_dir(1:config_dir_len) // 'input-aaindexdb.dat',
     & status='old')
      open(unit=6, 
     & file=output_dir(1:output_dir_len) // "output", 
     & status='unknown')
C
C     if output the information for each protein or not
C
      io6=1
C
C     Reading parameters controling the computation
C
      call readparam(ipp,iaac,idpc,iatd,nlag
     & ,imb,imoran,igeary,nid,iad,ictd,iqso,nphi,
     & ipaac,w1,lamda1,nset1,nacp1,idacp1,iapaac,w2,lamda2,
     & itop,ibcut,iapp,napp,iappn,methodpp,methodpl)
C
C
C    *************************Begining  for ligands***********************
C
      if(ipp.eq.1.or.ipp.eq.3)  goto 9
C
      open(unit=4, 
     & file=output_dir(1:output_dir_len) // 'output-ligand.dat',
     & status='unknown')
      open(unit=2, 
     & file=output_dir(1:output_dir_len) // 'output-ligand.nam', 
     & status='unknown')
      write(6,*)
      write(6,*) '*******************************'
      write(6,*) 'For ligands'
C
      open(unit=7,
     & file=config_dir(1:config_dir_len) // 'input-ligand.sdf',
     & status='old')
      nmol=0
1     continue
      call numbatom(natom,natoms,molname,iend,ierror)
      if(iend.eq.2) goto 9
      if(natom.eq.0) goto 1
C
      if(ierror.eq.0) then
            nmol=nmol+1
            write(6,*)
            write(6,*) 'nmol=',nmol,'    ', molname
C        write(6,*) 'number of atoms and heavy atoms', natom,natoms
C
        close(16)
        open(unit=16,
     & file=output_dir(1:output_dir_len) 
     & // 'output-temp/output-temp.dat',
     & status='unknown',err=10)
        goto 20
10      open(unit=16,
     & file=output_dir(1:output_dir_len)
     & //'output-temp.dat',
     & status='unknown')
20      continue
C
C       for ligands, ibcut=1
        ilbcut=1
C
        call topmol(natom,natoms,ilbcut,ierror)
        if(ierror.eq.1) goto 5
C
        if(io6.eq.0) rewind(6)
        call extractl(nmol,itop,xtop,ntop)
C    output on unit 9:output-des.dat
C
      write(4,4) molname
      write(4,*) ntop
      nn=0
2     nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.ntop) n2=ntop
      write(4,'(10E14.4)') (xtop(i),i=n1,n2)
      if(n2.lt.ntop) goto 2
      end if
C
5     continue
C
      if(iend.eq.2) goto 9
      if(iend.eq.0)  goto 1
C
4      format('>',a40)
9     continue
C
C ************End for ligands ***************************
C************Begining for protein sequences
C
      if(ipp.eq.2) goto 900
      open(unit=9, 
     & file=output_dir(1:output_dir_len) //'output-protein.dat', 
     & status='unknown')
      open(unit=10, 
     & file=output_dir(1:output_dir_len) //'output-protein.nam', 
     & status='unknown')
C
      write(6,*)
      write(6,*) '*******************************'
      write(6,*) 'For protein sequences'
C
      open(unit=5, file=protein_fasta_file, status='old')
      iend=0
      np=0
C
C     Reading the amino acid index provoded by the database 
C
        call readindex
C
C       reading sequence and calculating descriptors one by one
C
100   read(5,'(a101)',end=900) line
      if(line(1:1).ne.'>') goto 100  
C
200   name=line(2:101)
      np=np+1
      ln1=0
      ln2=0
300   read(5,'(a80)',end=700) line
      if(line(1:1).eq.';') goto 300
      if(line(1:1).eq.'>') goto 800
      call convert(line,ln)
      ln1=ln2+1
      ln2=ln2+ln
      seq(ln1:ln2)=line(1:ln)
      goto 300
700   iend=1
800   continue
      len=ln2
      if(io6.eq.0) rewind (6)
      write(6,*)
      write(6,*) '*************************************************'
      write(6,1001)np,name
      write(6,1002) len
      write(6,9001)seq(1:len)
      write(6,*)
C      nseg=iaac
      call des(seq,name,np,len,iaac,idpc,iatd,nlag,imb,imoran,igeary
     & ,nid,iad,ictd,iqso,nphi,ipaac,w1,lamda1,nset1,nacp1,idacp1
     & ,iapaac,w2,lamda2,itop,ibcut,iapp,napp,iappn,nv)
C
C
      if(iend.ne.1) goto 200
900   continue
C
C  ************    end for protein sequences*****************************
C
C ******************* begining for protein -protein interaction********************
C
      if(ipp.eq.3) then
            close(19)
            open(unit=19,
     & file=output_dir(1:output_dir_len) // 'output-ppi.dat',
     & status='unknown')
            close (8)
            open(unit=8,
     & file=config_dir(1:config_dir_len) // 'input-ppi.dat',
     & status='old')
C
            maxv=nv
C
            write(6,*)
            write(6,*) '*******************************'
            write(6,*) 'For protein-protein interaction'
C
            call ppides(maxv,methodpp)
      end if
C
C ******************* End  for protein -protein interaction*****************
C ******************* Begining for protein-ligand interaction********************
C
      if(ipp.eq.4) then
            close(15)
            open(unit=15,
     & file=output_dir(1:output_dir_len) // 'output-pli.dat',
     & status='unknown')
            close(14)
            open(unit=14,
     & file=config_dir(1:config_dir_len) // 'input-pli.dat',
     & status='old')
C
            maxp=nv
            maxl=ntop
C
            write(6,*)
            write(6,*) '*******************************'
            write(6,*) 'For protein-ligand interaction'
C
            call plides(maxp,maxl,methodpl)
C
      end if
C
C ******************* End  for protein-ligand interaction
C
1001  format('Sequence',1x,i5,2x,'>',a100)
1002  format('Length',1x,i6)
9001  format(1x,a)
      close(1)
      close(3)
      close(4)
      close(5)
      close(7)
      close(9)
      return
      end
C******************************************************
C    Convert character from lower case to upper case
C    and discard meaningless characters
C****************************************************
      subroutine convert(line,ln)
      implicit double precision(a-h,o-z)
      dimension ide(6)
      character line*80
C    
C     ASCII code for B,J,O,U,X,Z: the characters not 
C     appearing in the list of the standard amino acids' characters
C
      data ide/66,74,79,85,88,90/
C
C     convert lower case to upper case,eg. c to C
      ln=0
      do 100 i=1,80
      m=ichar(line(i:i))
      if(m.ge.97.and.m.le.122) then
        m=m-32
      end if
      if(m.ge.65.and.m.le.90) then
        is=1
        do k=1,6
        if(m.eq.ide(k)) is=0
        end do
        if(is.eq.1) then
          ln=ln+1
          line(ln:ln)=char(m)
        end if
      end if
100   continue
C
      return
      end
C***********************************************
C      Main programs for Descriptors calculation
C************************************************
      subroutine des(seq,name,np,len,iaac,idpc,iatd,nlag,imb
     &,imoran,igeary
     & ,nid,iad,ictd,iqso,nphi,ipaac,w1,lamda1,nset1,nacp1,idacp1
     & ,iapaac,w2,lamda2,itop,ibcut,iapp,napp,iappn,nv)
C
      implicit double precision(a-h,o-z)
      character seq*10000,name*100
      character str*20
      character line*80
      dimension v(10000),ip(len),n(20)
      dimension freq(20),dp(20,20,idpc)
      dimension xaacd(20,iaac)
      dimension xc(10000),xctd(1000,ictd)
      dimension xqso(1000,iqso),xqa(1000,iapaac)
      dimension xq(1000),xqs(1000,ipaac)
      dimension xtop(1000*itop)
      dimension iad(1000),iappn(1000)
      dimension xapp(napp,iapp)
      dimension nacp1(1000)
      dimension idacp1(1000,30)
      dimension idh(1000)
C
      str='ACDEFGHIKLMNPQRSTVWY'
C
      if(len.eq.0) then
       write(6,*) 'The sequence length is 0'
      end if
C     
C     Position in str for each amino in the sequence
C
      nv=0
      do 100 i=1,len
      ip(i)=index(str,seq(i:i))
100   continue      
      if(iaac.eq.0) goto 350
C
C     (1a) Amino acid composition*100 when iaac=1; (1b)  Amino acid composition distribution 
C         or sequence-segmented amino acid composition,when iaac>1, iaac is the number of segnemt
C
      call aacd(str,seq,len,iaac,ip,xaacd)
C
      do 360 i=1,20
       do 370 j=1,iaac
          nv=nv+1
          v(nv)=xaacd(i,j)
370   continue
360   continue
C     Feature names
      if(np.ne.1) goto 350  
      do 390  i=1,20
        do 391  j=1,iaac
          if(iaac.gt.1) then
            line(1:6)='AACD( '
          else
            line(1:6)=' AAC( '
          end if
          line(7:7)=str(i:i)
          line(8:9)=' )'
          if(iaac.gt.1) then 
             line(10:21)=' in segment '
             write(line(22:24),'(i3)') j
          else
             line(10:24)= ' '
          end if
          write(10,'(a24)') line(1:24)
391   continue
390   continue
C
350    continue
C
      if(idpc.eq.0) goto 550
C
C    (2)  Dipeptide composition*100
C
      call dipeptide(str,seq,len,ip,dp,idpc)
C
      do 405 k=1,idpc
      do 400 i=1,20
      do 500 j=1,20
      nv=nv+1
      v(nv)=dp(i,j,k)
500   continue
400   continue
405   continue
C
C     Feature names
      if(np.ne.1) goto 430
      do 440 k=1,idpc
      do 410 i=1,20
      do 420 j=1,20
      line(1:5)='DPC( '
      line(6:6)=str(i:i)
      line(7:7)=','
      line(8:8)=str(j:j)
      line(9:10)=' )'
      if(idpc.gt.1) then
        line(11:21)='in segment '
        write(line(22:24),'(i3)') k
      else
        line(11:24)=' '
      end if
      write(10,'(a24)') line(1:24)
420   continue
410   continue
440   continue
C
430   continue
C
C
 550  continue
C
      if(iatd.eq.0) goto 650
C
C     (3) Autocorrelation descriptors or Sequence-segmented autocorrelation
C
      nlag1=nlag
C
      call autodes(seq,len,np,nlag1,imb,imoran,igeary,nid,iad,xc,nx)
C
      do 600 i=1,nx
      nv=nv+1
      v(nv)=xc(i)
600   continue
C
 650  continue
C
      if(ictd.eq.0)  goto 750
C
C    (4)  composition,transition, distribution(CTD) or sequence-segmented CTD
C 
      call ctd(seq,len,np,ictd,xctd,nx)
C
      do 710 j=1,ictd
      do 700 i=1,nx
      nv=nv+1
      v(nv)=xctd(i,j)
700   continue
710   continue
C
 750  continue
C
      if(iqso.eq.0) goto 850
C
C     (5) Quasi-sequence-order(QSO) descriptors or sequence-segmented QSO
C
      nlag2=20+nphi
C
      call qsodes(seq,len,np,nlag2,iqso,xqso,nqso)
C
C
      do 810 j=1,iqso
      do 800 i=1,nqso
      nv=nv+1
      v(nv)=xqso(i,iqso)
800   continue
810   continue
C
 850  continue
C
      if(ipaac.eq.0) goto 950
C
C    (6) Pseudo-amino acid composition
C
      do 900 nac=1,nset1
      nid=nacp1(nac)
      do 901 k=1,nid
        idh(k)=idacp1(nac,k)
901   continue
C
C
      call spaac(seq,len,np,w1,lamda1,xqs,nq,ipaac,nid,idh)
C
      do 910 i=1,nq
       do 920  j=1,ipaac
      nv=nv+1
      v(nv)=xqs(i,j)
920   continue
910   continue
C
900   continue
C
950   continue
C
      if(iapaac.eq.0) goto 1050
C
C    (G7) amphiliphilic Pseudo-amino acid composition (APAAC) or sequence-segmented APAAC
C
      call apaac(seq,len,np,w2,lamda2,iapaac,xqa,nqa)
C
      do 1001 j=1,iapaac
      do 1000 i=1,nqa
      nv=nv+1
      v(nv)=xqa(i,j)
1000   continue
1001   continue
C
1050  continue
C
1150  continue
C
      if(itop.eq.0) goto 1350
C
C (G8)   Sequence-segmented topological descriptors
C
C    saved directly on unit 16
C
      close(16)
      open(unit=16,file='output-temp/output-temp.dat'
     & ,status='unknown',err=1151)
      goto 1152
1151  open(unit=16,file='output-temp.dat',status='unknown')
1152  continue
C
      call stopdes(seq,len,np,itop,ibcut)
C
C     Extract and output topological descriptors on unit 6 and descriptors name on unit 10
C
C
      call extractp(np,itop,xtop,ntop)
C
      do 1300 i=1,ntop
      nv=nv+1
      v(nv)=xtop(i)
1300  continue
C
1350   continue
C     
C
      if(iapp.eq.0) goto 1450
C   
C (G9) total amino acid properties
C
      call appdes(seq,len,np,iapp,napp,iappn,xapp)
C
      do 1410 i=1,napp
      do 1420 j=1,iapp
      nv=nv+1
      v(nv)=xapp(i,j)
1420   continue
1410   continue
C
1450  continue
C
C     End of the calculation of features for one sequence
C
C    output on unit 9:output-des.dat
C
C      write(9,10001) name
      write(9,*) ">" // name
      write(9,*) nv
      nn=0
10000 nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.nv) n2=nv
      write(9,'(10E14.4)') (v(i),i=n1,n2)
      if(n2.lt.nv) goto 10000
10001  format('>',a100)
      return
      end
C*********************************************************
C     Amino acid compositions
C------------------------------------
C     Refer to: Manoj Bhasin and Gajendra P.S.Raghava
C     "Classification of Nuclear Receptors Based on Amino Acid
C     Composition and dipide composition".
C     J.Bio.Chem. 2004,279,23262.
C************************************************************
      subroutine frequency(str,seq,len,ip,freq)
      implicit double precision(a-h,o-z)
      character seq*10000,str*20
      dimension ip(len),n(20),freq(20)
      do 100 i=1,20
      n(i)=0
100   continue
      do 200 i=1,len
      ii=ip(i)
      n(ii)=n(ii)+1
200   continue   
      do 300 i=1,20
      freq(i)=100.0d0*dble(n(i))/dble(len)
300   continue
      write(6,1001) 
      write(6,1002) (str(i:i),i=1,10)
      write(6,1003) (freq(i),i=1,10)
      write(6,*)
      write(6,1002) (str(i:i),i=11,20)
      write(6,1003) (freq(i),i=11,20)
C
1001  format('(G1) Amino acid composition(%)')
1002  format(10(4x,A1,5x))
1003  format(10f10.6)
      return
      end
C*********************************************************
C     Dipeptide compositions
C----------------------------------
C     Refer to: Manoj Bhasin and Gajendra P.S.Raghava
C     "Classification of Nuclear Receptors Based on Amino Acid
C     Composition and dipide composition".
C     J.Bio.Chem. 2004,279,23262.
C************************************************************
      subroutine dipeptide(str,seq,len,ip,dp,nseg)
      implicit double precision(a-h,o-z)
      character seq*10000,str*20
      dimension ip(len),n(20,20),dp(20,20,nseg)
      dimension lenseg(nseg)
C
      do 10 k=1,nseg
        do 20 i=1,20
           do 25 j=1,20
            dp(i,j,k)=0.0d0
25         continue
20      continue
10    continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
C      write(6,*) 'Length for segments'
C      write(6,*) (lenseg(i),i=1,nseg)
C
      write(6,*)
      if(nseg.gt.1) then 
         write(6,*) '(G2) Dipeptide Composition Distribution (DPCD)'
      else
          write(6,*) '(G2) Dipeptide Composition (DPC)'
      end if
C
      n2=0
      do 50 ns=1,nseg
      n1=n2+1
      n2=n1+lenseg(ns)-1
C
      do 100 i=1,20
      do 150 j=1,20
      n(i,j)=0
150   continue
100   continue
C
      do 200 i=n1,n2-1
      ii=ip(i)
      jj=ip(i+1)
      n(ii,jj)=n(ii,jj)+1
200   continue   
C
      do 300 i=1,20
      do 400 j=1,20
      dp(i,j,ns)=100.0d0*dble(n(i,j))/dble(lenseg(ns)-1)
400   continue
300   continue
      write(6,*)
      if(nseg.gt.1) then
         write(6,*) 'DPCD in segment ',ns
         write(6,*) seq(n1:n2)
         write(6,*)
      end if
      write(6,1002) (str(j:j),j=1,10)
      do 500 i=1,20
      write(6,1003) str(i:i),(dp(i,j,ns),j=1,10)
500   continue
      write(6,*)
      write(6,1002) (str(j:j),j=11,20)
      do 600 i=1,20
      write(6,1003) str(i:i),(dp(i,j,ns),j=11,20)
600   continue
50    continue
C
1002  format('Residue',3x,10(4x,A1,5x))
1003  format(4x,A1,2x,10f10.6)
      return
      end
C************************************************************************************
C     Autocorrelation descriptors
C----------------------------------
C     1. Autocorrelation
C     Definition: ac(n)=sum ( p(i)*p(i+n))/(N-n)
C     where p(i) is property of amino acid at position i 
C     in the sequence and this property must be standerized. 
C     All properties used here are amino acid induces from
C     AAinbdex1.
C     2. Moran Auotcorrelation
C
C       (1). Hydrophobicity scales: hs(i)
C            refer to: 
C            Cid,H.,Bunster,M.,Canales,M. and Gazitua,F
C            'Hydrophobicity and structural classes in proteins'
C            Protein Engeering 5, 373(1992)
C            source: AAINDEX1
C       (2). Average flexibility indices(Bhaskaran-Ponnuswamy,1988):xfi(i)
C            refer to:
C            Bhaskaran,R. and Ponnuswany,P.K
C            Positional flexibilities of amino acid residues in globual proteins
C            Int.J.Peptide Protein Res.1988,32,241-235.
C            source: AAINDEX1
C       (3). Polarizability parameter (Charton-Charton,1982):xpp(i)
C            refer to:
C            Charton,M.,and Charton, B.I.
C            The Structure Dependence of Amino Acid Hydrophobicity parameters
C            J.Theor.Biol.1982,99,629-644.
C       (4). Free energy of solution in water(Charton-Charton,1982):xfe(i)
C            refer to:
C            Charton M. and Charton,B.I.
C            The Structure dependence of amino acid hydrophobicity parameters
C            J.Theor.Biol.1982,99,629-644.
C       (5). Residue accessible surface area in trepeptide(Chothia,1976):xas(i)
C            refer to:
C            Chothia,C.
C            The Nature of the Accessible and buried surface in proteins
C            J.Mol.Biol.1976,105,1-14.
C        (6). Residue volume(bigelow,1967):xrv(i)
C             refer to:
C             Bigelow,C.C.
C             On the average hydrophobicity of proteins and the relation between
C             it and protein structure
C             J.Theor.biol.1967,16,187.
C        (7). Steric Parameter(Charton,1981):xsp(i)
C             Charton,M.
C             Protein folding and the genetic code:An alterative quantitative model
C             J.Theor.Biol.1981,91,115-123.
C        (8). Relative mutability (dayhoff et.al.,1978b):xrm(i)
C             Dayhoff, M.O.,, Schwartz,R.M., and Orcutt,B.C.
C             "Altas of Protein Sequence and Structure",Vol.5,Suppl.3(Dayhoff,M.O..,
C             ed.), National Biomedicinal Research Foundation, Washington,D.C.pp345,352. 
C
C*************************************************************************************
      subroutine autodes(seq,len,np,nlag,imb,imoran,igeary
     & ,nid,iad,xc,nx)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20,line*80
      character*80 title
      common/aaindex/aadb(1000,20)
      dimension xc(10000)
      dimension acmb(nlag,imb)
      dimension acmoran(nlag,imoran)
      dimension acgeary(nlag,igeary)
      dimension xindex(20)
      dimension iad(1000)
C
      aa='ARNDCQEGHILKMFPSTWYV'
C
      nx=0
C
C--------1. Autocorrelation------
C
      if(imb.eq.0) goto 200
      write(6,*)
      write(6,*) '(G3) Autocorrelation descriptors'
      write(6,*)
      do 101 k=1,nid
      iadk=iad(k)
      write(6,*) 'M-B autocorrelation by  amino acid index ',iadk
      do 102 j=1,20
       xindex(j)=aadb(iadk,j)
102   continue
C
      call automb(seq,len,imb,nlag,aa,xindex,acmb)
C
      do 104 m=1,imb
      do 103 i=1,nlag
      nx=nx+1
      xc(nx)=acmb(i,m)
      if(np.eq.1) write(10,*) 'MB-ATC by amino acid index ',iadk,
     & ' with lag ',i,' in segment ',m
103   continue
104   continue
101   continue
C
C   
200   continue
      if(imoran.eq.0) goto 300
C
C--------2. Moran Autocorrelation------
C
      write(6,*)
      write(6,*) 'Moran autocorrelation'
      write(6,*)
      do 201 k=1,nid
      iadk=iad(k)
      write(6,*) 'Moran autocorrelation by amino acid index ', iadk
      do 202 j=1,20
       xindex(j)=aadb(iadk,j)
202   continue
C
      call automoran(seq,len,imoran,nlag,aa,xindex,acmoran)
      do 204 m=1,imoran
      do 203 i=1,nlag
      nx=nx+1
      xc(nx)=acmoran(i,m)
      if(np.eq.1) write(10,*) 'Moran-ATC by amino acid index ',iadk,
     &  ' with lag ',i,' in segment ',m
203   continue
204   continue
201   continue
C
300   continue
      if(igeary.eq.0) goto 400
C
C--------3. Geary Autocorrelation------
C
      write(6,*)
      write(6,*) 'Geary autocorrelation descriptors'
      write(6,*)
      do 301 k=1,nid
      iadk=iad(k)
      write(6,*) 'Geary autocorrelation by  amino acid index ',iadk
      do 302 j=1,20
       xindex(j)=aadb(iadk,j)
302   continue
C
      call autogeary(seq,len,igeary,nlag,aa,xindex,acgeary)
      do 304 m=1,igeary
      do 303 i=1,nlag
      nx=nx+1
      xc(nx)=acgeary(i,m)
      if(np.eq.1) write(10,*) 'Geary-ATC by amino acid index ',iadk,
     &  ' with lag ',i,' in segment ',m
303   continue
304   continue
301   continue
C
400   continue
C
      return
      end
C*********************************************
C      M-B Aurocorrelation with property p
C**********************************************
      subroutine automb(seq,len,nseg,nlag,aa,p,ac)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20
      dimension ac(nlag,nseg),p(20)
      dimension ip(len)
      dimension lenseg(nseg)
C
      do 10 i=1,nlag
        do 20 j=1,nseg
        ac(i,j)=0.0d0
20      continue
10    continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
C     First standerized the properties
C
      av=0.0d0
      do 100 i=1,20
      av=av+p(i)
100   continue
      av=av/20.0d0
      s2=0.0d0
      do 200 i=1,20
      s2=s2+(p(i)-av)**2
200   continue
      s=dsqrt(s2/20.0d0)
      do 300 i=1,20
      p(i)=(p(i)-av)/s
300   continue
      do 400 i=1,len
      ip(i)=index(aa,seq(i:i))
400   continue      
C
      do 410 ns=1,nseg
       if(lenseg(ns).lt.nlag) then
          write(6,*) 'Error: Length of segment < maximum of the '
     &              ,'lag of the autocorrealtion'
         write(6,*) ' Decrease nlag or nseg' 
         return
       end if
410   continue
C
      n2=0
      do 50 ns=1,nseg
      n1=n2+1
      n2=n1+lenseg(ns)-1
C
C     Identification the position of each amino acid of the 
C     sequence in the sequence segement
C
C     Then autocorrelation
C
      do 600 i=1,nlag
      ac(i,ns)=0.0d0
      if(i.ge.lenseg(ns)) goto 600
      do 700 k=n1,n2-i
      i1=ip(k)
      i2=ip(k+i)
      ac(i,ns)=ac(i,ns)+p(i1)*p(i2)
700   continue
      ac(i,ns)=ac(i,ns)/dble(lenseg(ns)-i)
600   continue
C
      write(6,*) ' In segment ', ns
      write(6,*) seq(n1:n2)
      nn=0
800   nn=nn+1
      np1=(nn-1)*10+1
      np2=nn*10
      if(np2.gt.nlag) np2=nlag
      write(6,1001) (ac(i,ns),i=np1,np2)
      if(np2.lt.nlag) goto 800
50    continue
C
1001  format(10f10.6)
      return
      end
C*****************************************
C      Moran aurocorrelation with property p
C******************************************
      subroutine automoran(seq,len,nseg,nlag,aa,p,ac)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20
      dimension ac(nlag,nseg),p(20)
      dimension ip(len)
      dimension lenseg(nseg)
C
      do 10 i=1,nlag
        do 20 j=1,nseg
        ac(i,j)=0.0d0
20      continue
10    continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
      do 410 ns=1,nseg
C
       if(lenseg(ns).lt.nlag) then
          write(6,*) 'Error: Length of segment < maximum of the '
     &              ,'lag of the autocorrealtion'
         write(6,*) ' Decrease nlag or nseg' 
         return
       end if
410   continue
C
C     First standerized the properties
C
      av=0.0d0
      do 100 i=1,20
      av=av+p(i)
100   continue
      av=av/20.0d0
      s2=0.0d0
      do 200 i=1,20
      s2=s2+(p(i)-av)**2
200   continue
      s=dsqrt(s2/20.0d0)
      do 300 i=1,20
      p(i)=(p(i)-av)/s
300   continue
      do 400 i=1,len
      ip(i)=index(aa,seq(i:i))
400   continue      
C
C
      n2=0
      do 50 ns=1,nseg
      n1=n2+1
      n2=n1+lenseg(ns)-1
C
C     Then moran autocorrelation
C
C     average value of the property over the segment
C
      avc=0.0d0
      do 600 i=n1,n2
      ii=ip(i)
      avc=avc+p(ii)
600   continue
      avc=avc/dble(lenseg(ns))
C
      do 700 i=1,nlag
      if(i.ge.len) then
      ac(i,ns)=0.0d0
      goto 700
      end if
      ac1=0.0d0
      do 800 k=n1,n2-i
      i1=ip(k)
      i2=ip(k+i)
      ac1=ac1+(p(i1)-avc)*(p(i2)-avc)
800   continue
      ac1=ac1/dble(lenseg(ns)-i)
      ac2=0.0d0
      do 900 k=n1,n2
      ik=ip(k)
      ac2=ac2+(p(ik)-avc)*(p(ik)-avc)
900   continue
      ac2=ac2/dble(lenseg(ns))
      if(ac2.eq.0.0d0) then
       ac(i,ns)=1.0d0
      else
        ac(i,ns)=ac1/ac2
      end if
700   continue
C
      write(6,*) ' In segment ', ns
      write(6,*) seq(n1:n2)
      nn=0
1000  nn=nn+1
      np1=(nn-1)*10+1
      np2=nn*10
      if(np2.gt.nlag) np2=nlag
      write(6,1001) (ac(i,ns),i=np1,np2)
      if(np2.lt.nlag) goto 1000
50    continue
C
1001  format(10f10.6)
      return
      end
C********************************************
C      Geary aurocorrelation with property p
C*********************************************
      subroutine autogeary(seq,len,nseg,nlag,aa,p,ac)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20
      dimension ac(nlag,nseg),p(20)
      dimension ip(len)
      dimension lenseg(nseg)
C
      do 10 i=1,nlag
        do 20 j=1,nseg
        ac(i,j)=0.0d0
20      continue
10    continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
      do 410 ns=1,nseg
       if(lenseg(ns).lt.nlag) then
          write(6,*) 'Error: Length of segment < maximum of the '
     &              ,'lag of the autocorrealtion'
         write(6,*) ' Decrease nlag or nseg' 
         return
       end if
410   continue
C
C     First standerized the properties
C
      av=0.0d0
      do 100 i=1,20
      av=av+p(i)
100   continue
      av=av/20.0d0
      s2=0.0d0
      do 200 i=1,20
      s2=s2+(p(i)-av)**2
200   continue
      s=dsqrt(s2/20.0d0)
      do 300 i=1,20
      p(i)=(p(i)-av)/s
300   continue
      do 400 i=1,len
      ip(i)=index(aa,seq(i:i))
400   continue      
C
C
      n2=0
      do 50 ns=1,nseg
      n1=n2+1
      n2=n1+lenseg(ns)-1
      av=0.0d0
C
C     Then Geary autocorrelation
C
C     average value of the property over the segment
C
      avc=0.0d0
      do 600 i=n1,n2
      ii=ip(i)
      avc=avc+p(ii)
600   continue
      avc=avc/dble(lenseg(ns))
C
      do 700 i=1,nlag
      if(i.ge.len) then
      ac(i,ns)=0.0d0
      goto 700
      end if
      ac1=0.0d0
      do 800 k=n1,n2-i
      i1=ip(k)
      i2=ip(k+i)
      ac1=ac1+(p(i1)-p(i2))*(p(i1)-p(i2))
800   continue
      ac1=ac1/(2.0d0*dble(lenseg(ns)-i))
      ac2=0.0d0
      do 900 k=n1,n2
      ik=ip(k)
      ac2=ac2+(p(ik)-avc)*(p(ik)-avc)
900   continue
      ac2=ac2/dble(lenseg(ns)-1)
      if(ac2.eq.0.0d0) then
       ac(i,ns)=1.0d0
      else
        ac(i,ns)=ac1/ac2
      end if
700   continue
C
      write(6,*) ' In segment ', ns
      write(6,*) seq(n1:n2)
      nn=0
1000  nn=nn+1
      np1=(nn-1)*10+1
      np2=nn*10
      if(np2.gt.nlag) np2=nlag
      write(6,1001) (ac(i,ns),i=np1,np2)
      if(np2.lt.nlag) goto 1000
50    continue
C
1001  format(10f10.6)
      return
      end
C********************************************************************
C     Composition, transition, distribution
C-----------------------------------------
C     refer to:
C     (1).Inna Dubchak, Ilya Muchink,Stephen R.Holbrook and Sun-Hou Kim
C         Prediction of protein folding class using global description 
C         of amino acid sequence
C         Proc.Natl.Acad.Sci.USA,1995,92,8700-8704.
C     (2).Inna Dubchak, Ilya Muchink,Christopher Mayor,Igor Dralyuk and 
C         Sun-Hou Kim
C         Recognition of a Protein fold in the context of the SCOP
C         Classification
C         Proteins:1999,35,401-407.
C*********************************************************************
      subroutine ctd(seq,len,np,nseg,x,nx)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20,name*80,seqseg*10000
C
      dimension id(20)
      character namenew(100)*80
      dimension idnew(100,20)
C
C
      dimension x(1000,nseg),idhy(20),idv(20),idp1(20),idp2(20)
      dimension idch(20),idss(20),idsa(20)
C
C     adding a new set of 17 properties according to Dr. Chen Yuzong
C     2014,9,15 idnew(17,20)
C
C
      dimension xc(100)
      dimension lenseg(nseg)
C
      aa='ACDEFGHIKLMNPQRSTVWY'
C
C     Hydrophobicity class index
C
      data idhy /2,3,1,1,3,2,2,3,1,3,3,1,2,1,1,2,2,3,3,2/
C
C     Normalized van der Waals volume class index
C
      data idv /1,1,1,2,3,1,3,2,3,2,3,2,1,2,3,1,1,2,3,3/
C
C     Polarity class index
c
      data idp1/2,1,3,3,1,2,3,1,3,1,1,3,2,3,3,2,2,1,1,1/
c
C     Polarizability class index
C
      data idp2/1,2,1,2,3,1,3,2,3,2,3,2,2,2,3,1,1,2,3,3/
C
C     Charge class index from Cai Congzhong
C
      data idch/2,2,3,3,2,2,2,2,1,2,2,2,2,2,1,2,2,2,2,2/
C
C    Secondary Structure Class index from Cai Congzhong
C    
      data idss/1,2,3,1,2,3,1,2,1,1,1,3,3,1,1,3,2,2,2,2/
C
C    Solvent Accessible Class index from Cai Congzhong
C
      data idsa/1,1,2,2,1,1,3,1,2,1,3,2,3,2,2,3,3,1,1,3/
C
C     added a new set of 17 properties suggested by Dr.Chen Yuzong, 2014-9-14. 
C
C (1) Surface tension
      namenew(1)=' Surface tension'
      data (idnew(1,i),i=1,20)/1,2,1,2,3,1,1,3,2,3,3,1,3,1,1,2,2,3,3,3/
C (2)  Moleculaw weight
      namenew(2)=' Moleculaw weight'
      data (idnew(2,i),i=1,20)/1,2,2,2,3,1,2,2,2,2,2,2,2,2,3,1,2,2,3,3/
C  (3) solubility in water
      namenew(3)=' solubility in water'
      data (idnew(3,i),i=1,20)/1,1,3,2,2,1,2,2,1,2,2,2,2,2,1,2,1,2,2,3/
C (4) No of hydrogen bond donnor in side chain
      namenew(4)=' No.of hydrogen bond donnor in side chain'
      data (idnew(4,i),i=1,20)/3,3,2,2,3,3,1,3,1,3,3,1,3,1,1,2,2,3,2,2/
C (5) No of hydrogen bond acceptor in side chain
      namenew(5)=' No. of hydrogen bond acceptor in side chain'
      data (idnew(5,i),i=1,20)/3,3,1,1,3,3,1,3,2,3,3,1,3,1,1,2,2,3,2,2/
C (6) Clogp
      namenew(6)=' CLogP'
      data (idnew(6,i),i=1,20)/2,2,1,1,3,2,1,3,1,3,3,1,2,1,1,2,2,2,3,2/
C (7) Amino acid flexibility index 
      namenew(7)=' Amino acid flexibility index'
      data (idnew(7,i),i=1,20)/2,3,2,1,3,1,2,2,1,3,3,1,2,1,2,1,2,2,3,3/
C  (8) Protein-protein Interface hotspot propensity-Bogan
      namenew(8)=' Protein-protein Interface hotspot propensity-Bogan'
      data (idnew(8,i),i=1,20)/2,3,1,2,2,2,1,1,1,3,2,1,1,2,1,2,2,3,1,1/
C  (9) Protein-protein Interface (PPI) propensity-Ma
      namenew(9)=' Protein-protein Interface (PPI) propensity-Ma'
      data (idnew(9,i),i=1,20)/2,1,1,3,1,2,2,3,3,2,1,2,1,1,1,2,2,2,1,1/
C  (10) Protein-DNA Interface propensity-Schneider
      namenew(10)=' Protein-DNA Interface propensity-Schneider'
      data (idnew(10,i),i=1,20)/2,3,2,2,2,1,2,2,1,2,3,1,3,1,1,1,1,2,2,1/
C  (11) Protein-DNA Interface  propensity-Ahmad
      namenew(11)=' Protein-DNA Interface  propensity-Ahmad'
      data (idnew(11,i),i=1,20)/2,3,2,2,2,1,1,2,1,3,3,1,2,1,1,1,1,2,2,1/
C  (12) Protein-RNA Interface propensity-Kim
      namenew(12)=' Protein-RNA Interface propensity-Kim'
      data (idnew(12,i),i=1,20)/3,3,3,3,2,2,1,2,1,2,1,2,2,2,1,2,3,2,2,1/
C  (13) Protein-RNA Interface propensity-Kim
      namenew(13)=' Protein-RNA Interface propensity-Kim'
      data (idnew(13,i),i=1,20)/2,3,3,3,2,1,1,2,1,3,1,2,2,2,1,1,2,3,1,1/ 
C  (14) Protein-RNA Interface propensity-Phipps
      namenew(14)=' Protein-RNA Interface propensity-Phipps'
      data (idnew(14,i),i=1,20)/2,3,2,2,2,2,1,3,1,2,1,2,2,1,1,1,3,2,3,2/
C  (15) Protein-ligand binding site propensity-Khazanov
      namenew(15)=' Protein-ligand binding site propensity-Khazanov'
      data (idnew(15,i),i=1,20)/3,1,3,3,1,2,1,2,3,2,2,2,3,3,2,2,2,3,1,1/
C  (16) Protein-ligand valid binding site propensity-Khazanov
      namenew(16)=' Protein-ligand valid binding site propen-Khazanov'
      data (idnew(16,i),i=1,20)/3,1,2,3,1,2,1,2,3,2,1,2,3,3,3,2,2,2,1,1/
C  (17) propensity for Protein-ligand polar and aromatic non-bonded interactions-Imai
      namenew(17)=' propensity for Protein-ligand polar and arom-Imai'
      data (idnew(17,i),i=1,20)/3,2,1,1,2,3,1,3,2,3,2,2,3,2,1,2,2,3,2,1/
C
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
C
      n2=0
      do 50 ns=1,nseg
      nx=0
      n1=n2+1
      n2=n1+lenseg(ns)-1
      seqseg=seq(n1:n2)
      lens=lenseg(ns)
C
      write(6,*)
      if(nseg.gt.1) then 
        write(6,*) '(G4) Composition transition and distribution (CTD)'
        write(6,*) ' In segment ',ns
        write(6,*) seq(n1:n2)
      else
        write(6,*) '(G4) Composition transition and distribution (CTD)'
      end if
      write(6,*)
      write(6,*)'CTD according to hydrophobicity'
C
C
      name='  Hydrophobicity  '
      call ctdx(seqseg,lens,np,xc,nc,idhy,name,ns)
      do 100 i=1,nc
      nx=nx+1
      x(nx,ns)=xc(i)
100   continue
C
      write(6,*)
      write(6,*)'CTD according to normalized vdW volumes'
C
      name='   vdW Volumes    '
      call ctdx(seqseg,lens,np,xc,nc,idv,name,ns)
      do 200 i=1,nc
      nx=nx+1
      x(nx,ns)=xc(i)
200   continue
C
      write(6,*)
      write(6,*)'CTD according to polarity'
C
      name='Polarizability 1'
      call ctdx(seq,len,np,xc,nc,idp1,name,ns)
      do 300 i=1,nc
      nx=nx+1
      x(nx,ns)=xc(i)
300   continue
C
C
      write(6,*)
      write(6,*)'CTD according to polarizability'
C
      name='Polarizability 2'
      call ctdx(seqseg,lens,np,xc,nc,idp2,name,ns)
      do 400 i=1,nc
      nx=nx+1
      x(nx,ns)=xc(i)
400   continue
C
	write(6,*)
      write(6,*)'CTD according to charge'
C
        name='Charge'
      call ctdx(seqseg,lens,np,xc,nc,idch,name,ns)
	do 500 i=1,nc
        nx=nx+1
	x(nx,ns)=xc(i)	
500   continue
c
	write(6,*)
      write(6,*)'CTD according to secondary structure'
C
        name='Secondary Structure'
      call ctdx(seqseg,lens,np,xc,nc,idss,name,ns)
	do 600 i=1,nc
        nx=nx+1
	x(nx,ns)=xc(i)	
600   continue
c
	write(6,*)
      write(6,*)'CTD according to solvent accessibility'
C
        name='Solvent Accessible'
      call ctdx(seq,len,np,xc,nc,idsa,name,ns)
	do 700 i=1,nc
        nx=nx+1
	x(nx,ns)=xc(i)	
700   continue
c
C     added a new set of 17 properties suggested by Dr.Chen Yuzong, 2014-9-14. 
C
      do 800 iid=1,17
      write(6,*)
      write(6,*)'CTD according to'//namenew(iid)
C
      name=namenew(iid)
      do 810 j=1,20
      id(j)=idnew(iid,j)
810   con tinue
C
      call ctdx(seq,len,np,xc,nc,id,name,ns)
	do 900 i=1,nc
        nx=nx+1
	x(nx,ns)=xc(i)	
900   continue
800   continue
C
C
50    continue
c
      return
      end
C*****************************************************************************
C     Composition,transition, distribution according to property class id
C---------------------
C     Seq: sequence
C     len: length of the sequence
C     xc(i): the descriptors to be required
C     id(20): the amino acid type number for the 20 amino acid in aa
C     nc(i),i=1,2,3: number of the three types of amino acids in the sequence
C     nt(i,j): the number of transition from type i to type j
C     ix(i): the type number for amino acid in position i of the sequence
C     ip(i,j): the position in the sequence for the jth amino acid with type i 
C*******************************************************************************
      subroutine ctdx(seq,len,np,xc,nx,id,name,ns)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20,name*80
      dimension xc(100),id(20),nc(3),nt(3,3),ix(len)
      dimension ip(3,len),nd(3,5)
      aa='ACDEFGHIKLMNPQRSTVWY'
C
      do 100 i=1,len
      do 200 j=1,20
      if(seq(i:i).eq.aa(j:j)) then
      ix(i)=id(j)
      end if
200   continue
100   continue
C
C      write(6,*) 'Sequence coding'
C      write(6,'(10000i1)'),(ix(i),i=1,len)
C
C     composition
C
      do 300 i=1,3
      nc(i)=0
300   continue
C
      do 350 i=1,3
      do 350 j=1,len
      ip(i,j)=0
350   continue
C
      do 400 i=1,len
      ik=ix(i)
      nc(ik)=nc(ik)+1
      ip(ik,nc(ik))=i
400   continue
      do 500 i=1,3
      xc(i)=100.0d0*dble(nc(i))/dble(len)
500   continue
      write(6,*) 'Composition'
      write(6,'(3(2x,f10.6))') (xc(i),i=1,3)
      if(np.eq.1) write(10,*)'Composition for class 1 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)'Composition for class 2 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)'Composition for class 3 by '//name
     & ,' in segment ',ns
C
C     Transition
C
      do 600 i=1,3
      do 700 j=1,3
      nt(i,j)=0
700   continue
600   continue 
      do 800 i=1,len-1
      ix1=ix(i)
      ix2=ix(i+1)
      nt(ix1,ix2)=nt(ix1,ix2)+1
800   continue
      nt12=nt(1,2)+nt(2,1)
      nt13=nt(1,3)+nt(3,1)
      nt23=nt(2,3)+nt(3,2)
      xc(4)=100.0d0*dble(nt12)/dble(len-1)
      xc(5)=100.0d0*dble(nt13)/dble(len-1)
      xc(6)=100.0d0*dble(nt23)/dble(len-1)
      write(6,*) 'Transition'
      write(6,'(3(2x,f10.6))') (xc(i),i=4,6)
      if(np.eq.1) write(10,*)'Transition 12 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)'Transition 13 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)'Transition 23 by '//name
     & ,' in segment ',ns
C
C     Distribution
C
C
C     For first amino acid with the specifided type
C
C     write(6,*) 'Distribution'
      xc(7)=100.0d0*dble(ip(1,1))/dble(len)
      xc(8)=100.0d0*dble(ip(2,1))/dble(len)
      xc(9)=100.0d0*dble(ip(3,1))/dble(len)
      write(6,*) 'First amino acid type distribution'
      write(6,'(3(2x,f10.6))') (xc(i),i=7,9)
      if(np.eq.1) write(10,*)
     & 'First amino acid distribution for class 1 by '//name
     & ,'in segment ',ns
      if(np.eq.1) write(10,*)
     & 'First amino acid distribution for class 2 by '//name
     & ,'in segment ',ns
      if(np.eq.1) write(10,*)
     & 'First amino acid distribution for class 3 by '//name
     & ,'in segment ',ns
C
C     for the 25% amino acid with the specifided type
C
      n1=int(dble(nc(1))*0.25d0)
      n2=int(dble(nc(2))*0.25d0)
      n3=int(dble(nc(3))*0.25d0)
C
C     For too short peptides
C
      if(n1.lt.1) n1=1
      if(n2.lt.1) n2=1
      if(n3.lt.1) n3=1
C
      xc(10)=100*dble(ip(1,n1))/dble(len)
      xc(11)=100*dble(ip(2,n2))/dble(len)
      xc(12)=100*dble(ip(3,n3))/dble(len)
      write(6,*) '25% amino acid type distributioni'
      write(6,'(3(2x,f10.6))') (xc(i),i=10,12)
      if(np.eq.1) write(10,*)
     & '25% amino acid distribution for class 1 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '25% amino acid distribution for class 2 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '25% amino acid distribution for class 3 by '//name
     & ,' in segment ',ns
C
C     for the 50% amino acid with the specifided type
C
      n1=int(dble(nc(1))*0.50d0)
      n2=int(dble(nc(2))*0.50d0)
      n3=int(dble(nc(3))*0.50d0)
C
C     For too short peptides
C
      if(n1.lt.1)n1=1
      if(n2.lt.1)n2=1
      if(n3.lt.1)n3=1
C
      xc(13)=100*dble(ip(1,n1))/dble(len)
      xc(14)=100*dble(ip(2,n2))/dble(len)
      xc(15)=100*dble(ip(3,n3))/dble(len)
      write(6,*) '50% amino acid type distribution'
      write(6,'(3(2x,f10.6))') (xc(i),i=13,15)
      if(np.eq.1) write(10,*)
     & '50% amino acid distribution for class 1 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '50% amino acid distribution for class 2 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '50% amino acid distribution for class 3 by '//name
     & ,' in segment ',ns
C
C     for the 75% amino acid with the specifided type
C
      n1=int(dble(nc(1))*0.75d0)
      n2=int(dble(nc(2))*0.75d0)
      n3=int(dble(nc(3))*0.75d0)
C
      if(n1.lt.1)n1=1
      if(n2.lt.1)n2=1
      if(n3.lt.1)n3=1
C
      xc(16)=100*dble(ip(1,n1))/dble(len)
      xc(17)=100*dble(ip(2,n2))/dble(len)
      xc(18)=100*dble(ip(3,n3))/dble(len)
      write(6,*) '75% amino acid type distribution'
      write(6,'(3(2x,f10.6))') (xc(i),i=16,18)
      if(np.eq.1) write(10,*)
     & '75% amino acid distribution for class 1 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '75% amino acid distribution for class 2 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '75% amino acid distribution for class 3 by '//name
     & ,' in segment ',ns
C
C     for the 100% amino acid with the specifided type
C
      n1=nc(1)
      n2=nc(2)
      n3=nc(3)
C
      if(n1.lt.1)n1=1
      if(n2.lt.1)n2=1
      if(n3.lt.1)n3=1
C
      xc(19)=100*dble(ip(1,n1))/dble(len)
      xc(20)=100*dble(ip(2,n2))/dble(len)
      xc(21)=100*dble(ip(3,n3))/dble(len)
      write(6,*) '100% amino acid type distribution'
      write(6,'(3(2x,f10.6))') (xc(i),i=19,21)
      if(np.eq.1) write(10,*)
     & '100% amino acid distribution for class 1 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '100% amino acid distribution for class 2 by '//name
     & ,' in segment ',ns
      if(np.eq.1) write(10,*)
     & '100% amino acid distribution for class 3 by '//name
     & ,' in segment ',ns
C
      nx=21
C
      return
      end
C*****************************************************
C     Quasi-sequence-order features
C------------------------
C     refer to:
C     Kuo-Chen Chou
C     Prediction of protein Subcellar Locations by 
C     Incorpating Quasi-Sequence-Order Effect
C---------------
C     Amino Acid distance used here
C
C     (1).Normalized amino acid distance(Schneider and Wrede):xsw(i,j)
C         Gisbert Schneider and Paul Wrede
C        "The rational Design of Amino Acid Sequences by Artificial
C         Neural Networks and Simulated Molecular Evolution: De Novo
C         Design of an Idealized Leader Peptidase Cleavage site"
C         Biophysical Journal,1994,66-335
C`          
C     (2).Chemical distance(Grantham,1974):xgr(i,j)
C         Grantham,R.
C         Amino acid difference formula to help explain
C         protein evolution
C         Science,1974,185,862-864
C
C*****************************************************
      subroutine qsodes(seq,len,np,nlag,nseg,xq,nq)
      implicit double precision(a-h,o-z)
      character seq*10000,aasw*20,aagr*20,name*30
      character seqseg*10000
      dimension xq(1000,nseg),xd(1000)
      dimension xgr(20,20),xsw(20,20)
      dimension lenseg(nseg)
C     (1)
      data (xsw(1,i),i=1,20)/
     &     0.000,0.112,0.819,0.827,0.540,0.208,0.696,
     &     0.407,0.891,0.406,0.379,0.318,0.191,0.372,
     &     1.000,0.094,0.220,0.273,0.739,0.552/
      data (xsw(2,i),i=1,20)/
     &     0.114,0.000,0.847,0.838,0.437,0.320,0.660,
     &     0.304,0.887,0.301,0.277,0.324,0.157,0.341,
     & 1.000,0.176,0.233,0.167,0.639,0.457/
      data (xsw(3,i),i=1,20)/
     &     0.729,0.742,0.000,0.124,0.924,0.697,0.435,
     &     0.847,0.249,0.841,0.819,0.560,0.657,0.584,
     &     0.295,0.667,0.649,0.797,1.000,0.836/
      data (xsw(4,i),i=1,20)/
     &     0.790,0.788,0.133,0.000,0.932,0.779,0.406,
     &     0.860,0.143,0.854,0.830,0.599,0.688,0.598,
     &     0.234,0.726,0.682,0.824,1.000,0.837/
      data (xsw(5,i),i=1,20)/
     &     0.508,0.405,0.977,0.918,0.000,0.690,0.663,
     &     0.128,0.903,0.131,0.169,0.541,0.420,0.459,
     &     1.000,0.548,0.499,0.252,0.207,0.179/
      data (xsw(6,i),i=1,20)/
     &     0.206,0.312,0.776,0.807,0.727,0.000,0.769,
     &     0.592,0.894,0.591,0.557,0.381,0.323,0.467,
     &     1.000,0.158,0.272,0.464,0.923,0.728/
      data (xsw(7,i),i=1,20)/
     &     0.896,0.836,0.629,0.547,0.907,1.000,0.000,
     &     0.848,0.566,0.842,0.825,0.754,0.777,0.716,
     &     0.697,0.865,0.834,0.831,0.981,0.821/
      data (xsw(8,i),i=1,20)/
     &     0.403,0.296,0.942,0.891,0.134,0.592,0.652,
     &     0.000,0.892,0.013,0.057,0.457,0.311,0.383,
     &     1.000,0.443,0.396,0.133,0.339,0.213/
      data (xsw(9,i),i=1,20)/
     &     0.889,0.871,0.279,0.149,0.957,0.900,0.438,
     &     0.899,0.000,0.892,0.871,0.667,0.757,0.639,
     &     0.154,0.825,0.759,0.882,1.000,0.848/
      data (xsw(10,i),i=1,20)/
     &     0.405,0.296,0.944,0.892,0.139,0.596,0.653,
     &     0.013,0.893,0.000,0.062,0.452,0.309,0.376,
     &     1.000,0.443,0.397,0.133,0.341,0.205/
      data (xsw(11,i),i=1,20)/
     &     0.383,0.276,0.932,0.879,0.182,0.569,0.648,
     &     0.058,0.884,0.062,0.000,0.447,0.285,0.372,
     &     1.000,0.417,0.358,0.120,0.391,0.255/
      data (xsw(12,i),i=1,20)/
     &     0.424,0.425,0.838,0.835,0.766,0.512,0.780,
     &     0.615,0.891,0.603,0.588,0.000,0.266,0.175,
     &     1.000,0.361,0.368,0.503,0.945,0.641/
      data (xsw(13,i),i=1,20)/
     &     0.220,0.179,0.852,0.831,0.515,0.376,0.695,
     &     0.363,0.875,0.357,0.326,0.231,0.000,0.228,
     &     1.000,0.196,0.161,0.244,0.720,0.481/
      data (xsw(14,i),i=1,20)/
     &     0.512,0.462,0.903,0.861,0.671,0.648,0.765,
     &     0.532,0.881,0.518,0.505,0.181,0.272,0.000,
     &     1.000,0.461,0.389,0.464,0.831,0.522/
      data (xsw(15,i),i=1,20)/
     &     0.919,0.905,0.305,0.225,0.977,0.928,0.498,
     &     0.929,0.141,0.920,0.908,0.690,0.796,0.668,
     &     0.000,0.860,0.808,0.914,1.000,0.859/
      data (xsw(16,i),i=1,20)/
     &     0.100,0.185,0.801,0.812,0.622,0.170,0.718,
     &     0.478,0.883,0.474,0.440,0.289,0.181,0.358,
     &     1.000,0.000,0.174,0.342,0.827,0.615/
      data (xsw(17,i),i=1,20)/
     &     0.251,0.261,0.830,0.812,0.604,0.312,0.737,
     &     0.455,0.866,0.453,0.403,0.315,0.159,0.322,
     &     1.000,0.185,0.000,0.345,0.816,0.596/
      data (xsw(18,i),i=1,20)/
     &     0.275,0.165,0.900,0.867,0.269,0.471,0.649,
     &     0.135,0.889,0.134,0.120,0.380,0.212,0.339,
     &     1.000,0.322,0.305,0.000,0.472,0.310/
      data (xsw(19,i),i=1,20)/
     &     0.658,0.560,1.000,0.931,0.196,0.829,0.678,
     &     0.305,0.892,0.304,0.344,0.631,0.555,0.538,
     &     0.968,0.689,0.638,0.418,0.000,0.204/
      data (xsw(20,i),i=1,20)/
     &     0.587,0.478,1.000,0.932,0.202,0.782,0.678,
     &     0.230,0.904,0.219,0.268,0.512,0.444,0.404,
     &     0.995,0.612,0.557,0.328,0.244,0.000/
C     (2)
      data xgr(2,1)/112.0/
      data (xgr(3,i),i=1,2)  /111.0, 86.0/
      data (xgr(4,i),i=1,3)  /126.0, 96.0, 23.0/
      data (xgr(5,i),i=1,4)  /195.0,180.0,139.0,154.0/
      data (xgr(6,i),i=1,5)  / 91.0, 43.0, 46.0, 61.0,154.0/
      data (xgr(7,i),i=1,6)  /107.0, 54.0, 42.0, 45.0,170.0,
     &       29.0/
      data (xgr(8,i),i=1,7)  / 60.0,125.0, 80.0, 94.0,159.0,
     &       87.0, 98.0/
      data (xgr(9,i),i=1,8)  / 86.0, 29.0, 68.0, 81.0,174.0,
     &       24.0, 40.0, 98.0/
      data (xgr(10,i),i=1,9) / 94.0, 97.0,149.0,168.0,198.0,
     &      109.0,134.0,135.0, 94.0/
      data (xgr(11,i),i=1,10)/ 96.0,102.0,153.0,172.0,198.0,
     &      113.0,138.0,138.0, 99.0,  5.0/
      data (xgr(12,i),i=1,11)/106.0, 26.0, 94.0,101.0,202.0,
     &       53.0, 56.0,127.0, 32.0,102.0,107.0/
      data (xgr(13,i),i=1,12)/ 84.0, 91.0,142.0,160.0,196.0,
     &      101.0,126.0,127.0, 87.0, 10.0, 15.0, 95.0/
      data (xgr(14,i),i=1,13)/113.0, 97.0,158.0,177.0,205.0,
     &      116.0,140.0,153.0,100.0, 21.0, 22.0,102.0, 28.0/
      data (xgr(15,i),i=1,14)/ 27.0,103.0, 91.0,108.0,169.0,
     &       76.0, 93.0, 42.0, 77.0, 95.0, 98.0,103.0, 87.0,
     &      114.0/
      data (xgr(16,i),i=1,15)/ 99.0,110.0, 46.0, 65.0,112.0,
     &       68.0, 80.0, 56.0, 89.0,142.0,145.0,121.0,135.0,
     &      155.0, 74.0/
      data (xgr(17,i),i=1,16)/ 58.0, 71.0, 65.0, 85.0,149.0,
     &       42.0, 65.0, 59.0, 47.0, 89.0, 92.0, 78.0, 81.0,
     &      103.0, 38.0, 58.0/
      data (xgr(18,i),i=1,17)/148.0,101.0,174.0,181.0,215.0,
     &      130.0,152.0,184.0,115.0, 61.0, 61.0,110.0, 67.0,
     &       40.0,147.0,177.0,128.0/
      data (xgr(19,i),i=1,18)/112.0, 77.0,143.0,160.0,194.0,
     &       99.0,122.0,147.0, 83.0, 33.0, 36.0, 85.0, 36.0,
     &       22.0,110.0,144.0, 92.0, 37.0/
      data (xgr(20,i),i=1,19)/ 64.0, 96.0,133.0,152.0,192.0,
     &       96.0,121.0,109.0, 84.0, 29.0, 32.0, 97.0, 21.0,
     &       50.0, 68.0,124.0, 69.0, 88.0, 55.0/
C
      aagr='ARNDCQEGHILKMFPSTWYV'
      aasw='ACDEFGHIKLMNPQRSTVWY'
C
      do 150 i=1,20
      xgr(i,i)=0.0d0
150   continue
      do 200 i=1,19
      do 300 j=i+1,20
      xgr(i,j)=xgr(j,i)
300   continue
200   continue
      iprint=0
C
C     normalization of the Grantham Chemical distance matrix
C
      dmax=xgr(1,1)
      do 400 i=1,20
      do 500 j=1,20
      if(dmax.lt.xgr(i,j)) dmax=xgr(i,j)
500   continue
400   continue
      do 600 i=1,20
      do 700 j=1,20
      xgr(i,j)=xgr(i,j)/dmax
700   continue
600   continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
C
      n2=0
      do 50 ns=1,nseg
      nx=0
      n1=n2+1
      n2=n1+lenseg(ns)-1
      seqseg=seq(n1:n2)
      lens=lenseg(ns)
      nq=0
C
      write(6,*)
      write(6,*) '(5) Quasi-sequence-order descriptors'
      if(nseg.gt.1) then
        write(6,*) 'In segment',ns
        write(6,*) seq(n1:n2)
      end if
C
C     (1).Using Schneider-Wrede distance
C
      write(6,*)
      write(6,*) 'Quasi-sequence-order descriptors using ',
     & 'Schneider-Wrede distance'
C
      name='Scwneider-Wrede distance'
      call qsod(seqseg,lens,np,name,nlag,xsw,aasw,ns,xd,nd)
C
      do 100 i=1,nd
      xq(nq+i,ns)=xd(i)
100   continue
C   
      nq=nq+nd
C
C     (2).Using Normalized Grantham chemical distance
C
      write(6,*)
      write(6,*) 'Quasi-sequence-order descriptors using ',
     & 'normalized Grantham chemical distance'
C
C
      name='Grantham-Chemical distance'
      call qsod(seqseg,lens,np,name,nlag,xgr,aagr,ns,xd,nd)
C
      do 800 i=1,nd
      xq(nq+i,ns)=xd(i)
800   continue
C
      nq=nq+nd
C
50    continue
      return 
      end
C**********************************************************
C     Quasi-sequence-order features using dstance matrix d
C***********************************************************
      subroutine qsod(seq,len,np,name,nlag,d,aa,ns,xd,nd)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20,name*30
      dimension xd(1000),d(20,20)
      dimension ix(len)
      dimension t(1000),f(20),na(20)
C
      do 100 i=1,len
      ix(i)=index(aa,seq(i:i))
100   continue
C
      n=nlag-20
C
      do 300 j=1,n
      t(j)=0.0d0
      if(j.ge.len) goto 300
      do 400 i=1,len-j
      ii=ix(i)
      jj=ix(i+j)
      t(j)=t(j)+d(ii,jj)*d(ii,jj)
400   continue
300   continue
C
      do 500 i=1,20
      na(i)=0
500   continue
      do 600 i=1,len
      ii=ix(i)
      na(ii)=na(ii)+1
600   continue
      do 700 i=1,20
      f(i)=dble(na(i))/dble(len)
700   continue
C
      do 750 i=1,nlag-20
      xd(i)=t(i)
750   continue
C
      write(6,*) 'sequence-order-coupling numbers'
      nn=0
780   nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.(nlag-20)) n2=nlag-20
      write(6,'(10f15.6)') (xd(i),i=n1,n2)
      if(n2.lt.(nlag-20)) goto 780 
C
      do 760 i=1,nlag-20
      if(np.eq.1) write(10,*)
     & 'Sequence-order-coupling by '//name//'of order ',i
760   continue
C
      nd=nlag-20
C
      denom2=0.0d0
      do 800 j=1,n
      denom2=denom2+t(j)
800   continue
C
      w=0.10d0
C
      do 900 i=1,20
      xd(nd+i)=f(i)/(1.0d0+w*denom2)
900   continue
      do 1000 i=21,nlag
      xd(nd+i)=w*t(i-20)/(1.0d0+w*denom2)
1000  continue
C
      write(6,*) 'Quasi-sequence-order descriptors'
      np1=nd
      np2=nd+nlag
      nn=0
1100  nn=nn+1
      n1=(nn-1)*10+np1+1
      n2=nn*10+np1
      if(n2.gt.np2) n2=np2
      write(6,'(10f15.6)') (xd(i),i=n1,n2)
      if(n2.lt.np2) goto 1100
C
      do 1200 i=1,20
      if(np.eq.1) write(10,*)
     & 'QSO Xr by '//name//'for amino acid ',aa(i:i),' in segment ',ns
1200  continue
      do 1300 i=21,nlag
      if(np.eq.1) write(10,*)
     & 'QSO Xd by '//name//'with lag d= ', i, ' in segment ',ns
1300  continue
C
      nd=nd+nlag
C
      return
      end
C
      subroutine readparam(ipp,iaac,idpc,iatd,nlag
     & ,imb,imoran,igeary,nid,iad,ictd,iqso,nphi,
     & ipaac,w1,lamda1,nset1,nacp1,idacp1,iapaac,w2,lamda2,
     & itop,ibcut,iapp,napp,iappn,methodpp,methodpl)
C
      implicit double precision(a-h,o-z)
      character*80 title
      dimension iad(1000),iappn(1000)
      dimension nacp1(1000), nacp2(1000)
      dimension idacp1(1000,30),idacp2(10000,30)
C
C   .....(1)
C
C     ipp :=1: Descriptors for individual protein sequence; =2:  Descriptors for ligands;
C          =3: Descriptors for protein-protein interaction;=4: for protein-ligand intractions 
C
      read(1,'(a80)') title
      read(1,*) ipp
C
C......(2)
C
C iaac: Amino acid composition*100
      read(1,'(a80)') title
      read(1,*) iaac
C
C........(3) 
C
C idpc: Dipeptide composition*100
      read(1,'(a80)') title
      read(1,*) idpc
C
C........(4)
C
C iatd: Autocorrelation descriptors
      read(1,'(a80)') title
      read(1,*) iatd,nlag
      read(1,*) imb,imoran,igeary
      read(1,*,err=300) nid,(iad(i),i=1,nid)
      goto 310
300   write(6,*) 'Error in reading amino acid index for'
      write(6,*) 'autocorrelation,then, default values are used.'
      nid=8
      do 320 i=1,8
      iad(i)=i
320   continue
310   continue
C
C........(5)
C
c ictd:composition,transition, distribution
      read(1,'(a80)') title
      read(1,*) ictd
C
C......(6)
C
C iqso: Quasi-sequence-order descriptors
      read(1,'(a80)') title
      read(1,*) iqso,nphi
C
C.......(7)
C
C ipaac: pseudo-amino acid composition
      read(1,'(a80)') title
      read(1,*) ipaac
      read(1,*) w1,lamda1
      read(1,*) nst1
      if(nst1.eq.0) then
         nset1=1
         nacp1(1)=3
         idacp1(1,1)=115
         idacp1(1,2)=485
         idacp1(1,3)=486
      else
        nset1=nst1
      do i=1,nset1
        read(1,*) nacp1(i),(idacp1(i,j),j=1,nacp1(i))
      end do
      end if 
C
C.......(8)
C
C iapaac: amphilphilic pseudo amino acid composition
      read(1,'(a80)') title
      read(1,*) iapaac
      read(1,*) w2,lamda2
C
C.......(9)
C itop : sequence-segmented topological descriptors
      read(1,'(a80)') title
      read(1,*) itop,ibcut
C
C.........(10)
C
C iapp : total amino acid property or sequence-segmented sum of amino acid property
      read(1,'(a80)') title
      read(1,*) iapp
      read(1,*) napp,(iappn(i),i=1,napp)
C
C........(11)
C  methodpp: method for the construction of descriptors from two proteins
C 
      read(1,'(a80)') title
      read(1,*)methodpp 
C
C........(12)
C  methodpl: method for the construction of descriptors from one protein and one ligand
C 
      read(1,'(a80)') title
      read(1,*)methodpl 
C
      return
      end
C**********************************************************
C     pseudo amino acid composition
C  refer to:
C     Kuo-chen Chou, Proteins,2001,43-246.
C***********************************************************
      subroutine paac(seq,len,np,w,lamda,xq,nq,ns,nid,idh)
      implicit double precision(a-h,o-z)
      common/aaindex/aadb(1000,20)
      dimension xq(1000)
      dimension f(20),na(20)
      character seq*10000,aasw*20,aagr*20,name*30
      character aa*20
      dimension ip(len)
      dimension xh0(nid,20),xh(nid,20)
      dimension cf(20,20)
      dimension th(lamda)
      dimension idh(1000)
C
C     Amino acid indexes are in the following order
C
      aagr='ARNDCQEGHILKMFPSTWYV'
C
C
C     the length of the protein sequence must be greater than lamda
C     ,i.e, len >= lamda
C
       if(lamda.ge.len) then
         write(6,*)
         write(6,*) 'Error:lamda is not less than segment length '
         write(6,*) ',so set smaller value of lamda in parameter.dat'
         return
       end if
C
      do 5 i=1,nid
      ii=idh(i)
      do 8 j=1,20
       xh0(i,j)=aadb(ii,j) 
8     continue
5     continue
C
      do 10 i=1,len
      ip(i)=index(aagr,seq(i:i))
10    continue
C
C    occurence frequency of amnio acids
      do 100 i=1,20
      na(i)=0
100   continue
      do 200 i=1,len
      ii=ip(i)
      na(ii)=na(ii)+1
200   continue
      do 300 i=1,20
      f(i)=dble(na(i))/dble(len)
300   continue
C   
C     standard convesion of the the amino indices 
C
      do 350 k=1,nid
C
      avk=0.0d0
      do 400 i=1,20
      avk=av1+xh0(k,i)
400   continue
      avk=av1/20.0d0
      sk=0.0d0
      do 500 i=1,20
      sk=sk+(xh0(k,i)-avk)**2
500   continue
      sk=dsqrt(sk/20.0d0)
      do 600 i=1,20
      xh(k,i)=(xh0(k,i)-avk)/sk
600   continue
C
350   continue
C
C     The correlation function H(Ri,Rj): cf(i,j)
C
      do 800 i=1,20
      do 900 j=1,20
      hksum=0.0d0
      do 901 k=1,nid
      hksum=hksum+(xh(k,i)-xh(k,j))**2
901   continue
      cf(i,j)= hksum/dble(nid)
900   continue
800   continue
C
C      sequence order-correlated factor:th(i), i=1,2,lamda
C
      do 1000 i=1,lamda
         th(i)=0.0d0
         do 1100 j=1,len-i
         ii=ip(j)
         jj=ip(j+i)
         th(i)=th(i)+cf(ii,jj)
1100     continue
         th(i)=th(i)/dble(len-i)
1000  continue
C
C     Pseudo-amino acid composition
C
C      w=0.05d0
C
      sumf=0.0d0
      do 1200 i=1,20
      sumf=sumf+f(i)
1200  continue
      sumth=0.0d0
      do 1300 i=1,lamda
      sumth=sumth+th(i)
1300  continue
      denom=(sumf+w*sumth)
C
      do 1400 i=1,20
      xq(i)=f(i)/denom
1400  continue
      do 1500 i=1,lamda
      iu=20+i
      xq(iu)=w*th(i)/denom
1500  continue
C
      nq=20+lamda
C
      write(6,*)
      if(nid.gt.1) then
      write(6,*) '(G6) PAAC for amino acid index set ',(idh(k),k=1,nid)
      else
      write(6,*) '(G6) PAAC for amino acid index ',(idh(k),k=1,nid)
      end if
      nn=0
1600  nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.nq) n2=nq
      write(6,'(10f15.6)') (xq(i),i=n1,n2)
      if(n2.lt.nq) goto 1600
C
      do 1700 i=1,20
      if(np.eq.1) write(10,*)
     & 'Pseudo amino acid composition for amino acid ', aagr(i:i)
     & ,' in segment ',ns
1700  continue
      do 1800 i=21,nq
      if(np.eq.1) write(10,*)
     & 'Pseudo amino acid composition with lag u= ', i
     & ,' in segment ',ns
1800  continue
C
      return
      end
C**********************************************************
C      Amphiphilic pseudo amino acid composition
C  refer to:
C     Kuo-chen Chou, Using amphiphilic pseudo amino acid composition to predict enzyme subfamily classes
C      Bioinformatics,2005,21,10-19
C***********************************************************
      subroutine apaac(seq,len,np,w,lamda,nseg,xq,nq)
      implicit double precision(a-h,o-z)
      dimension xq(1000,nseg)
      dimension f(20),na(20)
      character seq*10000,aasw*20,aagr*20,name*30
      character aa*20
      dimension ip(len)
      dimension xh10(20),xh20(20),xm0(20)
      dimension xh1(20),xh2(20),xm(20)
      dimension h1(20,20),h2(20,20)
      dimension th(2*lamda)
      dimension lenseg(nseg)
C
C(1). Hydrophobicity index in order:aagr='ARNDCQEGHILKMFPSTWYV'
C
C     For old version of Profeat:
C
C     Hydrophobicity index H132 from : Jiri Damborsky, Protein Engeering ,1998,11,21-30
C      data  xh10/0.87,0.85,0.09,0.66,1.52,0.00,0.67,0.10,0.87,3.15
C     &        ,2.17,1.64,1.67,2.87,2.77,0.07,0.07,3.77,2.67,1.87/
C
C     For updated version of Profeat:
C
C     Hydrophobicity index xh10 from :
C     http://www.sjtu.edu.cn/bioinf/PseAAC/PseAAreadme.htm
C        
C
      data xh10/0.62,-2.53,-0.78,-0.90,0.29,-0.85,-0.74,0.48,-0.40,1.38
     &          ,1.06,-1.50,0.64,1.19,0.12,-0.18,-0.05,0.81,0.26,1.08/
C
C(2). Hydrophilicity value (Hopp-Woods, 1981)
C     Proc. Natl. Acad. Sci. USA 78, 3824-3828 (1981)
C    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
C    -0.5     3.0     0.2     3.0    -1.0     0.2     3.0     0.0    -0.5    -1.8
C    -1.8     3.0    -1.3    -2.5     0.0     0.3    -0.4    -3.4    -2.3    -1.5
      data xh20/-0.5,3.0,0.2,3.0,-1.0,0.2,3.0,0.0,-0.5,-1.8
     &         ,-1.8,3.0,-1.3,-2.5,0.0,0.3,-0.4,-3.4,-2.3,-1.5/
C(3). the mass of the amino acid side chain
      data xm0/15.027,100.136,58.052,59.037,47.093,72.079,73.064
     &        , 1.008, 81.089,57.108,57.108,72.118,57.147,91.125
     &        ,41.065, 31.026,45.053,130.161,107.124,43.081/
C
      aagr='ARNDCQEGHILKMFPSTWYV'
C
      do 10 i=1,len
      ip(i)=index(aagr,seq(i:i))
10    continue
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
C     standard convesion of the the amino indices 
C
      av1=0.0d0
      av2=0.0d0
      avm=0.0d0
      do 400 i=1,20
      av1=av1+xh10(i)
      av2=av2+xh20(i)
      avm=avm+xm0(i)
400   continue
      av1=av1/20.0d0
      av2=av2/20.0d0
      avm=avm/20.0d0
      s1=0.0d0
      s2=0.0d0
      sm=0.0d0
      do 500 i=1,20
      s1=s1+(xh10(i)-av1)**2
      s2=s2+(xh20(i)-av2)**2
      sm=sm+(xm0(i)-avm)**2
500   continue
      s1=dsqrt(s1/20.0d0)
      s2=dsqrt(s2/20.0d0)
      sm=dsqrt(sm/20.0d0)
      do 600 i=1,20
      xh1(i)=(xh10(i)-av1)/s1
      xh2(i)=(xh20(i)-av2)/s2
      xm(i)=(xm0(i)-avm)/sm
600   continue
C
C
C     The correlation function H1(Ri,Rj) and h2(i,j)
C
      do 800 i=1,20
      do 900 j=1,20
      h1(i,j)=xh1(i)*xh1(j)
      h2(i,j)=xh2(i)*xh2(j)
900   continue
800   continue
C
C
      n2=0
      do 50 ns=1,nseg
      n1=n2+1
      n2=n1+lenseg(ns)-1
C
       if(lamda.ge.lenseg(ns)) then
         write(6,*)
         write(6,*) 'Error: lamda is not less than segment length'
         write(6,*) ',so set smaller value of lamda in parameter.dat'
         write(6,*) ' or smaller value of number of segments.'
         return
       end if
C    occurence frequency of amnio acids
C
      do 100 i=1,20
      na(i)=0
100   continue
      do 200 i=n1,n2
      ii=ip(i)
      na(ii)=na(ii)+1
200   continue
      do 300 i=1,20
      f(i)=dble(na(i))/dble(lenseg(ns))
300   continue
C   
C      sequence order-correlated factor:th(2i),th(2i-1): i=1,2,lamda
C
      do 1000 i=1,lamda
         th(2*i)=0.0d0
         th(2*i-1)=0.0d0
         do 1100 j=n1,n2-i
         ii=ip(j)
         jj=ip(j+i)
         th(2*i)=th(2*i)+h1(ii,jj)
         th(2*i-1)=th(2*i-1)+h2(ii,jj)
1100     continue
         th(2*i)=th(2*i)/dble(lenseg(ns)-i)
         th(2*i-1)=th(2*i-1)/dble(lenseg(ns)-i)
1000  continue
C
C     Amphiphilic Pseudo-amino acid composition
C
C      w=0.5d0
C
      sumf=0.0d0
      do 1200 i=1,20
      sumf=sumf+f(i)
1200  continue
      sumth=0.0d0
      do 1300 i=1,lamda*2
      sumth=sumth+th(i)
1300  continue
      denom=(sumf+w*sumth)
C
      do 1400 i=1,20
      xq(i,ns)=f(i)/denom
1400  continue
      do 1500 i=1,lamda*2
      iu=20+i
      xq(iu,ns)=w*th(i)/denom
1500  continue
C
      nq=20+lamda*2
C
      write(6,*)
      write(6,*) ' (G7) Amphiphilic Pseudo amino acid composition'
      if(nseg.gt.1) then
         write(6,*) 'in segment ',ns
         write(6,*) seq(n1:n2)
      end if
      nn=0
1600  nn=nn+1
      np1=(nn-1)*10+1
      np2=nn*10
      if(np2.gt.nq) np2=nq
      write(6,'(10f15.6)') (xq(i,ns),i=np1,np2)
      if(np2.lt.nq) goto 1600
C
      do 1700 i=1,20
      if(np.eq.1) then 
        write(10,*)'APAAC for amino acid ',aagr(i:i), ' in segment',ns
      end if
1700  continue
      do 1800 i=21,nq
      if(np.eq.1) then 
        write(10,*)'APAAC with lag u= ', i,' in segment',ns
      end if
1800  continue
C
50    continue
C
      return
      end
C*************************************************************
C     Reading the amino acid index provoded 
C       by the database and by the user,respectively
C----------------------------------------------------------------
C    The index is provided in the form:
C    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
C    -0.5     3.0     0.2     3.0    -1.0     0.2     3.0     0.0    -0.5    -1.8
C    -1.8     3.0    -1.3    -2.5     0.0     0.3    -0.4    -3.4    -2.3    -1.5
C***************************************************************
      subroutine readindex
      implicit double precision(a-h,o-z)
      character*80 line
      character*20 aagr
C
      common/aaindex/aadb(1000,20)
C
      aagr='ARNDCQEGHILKMFPSTWYV'
C
        do 10 i=1,1000
          do 20 j=1,20
             aadb(i,j)=0.0d0
20        continue
10      continue
C
        do 40 i=1,1000
        do 50 j=1,20
          aadb(i,j)=0.0
50      continue
40    continue
      
C     Reading the amino acid index provoded by the database
200   continue 
      read(3,'(a80)',end=900) line
      if(line(1:9).eq.'[AAINDEX]') then 
          read(line(10:80),*) naa
      else 
          goto 200
      end if
C      write(6,*) naa
300   read(3,'(a80)',end=900) line
      if(line(1:2).eq.'//') goto 200
      if(line(1:1).eq.'I') then
        read(3,'(a80)') line
        read(line,*,err=400) (aadb(naa,i),i=1,10)
        read(3,'(a80)') line
        read(line,*,err=400) (aadb(naa,i),i=11,20)
        goto 500
400     write(6,*) 'Error in reading index values',
     & ' from database for index ',naa
        write(6,*) line
500     continue
      end if
      goto 300
900   continue
      return
      end
C*********************************************************
C     Amino acid composition distribution or sequence-segmented amino acid composition
C------------------------------------
C     Refer to: 
C      (1) Guang-Ya,Zhang, Bai-Shan Fang
C          "Predicting the cofactor of oxidoreductase based 
C           on amino acid composition distribution and Chou's 
C           amphilic pseudo-amino acid composition"
C         Journal of Theorteical Biology ,2008,253,310-315
C      (2) Jianyu Shi,Shaowu Zhang, Quan Pan, Guo-Ping Zhou
C          "Amino Acid Composition Distribution: A novel 
C           Sequence Representation for Prediction of Protein
C           Subcelleular Localization"
C           1st International Conference on Bioinformatics and
C           Biomedical Engineering (iCBBE 2007), vol.1, Wuhan, China. 
C      (3)  Shao-Wu Zhang, Wei Chen, Feng Yang
C           "Using Chou's pseudo amino acid composition to predict protein quaternary structure:
C            a sequence-segmented PseAAC approach"
C            Amino Acids,2008,35,591
C
C************************************************************
      subroutine aacd(str,seq,len,nseg,ip,xaacd)
      implicit double precision(a-h,o-z)
      character seq*10000,str*20
      dimension ip(len),n(20),xaacd(20,nseg)
      dimension lenseg(nseg)
C
      do 10 i=1,20
        do 20 j=1,nseg
        xaacd(i,j)=0.0d0
20      continue
10    continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
C      write(6,*) 'Length for segments'
C      write(6,*) (lenseg(i),i=1,nseg)
C
      write(6,*)
      if(nseg.gt.1) then 
         write(6,*) '(1) Amino Acid Composition Distribution (AACD)'
      else
          write(6,*) '(1) Amino Acid Composition (AAC)'
      end if
      n2=0
C
      do 50 ns=1,nseg
      n1=n2+1
      n2=n1+lenseg(ns)-1
C
      do 100 i=1,20
      n(i)=0
100   continue
C
      do 200 j=n1,n2
      jj=ip(j)
      n(jj)=n(jj)+1
200   continue   
C
      do 300 i=1,20
      xaacd(i,ns)=100.0d0*dble(n(i))/dble(lenseg(ns))
300   continue
C
      write(6,*)
      if(nseg.gt.1) then
         write(6,*) 'AACD in segment ', ns
         write(6,*) seq(n1:n2)
         write(6,*)
      end if
      write(6,1002) (str(i:i),i=1,10)
      write(6,1003) (xaacd(i,ns),i=1,10)
      write(6,1002) (str(i:i),i=11,20)
      write(6,1003) (xaacd(i,ns),i=11,20)
C
50    continue
C
1002  format(10(4x,A1,5x))
1003  format(10f10.6)
      return
      end
C*********************************************************************
C       sequence-segmented PseAAC 
C------------------------------------------------------------
C       refer to:
C
C      (1)  Shao-Wu Zhang, Wei Chen, Feng Yang
C           "Using Chou's pseudo amino acid composition to predict protein quaternary structure:
C            a sequence-segmented PseAAC approach"
C            Amino Acids,2008,35,591
C
      subroutine spaac(seq,len,np,wpseaa,lamda,xqs,nq,nseg,nid,idh)
      implicit double precision(a-h,o-z)
      common/aaindex/aadb(1000,20)
      character seq*10000,str*20
      character seqseg*10000
      dimension ip(len),n(20)
      dimension xq(1000),xqs(1000,nseg)
      dimension idh(1000)
      dimension lenseg(nseg)
C
      do 10 i=1,20
        do 20 j=1,nseg
        xqs(i,j)=0.0d0
20      continue
10    continue
C
      if(len.lt.nseg) then
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C      write(6,*) 'Length for segments'
C      write(6,*) (lenseg(i),i=1,nseg)
C
C      write(6,*) '(G6) Sequence-segmented pseudo amino acid composition'
C      write(6,*)
      n2=0
      do 100 ns=1,nseg
C         write(6,*) 'for segment ',ns
         n1=n2+1
         n2=n1+lenseg(ns)-1
C         write(6,*) seq(n1:n2)
         lens=lenseg(ns)
         seqseg(1:lens)=seq(n1:n2)
         call paac(seqseg,lens,np,wpseaa,lamda,xq,nq,ns,nid,idh)
         do 200 j=1,nq
         xqs(j,ns)=xq(j)
200      continue
C
100   continue
      return
      end
C*********************************************************************
C       sequence-segmented  topological descriptors
C------------------------------------------------------------
C       refer to:
C
C      (1)  Philip D.Mosier, Anne E. Counterman and Peter C. Jurs
C           "Prediction of peptide icon collision cross sections from
C            topological molecular structure adn amino acid parameters'
C            Anal Chem,2002,74,1360-1370
C      (2) Shuo Mao,Huo Dan-Quan, Mei Hu,Liang Guo-Zhao,Zhang Mei, Li Zhi-Liang
C          " New descriptors of amino acids and its applications to peptide quantitative 
C            structure-activity relationships"
C           Chinese J.Struct.Chem. 2008,27,1375-1383
C
      subroutine stopdes(seq,len,np,nseg,ibcut)
      implicit double precision(a-h,o-z)
      character seq*10000,str*20
      character seqseg*10000
      dimension ip(len),n(20)
      dimension lenseg(nseg)
C
C
      if(len.lt.nseg) then
      write(6,*)
      write(6,*) 'error: len < nseg'
      write(6,*) 'Decrease nseg in input-param.dat so that len > nseg'
      return
      end if
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C      write(6,*) 'Length for segments'
C      write(6,*) (lenseg(i),i=1,nseg)
      if(nseg.gt.1) then 
        write(16,*) '(G8) Sequence-segmented topological descriptors'
      else
        write(16,*) '(G8) Topological descriptors'
      end if
      write(16,*)
      n2=0
      do 100 ns=1,nseg
         n1=n2+1
         n2=n1+lenseg(ns)-1
         lens=lenseg(ns)
         seqseg(1:lens)=seq(n1:n2)
         if(ns.eq.1) then
           ibegin=1
         else
           ibegin=0
         end if
         if(ns.eq.nseg) then
           iend=1
         else
           iend=0
         end if
C
C        number of heavy atoms for segment, calculated for defining arrays
C
         call seqnatom(seqseg,ibegin,iend,lens,natom)
C
C        topological descriptors for segment
C
         call topseq(seqseg,ibegin,iend,lens,np,ns,natom,ibcut)
C
100   continue
C
      return
      end
C*********************************************************************
C     Number of heavy atoms for segment
C*********************************************************************
      subroutine seqnatom(seq,ibegin,iend,len,natom)
      implicit double precision(a-h,o-z)
      character seq*10000
      character aa*20
      dimension nsd(20)
C     number of heavy atoms for side chain for amino acids in order :
C                     aa='ACDEFGHIKLMNPQRSTVWY'
      data nsd/1,2,4,5,7,0,6,4,5,4,4,4,3,5,7,2,3,3,10,8/
C
      aa='ACDEFGHIKLMNPQRSTVWY'
      natom=0
C     for the end amino acid,a OH group exists.
      if(iend.eq.1) natom=1
      do 100 i=1,len
      ixi=index(aa,seq(i:i))
      if(ixi.eq.0) then 
        write(16,*) 'A non-natural amino acid ', seq(i:i) 
     & , ' in this sequence'
        return
      end if
      na=nsd(ixi)
      natom=natom+na+4
100   continue
      return
      end
C*********************************************************************
C      Topological descriptors for sequence seq with length len
C*********************************************************************
      subroutine topseq(seq,ibegin,iend,len,np,ns,natom,ibcut)
      implicit double precision(a-h,o-z)
      character seq*10000
      character aa*20
      dimension ic(natom,natom),number(natom),numh(natom)
      dimension is(natom,natom),nbond(natom,natom)
C
C     Toplogy from sequence segment
C
      call topgeo(seq,ibegin,iend,len,natom,ic,number,numh,nbond)
C      do i=1,natom
C         do j=1,natom
C         if(ic(i,j).ne.0) write(6,*)i,j,'....', ic(i,j)
C         end do
C      end do
C
       call topdis(natom,ic,is)
C
       write(16,*) 'For segment ',ns
       call top(natom,ic,is,number,numh,nbond,ibcut)
C
      return
      end
C
C*********************************************************************
C     Topology geomrtey for peptide determined from its sequence
C*********************************************************************
      subroutine topgeo(seq,ibegin,iend,len,natom,ic,number,numh,nbond)
      implicit double precision(a-h,o-z)
      character seq*10000
      character aa*20
      dimension ic(natom,natom),number(natom),numh(natom)
      dimension iaa(natom)
      dimension lcon(natom),nbond(natom,natom)
C
      aa='ACDEFGHIKLMNPQRSTVWY'
      do 10 i=1,natom
        do 20 j=1,natom
         ic(i,j)=0
         nbond(i,j)=0
20      continue
10    continue
C
C     Adding Amino acids to form peptide
C
      na=0
      n3=0
      do 150 i=1,len
      na0=na+1
      n3a=n3
      n1=na+1
      n2=na+2
      n3=na+3
      n4=na+4
      n5=na+5
      n6=na+6
      n7=na+7
      n8=na+8
      n9=na+9
      n10=na+10
      n11=na+11
      n12=na+12
      n13=na+13
      n14=na+14
C
      nbond(n3,n4)=2
C
      if(seq(i:i).eq.'G') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=2
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       na=na+4
       goto 100
      end if
      if(seq(i:i).eq.'A') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=3
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       na=na+5
       goto 100
      end if
      if(seq(i:i).eq.'V') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=1
       number(n6)=6
       numh(n6)=3
       number(n7)=6
       numh(n7)=3
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n5,n7)=1
       na=na+7
       goto 100
      end if
      if(seq(i:i).eq.'L') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=1
       number(n7)=6
       numh(n7)=3
       number(n8)=6
       numh(n8)=3
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n6,n8)=1
       na=na+8
       goto 100
      end if
      if(seq(i:i).eq.'I') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=1
       number(n6)=6
       numh(n6)=3
       number(n7)=6
       numh(n7)=2
       number(n8)=6
       numh(n8)=3
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n5,n7)=1
       ic(n7,n8)=1
       na=na+8
       goto 100
      end if
      if(seq(i:i).eq.'P') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=2
       number(n7)=6
       numh(n7)=2
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n1,n7)=1
       na=na+7
       goto 100
      end if
      if(seq(i:i).eq.'F') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=0
       number(n7)=6
       numh(n7)=1
       number(n8)=6
       numh(n8)=1
       number(n9)=6
       numh(n9)=1
       number(n10)=6
       numh(n10)=1
       number(n11)=6
       numh(n11)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n8,n9)=1
       ic(n9,n10)=1
       ic(n10,n11)=1
       ic(n6,n11)=1
C
       nbond(n6,n7)=4
       nbond(n7,n8)=4
       nbond(n8,n9)=4
       nbond(n9,n10)=4
       nbond(n10,n11)=4
C
       na=na+11
       goto 100
      end if
      if(seq(i:i).eq.'Y') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=0
       number(n7)=6
       numh(n7)=1
       number(n8)=6
       numh(n8)=1
       number(n9)=6
       numh(n9)=0
       number(n10)=6
       numh(n10)=1
       number(n11)=6
       numh(n11)=1
       number(n12)=8
       numh(n12)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n8,n9)=1
       ic(n9,n10)=1
       ic(n10,n11)=1
       ic(n6,n11)=1
       ic(n9,n12)=1
C
       nbond(n6,n7)=4
       nbond(n7,n8)=4
       nbond(n8,n9)=4
       nbond(n9,n10)=4
       nbond(n10,n11)=4
C
       na=na+12
       goto 100
      end if
      if(seq(i:i).eq.'W') then
       number(n1)=7
       numh(n1)=1
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=0
       number(n7)=6
       numh(n7)=1
       number(n8)=7
       numh(n8)=1
       number(n9)=6
       numh(n9)=0
       number(n10)=6
       numh(n10)=0
       number(n11)=6
       numh(n11)=1
       number(n12)=6
       numh(n12)=1
       number(n13)=6
       numh(n13)=1
       number(n14)=6
       numh(n14)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n8,n9)=1
       ic(n9,n10)=1
       ic(n6,n10)=1
       ic(n10,n11)=1
       ic(n11,n12)=1
       ic(n12,n13)=1
       ic(n13,n14)=1
       ic(n9,n14)=1
C
       nbond(n6,n7)=4
       nbond(n7,n8)=4
       nbond(n8,n9)=4
       nbond(n9,n10)=4
       nbond(n10,n11)=4
       nbond(n11,n12)=4
       nbond(n12,n13)=4
       nbond(n13,n14)=4
       nbond(n9,n14)=4
C
       na=na+14
       goto 100
      end if
      if(seq(i:i).eq.'D') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=0
       number(n7)=8
       numh(n7)=0
       number(n8)=8
       numh(n8)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n6,n8)=1
       na=na+8
       goto 100
      end if
      if(seq(i:i).eq.'E') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=2
       number(n7)=6
       numh(n7)=0
       number(n8)=8
       numh(n8)=0
       number(n9)=8
       numh(n9)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n7,n9)=1
C
       nbond(n7,n8)=2
C
       na=na+9
       goto 100
      end if
      if(seq(i:i).eq.'H') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=0
       number(n7)=6
       numh(n7)=1
       number(n8)=7
       numh(n8)=1
       number(n9)=6
       numh(n9)=1
       number(n10)=7
       numh(n10)=0
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n8,n9)=1
       ic(n9,n10)=1
       ic(n6,n10)=1
C
       nbond(n6,n7)=2
       nbond(n9,n10)=2
C
       na=na+10
       goto 100
      end if
      if(seq(i:i).eq.'C') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=16
       numh(n6)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       na=na+6
       goto 100
      end if
      if(seq(i:i).eq.'M') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=2
       number(n7)=16
       numh(n7)=0
       number(n8)=6
       numh(n8)=3
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       na=na+8
       goto 100
      end if
      if(seq(i:i).eq.'S') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=8
       numh(n6)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       na=na+6
       goto 100
      end if
      if(seq(i:i).eq.'T') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=1
       number(n6)=8
       numh(n6)=1
       number(n7)=6
       numh(n7)=3
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n5,n7)=1
       na=na+7
       goto 100
      end if
      if(seq(i:i).eq.'N') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=0
       number(n7)=8
       numh(n7)=0
       number(n8)=7
       numh(n8)=2
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n6,n8)=1
C
       nbond(n6,n7)=2
C
       na=na+8
       goto 100
      end if
      if(seq(i:i).eq.'Q') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=2
       number(n7)=6
       numh(n7)=0
       number(n8)=8
       numh(n8)=0
       number(n9)=7
       numh(n9)=2
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n7,n9)=1
C
       nbond(n7,n8)=2
C
       na=na+9
       goto 100
      end if
      if(seq(i:i).eq.'K') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=2
       number(n7)=6
       numh(n7)=2
       number(n8)=6
       numh(n8)=2
       number(n9)=7
       numh(n9)=2
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n8,n9)=1
       na=na+9
       goto 100
      end if
      if(seq(i:i).eq.'R') then
       number(n1)=7
       numh(n1)=0
       number(n2)=6
       numh(n2)=1
       number(n3)=6
       numh(n3)=0
       number(n4)=8
       numh(n4)=0
       number(n5)=6
       numh(n5)=2
       number(n6)=6
       numh(n6)=2
       number(n7)=6
       numh(n7)=2
       number(n8)=7
       numh(n8)=1
       number(n9)=6
       numh(n9)=0
       number(n10)=7
       numh(n10)=2
       number(n11)=7
       numh(n11)=1
       if(n3a.ne.0) ic(n1,n3a)=1
       ic(n1,n2)=1
       ic(n2,n3)=1
       ic(n3,n4)=1
       ic(n2,n5)=1
       ic(n5,n6)=1
       ic(n6,n7)=1
       ic(n7,n8)=1
       ic(n8,n9)=1
       ic(n9,n10)=1
       ic(n9,n11)=1
C
       nbond(n9,n11)=2
C
       na=na+11
       goto 100
      end if
      write(6,*) 'No topological Geometry for ',seq(i:i)
C
100   continue
C
C     Which amino acid the atom  belongs to
C
      do 140 k=na0,na
       iaa(k)=i
140   continue
C
      if(ibegin.eq.1.and.i.eq.1) then
         if(seq(i:i).ne.'P') then
            numh(n1)=2
         else 
            numh(n1)=1
         end if
      end if
C
150   continue
C     For the end amino acid, a OH group exists.
      if(iend.eq.1) then
         na=na+1
         number(na)=8
         numh(na)=1
         ic(n3,na)=1
         iaa(na)=len
      end if
C      write(6,*) 'Number of heavy atoms in amino acids ',na,natom
C
      do 200 i=1,na
       do 300 j=1,na
          if(ic(i,j).ne.0) ic(j,i)=ic(i,j)
          if(nbond(i,j).gt.1) nbond(j,i)=nbond(i,j)
300    continue
200   continue
C
      do 210 i=1,na
       do 310 j=1,na
         if(nbond(i,j).eq.0) nbond(i,j)=ic(i,j)
310    continue
210   continue
C
      icheck=0
      if(icheck.ne.1) goto 600
      write(6,*) 'In segment:',seq(1:len)
      do 400 i=1,natom
      ncon=0
       do 500 j=1,natom
          if(ic(i,j).eq.1) then
             ncon=ncon+1
             lcon(ncon)=j
          end if
500    continue
      write(6,*) i, number(i),numh(i),' residue:'
     &,seq(iaa(i):iaa(i)),' bonded atoms:'
     & ,(lcon(j),j=1,ncon),' bond order:',(nbond(i,lcon(j)),j=1,ncon)
C
      if((numh(i)+ncon).gt.4) write(6,*) i,' > 4 atoms connected'
400   continue
600   continue
C
      return
      end
C****************************************************************
c     Topological Indices Calculation
C     Note: delete-H graph should be used for topological indices here
C*****************************************************************
      subroutine top(natom,ic,is,number,numh,nbond,ibcut)
      implicit  double precision(a-h,o-z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Unit-model: H-deleted model
C
C     input Variables description
C     natom:      number of atoms
C     number(i):  atomic number of atom i
C     name(i):    atomic symbol of atom i
C     numh(i):    number of H atoms bonded to atom i
C     sand(i):    sanderson electronegativity for atom i
C     rvd(i):     radius of van der Waals
C     vdw(i):     van der Waals volume for atom
C     alpha(i):   polarizability for atom i
C
C     hybrid(i):  hybrid type of atom i
C     ic(i,j):    adjacent matrix:
C                 ic(i,j).ne.0: if atom and j is bonded
C                 ic(i,j)=0:    if atom i and j is not directly bonded
C     iarom(i):   =1: if atom i is ia on an aromatic ring
C                 =0: if atom i is not on an aromatic ring
C     natom(i,j):  formal bond order between atom i,j
C     Used Variables
C
C     is(i,j):   topological distance matrix
C     nvd(i):    vertex degree of atom i
C     vvd(i):    valence vertex degree of atom i
C     vdm(i):    modified vertex degree of atom i
C     nve(i):    number of valence electrons for atom with atomic number i
C     nvds(i):   vertex distance degree, which is calculated as a sum
C                over the rows or column of topological distance matrix
C     rdm(i,j):  reciprocal distance matrix
C     r2dm(i,j): reciprocal square distance matrix
C     gm(i,j):   Galvez matrix
C     ctm(i,j):  charge term matrix
C     es(i):     estate index for atom i
C     esh(i):    estate index for Hydrogen bonded to heavy atom i 
C------------------------------------------------------------------
      character name(natom)*2,hybrid(natom)*3
      dimension number(natom),ic(natom,natom)
      dimension is(natom,natom),nvd(natom),numh(natom)
      dimension vvd(natom),iarom(natom)
      dimension nvds(natom)
      dimension sand(natom),rvd(natom),alpha(natom)
      dimension nve(100),vdw(natom)
      dimension vdm(natom),es(natom),esh(natom)
      dimension rdm(natom,natom),r2dm(natom,natom)
      dimension gm(natom,natom),ctm(natom,natom)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension nringarom(natom),nbond(natom,natom)
      dimension bond(natom,natom)
C
C      constitutational descriptors
C
      call constitution(natom,number,ic,numh)
C
C      Ring information
C
      call ringinf(natom,number,ic,is,iarom,
     & nring,nringatm,nringsize,nringarom)
C
C     reassign  bond order for bonds on aromatic rings(=4)
C
      do 100 i=1,nring
      if(nringarom(i).eq.0) goto 100
      nsi=nringsize(i)
      do 200 m=1,nsi-1
      do 300 n=m+1,nsi
      ma=nringatm(i,m)
      na=nringatm(i,n)
      if(nbond(ma,na).ge.1) nbond(ma,na)=4
 300  continue
 200  continue
 100  continue
C
C     fingerprint: according to function groups
C
      call groups(natom,number,nring,nbond
     & ,nringatm,nringsize,nringarom,ic,numh,iarom)
C
C     Number of rotable bonds
C
      call numbrot(natom,number,nring,nbond
     & ,nringatm,nringsize,nringarom,ic,is,numh,iarom)
C
C    Reciprocal distance matrix rdm(i,j) and reciprocal square distance matrix r2dm(i,j)
C    Refer to: Handbook of MOlecular Descriptors, pp116 - pp117
C
      do  i=1,natom
        do j=1,natom
          if(i.eq.j) then
            rdm(i,j)=0.0d0
            r2dm(i,j)=0.0d0
          else
            rdm(i,j)=1.0d0/dble(is(i,j))
            r2dm(i,j)=1.0d0/dble(is(i,j)*is(i,j))
          end if
        end do
      end do
C
C     Atomic property parameters
C
      call atompt(natom,number,sand,alpha,vdw,ierror)
C
      call radiusvdw(natom,number,rvd)
C
C     Number of Valence Electron
C
      call numbve(nve)
C
C     Vertex degrees nvd(i) determination
C
      call  vdd(natom,ic,nvd)
C
C     Valence vertex degree
C
      call nvvd(natom,number,numh,nve,vvd)
C
C     Vertex distance degrees nvds(i) determination
C
      call  vddd(natom,is,nvds)
C
C     Schultz molecular topological index
C
      call schultz(natom,is,nvd,nvds)
C
C     Galvez matrix gm=M
C         refter to: Handbook of molecular descriptors,pp445
C    
       do i=1,natom
         do j =1,natom
            gm(i,j)=0.0d0
            do k=1,natom
              gm(i,j)=gm(i,j)+ic(i,k)*r2dm(k,j)
            end do
          end do
       end do
C
C      Charge Trem Matrix = CT
C      refter to: Handbook of molecular matrix,pp445
C 
      do i=1,natom
         do j=1,natom
           if(i.eq.j) then
              ctm(i,j)=nvd(i)
            else
              ctm(i,j)=gm(i,j)-gm(j,i)
            end if
          end do
      end do
C
C     topological charge index
C
      call tci(natom,is,ctm)
C
C     wiener index
C
      call wiener(natom,is)
C
C
C     Harry Number H (also called Harary index)
C
      call harary(natom,is)
C
C     Gravitational topological index
C
      call gravittop(natom,number,is)
C
C     Total path Count
C
      call pathcount(natom,is)     
C
C     Modified vertex degree
C
      call mvd(natom,number,numh,nve,nvd,vdm)
C
C     Estate (Electrotopological state index) for atoms
C
      call estateatom(natom,number,numh,nve,nvd,ic,is,es,esh)
C
C    Xu index and modified Xu index
C
      call xuindex(natom,nve,nvd,nvds,vdm)
C
C     Balaban Index
C     nb: returned value of number of bonds
C
      call balaban(natom,ic,nvds,nb)
C
C     Platto Number
C
      call platt(natom,nvd,ic)
C
C     Zagreb index
C
      call zagreb(natom,ic,nvd)
C
C     Edge connective index
C
C     nb:number of edges=number of bonds
      neg=natom*(natom-1)/2
      call edge(natom,ic,is,nvd,nb)
C
C     Kier-Hall Connectivity indices
C     i.e.,Molecular connectivity chi index
C     Radic indices of different orders
C
      call  kierhall(natom,ic,nvd,x0,x1,x2)
C
C     Kier-Hall Valence connectivity index
C
      call valencecon(natom,number,numh,ic,nve,vvd,x0,x1,x2)
C
C     Solvation Connection index of different orders
C
      call sci(natom,number,ic,nvd,x0,x1,x2)
C
C     Path count,i.e.,
C     Number of path of length 1,2,3:  np1,np2,np3
C
      call npl(natom,is,np1,np2,np3)
C
C     Kappa-shape index
C
      call kappashape(natom,ic,np1,np2,np3)
C
C     Kappa-alpha-indices
C
      call kappaalpha(natom,number,numh,ic,np1,np2,np3)
C
C     Topological distance related descriptors
C  
      call dmrd(natom,is)
C
C     LogP from topological structure
C
      call logptop(natom,number,numh,ic,vvd,iarom)
C
C     Topological BCUT descriptors
C 
C     Bonds order in double precision  for bcut
C
      do i=1,natom
        do j=1,natom
        if(nbond(i,j).le.3) bond(i,j)=dble(nbond(i,j)) 
        if(nbond(i,j).eq.4) bond(i,j)=1.50d0
        if(nbond(i,j).ge.5) bond(i,j)=1.0d0
        end do
      end do
C
      if(ibcut.eq.1) then
        call bcuttop(natom,number,ic,sand,bond,es,rvd,alpha,vdw)
      end if
C
C    Autocorrelation descriptors
C
C     Moreau-Broto autocorrelation of a topological  structure
C
      call moreau(natom,is,number,sand,es,rvd,alpha,vdw)
C
C     Moran autocorrelation of a topological  structure
C
      call moran(natom,is,number,sand,es,rvd,alpha,vdw)
C
C     Geary autocorrelation of a topological  structure
C
      call geary(natom,is,number,sand,es,rvd,alpha,vdw)
C
C       Topological polar surface area(TPSA)
C
        call tpsa(natom,number,ic,nbond,numh,iarom 
     & ,nring,nringatm,nringsize)
C
      return
      end
C******************************************************************
C     Topological distance is(i,j) calculation 
C     output: the shortest path from atom i to atom j: is(i,j)
C     Input:
C     natom: number of atoms
C     ic: adjacent matrix:
C     ic(i,j)=1: if atom and j is bonded
C     ic(i,j)=0: if atom i and j is not directly bonded
C     ic(i,i)=0
C*****************************************************************
      Subroutine topdis(natom,ic,is)
      implicit  double precision(a-h,o-z)
      dimension ic(natom,natom),is(natom,natom),kill(natom) 
      dimension node(natom),nw(natom,natom)
C      do i=1,natom
C         do j=1,natom
C         if(ic(i,j).ne.0) write(6,*) i,j,'....', ic(i,j)
C         end do
C      end do
C
      do 100 i=1,natom
      do 200 j=1,natom
      is(i,j)=0
 200  continue
 100  continue
      do 300 i=1,natom
      do 400 j=1,natom
      kill(j)=0
 400  continue
      kill(i)=1
      norder=1
      node(1)=1
      nw(1,1)=i
 450  continue
      num=0
      do 500 j=1,node(norder)
      jc=nw(norder,j)
      do 600 k=1,natom
      if(kill(k).eq.1) goto 600
      if(ic(jc,k).eq.1) then
         kill(k)=1
         is(i,k)=norder
         num=num+1
         nw(norder+1,num)=k
      endif
 600  continue
 500  continue
      if(num.ne.0) then
      node(norder+1)=num
      norder=norder+1
      goto 450
      endif
 300  continue
C     
      iprint=0
      if(iprint.eq.0) goto 1000
      write(6,*) 'Topological Distance Matrix'
      do 900 i=1,natom
      do 950 j=1,natom
      write(6,'(3i5)') i,j,is(i,j)
950   continue
900   continue
1000  continue
      return
      end
C----------------------------------------------------
C      Kier-Hall Connectivity indices
C-----------------------------------------------
C     Refer to: Handbook of molecular descriptors,
C     pp85
C
C---------------------------------------------------
      Subroutine kierhall(natom,ic,nvd,x0,x1,x2)
      implicit  double precision(a-h,o-z)
      dimension ic(natom,natom),nvd(natom)
      write(16,*)
C     0th: X0
      X0=0.0d0
      do 500 i=1,natom
      X0=X0 + 1.0d0/dsqrt(dble(nvd(i)))
 500  continue
      write(16,1100) X0
C     1th order : X1
C     Number of edge: nb
      nb=0
      X1=0.0d0
      do 600 i=1,natom+1
      do 700 j=i+1,natom
      if(ic(i,j).ne.1) goto 700
      nb=nb+1
      X1=X1 + 1.0d0/dsqrt(dble(nvd(i)*nvd(j)))
 700  continue
 600  continue
C     mean Randic connectivity index
      x1bar=x1/dble(nb)
      write(16,1200) X1
      write(16,1210) X1bar
C
C     2th order: X2
C note: it is not required that i and j are adjacent
C here in this algorithm
      X2=0.0d0
      do 800 i=1,natom-1
      do 900 j=i+1,natom
      do 1000 k=1,natom
      if(ic(k,i).eq.1.and.ic(k,j).eq.1) then
      X2=X2+1.0d0/dsqrt(dble(nvd(i)*nvd(j)*nvd(k)))
      endif
1000  continue
 900  continue
 800  continue
      write(16,1300) X2
C
C Log of Simple topological index by Narumi(S)=sti
C: refer to: handbook,page 476
C and Narumi,H.(1987).New Topological Indices for Finite and 
C Infinite Systems.MATCH(Commm.Math.Comp.Chem.),22,195-207
C
C and Harmonic topological index(H)=hti and Geometric topologica index (G)=gti
C and  arithmetic topological index A=ati
C
      sti=0.0d0
      hti=0.0d0
      ati=0.0d0
      do 2000 i=1,natom 
      sti=sti+dble(nvd(i))
      hti=hti+1.0d0/dble(nvd(i))
      ati=ati+dble(nvd(i))
 2000 continue
      sti=log10(sti)
      gti=sti**(1.0d0/dble(natom))
      ati=ati/dble(natom)
      write(16,2010) sti
      write(16,2020) hti         
      write(16,2030) gti
      write(16,2040) ati
1100  Format('0th Kier-Hall connectivity index:', e12.4)
1200  format('1th Kier-Hall connectivity index:', e12.4)
1300  format('2th Kier-Hall connectivity index:', e12.4)
1210  format('Mean Randic Connectivity index : ', e12.4)
2010  format('Simple topological index by Narumi: ', e12.4)
2020  format('Harmonic topological index by Narumi:',e12.4)
2030  format('Geometric topological index by Narumi:',e12.4)
2040  format('Arithmetic topological index by Narumi:',e12.4)
      return
      end
C
C------------------------------------------------------------
      subroutine valencecon(natom,number,numh,ic,
     & nve,vvd,x0,x1,x2)
      implicit  double precision(a-h,o-z)
C
C      n(i): principal quantum number of atom with atomic number of i
C      nv(i): number of valence electrons for atom  i
C      nve(i): number of valence electrons for atom of atomic number i
C      Valence-connectivity indices
C
C     Refer to: 
C               (1).Lowell H.Hall, Lemont B.Kier,J.Mol.Graph.Model.20(2001)4-18.
C
      dimension number(natom),ic(natom,natom)
      dimension numh(natom),nv(natom)
      dimension nve(100),deltav(natom)
      dimension vvd(natom)
C n(i): Principal quantum number for atom with atomic number i
      dimension n(100)
      data n/1,1
     &      ,2,2,2,2,2,2,2,2
     &      ,3,3,3,3,3,3,3,3
     &      ,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
     &      ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
     &      ,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6
     &              ,6,6,6,6,6,6,6,6,6,6,6,6,6,6
     &      ,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
C nv(i):number of valence electrons for atom i
      do 100 i=1,natom
      k=number(i)
      nv(i)=nve(k)
 100  continue
C
C     write(6,*)'Number of valence electrons',(nv(i),i=1,natom)
C     write(6,*)'Number of attached H',(numh(i),i=1,natom)
C
C Valence Vertex Degree
      do 200 i=1,natom
        deltav(i)=vvd(i)
 200  continue
C
C     write(6,*)'dettav',(deltav(i),i=1,natom)
C     0th: Xv0
      Xv0=0.0d0
      do 500 i=1,natom
      Xv0=Xv0 + 1.0d0/dsqrt(dble(deltav(i)))
 500  continue
      write(16,1100) Xv0
C     1th order : Xv1
      Xv1=0.0d0
      do 600 i=1,natom+1
      do 700 j=i+1,natom
      if(ic(i,j).ne.1) goto 700
      Xv1=Xv1 + 1.0d0/dsqrt(dble(deltav(i)*deltav(j)))
 700  continue
 600  continue
      write(16,1200) Xv1
C     2th order: Xv2
      Xv2=0.0d0
      do 800 i=1,natom-1
      do 900 j=i+1,natom
      do 1000 k=1,natom
      if(ic(k,i).eq.1.and.ic(k,j).eq.1) then
      Xv2=Xv2+1.0d0/dsqrt(dble(deltav(i)*deltav(j)*deltav(k)))
      endif
1000  continue
 900  continue
 800  continue
      write(16,1300) Xv2
C
1100  Format('0th valence connectivity index:',e14.4)
1200  format('1th valence connectivity index:',e14.4)
1300  format('2th valence connectivity index:',e14.4)
C
C     delta chi index
C     refer to: Loewll H.Hall,Lemont B.Kier,J.Mol.Graph.Model.
C               20(2001)1-18
C
       dx0=x0-xv0
       dx1=x1-xv1
       dx2=x2-xv2
       write(16,1400) dx0
       write(16,1500) dx1
       write(16,1600) dx2
C
C Poligani index: refer to: Handbook,p476
C
      dz=0.0d0
      do 2000 i=1,natom
      k=number(i)
      li=n(k)
      dz=dz+dble(nv(i))/dble(li)
 20000 continue
      write(16,2010) dz
C
 1400 format('0th order delta chi index:',e12.4)
 1500 format('1th order delta chi index:',e12.4)
 1600 format('2th order delta chi index:',e12.4)
 2010 format('Pogliani index:',e12.4)
      return
      end
C-------------------------------------------------------------
C      Kier shape index of different orders
C
C     Refer to: (1).Reviews of computational chemistry,vol.2,1990,pp394
C                (2).Handbook of molecular descriptors,page 248-249
C--------------------------------------------------------------
      Subroutine kappashape(natom,ic,ip1,ip2,ip3)
      implicit  double precision(a-h,o-z)
      dimension ic(natom,natom)
C     assumed of an value for non-existing kappa values
C
      ppa1=1.0d0
      if(natom.eq.1) then
        ppa2=0.0d0
      else
        ppa2=1.0d0
      endif
      ppa3=0.0d0
      if(natom.eq.2) ppa3=1.450d0
      if(natom.eq.3) ppa3=2.0d0
      if(natom.eq.4) ppa3=3.378d0 
C
C
      a=dble(natom)
C     Kappa 1
      if(ip1.ne.0) then
      ppa1=a*(a-1.0d0)*(a-1.0d0)/dble(ip1*ip1)
      endif
C     Kappa 2
      if(ip2.ne.0) then
      ppa2=(a-1.0d0)*(a-2.0d0)*(a-2.0d0)/dble(ip2*ip2)
      endif
C     Kappa 3
      neven=2*int(natom/2)
      if(ip3.ne.0) then
        if(neven.eq.natom) then
          ppa3=(a-2.0d0)*(a-2.0d0)*(a-3.0d0)/dble(ip3*ip3)
        else
          ppa3=(a-1.0d0)*(a-3.0d0)*(a-3.0d0)/dble(ip3*ip3)
        endif
      endif
C
      write(16,1100) ppa1
      write(16,1200) ppa2
      write(16,1300) ppa3
C
1100  Format('1th order Kier shape index:', e12.4)
1200  format('2th order Kier shape index:', e12.4)
1300  format('3th order Kier shape index:', e12.4)
      return
      end
C**********************************************************
C      Kier alpha-modified shape index of different orders
C************************************************************
C     Refer to: (1).Reviews of computational chemistry,vol.2,1990,pp401
C                (2). Handbook of molecular descriptors,page 249-250
C
C     Note: (1).only neutral molecules can be used here
C               because atom type is judged from its connected atom number
C            (2).Ref1 and 2 have parameters for atoms: C,N,O,S,P,F,Cl,Br,I
C                with its different hybrid type.For other atoms, the alpha 
C                is calculated according to the formula ref1:
C                     alpha(x)=r(x)/r(Csp3)-1
C                where r(x) and r(csp3) are the radius of atom x and Csp3
C                respectively.( r(Csp3=0.77) ).
C
C-------------------------------------------------------------
      subroutine  kappaalpha(natom,number,numh,ic,ip1,ip2,ip3)
      implicit  double precision(a-h,o-z)
      dimension ic(natom,natom),number(natom)
      dimension nc(natom),alpha(natom),numh(natom)
      dimension covar(100)
      character*3 hybrid(natom)
C
C     assumed of an value for non-existing kappa values
C
      ppa1=1.0d0
      if(natom.eq.1) then
        ppa2=0.0d0
      else
        ppa2=1.0d0
      endif
      ppa3=0.0d0
      if(natom.eq.2) ppa3=1.450d0
      if(natom.eq.3) ppa3=2.0d0
      if(natom.eq.4) ppa3=3.378d0 
C
C Number of atoms
C
      a=dble(natom)
C
C Hybrid state by connection number including H atom    
C
      do 500 i=1,natom
      nc(i)=numh(i)
      do 600 j=1,natom
      if(ic(i,j).eq.1) nc(i)=nc(i)+1
 600  continue
 500  continue
C
C     write(6,*)'nc',(nc(i),i=1,natom)
c     atom type
C
      do 700 i=1,natom
C
C     for B
C
      if(number(i).eq.5) then
        alpha(i)=0.052
        goto 700
      end if
C
C     For Si
C
      if(number(i).eq.14)then
       alpha(i)=0.53
       goto 700
      end if
C
C     For C
C
      if(number(i).ne.6) goto 800
C C(SP3)
C     for the case that 3D structure has error
      if(nc(i).ge.5) then
        write(6,*) 'Warning** too many atoms(including H) are'
        write(6,*) '         connected to atom ',i
        write(6,*) 'so the results may have error'
        write(6,*) 'numh(i)=',numh(i)
        do j=1,natom
         if(ic(i,j).eq.1) write(6,*) 'j=',j
        end do
        goto 750
      endif
C
      if(nc(i).eq.4) then
        alpha(i)=0.0d0
        goto 700
      endif
C C(SP2)
      if(nc(i).eq.3) then
        alpha(i)=-0.13d0
        goto 700
      endif
C C(SP)
      if(nc(i).eq.2) then
        alpha(i)=-0.22d0
        goto 700
      endif
      write(6,*) 'No alpha value for C,atom type error'
      goto 750
 800  continue
C
C     For N
C
      if(number(i).ne.7) goto 900
C N(sp3)
      if(nc(i).eq.3) then
        alpha(i)=-0.04d0
        goto 700
      endif
C N(SP2)
      if(nc(i).eq.2) then
        alpha(i)=-0.20d0
        goto 700
      endif
C N(SP)
      if(nc(i).eq.1) then
        alpha(i)=-0.29d0
        goto 700
      endif
      write(6,*) 'No alpha value for N,atom type error'
      goto 750
 900  continue
C
C     For O
C
      if(number(i).ne.8) goto 1000
C O(SP3)
      if(nc(i).eq.2) then
        alpha(i)=-0.04d0
        goto 700
      endif
C O(SP2)
      if(nc(i).eq.1) then
        alpha(i)=-0.20d0
        goto 700
      endif
C
 1000 continue
C
C     For F
C
      if(number(i).eq.9) then
         alpha(i)=-0.07d0
         goto 700
      endif
C
C    For P 
C
      if(number(i).ne.15) goto 1010
         if(nc(i).ge.3) then
           alpha(i)=0.430d0 
           goto 700
         end if
         if(nc(i).eq.2) then
           alpha(i)=0.30
           goto 700
         end if
         if(nc(i).eq.1) then
           alpha(i)=0.20
           goto 700
         end if
 1010 continue
C
C     S
C
      if(number(i).ne.16) goto 1020
C S(SP3)
      if(nc(i).eq.2) then
        alpha(i)=0.350
        goto 700
      endif
C S(SP2)
      if(nc(i).eq.1) then
        alpha(i)=0.22d0
        goto 700
      endif
 1020 continue
C
C     Cl
C
      if(number(i).eq.17) then
      alpha(i)=0.29d0
      goto 700
      endif
C
C      Br
C
      if(number(i).eq.35) then
      alpha(i)=0.48d0
      goto 700
      endif
C
C     I 
C
      if(number(i).eq.53) then
      alpha(i)=0.73d0
      goto 700
      endif
C
C     For atoms without alpha parameters in ref1
C     or with atom type determination error in the above
C     alpha is calculated now
C
 750  continue
      call covalent(covar)
      rcsp3=0.77d0
      rx=covar(number(i))
      alpha(i)=rx/rcsp3-1.0d0
 700  continue
C
C     alpha
C
      alph=0.0d0
      do 1080 i=1,natom
      alph=alph+alpha(i)
 1080 continue
C     write(6,1001)(alpha(i),i=1,natom)
C     write(6,*) 'alph',alph
 1001 format('alpha',10f10.4)
C
C     Kappa 1
      if(ip1.ne.0) then
      ppa1=(a+alph)*(a+alph-1.0d0)*(a+alph-1.0d0)
      ppa1=ppa1/((dble(ip1)+alph)*(dble(ip1)+alph))
      endif
C     Kappa 2
      if(ip2.ne.0) then
      ppa2=(a+alph-1.0d0)*(a+alph-2.0d0)*(a+alph-2.0d0)
      ppa2=ppa2/((dble(ip2)+alph)*(dble(ip2)+alph))
      endif
C     Kappa 3
      neven=2*int(natom/2)
      if(ip3.ne.0) then
        if(neven.eq.natom) then
          ppa3=(a+alph-2.0d0)*(a+alph-2.0d0)*(a+alph-3.0d0)
        else
          ppa3=(a+alph-1.0d0)*(a+alph-3.0d0)*(a+alph-3.0d0)
        endif
        ppa3=ppa3/((dble(ip3)+alph)*(dble(ip3)+alph))    
      endif
C
C     Kier Molecular Flexibility Index
C     refer to (1): Handbook of molecular descriptors,pp178
C              (2): Kier,L.B.,QSAR,1989,8,221-224
      phi=ppa1*ppa2/a
C
      write(16,1100) ppa1
      write(16,1200) ppa2
      write(16,1300) ppa3
      write(16,1400) phi
C
1100  Format('1th order Kappa alpha shape index:', e12.4)
1200  format('2th order Kappa alpha shape index:', e12.4)
1300  format('3th order Kappa alpha shape index:', e12.4)
1400  format('Kier Molecular Flexibility Index:', e12.4)
      return
      end
C*********************************************************
C      Atomic properties:
C*********************************************************
C      refer to:
C      Perspective in Drug Discovery and Design,1998,9-11,355-380
C
C      Note: all properties listed in data block are normalized to C
C--------------------------------------------------------------------
      subroutine  atompt(natom,number,sand,alpha,vdw,ierror)
      implicit double precision(a-h,o-z)
      dimension number(natom)
      dimension sand(natom),en(100),nb(100)
      dimension vdw(natom),v(100)
      dimension alpha(natom),a(100)
C
C     Corresponding atomic numbers
C
      data(nb(i),i=1,19)/ 1,5,6,7,8,9,13,14,15,16,17,26,
     & 27,28,29,30,35,50,53/
C
C     Corresponding vDW volume
C
      data(v(i),i=1,19)/0.299,0.796,1.000,0.695,0.512,
     & 0.410,1.626,1.424,1.181,1.088,1.035,1.829,1.561,
     & 0.764,0.512,1.708,1.384,2.042,1.728/
C
C     Corresponding Sanderson Electronegativity
C
      data(en(i),i=1,19)/0.944,0.828,1.000,1.163,1.331,
     & 1.457,0.624,0.779,0.916,1.077,1.265,0.728,0.728,
     & 0.728,0.740,0.810,1.172,0.837,1.012/
C
C     Corresponding Polarizability
C
      data(a(i),i=1,19)/0.379,1.722,1.000,0.625,0.456,
     & 0.316,3.864,3.057,2.063,1.648,1.239,4.773,4.261,
     & 3.864,3.466,4.034,1.733,4.375,3.040/
      do 100 i=1,natom
      numb=number(i)
      do 200 k=1,19
      if(numb.eq.nb(k)) then
         sand(i)=en(k)*2.746
         alpha(i)=a(k)*1.760
         vdw(i)=v(k)*22.449
         goto 100
      end if
 200  continue
      write(16,*) 'No atomic properties for atom ',i
      write(16,*) 'Its atomic number is',numb
      ierror=1 
      return
100   continue
      ierror=0
      return
      end
C-----------------------------------------------------------
C     covalent radii for atoms
C     (1).for nonexist covalent radii data,set to be 1.5
C     most probally it is for transition metal atoms
C     (2). References are from different resources on
C          the internet web
C          http://EnvironmentalChemistry.com
C
C-----------------------------------------------------------
      subroutine covalent(covar)
      implicit double precision(a-h,o-z)
      dimension covar(100)
      do 100  i=1,100
 100  covar(i)=1.5
C
      covar(1)=0.32
      covar(2)=0.93
      covar(3)=1.23
      covar(4)=0.90
      covar(5)=0.82
      covar(6)=0.77
      covar(7)=0.75
      covar(8)=0.73
      covar(9)=0.72
      covar(10)=0.71
      covar(11)=1.54
      covar(12)=1.36
      covar(13)=1.18
      covar(14)=1.11
      covar(15)=1.06
      covar(16)=1.02
      covar(17)=0.99
      covar(18)=0.98
      covar(19)=2.03
      covar(20)=1.74
      covar(21)=1.44
      covar(22)=1.32
      covar(23)=1.22
      covar(24)=1.18
      covar(25)=1.17
      covar(26)=1.17
      covar(27)=1.16
      covar(28)=1.15
      covar(29)=1.17
      covar(30)=1.25
      covar(31)=1.26
      covar(32)=1.22
      covar(33)=1.20
      covar(34)=1.16
      covar(35)=1.14
      covar(36)=1.12
      covar(37)=2.16
      covar(38)=1.91
      covar(39)=1.62
      covar(40)=1.45
      covar(41)=1.34
      covar(42)=1.30
      covar(43)=1.27
      covar(44)=1.25
      covar(45)=1.25
      covar(46)=1.28
      covar(47)=1.34
      covar(48)=1.48
      covar(49)=1.44
      covar(50)=1.41
      covar(51)=1.40
      covar(52)=1.36
      covar(53)=1.33
      covar(54)=1.31
      covar(55)=2.35
      covar(56)=1.98
      covar(57)=1.69
      covar(58)=1.65
      covar(59)=1.65
      covar(60)=1.64
      covar(61)=1.63
      covar(62)=1.62
      covar(63)=1.85
      covar(64)=1.61
      covar(65)=1.59
      covar(66)=1.59
      covar(67)=1.58
      covar(68)=1.56
      covar(69)=1.56
      covar(70)=1.74
      covar(71)=1.56
      covar(72)=1.44
      covar(73)=1.34
      covar(74)=1.30
      covar(75)=1.28
      covar(76)=1.26
      covar(77)=1.27
      covar(78)=1.30
      covar(79)=1.34
      covar(80)=1.49
      covar(81)=1.48
      covar(82)=1.47
      covar(83)=1.46
      covar(84)=1.46
      covar(85)=1.45
      covar(90)=1.65
      covar(92)=1.42
      return
      end 
C-------------------------------------------------------
C First and Second Zagreb index
C Refer to: Gutman,I.et al.,J.Chem.Phys.1975,62,3399-3405
C           Handbook of molecular descriptors,pp509
C-------------------------------------------------------
      subroutine zagreb(natom,ic,nvd)
      implicit double precision(a-h,o-z)
      dimension ic(natom,natom),nvd(natom)
      m1=0
      m2=0
      rm1=0.0d0
      rm2=0.0d0
      do 100 i=1,natom
      m1=m1+nvd(i)*nvd(i)
      rm1=rm1+1.0d0/(nvd(i)*nvd(i))
      do 200 j=i+1,natom
      if(ic(i,j).eq.1) then
      m2=m2+nvd(i)*nvd(j)
      rm2=rm2+1.0d0/(nvd(i)*nvd(j))
      endif
 200  continue
 100  continue
C Quadratic index
      q=3.0-2.0*dble(natom)+dble(m1)/2.0d0
C
      write(16,101) m1
      write(16,102) m2
      write(16,103) rm1
      write(16,104) rm2
      write(16,105) q
 101  format('First  Zagreb Index(M1): ',i8)
 102  format('Second Zagreb Index(M2):',i8)
 103  format('First  Modified Zagreb Index:',e12.4)
 104  format('Second Modified Zagreb Index:',e12.4)
 105  format('Quadratic index(Q):          ',e12.4)
      return
      end
C-------------------------------------------------
C wiener index
C refer to: Handbook of molecular descriptors,p497
C--------------------------------------------------
      subroutine wiener(natom,is)
      implicit double precision(a-h,o-z)
      dimension is(natom,natom)
      iw=0
      do 100 i=1,natom-1
      do 200 j=i+1,natom
      iw=iw+is(i,j)
 200  continue
 100  continue
C     Average Wiener index
       wbar=2.0d0*dble(iw)/dble(natom*(natom-1))
       write(16,*)
       write(16,1100)iw
       write(16,1200) wbar
 1100 format('Wiener index:           ',i8)
 1200 format('Mean Wiener index:   ',e14.4)
       return
       end
     	
C-------------------------------------------------
C     Harary index
C     refer to: Handbook of molecular descriptors,p209-210
C
C--------------------------------------------------
      subroutine harary(natom,is)
      implicit double precision(a-h,o-z)
      dimension is(natom,natom)
      har=0.0d0
      do 100 i=1,natom-1
      do 200 j=i+1,natom
      har=har+1.0d0/dble(is(i,j))
 200  continue
 100  continue
       write(16,1100)har
 1100 format('Harary index:           ',e12.4)
       return
       end
     	
C-----------------------------------------------------------
C      Solvation  Connectivity indices
C
C     Refer to:
C    (1).Hand book of molecular descriptors, pp88
C    (2). Zefirov,N.S.,Palyulin,V.A.,JCICS,2001,41,1022-1027
C------------------------------------------------------------
      subroutine sci(natom,number,ic,nvd,x0,x1,x2)
      implicit  double precision(a-h,o-z)
      dimension ic(natom,natom),nvd(natom)
      dimension number(natom),n(100)
      dimension pr(natom)
C     n(i): principal number of atom  with atomic number(i)
C
      data n/1,1
     &      ,2,2,2,2,2,2,2,2
     &      ,3,3,3,3,3,3,3,3
     &      ,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
     &      ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
     &      ,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6
     &              ,6,6,6,6,6,6,6,6,6,6,6,6,6,6
     &      ,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
C
C      write(6,*)
C      write(6,*)    '    Solvation connectivity index'
C
C Principal number determination
      do 100 i=1,natom
      nb=number(i)
      pr(i)=dble(nb)
 100  continue
C 0th: Xs0
      Xs0=0.0d0
      do 200 i=1,natom
      Xs0=Xs0 + pr(i)/dsqrt(dble(nvd(i)))
 200  continue
      xs0=0.5d0*xs0
      write(16,1100) Xs0
C1th order : Xs1
C Number of edge: nb
      nb=0
      Xs1=0.0d0
      do 600 i=1,natom+1
      do 700 j=i+1,natom
      if(ic(i,j).ne.1) goto 700
      nb=nb+1
      Xs1=Xs1 + pr(i)*pr(j)/dsqrt(dble(nvd(i)*nvd(j)))
 700  continue
 600  continue
      xs1=0.25d0*Xs1
      write(16,1200) xs1
C     2th order: Xs2
      Xs2=0.0d0
      do 800 i=1,natom-1
      do 900 j=i+1,natom
      do 1000 k=1,natom
      if(ic(k,i).eq.1.and.ic(k,j).eq.1) then
      Xs2=Xs2+pr(i)*pr(j)*pr(k)/dsqrt(dble(nvd(i)*nvd(j)*nvd(k)))
      endif
1000  continue
 900  continue
 800  continue
      xs2=xs2/8.0d0
      write(16,1300) Xs2
C
1100  Format('0th Solvation connectivity index: ',e12.4)
1200  format('1th Solvation connectivity index: ',e12.4)
1300  format('2th Solvation connectivity index: ',e12.4)
      return
      end
C-------------------------------------------------------
C Platt Number F
C refer to:(1).Handbook of molecular descriptors,page125
C          (2).Platt,J.R.J.Chem.Phys.1947,15,419-420;
C                        J.Phys.Chem.1952,56,328
C-------------------------------------------------------
       subroutine platt(natom,nvd,ic)
       dimension ic(natom,natom),nvd(natom)
       nf=0
       do 100 i=1,natom-1
       do 200 j=i+1,natom
       if(ic(i,j).eq.1) then
       nf=nf+nvd(i)+nvd(j)-2
       endif
 200   continue
 100   continue
       write(16,101) nf
 101   format('Platt Number:',9x,i5)
       return
       end
C-------------------------------------
C Number of paths of length 1,2,3
C------------------------------------- 
      Subroutine npl(natom,is,ip1,ip2,ip3)
      implicit  double precision(a-h,o-z)
      dimension is(natom,natom)
C     path counts of order 1,2,3
      ip1=0
      ip2=0
      ip3=0
      do 100 i=1,natom-1
      do 200 j=i+1,natom
      if(is(i,j).eq.1) ip1=ip1+1
      if(is(i,j).eq.2) ip2=ip2+1
      if(is(i,j).eq.3) ip3=ip3+1
 200  continue
 100  continue
C     write(6,*) 'p1,p2,p3',ip1,ip2,ip3
      return
      end
C-----------------------------------------
C      Vertex degrees nvd(i) determination
C-------------------------------------------
       subroutine vdd(natom,ic,nvd)
       dimension ic(natom,natom),nvd(natom)
       do 100  i=1,natom
       nvd(i)=0
       do 200 j=1,natom
       if(ic(i,j).eq.1) then
        nvd(i)=nvd(i)+1
       end if
 200   continue
       if(nvd(i).eq.0) write(6,*) 'Warning** atom is ',i,'isolated'
 100   continue
       return
       end
C--------------------------------------------------
C      Vertex distance degrees nvd(i) determination
C--------------------------------------------------
       subroutine vddd(natom,is,nvds)
       dimension is(natom,natom),nvds(natom)
       do 100  i=1,natom
       nvds(i)=0
       do 200 j=1,natom
       nvds(i)=nvds(i)+is(i,j)
 200   continue
 100   continue
       return
       end
C---------------------------------------------------
C        Balaban Index
C
C Refre to: handbook of molecular descriptors,page 21
C---------------------------------------------------
       subroutine balaban(natom,ic,nvds,nb)
       implicit double precision (a-h,o-z)
       dimension ic(natom,natom),nvds(natom)
C      Number of bonds
       nb=0
       do 100 i=1,natom-1
       do 200 j=i+1,natom
       if (ic(i,j).eq.1) nb=nb+1
 200   continue
 100   continue
C number of cycles
       nc=nb-natom+1
       bj=0.0d0
       do 300 i=1,natom-1
       do 400 j=i+1,natom
       if(ic(i,j).ne.1) goto 400
       bj=bj+1.0d0/dsqrt(dble(nvds(i)))*1.0/dsqrt(dble(nvds(j)))
 400   continue
 300   continue
       bj=bj*dble(nb)/dble(nc+1)
       write(16,101) bj
101    format('Balaban Index J:  ',e14.4)
       return
       end
C---------------------------------------------------------
C     Path counts       
C
C     Refer to:Handbook of Molecular Descriptors,page 344-346
C-----------------------------------------------------------
      subroutine pathcount(natom,is)
      implicit double precision(a-h,o-z)
      dimension is(natom,natom),mpath(natom-1)
C
C     path count of length from 1 to natom-1
C
C      write(6,*) '   Molecular path count'
      do 100 k=1,natom-1
      mpath(k)=0
      do 200 i=1,natom-1
      do 300 j=i+1,natom
      if(is(i,j).eq.k) mpath(k)=mpath(k)+1
 300  continue
 200  continue
 100  continue
C Total path count
      mtpc=0
      do 400 i=1,natom-1
      mtpc=mtpc+mpath(i)
 400  continue      
      n0=0
      if(natom.gt.1) then
         write(16,101) mpath(1)
      else
         write(16,101) n0
      end if
      if(natom.gt.2) then
         write(16,102) mpath(2)
      else
         write(16,102) n0
      end if
      if(natom.gt.3) then
         write(16,103) mpath(3)
      else
         write(16,103) n0
      end if
      if(natom.gt.4) then
         write(16,104) mpath(4)
      else
         write(16,104) n0
      end if
      if(natom.gt.5) then
         write(16,105) mpath(5)
      else
         write(16,105) n0
      end if
      if(natom.gt.6) then
         write(16,106) mpath(6)
      else
         write(16,106)  n0
      end if
C
      write(16,*)
      write(16,107) mtpc
C
 101  format('Molecular path count of length 1:',i9)
 102  format('Molecular path count of length 2:',i9)
 103  format('Molecular path count of length 3:',i9)
 104  format('Molecular path count of length 4:',i9)
 105  format('Molecular path count of length 5:',i9)
 106  format('Molecular path count of length 6:',i9)
 107  format('Total path count:                ',i9)
      return
      end
C-------------------------------------------------
C Gravitational index using topological distance
C
C refer to: Handbook of molecular descriptors,p209-210
C
C--------------------------------------------------
      subroutine gravittop(natom,number,is)
      implicit double precision(a-h,o-z)
      dimension number(natom)
      dimension is(natom,natom),aw(100)
C
C Atomic weight
C
      call atomw(aw)
C
      gravt=0.0d0
      do 100 i=1,natom-1
      do 200 j=i+1,natom
      nbi=number(i)
      nbj=number(j)
      wi=aw(nbi)/dble(natom)
      wj=aw(nbj)/dble(natom)
      gravt=gravt+wi*wj/dble(is(i,j)*is(i,j))
 200  continue
 100  continue
       write(16,1100)gravt
 1100 format('Gravitational topological index: ',e12.4)
       return
       end
C----------------------------------------------
C    Atomic weight for atom with atomic number of i
C
C     Atomic weight:refer to :http://dayah.com/periodic/
C----------------------------------------------
      subroutine  atomw(aw)
      implicit  double precision(a-h,o-z)
      dimension atw(100),aw(100)
      data (atw(i),i=1,53)/
     &         1.0079, 4.0026, 6.9410, 9.0122,10.8110,12.0107,14.0067
     &       ,15.9994,18.9984,20.1797,22.9898,24.3050,26.9815,28.0855   
     &       ,30.9738,32.0660,35.4527,39.9480,39.0983,40.0780,44.9559
     &       ,47.8670,50.9415,51.9961,54.9380,55.8457,58.9332,58.6934
     &       ,63.5460,65.3900,69.7230,72.6100,74.9216,78.9600,79.9040
     &       ,83.8000,85.4678,87.6200,88.9058,91.2240,92.9064,95.9400
     &       ,98.0000,101.070,102.905,106.420,107.868,112.411,114.810
     &       ,118.710,121.760,127.600,126.904/
      data (atw(i),i=54,100)/
     &       131.290,132.905,137.327,139.905,140.116,140.908,144.240
     &      ,145.000,150.360,151.964,157.250,158.925,162.500,164.930
     &      ,167.260,168.934,173.040,174.967,178.490,180.948,183.840
     &      ,186.207,190.230,192.217,195.078,196.966,200.590,204.383
     &      ,207.200,208.980,209.000,210.000,222.000,223.000,226.000
     &      ,227.000,232.038,231.036,238.029,237.000,244.000,243.000
     &      ,247.000,247.000,251.000,252.000,257.000/
       do 100 i=1,100
       aw(i)=atw(i)
 100   continue
      return
      end
C--------------------------------
C       Number of Valence Electrons
C-----------------------------------
      subroutine numbve(nve)
      implicit double precision(a-h,o-z)
      dimension nve(100)
C First row (1-2)
      nve(1)=1
      nve(2)=0
C Second row (3-8)   
      nve(3)=1
      nve(4)=2
      nve(5)=3
      nve(6)=4
      nve(7)=5
      nve(8)=6
      nve(9)=7
      nve(10)=0
C Third row(11-18)
      do 10 i=3,10
      nve(i+8)=nve(i)
 10   continue
C Fourth row main elements
      do  20  i=11,12
      nve(i+8)=nve(i)
 20   continue
      do  30 i=13,18
      nve(i+18)=nve(i)
 30   continue
C Fifth row main elements
      do  40 i=3,4
      nve(i+34)=nve(i)
 40   continue
      do  50 i=5,10
      nve(i+44)=nve(i)
 50   continue
C Sixth row main elements
      do  60 i=3,4
      nve(i+52)=nve(i)
 60   continue
      do  70 i=5,10
      nve(i+76)=nve(i)
 70   continue
C  Seventh row main elements
      nve(87)=1
      nve(88)=2
C  Third row: Sc to Zn
      do  80 i=21,26
      nve(i)=3+(i-21)
80    continue 
      nve(27)=8
      nve(28)=8
      nve(29)=1
      nve(30)=2
C  Fouth row:Y-Cd
      do  91 i=21,30
      nve(i+18)=nve(i)
  91  continue
C Fifth row:Lu to Hg
      do  92 i=39,48
      nve(i+32)=nve(i)
 92   continue
C Lanthanoids
      do  93 i=57,70
      nve(i)=3
 93   continue
C  Actinoids
      do  94 i=89,100
      nve(i)=3
 94   continue
C    
      return
      end
C-----------------------------------
C    Modified Vertex Degree
C------------------------------------
      subroutine mvd(natom,number,numh,nve,nvd,vdm)
      implicit double precision(a-h,o-z)
      dimension number(natom),numh(natom),nve(natom),nvd(natom)
      dimension vdm(natom)
      a=2.0d0/dble(natom)
      a=a*a
      do 100 i=1,natom
      nb=number(i)
      bi=nve(nb)-numh(i)
      ci=nb-nve(nb)
      di=dble(bi)/dble(ci)
      di=di*a+1.0d0
      vdm(i)=1.0d0/di
 100  continue
      return
      end
C------------------------------------------------------------
C     Xu index and Modified Xu index
C
C Refer to:(1).Handbook of Molecular Descriptors,page 507
C         :(2). Ren,B.,J.Chem.Inform.Comput.Sci.1999,39,139.
C------------------------------------------------------------
      subroutine xuindex(natom,nve,nvd,nvds,vdm)
      implicit double precision(a-h,o-z)
      dimension nve(natom),nvds(natom),vdm(natom)
      dimension nvd(natom)
      a=0.0d0
      b=0.0d0
      c=0.0d0
      d=0.0d0
      sqa=dsqrt(dble(natom))
      do 100 i=1,natom
      a=a+dble(nvd(i)*nvds(i)*nvds(i))
      b=b+dble(nvd(i)*nvds(i))
      c=c+vdm(i)*dble(nvds(i)*nvds(i))
      d=d+vdm(i)*dble(nvds(i))
 100  continue
      xu=sqa*log10(a/b)
      xum=sqa*log10(c/d)
      write(16,101) xu
      write(16,102) xum
  101 format('Xu index:',e12.4)
  102 format('Modified Xu Index:',e12.4)
      return
      end
Ca--------------------------------------------------------
C Topological charge index = Gk
C      refer to: Hand book of molecular descriptors,pp445
C---------------------------------------------------------
C  Natom: number of atoms (not includeng H atoms)
C  is(i,j): topological distance between atom i and j
C  ctm: charge trem matrix
C
C  g(k): topological charge index of order k
C  gm(k): mean topological charge index:Jk
C  J: global topological charge index
C------------------------------------------------------------
      subroutine tci(natom,is,ctm)
      implicit double precision(a-h,o-z)
      dimension is(natom,natom),ctm(natom,natom),g(5)
      dimension gm(5)
C
C      write(16,*) '     Topological charge index'
C
      do 100  k=1,5
         g(k)=0.0d0
 100  continue
      do 200 k=1,5
          g(k)=0.0d0
      do  300 i=1,natom
      do 400 j=1,natom
            if (is(i,j).eq.k) then
              delta=1.0d0
            else
               delta=0.0d0
             end if
            g(k)=g(k)+dabs(ctm(i,j))*delta
  400 continue
  300 continue
          g(k)=0.5d0*g(k)
  200 continue
      write(16,*)
      do 500 k=1,5
         write(16,1010)k,g(k)
  500 continue
 1010 format('Topological charge index G',i1,':',e12.4)
C     
C Mean topological charge index Jk
C
      write(16,*) '      Mean topological charge index'
C
      do 600 k=1,5
        gm(k)=g(k)/dble(natom-1)
 600  continue
      
      do 700 k=1,5
         write(16,1020)k,gm(k)
  700 continue
 1020 format('Mean topological charge index J',i1,':',e12.4)
C     
C Global topological charge index J
C
      gj=0.0d0
      do 800 k=1,5 
         gj=gj+gm(k)
 800  continue
      write(16,1030) gj
 1030 format('Global topological charge index J:',e12.4)
C
      return
      end
C--------------------------------------------------------------------------------------
C Schultz molecular topological index(MTI)
C-------------------------
C  refer to: Handbook,pp381-pp383
C  and Schultz,H.P.(1989). Topological Organic Chemistry.i
C      1.Graph Theory  and Topologica indices of Alkanes.J.Chem.
C        Inf.Comput.Sci.1989,29,227-228.
C------------------------
C  is(i,j): topological distance between i and atom j
C  nvd(i): vertex degree of atom i
C  nvds(i):vertex distance degree of atom i
C-------------------------------------------------------------------------------------------
      subroutine schultz(natom,is,nvd,nvds)
      implicit double precision(a-h,o-z)
      dimension nvd(natom),nvds(natom)
      dimension is(natom,natom)
C
C  Schultz molecular topological index
C  The Second Zagreb index=M2:szi2
C  S index: si
      szi2=0.0d0
      si=0.0d0
      do 100 i=1,natom
      szi2=szi2+dble(nvd(i)*nvd(i))
      si=si+dble(nvd(i)*nvds(i))
 100  continue
      dmti=szi2+si
      write(16,1010) dmti
C    
C Gutman molecular topological index SG
      nsg=0
      do 200 i=1,natom
      do 300 j=1,natom
      nsg=nsg+nvd(i)*nvd(j)*is(i,j)
  300 continue
  200 continue
      sg=dble(nsg)
      write(16,1020) sg
C
 1010 format('Schultz molecular topological index: ',e14.4)
 1020 format('Gutman molecular topological index : ',e14.4)
      return
      end
C------------------------------------------------------------
C valence vertex degree
C----------------------------------
C  Refer to: Handbook,page 474-475
C------------------------------------------------------------
      subroutine nvvd(natom,number,numh,nve,vvd)
      implicit  double precision(a-h,o-z)
C
C      nv(i): number of valence electrons for atom  i
C      nve(i): number of valence electrons for atom of atomic number i
C
C
      dimension number(natom)
      dimension numh(natom),nv(natom)
      dimension nve(100)
      dimension vvd(natom)
C nv(i):number of valence electrons for atom i
      do 100 i=1,natom
      k=number(i)
      nv(i)=nve(k)
 100  continue
C
C     write(16,*)'Number of valence electrons',(nv(i),i=1,natom)
C     write(16,*)'Number of attached H',(numh(i),i=1,natom)
C
C Valence Vertex Degree
      do 200 i=1,natom
      if(number(i).le.9) then
        vvd(i)=dble(nv(i)-numh(i))
      else
        vvd(i)=dble(nv(i)-numh(i))/dble(number(i)-nv(i)-1)
      endif
 200  continue
C      write(6,*) 'vvd=',(vvd(i),i=1,natom)
C
      return
      end
C------------------------------------------------------------
C  Edge connectivity index
C------------------------
C  refer to: (1).Handbook,pp124-130
C            (2).Estrada,E. "Edge adjacentcy relationshipas and 
C                a novel topological index required to molecule 
C                volume",J.Chem.Inf.Comput.Sci.,1995,35.31.
C-------------------------
C  Spectral moment :refer to:
C--------------------------
C  Ernesto Estrada,Santa Clara,"Spectral Moments of Edge adjacent 
C  Matrix in Molecular Graphs.1.Definition and Applications to thei
C  Prediction of Physical Properties of Alkanes.". J.CHem.Inf.Comput
C .Sci.,1996,36,844
C-------------------------------------
C neg: number of edges  =   nb: number of bonds
C ned(i,j): edge adjacency matrix
C negd(i):  edge degree of edge i
C n1(i) and n2(i): edge i  is connected bt vertex n1(i) and n2(i)
C nw(i,j): a working matrix)
C-----------------------------------------------------------------
      subroutine edge(natom,ic,is,nvd,nb)
      implicit double precision(a-h,o-z)
      dimension ic(natom,natom),is(natom,natom),nvd(natom)
      dimension ned(nb,nb),negd(nb)
      dimension n1(nb),n2(nb)
      dimension nw2(nb,nb),nw3(nb,nb),nw4(nb,nb),nw5(nb,nb),nw6(nb,nb)
      dimension nw7(nb,nb),nw8(nb,nb),nw9(nb,nb),nw10(nb,nb)
C
      do 100 i=1,nb
      negd(i)=0
      do 200 j= 1,nb
      ned(i,j)=0
  200 continue
  100 continue
C
C sorting edges
C
      neg=0
      do 300 i=1,natom-1
      do 400 j=i+1,natom
      if (ic(i,j).eq.0) goto 400
      neg=neg+1
      n1(neg)=i
      n2(neg)=j
 400  continue
 300  continue
C
C Edge adjacency matrix E: ned(i,j)
C
       do 500 i=1,neg-1
       do 600 j=i+1,neg
       if(n1(i).eq.n1(j).or.n1(i).eq.n2(j).or
     & .n2(i).eq.n1(j).or.n2(i).eq.n2(j)) then
       ned(i,j)=1
       ned(j,i)=1
       end if
 600   continue
 500   continue
C      write(6,*) 'Edge adjacency matrix'
C      do 700 i=1,neg
C      write(6,*) (ned(i,j),j=1,i)
C700   continue
C
C Edge Degree
C
      do 800 i=1,neg
      negd(i)=0
      do 900 j=1,neg
         if(ned(i,j).eq.1) then
            negd(i)=negd(i)+1
         end if
 900  continue
 800  continue
C
C  Edge connectivity index 
C
      epson1=0.0d0
      if(neg.eq.1) goto 915
      do 910 i=1,neg
      epson1=epson1+1.0d0/dsqrt(dble(negd(i)))
 910  continue
 915  write(16,1002) epson1
C
      epson2=0.0d0
      if(neg.eq.1) goto 1115
      do 1000 i=1,neg-1
      do 1100 j=i+1,neg
      if(ned(i,j).eq.1) then
         epson2=epson2+1.0d0/dsqrt(dble(negd(i)*negd(j)))
      end if
 1100 continue
 1000 continue
 1115 write(16,1003) epson2
C      
C    extented edge connectivity indices 
C
      epson3=0.0d0
      if(neq.eq.1) goto 1225
      do 1200 i=1,neg-1
      do 1300 j=i+1,neg
      do 1400 k=1,neg
      if(ned(k,i).eq.1.and.ned(k,j).eq.1) then
      epson3=epson3+1.0d0/dsqrt(dble(negd(i)*negd(j)*negd(k)))
      endif
1400  continue
1300  continue
1200  continue
1225  write(16,1004) epson3
C
C
        ispec=0
        if(ispec.eq.0) goto 9000
C
C Spectral moments of the edge adjacency matrix
C
C E*E
       call  matrixp(neg,ned,ned,nw2)
       nu2=0
       do 1500 i=1,neg
       nu2=nu2+nw2(i,i)
 1500  continue
       write(16,1005) nu2    
C    
C E*E*E
       call  matrixp(neg,ned,nw2,nw3)
       nu3=0
       do 1600 i=1,neg
       nu3=nu3+nw3(i,i)
 1600  continue
       write(16,1006) nu3
C E**4
       call  matrixp(neg,ned,nw3,nw4)
       nu4=0
       do 1700 i=1,neg
       nu4=nu4+nw4(i,i)
 1700  continue      
       write(16,1007)  nu4
C E**5
       call  matrixp(neg,ned,nw4,nw5)
       nu5=0
       do 1800 i=1,neg
       nu5=nu5+nw5(i,i)
 1800  continue      
       write(16,1008)  nu5
C E**6
       call  matrixp(neg,ned,nw5,nw6)
       nu6=0
       do 1900 i=1,neg
       nu6=nu6+nw6(i,i)
 1900  continue      
       write(16,1009) nu6
C E**7
       call  matrixp(neg,ned,nw6,nw7)
       nu7=0
       do 2000 i=1,neg
       nu7=nu7+nw7(i,i)
 2000  continue      
       write(16,1010)  nu7
C E**8
       call  matrixp(neg,ned,nw7,nw8)
       nu8=0
       do 2100 i=1,neg
       nu8=nu8+nw8(i,i)
 2100  continue      
       write(16,1020) nu8
C E**9
       call  matrixp(neg,ned,nw8,nw9)
       nu9=0
       do 2200 i=1,neg
       nu9=nu9+nw9(i,i)
 2200  continue      
       write(16,1030) nu9
C E**10
       call  matrixp(neg,ned,nw9,nw10)
       nu10=0
       do 2300 i=1,neg
       nu10=nu10+nw10(i,i)
 2300  continue      
       write(16,1040) nu10
C
 9000  continue
C
 1001 format('Error in edge: number of edges is not self-consistent')
 1002 format('0th edge connectivity index:',e14.4)
 1003 format('Edge connectivity index:      ',e14.4)
 1004 format('Extened edge connectivity inndex:',e14.4)
 1005 format('2th  spectral moment:',i10)
 1006 format('3th  spectral moment:',i10)
 1007 format('4th  spectral moment:',i10)
 1008 format('5th  spectral moment:',i10)
 1009 format('6th  spectral moment:',i10)
 1010 format('7th  spectral moment:',i10)
 1020 format('8th  spectral moment:',i10)
 1030 format('9th  spectral moment:',i10)
 1040 format('10th spectral moment:',i10)
      return
      end
C---------------------------------
C Product of square matrix mat1 and mat2
C mat=mat1*mat2
C where the matrice are integeral type
C-----------------------------------------
      subroutine matrixp(n,mat1,mat2,mat)
      implicit double precision(a-h,o-z)
      dimension mat1(n,n),mat2(n,n),mat(n,n)
      do 100 i=1,n
      do 200 j=1,n
      mat(i,j)=0
      do 300 k=1,n
      mat(i,j)=mat(i,j)+mat1(i,k)*mat2(k,j)
 300  continue
 200  continue
 100  continue
C     write(6,*) 'Matrix A*B'
C     do 400 i=1,n
C     write(6,*) (mat(i,j),j=1,n)
C400  continue
      return
      end
C
C------------------------------------------------
C     Number of H atoms that is attached to heavy atoms
C----------------------------------------------
      subroutine numbh(natom,natoms,number,intercon,numh)
      implicit double precision(a-h,o-z)
      dimension number(natom),intercon(natom,natom)
      dimension numh(natoms)
C
C     numh(i) (number of H attached to heavy atom i) 
C
      do 70 i=1,natoms 
      numh(i) = 0
 70   continue
      num=0
      do 80  i=1,natom
      if(number(i).eq.1) goto 80
      num=num+1
      do 90 j=1,natom
      if(j.eq.i) goto 90
      if(intercon(i,j).eq.0) goto 90
      if(number(j).ne.1) goto 90
      numh(num)=numh(num)+1
 90   continue
 80   continue
      return
      end
C***************************************************************************
C Elestrotopological state (Estate) indices for atoms
C****************************************************
C Refer to: (1).J.Chem.Inf.Comput.Sci.1995,35,1039-1045;1074-1080.
C           (2).Molconn-Z manual: chapter 2
C           (3).'Molecular structure descriptio.The Electropological
C                state.",L.B.Kier and L.H.Hall 
C H-depleted Model 
C     number(i): atomic number for atom  i
C     numh(i):   number of H atom bonded to heavy atom i
C     nve(i):    number of valence electrons of atom  with atomic number i
C     nvd(i):    number of connected edges to atom i in H-depleted graph
C     ic(i,j):   adjacent matrix
C     is(i,j):   topological distance matrix
C     es(i):     estate index for atom i
C     esh(i):    estate index for H atoms bonded to atom i
C     ei(i):     intric state  of atom i
C     eih(i):    intrix state for atom H atoms bonded to atom i 
C     n(i):      Principal quantum number for atom with atomic number i
C estate for H : using the X-H model
C      see: handbook of molecular descriptors, the frist equation
C      on page 162
C****************************************************************************
      subroutine  estateatom(natom,number,numh,nve,nvd,ic,is,es,esh)
      implicit double precision(a-h,o-z)
      dimension number(natom),ic(natom,natom),is(natom,natom)
      dimension numh(natom),nve(100),nvd(natom),es(natom),esh(natom)
      dimension ei(natom),eih(natom),n(100)
      data n/1,1
     &      ,2,2,2,2,2,2,2,2
     &      ,3,3,3,3,3,3,3,3
     &      ,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
     &      ,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
     &      ,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6
     &              ,6,6,6,6,6,6,6,6,6,6,6,6,6,6
     &      ,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
C     Intric state for heavy atoms
      do i=1,natom
        ei(i)=0.0d0 
        nb=number(i)
        delta=dble(nvd(i))
        deltav=dble(nve(nb)-numh(i))
        if((nve(nb)-numh(i)).ne.0) then
          ei(i)=(2.0d0/dble(n(nb)))**2
          ei(i)=(ei(i)*deltav+1.0d0)/delta
        end if
      end do 
C      write(6,*) 'Intric state for heavy atoms'
C      do i=1,natom
C       write(6,'(i4,e12.4)') i,ei(i)
C      end do
C     Pertrubation term for heavy atoms:pi
      do i=1,natom
         pi=0.0d0
       do j=1,natom
         if(i.ne.j)  then
           rij=dble(is(i,j)+1.0d0)
           rij2=rij*rij
           pi=pi+(ei(i)-ei(j))/rij2
         end if
       end do
       es(i)=ei(i)+pi
      end do  
C
C      write(6,*) 'Estate for  heavy atoms'
C      do i=1,natom
C        write(6,'(i4,e12.4)') i,es(i)
C      end do
C     Intric state for H atoms bonded to heavy atom i
      do i=1,natom
      nb=number(i)
      esh(i)=0.0d0
      eih(i)=0.0d0
      if(numh(i).ge.1) then
        delta=dble(nvd(i))
        deltav=dble(nve(nb)-numh(i))
        if((nve(nb)-numh(i)).ne.0) then
          eih(i)=(deltav-delta)**2/delta
        end if
      end if
      end do
C      write(6,*) 'Intrinc state  for H atoms bonded to atom i'
C      do i=1,natom
C        write(6,'(i4,f10.4)') i,eih(i)
C      end do
      do i=1,natom
        pe=0.0d0
C       Perturbation term for H atoms bonded to heavy atom i
        if(numh(i).ge.1) then 
          do j=1,natom
             if(i.ne.j) then
               rij=dble(is(i,j)+1.0d0)
               rij2=rij*rij
               pe=pe+(ei(i)-ei(j))/rij2
             end if
          end do
          esh(i)=eih(i)+pe
        end if
      end do
C
C      write(6,*) 'Estate for H atoms bonded to atom i'
C      do i=1,natom
C      write(6,'(i4,f10.4)') i,esh(i)
C      end do
      return
      end
C*********************************************************
C     Atomic number assignment from atom name in mol file
C**********************************************************
      subroutine  atomnumb(name,numb,ierror)
      implicit double precision(a-h,o-z)
      character name*2 
      character element(100)*2,element1(100)*2
      data (element(i),i=1,50)/
     &            ' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
     &            'NA','Mg','Al','SI',' P',' S','CL','AR',' K','CA',
     &            'SC','TI',' V','CR','MN','FE','CO','NI','CU','ZN',
     &            'GA','GE','AS','SE','BR','KR','RB','SR',' Y','ZR',
     &            'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN'/
      data (element(i),i=51,100)/
     &            'SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',
     &            'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     &            'LU','HF','TA',' W','RE','OS','IR','PT','AU','HG',
     &            'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     &            'PA',' U','NP','PU','AM','CM','BK','CF','ES','FM'/
C
      data (element1(i),i=1,50)/
     &            'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     &            'NA','Mg','Al','SI','P ','S ','CL','AR','K ','CA',
     &            'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     &            'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     &            'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN'/
      data (element1(i),i=51,100)/
     &            'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
     &            'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     &            'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
     &            'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     &            'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM'/
C
C     convert lower case to upper case,ex. Cl to CL,Br to BR
C
      m1=ichar(name(1:1))
      m2=ichar(name(2:2))
      if(m1.ge.97.and.m1.le.122) then
          m1=m1-32
          name(1:1)=char(m1)
      endif
      if(m2.ge.97.and.m2.le.122) then
          m2=m2-32
          name(2:2)=char(m2)
      endif
C
      do 20 j=1,100
      if(name.eq.element1(j).or.name.eq.element(j)) goto 30
 20   continue
C     write(6,*) 'Atom symbol not recognized ',name
      ierror=1
      numb=0
      return
 30   numb=j
      ierror=0
      return
      end 
C*********************************************************
C     extracting topological descriptors from  unit 16 for peptides
C**********************************************************
      subroutine extractp(np,itop,xtop,ntop)
      implicit double precision(a-h,o-z)
      dimension xtop(1000*itop)
      character molname*40,line*120,idname*30,sgname*15
      ierror=0
      rewind(16)
      ntop=0
      nsg=1
 100  continue
      read(16,'(a120)',end=300) line
      write(6,*) line
      ind=index(line,':')
      inds=index(line,'For segment')
      if(inds.gt.0) then 
         read(line(13:120),*,err=130) nsg
         goto 140
130      write(6,*) 'Error in reading segment number'
         write(6,*) line
140      continue
      end if
      if(ind.gt.0) read(line(ind+1:120),*,err=120) v
      goto 110
120   write(6,*) 'Error in reading topological descriptors in'
      write(6,*) line
110   continue
      if(ind.gt.0) then
      sgname(1:12)=' in segment '
      write(sgname(13:15),'(i3)') nsg
      if(np.eq.1)  write(10,'(A)') line(1:ind-1)//sgname(1:15)  
      ntop=ntop+1
      xtop(ntop)=v
      end if
      goto 100
C
 300  continue
      return
      end 
C*********************************************************
C     extracting topological descriptors from  unit 16
C**********************************************************
      subroutine extractl(np,itop,xtop,ntop)
      implicit double precision(a-h,o-z)
      dimension xtop(1000*itop)
      character molname*40,line*120,idname*30,sgname*15
      ierror=0
      rewind(16)
      ntop=0
      nsg=1
 100  continue
      read(16,'(a120)',end=300) line
      write(6,*) line
      ind=index(line,':')
      inds=index(line,'For segment')
      if(inds.gt.0) then 
         read(line(13:120),*,err=130) nsg
         goto 140
130      write(6,*) 'Error in reading segment number'
         write(6,*) line
140      continue
      end if
      if(ind.gt.0) read(line(ind+1:120),*,err=120) v
      goto 110
120   write(6,*) 'Error in reading topological descriptors in'
      write(6,*) line
110   continue
      if(ind.gt.0) then
      sgname(1:15)='  '
C      write(sgname(13:15),'(i3)') nsg
      if(np.eq.1)  write(2,'(A)') line(1:ind-1)//sgname(1:15)  
      ntop=ntop+1
      xtop(ntop)=v
      end if
      goto 100
C
 300  continue
      return
      end 
C*****************************************************************
C     van der Waals Radii of Elements
C     for all-atom molecule modles (including H atoms)
C
C     data are taken from:
C     Inorganic Materials,2001,vol.37(no.9),pp1031-1046
C     the recommended crystallographic data in table 9  
C     except H and inert atoms( radius for inert atom is useless)
C*******************************************************************
      subroutine radiusvdw(natom,number,rvd)
      implicit double precision(a-h,o-z)
      dimension rvdw(83),rvd(natom),number(natom)
C
C     Explicit H model(all-atom model)--rasmol radii(model=1)
C
C H

       rvdw(1)= 1.10
C Li-F
       data (rvdw(i),i=3,9) / 
     + 2.20,1.90,1.80,1.70,1.60,1.55,1.50/
C Na-Cl
       data (rvdw(i),i=11,17)/
     + 2.40,2.20,2.10,2.10,1.95,1.80,1.80/
C K-Br
       data (rvdw(i),i=19,35)/
     + 2.80,2.40,2.30,2.15,2.05,2.05,2.05,2.05,2.00,2.00,
     + 2.00,2.10,2.10,2.10,2.05,1.90,1.90/
C Rb-I
       data (rvdw(i),i=37,53)/
     + 2.90,2.55,2.40,2.30,2.15,2.10,2.05,2.05,2.00,2.05,
     + 2.10,2.20,2.20,2.25,2.20,2.10,2.10/
C Cs-Bi
       data (rvdw(i),i=55,56)/
     + 3.00,2.70/
C
       do 100 i=57,71
       rvdw(i)=2.50
 100   continue 
C
       data (rvdw(i),i=72,83)/
     + 2.25,2.20,2.10,2.05,2.00,2.00,2.05,
     + 2.10,2.05,2.20,2.30,2.30/
C
      do 200 i=1,natom
      nbi=number(i)
      if (nbi.le.83) then
         rvd(i)=rvdw(nbi)
      else 
         rvd(i)=2.2d0
      endif
 200  continue
      return
      End
C***********************************************************
C     Topological Distance Matrix related descriptors
C     refer to: Handbook of Molecular descriptors, pp112    
C*********************************************************
C     is(i,j): topological disance
C     ne(i): atom eccentricity 
C     nth(i): vertex distance degree
C     rdm(i,j): reciprocal distance matrix element
C     rds(i): reciprocal distance sum
C------------------------------------------------------
      subroutine dmrd(natom,is)
      implicit  double precision(a-h,o-z)
      dimension is(natom,natom),ne(natom),nth(natom)
      dimension rdm(natom,natom),rds(natom)
C      write(6,*)
C      write(6,*)'   Topological distance related'
      do 100 i=1,natom
      ne(i)=0
      nth(i)=0
      do 200 j=1,natom
      if(is(i,j).gt.ne(i)) ne(i)=is(i,j)
      nth(i)=nth(i)+is(i,j)
 200  continue
 100  continue
C     topological radius R and topological diameter D
      nd=ne(1)
      nr=ne(1)
      do 300 i=2,natom
      if(ne(i).lt.nr) nr=ne(i)
      if(ne(i).gt.nd) nd=ne(i)
 300  continue
      write(16,1001) nr
      write(16,1002) nd
C     Graph-theoretical shape coefficient I2: 
C     refer to: handbook of molecular descriptors,pp390
      pi2=dble(nd-nr)/dble(nr)
      write(16,2007) pi2
C     molecular eccentricity,average eccentricity,eccendirc difference
      avgn=0.0d0
      do 400 i=1,natom
      avgn=avgn+dble(ne(i))
  400 continue
      avgn=avgn/dble(natom)
      nmol=0
      ed=0.0d0
      do 500 i=1,natom
      nmol=nmol+ne(i)
      ed=ed+dabs(dble(ne(i))-avgn  )
 500  continue
      aae=dble(nmol)/dble(natom)
      ed=ed/dble(natom)
      write(16,1003) nmol
      write(16,1004) aae
      write(16,1005) ed
C     the average distance  degree and mean distance degree deviation
      thavg=0.0d0
      do 600 i=1,natom
      thavg=thavg+dble(nth(i))
  600 continue
      thavg=thavg/dble(natom)
      ddd=0.0d0
      do 700 i=1,natom
      ddd=ddd+dabs(dble(nth(i))-thavg)
  700 continue
       ddd=ddd/dble(natom)
C     unipolarity
      nthmin=nth(1)
      do 800 i=2,natom
      if(nth(i).lt.nthmin) nthmin=nth(i)
  800 continue
C     Rouvary index:Irouv
      irouv=0
      do 900 i=1,natom
      irouv=irouv+nth(i)
 900  continue
C     Centralization,variation, dispersion
      ndth=irouv-natom*nthmin
      ndthp=nth(1)-nthmin
      do 1000 i=2,natom
      mm=nth(i)-nthmin
      if(mm.gt.ndthp) ndthp=mm
 1000 continue
      disp1=0.0d0
      do 1100 i=1,natom
      disp1=disp1+dble(is(1,i)*is(1,i))
 1100 continue
      disp1=disp1/dble(natom)
      dispmin=disp1
      do 1200 i=2,natom
      dispi=0.0d0
      do 1300 j=1,natom
      dispi=dispi+dble(is(i,j)*is(i,j))
 1300 continue
      if(dispi.lt.dispmin) dispmin=dispi
 1200 continue
C     PRS index
      prs=0.0d0
      do 1400 i=1,natom
      prs=prs+log(dble(nth(i)))/log(10.0d0)
 1400 continue
      write(16,1006) thavg
      write(16,1007) ddd
      write(16,1008) nthmin
      write(16,1009) irouv
      write(16,2001) ndth
      write(16,2002) ndthp
      write(16,2003) dispmin
      write(16,2004) prs
C
C     Reciprocal distance matrix
C    
      do 1500 i=1,natom
      do 1600 j=1,natom
      if(i.eq.j) then
        rdm(i,j)=0.0d0
      else
        rdm(i,j)=1.0d0/dble(is(i,j))
      end if
 1600 continue
 1500 continue
C     reciprocal distance sum
      do 1700 i=1,natom
      rds(i)=0.0d0
      do 1800 j=1,natom
      rds(i)=rds(i)+rdm(i,j)
 1800 continue
 1700 continue
C     RDSQ index and RDCHI index
      rdsq=0.0d0
      rdchi=0.0d0
      do 1900 i=1,natom-1
      do 2000 j=i+1,natom
      if(is(i,j).ne.1) goto 2000
      aa=dsqrt(rds(i)*rds(j))
      rdsq=rdsq+aa
      rdchi=rdchi+1.0d0/aa
2000  continue 
1900  continue
      write(16,2005)rdsq
      write(16,2006) rdchi
 1001 format('Topological radius:  ',i9)
 1002 format('Topological diameter:',i9)
 1003 format('Eccentricity        :',i9)
 1004 format('Average atom eccentricity: ',e12.4)
 1005 format('Mean eccentricity deviation:',e12.4)
 1006 format('Average distance degree:    ',e12.4)
 1007 format('Mean distance degree deviation:',e12.4)
 1008 format('Unipolarity:        ',i9)
 1009 format('Rouvary index Irouv:',i9)
 2001 format('Centralization:     ',i9)
 2002 format('Variation :         ',i9)
 2003 format('Dispersion:      ',e12.4)
 2004 format('Log of PRS INDEX:',e12.4)
 2005 format('RDSQ ondex:      ',e12.4)
 2006 format('RDCHI index:     ',e12.4)
 2007 format('Graph-theoretical shape coefficient:',e12.4)
      return
      end
C*************************************
C     LogP from topological structure
C     refer to:
C       Milan Soskie, Dejan Plavsic
C       Moldeling the Octanol-Water Partion Coefficients by an
C       optimized molecular connectivity index.
C       J.Chem.Inf.Model.2005,45,930-938
C**************************************
      subroutine logptop(natom,number,numh,ic,vvd,iarom)
      implicit  double precision(a-h,o-z)
      dimension number(natom),ic(natom,natom)
      dimension numh(natom),iarom(natom)
      dimension dopt(natom)
      dimension vvd(natom)
C
      call deltaopt(natom,number,numh,ic,vvd,iarom,xopt)
C
C     neq=i: using equation i in the reference
C     (neq may be : 11,12,13,14,15,16,17,18)
C
       neq=18
C
      if(neq.eq.11) then
         plog=0.816+1.061*xopt
         goto 100
      else
         call indicator(natom,number,numh,ic,iarom,isc,itc,imet,
     & iphyd,ialrin,iconjug,ihg2,ihg3,ihvic,ipg2,ipvic)
      end if
      if(neq.eq.12) then
       plog=0.818 + 1.060*xopt + 0.553*dble(isc)-0.837*dble(itc)
       goto 100
      end if
      if(neq.eq.18) then
       plog=0.829 + 1.055*xopt + 0.58*dble(imet)
       plog=plog + 0.367*dble(iphyd) - 0.627*dble(ialrin)
       plog=plog + 0.454*dble(iconjug) + 0.658*dble(ihg2)
       plog=plog + 1.726*dble(ihg3) + 0.381*dble(ihvic)
       plog=plog + 1.271*dble(ipg2) + 0.605*dble(ipvic)
       goto 100
      end if
C
100   continue
      write(16,1001) xopt
      write(16,1002) plog
1001  format('Optimized 1th connectivity index:',e12.4)
1002  format('Logp from connectivity:',e12.4)
      return 
      end
C*******************************************************
C     Optimized first-order molecular connectivity index
C********************************************************
      subroutine deltaopt(natom,number,numh,ic,vvd,iarom,xopt)
      implicit  double precision(a-h,o-z)
      dimension number(natom),ic(natom,natom)
      dimension numh(natom),iarom(natom)
      dimension dopt(natom)
      dimension vvd(natom)
      dimension nc(natom),linc(natom,8)
      dimension itypec(natom)
C     default values for atoms without in the list of 
C     the table2 in the reference article
C
      do 100 i=1,natom
      dopt(i)=vvd(i)
100   continue
C
C     nc(i): number of non-H atoms bonded to atoms
C     linc(i,j): the jth non-h atom bonded to atom i is atom linc(i,j)
C
      do 200 i=1,natom
        nu=0
        do 300 j=1,natom
         if(ic(i,j).eq.1) then
          nu=nu+1
          linc(i,nu)=j
         end if
         nc(i)=nu
300   continue
200   continue
C
C     atom type classification according to the table 3
C
      do 400 i=1,natom
C C atom
      if(number(i).ne.6) goto 401
C     for other cases that are not included in the following
      itypec(i)=9
C
C(1). C in -C-trbond  triple bond
      if(nc(i).eq.1.and.numh(i).eq.0) itypec(i)=1        
C(2). -CH2-
      if(nc(i).eq.2.and.numh(i).eq.2) itypec(i)=2        
C(3). -CH--
      if(nc(i).eq.3.and.numh(i).eq.1) itypec(i)=3        
C(4). --C--
      if(nc(i).eq.4.and.numh(i).eq.0) itypec(i)=4        
C(5)  C=O : bonded to =O
      do 403 j=1,nc(i)
      nj=linc(i,j)
      if(number(nj).eq.8.and.numh(nj).eq.0.and.nc(nj).eq.1) itypec(i)=5
403   continue    
C(0). C atom in aromatic ring
      if(iarom(i).eq.1) itypec(i)=0
401   continue
C
400   continue
C   
C     Assignment optimized delta values
C
      do 900 i=1,natom
C  C atom
      if(number(i).ne.6) goto 910
C   -CH contected to triple c-c bond
      if(nc(i).eq.1.and.numh(i).eq.1) dopt(i)=-3.00        
      if(iarom(i).eq.1) dopt(i)=4.30
      if(iarom(i).eq.1.and.nc(i).gt.2) dopt(i)=5.30
      goto 900
910   continue
C  O atom
      if(number(i).ne.8) goto 920
C  O atom in -OH
      if(numh(i).ne.1) goto 921
      nj=linc(i,1)
      if(number(nj).ne.6) goto 921
C (0) ar-OH
      if(itypec(nj).eq.0) dopt(i)=-0.35
C  -CH2-OH
      if(itypec(nj).eq.2) dopt(i)=-0.20
C  -CH--OH
      if(itypec(nj).eq.3) dopt(i)=-0.10
C  -C--OH
      if(itypec(nj).eq.4) dopt(i)=-0.07
C  CO-OH
      if(itypec(nj).eq.5) dopt(i)=-0.25
      goto 900
921   continue
C  O atom in  -O-
      if(nc(i).ne.2) goto 922
      dopt(i)=-1.20
      do 923 j=1,2
      nj=linc(i,j)
      if(number(nj).eq.6.and.itypec(nj).eq.5) dopt(i)=-1.40
923   continue
      goto 900
922   continue
C  O atom in C=O
      if(nc(i).ne.1.and.numh(i).ne.0) goto 924
      dopt(i)=-0.06
      j1=linc(i,1)
      do 925 j=1,nc(j1)
      nj=linc(j1,j)
      if(nj.eq.i) goto 925
      if(number(nj).eq.8.and.numh(nj).eq.1) dopt(i)=-1.25
C      if(number(nj).eq.7.and.numh(nj).eq.2) dopt(i)=-1.25
      if(number(nj).eq.7) dopt(i)=-1.25
      if(number(nj).eq.8.and.numh(nj).eq.0) dopt(i)=-2.40
925   continue
924   continue
      goto 900
920   continue
C  N atom
      if(number(i).ne.7) goto 930
C  -NH2
      if(numh(i).ne.2) goto 931
      nj=linc(i,1)
      if(number(nj).eq.6.and.itypec(nj).eq.2) dopt(i)=-0.20
      if(number(nj).eq.6.and.itypec(nj).eq.3) dopt(i)=-0.11
      if(number(nj).eq.6.and.itypec(nj).eq.4) dopt(i)=-0.07
      if(number(nj).eq.6.and.itypec(nj).eq.5) dopt(i)=-0.06
931   continue
C -NH
      if(numh(i).ne.1) goto 932
      if(nc(i).ne.2) goto 932
      j1=linc(i,1)
      j2=linc(i,2)
      if(itypec(j1).eq.2.and.itypec(j2).eq.2) dopt(i)=-0.70
      if(itypec(j1).eq.3.and.itypec(j2).eq.3) dopt(i)=-0.34
      if(itypec(j1).eq.5.or.number(j2).eq.6) dopt(i)=-0.49
      if(number(j1).eq.5.and.itypec(j2).eq.5) dopt(i)=-0.49
932   continue
C  -N
      if(numh(i).ne.0) goto 933
      if(nc(i).eq.1.and.number(linc(i,1)).eq.6) dopt(i)=-0.09
      if(nc(i).eq.3) then
        dopt(i)=-2.00
        do 934 j=1,3
        nj=linc(i,j)
        if(number(nj).eq.6.and.itypec(nj).eq.5) dopt(i)=-0.97
934     continue
      end if
933   continue
C
       goto 900
930   continue
C   F
      if(number(i).ne.9)  goto 940 
      dopt(i)=-1.00
      goto 900
940   continue
C  Cl
      if(number(i).ne.17) goto 950
      nj=linc(i,1)
      if(number(nj).eq.6.and.itypec(nj).eq.2) dopt(i)=-100.0
      if(number(nj).eq.6.and.itypec(nj).eq.3) dopt(i)=-20.0
      if(number(nj).eq.6.and.itypec(nj).eq.4) dopt(i)=-11.0
      if(number(nj).eq.6.and.itypec(nj).eq.0) dopt(i)=0.55
      goto 900
950   continue
C Br
      if(number(i).ne.35) goto 960
      nj=linc(i,1)
      if(number(nj).eq.6.and.itypec(nj).eq.0) then
        dopt(i)=0.41
      else
        dopt(i)=50.0
      end if
      goto 900
960   continue
      if(number(i).eq.53) dopt(i)=3.00
900   continue    
C
C   The optimized first-order molecular connectivity index
C
      xopt=0.0d0
      xv=0.0d0
      do 1100 i=2,natom
      do 1200 j=1,i-1
      if(ic(i,j).eq.0) goto 1200
      if(dopt(i).lt.0.0.or.dopt(j).lt.0.0) then
         sign=-1.0d0
      else
         sign=1.0d0
      end if
C
C
      a1=dabs(dopt(i)*dopt(j))
      a2=dsqrt(a1)
      a3=1.0/a2
      xopt=xopt+sign*a3
      xv=xv+1.0d0/dsqrt(vvd(i)*vvd(j))
1200  continue
1100  continue
      return 
      end
      subroutine indicator(natom,number,numh,ic,iarom,isc,itc,imet,
     & iphyd,ialrin,iconjug,ihg2,ihg3,ihvic,ipg2,ipvic)
      implicit  double precision(a-h,o-z)
      dimension number(natom),ic(natom,natom)
      dimension numh(natom),iarom(natom)
      dimension nc(natom),linc(natom,8)
      dimension itypec(natom)
C
C     Initialization of the indicator variables
C
      isc=0
      itc=0
      imet=0
      iphyd=0
      ialrin=0
      iconjug=0
      ihg2=0
      ihg3=0
      ihvic=0
      ipg2=0
      ipvic=0
C
C     nc(i): number of non-H atoms bonded to atoms
C     linc(i,j): the jth bonded to atom i is atom linc(i,j)
C
      do 200 i=1,natom
        nu=0
        do 300 j=1,natom
         if(ic(i,j).eq.1) then
          nu=nu+1
          linc(i,nu)=j
         end if
         nc(i)=nu
300   continue
200   continue
C
C     atom type classification according to the table 3
C
      do 400 i=1,natom
C C atom
      if(number(i).ne.6) goto 401
C     for other cases that are not included in the following
      itypec(i)=9
C
C(1). C in -C-trbond  triple bond
      if(nc(i).eq.1.and.numh(i).eq.0) itypec(i)=1        
C(2). -CH2-
      if(nc(i).eq.2.and.numh(i).eq.2) itypec(i)=2        
C(3). -CH--
      if(nc(i).eq.3.and.numh(i).eq.1) itypec(i)=3        
C(4). --C--
      if(nc(i).eq.4.and.numh(i).eq.0) itypec(i)=4        
C(5)  C=O : bonded to =O
      do 403 j=1,nc(i)
      nj=linc(i,j)
      if(number(nj).eq.8.and.numh(nj).eq.0.and.nc(nj).eq.1) itypec(i)=5
403   continue    
C(0). C atom in aromatic ring
      if(iarom(i).eq.1) itypec(i)=0
401   continue
C
400   continue
C   
C     judge the indicator
C
C     ISC and ITC
C           
      do 500 i=1,natom
C     R-OH
      if(number(i).eq.8.and.numh(i).eq.1) then
      nj=linc(i,1)
        if(number(nj).eq.6.and.nc(nj).eq.2) then
          isc=1
        end if 
        if(number(nj).eq.6.and.nc(nj).eq.3) then
          itc=1
        end if 
      end if
C     R-NH2
      if(number(i).eq.7.and.numh(i).eq.2) then
      nj=linc(i,1)
        if(number(nj).eq.6.and.nc(nj).eq.2) isc=1
        if(number(nj).eq.6.and.nc(nj).eq.3) itc=1
      end if
500   continue
C IMET
       do 600 i=1,natom
       if(number(i).eq.6.and.numh(i).eq.3) then
       nj=linc(i,1)
       nb=number(nj)
C     if(nb.eq.7.or.nb.eq.8.or.nb.eq.9) imet=1
       if(nb.eq.7.or.nb.eq.8.or.nb.eq.9) imet=imet+1
       end if
600    continue
C
C  IPHYD
C
       do 650 i=1,natom
       if(number(i).ne.6) goto 650
       if(numh(i).lt.1) goto 650
       if(itypec(i).ne.5) goto 650 
       nx=0
       do 660 j=1,nc(i)
       nj=linc(i,j)
       if(number(nj).ne.6) nx=nx+1
660    continue
       if(nx.eq.2) iphyd=1
650    continue
C
C  IALRIN
C
C     number of bonds
      nbd=0
      do 710 i=1,natom-1
      do 720 j=i+1,natom
      if(ic(i,j).eq.1) nbd=nbd+1
720   continue
710   continue
C     number of rings according to the  Euler rule
C
      nring=nbd-natom+1
      if(nring.eq.0) goto 700
      ialrin=1
      do 730 i=1,natom
      if(iarom(i).eq.1) ialrin=0
730   continue
700   continue
C
C  ICONJUG
C
C   O=C-C(sp2)
C
       do 800 i=1,natom
       if(number(i).ne.6) goto 800
       if(itypec(i).ne.5) goto 800
       do 810 j=1,nc(i)
       nj=linc(i,j)
       if(number(nj).ne.6) goto 810
       ncon=nc(nj)+numh(nj)
       if(ncon.lt.4) iconjug=1
810    continue
800    continue
C      IHG2 and IHG3
C
      do 900 i=1,natom
      nxi=0
      do 910 j=1,nc(i)
      nj=linc(i,j)
      nbj=number(nj)
      nxp=(nbj-9)*(nbj-17)*(nbj-35)*(nbj-53)
      if(nxp.eq.0) nxi=nxi+1
910   continue
C      if(nxi.eq.2) ihg2=1
C      if(nxi.eq.3) ihg3=1
      if(nxi.eq.2) ihg2=ihg2+1
      if(nxi.eq.3) ihg3=ihg3+1
900   continue
C
C     IHVIC
C
      do 1000 i=1,natom
         nxi=0
         do 1010 ii=1,nc(i)
         ni=linc(i,ii)
         nbi=number(ni)
         nxpi=(nbi-9)*(nbi-17)*(nbi-35)*(nbi-53)
         if(nxpi.eq.0) nxi=nxi+1
1010     continue 
      if(nxi.eq.0) goto 1000
      do 1020 j=1,nc(i)
      nj=linc(i,j)
        nxj=0   
        do 1030 jj=1,nc(nj)
         njj=linc(nj,jj)
         nbj=number(njj)
         nxpj=(nbj-9)*(nbj-17)*(nbj-35)*(nbj-53)
         if(nxpj.eq.0) nxj=nxj+1
1030     continue 
         if(nxj.ne.0) ihvic=1 
1020  continue
1000  continue
C
C     IPG2 
C
      do 1100 i=1,natom
      if(number(i).ne.6) goto 1100
      idon=0
      do 1110 j=1,nc(i)
      nj=linc(i,j)
      ncon=nc(nj)+numh(nj)
      if(number(nj).eq.8.and.ncon.eq.2) idon=idon+1
      if(number(nj).eq.9) idon=idon+1
1110  continue
      if(idon.ge.2) ipg2=1
1100  continue
C
C     IPVIC 
C
      do 1200 i=1,natom
      if(number(i).ne.6) goto 1200
      idon=0
      do 1210 j=1,nc(i)
      nj=linc(i,j)
      jcon=nc(nj)+numh(nj)
      if(number(nj).eq.8.and.jcon.eq.2) idon=idon+1
      if(number(nj).eq.9) idon=idon+1
1210  continue
      if(idon.eq.0) goto 1200
      do 1220 k=1,nc(i)
        nk=linc(i,k)
        if(number(nk).ne.6) goto 1220
        kdon=0
        do 1230 m=1,nc(nk)
        nm=linc(nk,m)
        mcon=nc(nm)+numh(nm)
        if(number(nm).eq.8.and.mcon.eq.2) kdon=kdon+1
        if(number(nm).eq.9) kdon=kdon+1
1230    continue
      if(kdon.gt.0) ipvic=1
1220  continue
C
1200  continue
      if(idon.ge.2) ipg2=1
C
1300  continue
C
      return
      end
C**************************************
C     Protein-protein interaction descriptors constructed from  two protein descriptors
C     methods:
C    methodpp=1: out put two vectors for protein A and B as V=(Va,Vb)  for A+B and V=(Vb,Va) for B+A
C    methodpp=2: V=(Vai+Vbi,Vai*Vbi)
C**************************************
      subroutine ppides(maxv,methodpp)
      implicit double precision(a-h,o-z)
      character name*100,namea*100,nameb*200,line*201,outfile*200
      dimension v(maxv),va(maxv),vb(maxv)
C
      rewind (9)
 100  read(9,'(a101)',end=900) line
      name=line(2:101)
      idm=100
      call trans(name,idm,outfile,ii)
      close (18)
      open(18,file='output-temp/'//outfile(1:ii)//'.dat'
     & ,status='unknown',err=110)

      read(9,*,err=200) nv
      write(18,*) nv
      read(9,*,err=400) (v(i),i=1,nv)
      write(18,*) (v(i),i=1,nv)
C
      goto 100
C
110   write(6,*) 'Error in opening file '
     &,'output-temp/'//outfile(1:ii)//'dat'
      write(6,*) 'Make directory: output-temp?'
      return
200   write(6,*)'Error in reading number of descriptors in ',name
      return
400   write(6,*) 'Error in reading descriptors in ',name
      write(6,*) 'Stop: because of error in output-des.dat'
300   return 
C
900   continue
C
C     Reading names in interaction protein sequences
C
700   read(8,'(a201)',end=1900) line
      ind=index(line,'+')
      if(ind.lt.2) then
        write(6,*) 'Error in input-ppi.dat: no sign + found in'
        write(6,*) line
        goto 700
      end if
C
      namea=' '
      nameb=' '
      namea=line(1:ind-1)
      nameb=line(ind+1:201)
C
      close(18)
      idm=ind-1
      call trans(namea,idm,outfile,iia)
      open(18,file='output-temp/'//outfile(1:iia)//'.dat'
     & ,status='old',err=710)
      namea=outfile(1:iia)
      read(18,*) na
      read(18,*) (va(i),i=1,na)
      goto 720
710   write(6,*) 'error in opening file '
     & ,'output-temp/'//outfile(1:ii)//'.dat'
      write(6,*) 'Abnormal termination for ',line
      goto 700
720   continue
C
      close(17)
      idm=201-ind
      call trans(nameb,idm,outfile,iib)
      open(17,file='output-temp/'//outfile(1:iib)//'.dat'
     &,status='old',err=730)
      nameb=outfile(1:iib)
      read(17,*) nb
      read(17,*) (vb(i),i=1,nb)
      goto 740
730   write(6,*) 'error in opening file '
     & ,'output-temp/'//outfile(1:ii)//'.dat'
      write(6,*) 'Abnormal termination for ',line
      goto 700
740   continue
C
      nunit=19
C
      if(methodpp.ne.1) goto 1300
C
      nab=na+nb
C
      write(6,*) namea(1:iia),'  +    ',nameb(1:iib)
C
      write(19,*)'>'//namea(1:iia),'  +    ',nameb(1:iib)
      call methodpp1(na,nb,nab,va,vb,nunit)
      write(19,*) '>'//nameb(1:iib),'  +    ',namea(1:iia)
      call methodpp1(nb,na,nab,vb,va,nunit)
C
      write(6,*) 'Normal Termination'
      write(6,*)
C     
1300  continue
C
      if(methodpp.ne.2) goto 1400
C
      nab=2*(na+nb)
C
      write(6,*) namea(1:iia),'  +  ',nameb(1:iib)
C
      write(19,*)'>'//namea(1:iia),'  +   ',nameb(1:iib)
      call methodpp2(na,nb,nab,va,vb,nunit)
C
      write(6,*) 'Normal Termination'
      write(6,*)
C     
1400  continue
C
      if(methodpp.ne.3) goto 1500
C
      nab=na*nb
C
      write(6,*) namea(1:iia),'  +  ',nameb(1:iib)
C
      write(19,*)'>'//namea(1:iia),'  +  ',nameb(1:iib)
      call methodpp3(na,nb,nab,va,vb,nunit)
C
      write(6,*) 'Normal Termination'
      write(6,*)
C
1500  continue
C
       goto 700
C
1900  continue
C
      return
      end
C***********************************
C      construction of vector from two protein vectors VA and VB
C      methodpp=1: V=(VA,VB)
C**********************************************
      subroutine methodpp1(na,nb,nab,va,vb,nunit)
      implicit double precision(a-h,o-z)
      dimension va(na),vb(nb),vab(nab)
C
      do  100 i=1,na
        vab(i)=va(i)
100   continue
      do 200 i=1,nb
        vab(i+na)=vb(i)
200   continue
C
      write(nunit,*) na+nb
C
      nn=0
300   nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.(na+nb)) n2=na+nb
      write(nunit,'(10E14.4)') (vab(i),i=n1,n2)
      if(n2.lt.(na+nb)) goto 300   
C
      return
      end 
C***********************************
C      construction of vector from two protein vectors VA and VB with same dimension: na=nb
C      methodpp=2: V=(Vai*Vbi,Vai+Vbi)
C**********************************************
      subroutine methodpp2(na,nb,nab,va,vb,nunit)
      implicit double precision(a-h,o-z)
      dimension va(na),vb(nb),vab(nab)
      if(na.ne.nb) then
         write(6,*) 'The two vectors have not the same dimension'
         write(6,*) 'Abnormal termination'
         return
      end if
C
      do  100 i=1,na
        vab(i)=va(i)*vb(i)
100   continue
      do 200 i=1,nb
        vab(i+na)=va(i)+vb(i)
200   continue
C
      write(nunit,*) na+nb
C
      nn=0
300   nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.(na+nb)) n2=na+nb
      write(nunit,'(10E14.4)') (vab(i),i=n1,n2)
      if(n2.lt.(na+nb)) goto 300   
C
      return
      end 
C***********************************
C      construction of vector from two protein vectors VA and VB
C      methodpp=3: V={Vij}={(Vai*Vbj,i=1,na,j=1,nb)}
C**********************************************
      subroutine methodpp3(na,nb,nab,va,vb,nunit)
      implicit double precision(a-h,o-z)
      dimension va(na),vb(nb),vab(nab)
C
      nv=0
      do  100 i=1,na
       do 200 j=1,nb
       nv=nv+1
        vab(nv)=va(i)*vb(j)
200   continue
100   continue
C
      write(nunit,*) nab
C
      nn=0
300   nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.(nab)) n2=nab
      write(nunit,'(10E14.4)') (vab(i),i=n1,n2)
      if(n2.lt.(nab)) goto 300   
C
      return
      end 
C**************************************
C     Protein-ligand interaction descriptors constructed from protein descriptor Vector Va and ligand descriptor vector Vb
C     methods:
C    methodpl=1: output one vector for protein A and ligand B  V=as (Va,Vb)
C    methodpl=2: V={Vij}={Vai*Vbj,i=1,na, j=1,nb)
C**************************************
      subroutine plides(maxp,maxl,methodpl)
      implicit double precision(a-h,o-z)
      character name*100,namea*100,nameb*200,line*201,outfile*200
      dimension vp(maxp),vl(maxl),va(maxp),vb(maxl)
C
C     reading data for protein from unit 9:output-protein.dat and saving data for each sequence in a individual file
C
      rewind (9)
 100  read(9,'(a101)',end=109) line
      name=line(2:101)
      idm=100
      call trans(name,idm,outfile,ii)
      close (18)
      open(18,file='output-temp/'//outfile(1:ii)//'.dat'
     & ,status='unknown',err=101)

      read(9,*,err=102) np
      write(18,*) np
      read(9,*,err=103) (vp(i),i=1,np)
      write(18,*) (vp(i),i=1,np)
C
      goto 100
C
101   write(6,*) 'Error in opening file '
     &,'output-temp/'//outfile(1:ii)//'dat'
      write(6,*) 'Make directory: output-temp?'
      return
102   write(6,*)'Error in reading number of descriptors in ',name
      return
103   write(6,*) 'Error in reading descriptors in ',name
      write(6,*) 'Stop: because of error in output-des.dat'
      return 
C
109   continue
C
C     reading data for ligand from unit 4:output-ligand.dat and saving data for each ligand in a individual file
C
      rewind (4)
 200  read(4,'(a41)',end=209) line
      name=line(2:41)
      idm=40
      call trans(name,idm,outfile,ii)
      close (18)
      open(18,file='output-temp/'//outfile(1:ii)//'.dat'
     & ,status='unknown',err=201)

      read(4,*,err=202) nl
      write(18,*) nl
      read(4,*,err=203) (vl(i),i=1,nl)
      write(18,*) (vl(i),i=1,nl)
C
      goto 200
C
201   write(6,*) 'Error in opening file '
     &,'output-temp/'//outfile(1:ii)//'dat'
      write(6,*) 'Make directory: output-temp?'
      return
202   write(6,*)'Error in reading number of descriptors in ',name
      return
203   write(6,*) 'Error in reading descriptors in ',name
      write(6,*) 'Stop: because of error in output-des.dat'
      return 
C
209   continue
C
C     Reading names in protein-ligand interaction input: input-
C
700   read(14,'(a201)',end=1900) line
      ind=index(line,'+')
      if(ind.lt.2) then
        write(6,*) 'Error in input-ppi.dat: no sign + found in'
        write(6,*) line
        goto 700
      end if
C
      namea=' '
      nameb=' '
      namea=line(1:ind-1)
      nameb=line(ind+1:201)
C
      close(18)
      idm=ind-1
      call trans(namea,idm,outfile,iia)
      open(18,file='output-temp/'//outfile(1:iia)//'.dat'
     & ,status='old',err=710)
      read(18,*) na
      read(18,*) (va(i),i=1,na)
      goto 720
710   write(6,*) 'error in opening file for protein '
     & ,'output-temp/'//outfile(1:iia)//'.dat'
      write(6,*) 'Abnormal termination for ',line
      goto 700
720   continue
      namea=outfile(1:iia)
C
      close(17)
      idm=201-ind
      call trans(nameb,idm,outfile,iib)
      open(17,file='output-temp/'//outfile(1:iib)//'.dat'
     &,status='old',err=730)
      read(17,*) nb
      read(17,*) (vb(i),i=1,nb)
      goto 740
730   write(6,*) 'error in opening file for ligand '
     & ,'output-temp/'//outfile(1:iib)//'.dat'
      write(6,*) 'Abnormal termination for ',line
      goto 700
740   continue
      nameb=outfile(1:iib)
C
      nunit=15
C
      if(methodpl.ne.1) goto 1300
C
      nab=na+nb
C
      write(6,*) namea(1:iia),'  +  ',nameb(1:iib)
C
      write(15,*) '>'//namea(1:iia),'  +  ',nameb(1:iib)
      call methodpl1(na,nb,nab,va,vb,nunit)
      write(6,*) 'Normal Termination'
      write(6,*)
C     
1300  continue
C
      if(methodpl.ne.2) goto 1400
C
      nab=na*nb
C
      write(6,*) namea(1:iia),'  +  ',nameb(1:iib)
C
      write(15,*) '>'//namea(1:iia),'  +  ',nameb(1:iib)
      call methodpl2(na,nb,nab,va,vb,nunit)
      write(6,*) 'Normal Termination'
      write(6,*)
C     
1400  continue
C
       goto 700
C
1900  continue
C
      return
      end
C***********************************
C      construction of vector from protein vector VA and  and ligand vector VB
C      method1: V=(VA,VB)
C**********************************************
      subroutine methodpl1(na,nb,nab,va,vb,nunit)
      implicit double precision(a-h,o-z)
      dimension va(na),vb(nb),vab(nab)
C
      do  100 i=1,na
        vab(i)=va(i)
100   continue
      do 200 i=1,nb
        vab(i+na)=vb(i)
200   continue
C
      write(nunit,*) na+nb
C
      nn=0
300   nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.(na+nb)) n2=na+nb
      write(nunit,'(10E14.4)') (vab(i),i=n1,n2)
      if(n2.lt.(na+nb)) goto 300   
C
      return
      end 
C***********************************
C      construction of vector from protein vector VA and ligand vector VB
C      methodpl=2: V={Vij}={(Vai*Vbj,i=1,na,j=1,nb)}
C**********************************************
      subroutine methodpl2(na,nb,nab,va,vb,nunit)
      implicit double precision(a-h,o-z)
      dimension va(na),vb(nb),vab(nab)
C
      nv=0
      do  100 i=1,na
       do 200 j=1,nb
       nv=nv+1
        vab(nv)=va(i)*vb(j)
200   continue
100   continue
C
      write(nunit,*) nab
C
      nn=0
300   nn=nn+1
      n1=(nn-1)*10+1
      n2=nn*10
      if(n2.gt.(nab)) n2=nab
      write(nunit,'(10E14.4)') (vab(i),i=n1,n2)
      if(n2.lt.(nab)) goto 300   
C
      return
      end 
C*****************************************************************************************************
C     Total amino acid properties
C--------------------------------------
C     refer to:
C     M.Michael Gromiha, Makiko Suwa.
C     Influence of amino acid properties for discriminating outer membrane proteins at better accuracy
C     Biochimica of Biophysica Acta, 2006,1764,1493-1497.
C--------------------------------------
C      seq: original protein sequence
C      len: length of the seqence
C      np: if or not to print the feature name on unit 10
C      nseg: number of the segment that the sequence will be cutted in to
C      napp: number of the amino acid properties to be calculated
C      iappn: the serial number of the amino acid index in index-database.dat
C      aadb: all the properties in index-database.dat
C      xpp(i,j): the sum of ith property for segment j
C******************************************************************************************************
      subroutine appdes(seq,len,np,nseg,napp,iappn,xapp)
      implicit double precision(a-h,o-z)
      character seq*10000,aa*20,line*80
      character*80 title
      common/aaindex/aadb(1000,20)
      dimension p(20)
      dimension iappn(1000)
      dimension ip(len)
      dimension lenseg(nseg)
      dimension xapp(napp,nseg)
C
      aa='ARNDCQEGHILKMFPSTWYV'
C
C     poisition of the amino acid of the sequence in sring 
C     aa='ARNDCQEGHILKMFPSTWYV'
C
      do 10 i=1,len
      ip(i)=index(aa,seq(i:i))
10    continue      
C
      write(6,*)
      if(nseg.eq.1) then
        write(6,*)'(G9) Total amino acid property'
      else
        write(6,*)'(G9) Total amino acid property for sequence segments'
      end if
C
C     dividing the sequences into segments
C
      nav=len/nseg
      do 30 i=1,nseg
        lenseg(i)=nav
30    continue
C
      nr=len-nav*nseg
      if(nr.gt.0) then
        do 40 i=1,nr
          lenseg(i)=lenseg(i)+1
40      continue
      end if
C
      write(6,*)
      do 90 kk=1,napp
      iadk=iappn(kk)
      do 102 j=1,20
       p(j)=aadb(iadk,j)
102   continue
C
C
C     Normalization of  the properties
C
      pmin=p(1)
      pmax=p(1)
      do 50  i=2,20
      if(p(i).gt.pmax) pmax=p(i)
      if(p(i).lt.pmin) pmin=p(i)
50    continue
      do 60 i=1,20
      p(i)=(p(i)-pmin)/(pmax-pmin)
60    continue
C
      n2=0
      do 70 ns=1,nseg
      xapp(kk,ns)=0.0d0
      n1=n2+1
      n2=n1+lenseg(ns)-1
C
      if(nseg.gt.1) then
C         write(6,*)
C         write(6,*) ' In segment ', ns
C         write(6,*) seq(n1:n2)
      end if
C
C     Total properties in the segment
C
      do 80 i=n1,n2
      ipi=ip(i)
      xapp(kk,ns)=xapp(kk,ns)+p(ipi)
80    continue
      write(6,*) 'Total property by aaindex in segment ',ns
     & , 'for property ',iadk,':', xapp(kk,ns)
      if(nseg.gt.1) then 
      if(np.eq.1) write(10,*) ' total property of aaindex '
     & ,iadk,' for segment ' ,ns
      else
      if(np.eq.1) write(10,*) ' total property of aaindex ',iadk
      end if
C
70    continue
C
90    continue
C
      return
      end 
C*****************************************************
C     natom: number of atoms including H and Heavy atoms
C     natoms: number of atoms including only heavy atoms
C*********************************************************
      subroutine numbatom(natom,natoms,molname,iend,ierror)
      implicit  double precision(a-h,o-z)
      character line*80,molname*40
      character name*2,name3*3 
      ierror=0
      close(11)
      open(11,file='output-temp/output-temp.mol',status='unknown',err=1)
      goto 2
1     open(11,file='output-temp.mol',status='unknown')
2     continue
      read(7,'(A40)',end=90) molname
      read(7,'(A80)') line
      read(7,'(A80)') line
      read(7,'(2i3)') natom,numbond
      write(11,*) natom,numbond
      if(natom.eq.0) goto 30
C
      natoms=0
      do 10 i=1,natom
      read(7,1000) x,y,z,name3
       if(name3(1:1).eq.' ') then
             name(1:2)=name3(2:3)
       else
             name(1:2)=name3(1:2)
       end if
      call atomnumb(name,numb,ierror)
      if(numb.ne.1) natoms=natoms+1
      if(ierror.eq.1) goto 30
      write(11,*) numb
 10   continue
C
      do 20 i=1,numbond
      read(7,'(3i3)') k1,k2,k3
      write(11,*) k1,k2,k3
 20   continue
C
C
 30   read(7,'(A)',end=80) line
      if(line(1:4).ne.'$$$$') goto 30
      iend=0
      return
  80  iend=1
      return
C
  90  iend=2
C
1000  format(3f10.4,1x,a3)
      return
      end
C****************************************
C     Natom: number of atoms
C     natoms: number of heavy atoms
C****************************************
      subroutine  topmol(natom,natoms,ibcut,ierror)
      implicit  double precision(a-h,o-z)
      dimension ic(natoms,natoms),numh(natoms),number(natoms)
      dimension is(natoms,natoms),nbond(natoms,natoms)
      
C
      call readmol(natom,natoms, number,ic,numh,nbond)
C
C      Topological distance between atoms for H-deleted model
C
      call topdis(natoms,ic,is)
C
C     check if there are isolated atoms or if the molecules has more than one fragments
C
      call check(natoms,ic,is,ierror)
      if(ierror.eq.1) then
      write(6,*) 'Stop: abnormal termination because of input steucture'
        return
      end if
C
       call top(natoms,ic,is,number,numh,nbond,ibcut)
C
      return
      end

      subroutine readmol(natom,natoms, numbers,intercons,numhs,nbonds)
      implicit double precision(a-h,o-z)
      dimension number(natom),numbers(natoms)
      dimension intercon(natom,natom),intercons(natoms,natoms)
      dimension ind(natoms)
      dimension numh(natom), numhs(natoms)
      dimension nbond(natom,natom),nbonds(natoms,natoms)
      character*80 title,line
C
      rewind(11)
      read(11,*) n1,numbond
      do 10 i=1,natom
      read(11,*) number(i)
 10   continue
C
      do 40 i=1,natom
      do 50 j=1,natom
      intercon(i,j)=0
      nbond(i,j)=0
 50   continue
 40   continue
      do 60 i=1,numbond
      read(11,*,end=70,err=70) k1,k2,k3
      intercon(k1,k2)=1
      intercon(k2,k1)=1
      nbond(k1,k2)=k3
      nbond(k2,k1)=k3
 60   continue
      goto 75
 70   ierror=1
 75   continue
C    
      num=0
      do 300 i=1,natom
      if(number(i).ne.1) then
        num=num+1
        numbers(num)=number(i)
        ind(num)=i
       endif
 300  continue
C
      if(num.ne.natoms) then
       write(6,*) 'Error in reading heavy atoms'
       ierror=1
      end if
C
      do 80 i=1,natoms
      do 90 j=1,natoms
      intercons(i,j)=0
      nbonds(i,j)=0
 90   continue
 80   continue
C
      do 110 i=1,natoms
      do 120 j=1,natoms
      k1=ind(i)
      k2=ind(j)
      intercons(i,j)=intercon(k1,k2)  
      nbonds(i,j)=nbond(k1,k2)
 120  continue
 110  continue
C     
      do 130 i=1,natom
       numh(i)=0
       do 140 j=1,natom
        if(intercon(i,j).eq.1.and.number(j).eq.1) numh(i)=numh(i)+1
140    continue
130    continue
C
        do 150 i=1,natoms
          k=ind(i)
          numhs(i)=numh(k)
150     continue
C
      return
      end
C*******************************************
C     check if there are isolated atoms
C*******************************************
C     ic(i,j): adjacent matrix
C     is(i,j): topological distance matrix
C     nc(i): connectivity
C------------------------------------------
      subroutine check(natom,ic,is,ierror)
      implicit  double precision(a-h,o-z)
      dimension ic(natom,natom),nc(natom)
      dimension is(natom,natom)
      ierror=0
      do 100 i=1,natom
      nc(i)=0
      do 200 j=1,natom
      if(ic(i,j).eq.1) nc(i)=nc(i)+1
 200  continue
 100  continue
C
      do 300 i=1,natom
      if(nc(i).eq.0) then
        ierror=1
        write(6,*) 'atom ',i, 'is isolated.'
      end if
 300  continue
      do 400 i=1,natom-1
       do 500 j=i+1,natom
        if(is(i,j).eq.0) then
           ierror=1
           write(6,*) 'more than one fragment'
          goto 600
        end if
500   continue
400   continue
600   continue
      return
      end 
C******************************************************
C    constitutational descriptors
C    H-deleted model
C------------------------------------------------
      subroutine constitution(natom,number,intercon,numh)
C
      implicit double precision(a-h,o-z)
      dimension intercon(natom,natom)
      dimension number(natom),numh(natom)
C
      write(16,*) '  Constitutational Descriptors'
C
      nh=0
      do 10 i=1,natom
      nh=nh+numh(i)
10    continue
      write(16,101) natom+nh
      write(16,102) natom
101   format('Number of Atoms:',i4)
102   format('Number of Heavy atoms:',i4)
C
C  number of specified elements
      nb=0
      nc=0
      nn=0
      no=0
      nf=0
      np=0
      ns=0
      ncl=0
      nbr=0
      ni=0
      do 100 i=1,natom
      numb=number(i)
      if(numb.eq.5) nb=nb+1
      if(numb.eq.6) nc=nc+1
      if(numb.eq.7) nn=nn+1
      if(numb.eq.8) no=no+1
      if(numb.eq.9) nf=nf+1
      if(numb.eq.15) np=np+1
      if(numb.eq.16) ns=ns+1
      if(numb.eq.17) ncl=ncl+1
      if(numb.eq.35) nbr=nbr+1
      if(numb.eq.53) ni=ni+1
 100  continue
      write(16,103) nh
103   format('Number of H atoms: ',i4) 
      write(16,104) nb
104   format('Number of B atoms: ',i4) 
      write(16,105) nc
105   format('Number of C atoms: ',i4) 
      write(16,106) nn
106   format('Number of N atoms: ',i4) 
      write(16,107) no
107   format('Number of O atoms: ',i4) 
      write(16,108) nf
108   format('Number of F atoms: ',i4) 
      write(16,109) np
109   format('Number of P atoms: ',i4) 
      write(16,110) ns
110   format('Number of S atoms: ',i4) 
      write(16,120) ncl
120   format('Number of Cl atoms: ',i4) 
      write(16,130) nbr
130   format('Number of Br atoms: ',i4) 
      write(16,140) ni
140   format('Number of I atoms: ',i4) 
C
C number of bonds (including H)
      nb=0
      do 200 i=1,natom-1
      do 300 j=i+1,natom
      if(intercon(i,j).ne.0) nb=nb+1
 300  continue
 200  continue
      write(16,150) nb+nh
 150  format('Number of Bonds: ',i4)
      write(16,160) nb
 160  format('Number of non-H Bonds: ',i4)
C
C Number of rings
      nr=nb-(natom-1)
      write(16,170) nr
 170  format('Number of rings: ',i4)
C  Molecular weight
      call weight(natom,number,numh)
C
C Number of H-bond donors(any N O F) and acceptors (any N-H,O-H,F-H)
      call hbond(natom,number,intercon,numh)
C
       return
       end
C**************************************************
C      Number of H-bond donors and H-acceptors
C      H-deleted model
C**************************************************
      subroutine  hbond(natom,number,intercon,numh)
      implicit double precision(a-h,o-z)
      dimension intercon(natom,natom),number(natom)
      dimension numh(natom)
      ndonnor=0
      naccept=0
      do 100 i=1,natom
      nb=number(i)
      if(nb.eq.7.or.nb.eq.8.or.nb.eq.9) then
         naccept=naccept+1
         if(numh(i).ge.1) ndonnor=ndonnor+1
      endif
 100  continue
      write(16,1010) ndonnor        
      write(16,1020) naccept      
1010  format('Number of H-bond donnor :  ',I5)
1020  format('Number of H-bond acceptor: ',I5)
      return
      end
C----------------------------------------------
C  Molecule Weight for all-atom model
C     H-deleted model
C-----------------------------------------------
      subroutine  weight(natom,number,numh)
      implicit  double precision(a-h,o-z)
      dimension aw(100),number(natom),numh(natom)
      call  atomw(aw)
      w=0.0d0
      do 100 i=1,natom
      ni = number(i)
      w = w + aw(ni)+dble(numh(i))*1.0079
 100  continue      
C
C
C Average molecular weight
C
      amw=w/dble(natom)
      write(16,101) w
      write(16,102) amw
 101  format('Molecular weight(MW) :       ',e12.4)
 102  format('Average molecular weight(AMW):',e12.4)
      return
      end
C*****************************************************
C     identification ringinformation
C
C*******************************************************
      subroutine ringinf(natom,number,ic,is,iarom
     &,nring,nringatm,nringsize,nringarom)
      implicit double precision(a-h,o-z)
      dimension iarom(natom),number(natom),ic(natom,natom)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension nringarom(natom)
      dimension nc(natom),linc(natom,8)
C     
      do 100 i=1,natom
      nt=0
      do 200 j=1,natom
      if(ic(i,j).eq.1) then
      nt=nt+1
      linc(i,nt)=j
      end if
  200 continue
      nc(i)=nt
  100 continue  
C
C     identification rings and atoms on rings       
C
      call ringpercept(natom,ic,is,nc,linc,nring,nringatm,nringsize)
C
      call printring(natom,number,nring,nringatm,nringsize,nringarom)
C     Judge if an ring is aromatic
C
      call ringarom(natom,number,nring,nringatm,nringsize,nringarom)
C
      do 300 i=1,natom
      iarom(i)=0
300   continue
      if(nring.ge.1) then
        do 400 i=1,nring
         do 500 j=1,nringsize(i)
         ija=nringatm(i,j)
         if(nringarom(i).eq.1) iarom(ija)=1
500      continue
400     continue
      end if
C
      return
      end
C*******************************************
C     print ring-derived descriptors
C******************************************
      subroutine  printring(natom,number,nring
     & ,nringatm,nringsize,nringarom)
      implicit double precision(a-h,o-z)
      character*3 hybrid(natom)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension nringarom(natom),number(natom)
C     number of n-member nonaromatic rings
      n3=0
      n4=0
      n5non=0
      n6non=0
      n7=0
C     number of n-member aromatic rings
      n5ar=0
      n6ar=0
C
      do 100 i=1,nring
      if(nsi.eq.3) n3=n3+1
      if(nsi.eq.4) n4=n4+1
      if(nsi.eq.7) n7=n7+1
      nsi=nringsize(i)
      if( nringarom(i).eq.0) then
        if(nsi.eq.5) n5non=n5non+1
        if(nsi.eq.6) n6non=n6non+1
      else
        if(nsi.eq.5) n5ar=n5ar+1
        if(nsi.eq.6) n6ar=n6ar+1
      end if
100   continue
      write(16,1001) n3
      write(16,1002) n4
      write(16,1003) n7
      write(16,1004) n5non
      write(16,1005) n6non
      write(16,1006) n5ar
      write(16,1007) n6ar
C
C     heterocyclic rings
C     nhr:number of heterocyclic rings
C     nhrn:number of heterocyclic rings containing at least one N atom
C     nhro:number of heterocyclic rings containing at least one O atom
C     nhrs:number of heterocyclic rings containing at least one S atom
C
      nhr=0
      nhrn=0
      nhro=0
      nhrs=0
      do 200 i=1,nring
      do 300 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.eq.7) then
       nhrn=nhrn+1
       goto 200
      end if
300   continue
200   continue
      do 400 i=1,nring
      do 500 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.eq.8) then
       nhro=nhro+1
       goto 400
      end if
500   continue
400   continue
      do 600 i=1,nring
      do 700 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.eq.16) then
       nhrs=nhrs+1
       goto 600
      end if
700   continue
600   continue
      do 800 i=1,nring
      do 900 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.ne.6) then
       nhr=nhr+1
       goto 800
      end if
900   continue
800   continue
      write(16,1101) nhr
      write(16,1102) nhrn
      write(16,1103) nhro
      write(16,1104) nhrs
1001  format('Number of 3-member rings:',i5)
1002  format('Number of 4-member ings:',i5)
1003  format('Number of 7-member rings:',i5)
1004  format('Number of 5-member non-aromatic rings:',i5)
1005  format('Number of 6-member non-aromatic rings:',i5)
1006  format('Number of 5-member aromatic rings:',i5)
1007  format('Number of 6-member aromatic rings:',i5)
1101  format('Number of heterocyclic rings:',i5)
1102  format('Number of N heterocyclic rings:',i5) 
1103  format('Number of O heterocyclic rings:',i5) 
1104  format('Number of S heterocyclic rings:',i5) 
      return
      end
C*************************************************************************
C       Program  for Ring Perception from Molecular Topological Structures
C**************************************************************************
C     Variables Description 
C     cor(j,i),j=1,3):x,y,z cordinates of atom i 
C     number(i): atomic number of atom i
C     nc(i): number of covalently bonded atoms  to atom (Vertex degree)
C     linc(i,j): the jth atom that is  covalently bonded to atom i is linc(i,j)
C-------------------------------------------------------------------------------
      subroutine ringpercept(natom,ic,is,nc,
     &  linc,nring,nringatm,nringsize)
      implicit double precision (a-h,o-z)
      dimension number(natom)
      dimension nc(natom),linc(natom,8),ic(natom,natom)
      dimension is(natom,natom)
      dimension nringatm(natom,natom),nringsize(natom)
C
      call ring(natom,ic,is,nc,linc,nring,nringatm,nringsize)
C
       iprint=0
       if(iprint.eq.0) goto 200
       write(16,*) 'Ring information'
       do 100 i=1,nring
            ns=nringsize(i)
            write(16,*)
            write(16,1001) i,ns
            write(16,1002) (nringatm(i,j),j=1,ns)
 100  continue
 200  continue
 1001 format('The size of ring ',i2,' is',i2)
 1002 format(10i4)
      return
      end
C*************************************************
C     Perception of rings from toplogical structures
C************************************************
      subroutine ring(natom,ic,is,nc,linc,nring,nringatm,nringsize)
      implicit  double precision(a-h,o-z)
C-----------------------------------------
C     natom: number of atoms
C     ic: adjacent matrix:
C     nring: number of ring
C------------------------------------------
      dimension ic(natom,natom),is(natom,natom)
      dimension nc(natom)
      dimension kill(natom),linc(natom,8)
      dimension ncw(natom),killw(natom)
      dimension nringatm(natom,natom),nringsize(natom)
C
       do 100 i=1,natom
       kill(i)=0
       ncw(i)=nc(i)
 100   continue
C
C     Removing rerecursively the atoms with vertex degree=1
C     then,only atoms on rings or between rings are left. 
C     (removed atoms are flagged by kill(i)=1)
C
 200   num=0
       do 300 i=1,natom
       if(kill(i).eq.1) goto 300
       if(ncw(i).eq.1) then
          num=num+1
          kill(i)=1
          ncw(i)=ncw(i)-1 
          do 400 k=1,nc(i)
          ji=linc(i,k)
          if(kill(ji).ne.1) ncw(ji)=ncw(ji)-1
 400      continue
      end if
 300  continue
      if(num.gt.0) goto 200
C
      do 500 i=1,natom
      killw(i)=kill(i)
 500  continue
C
C     then judge if an edge is inside on inside cycle or between two cycles
C     i,e.,ring linkage): if an edge  is inside one cycle,there must be one 
C     of the two cases:
C         (1).ithere is one vertex  that the topological distances 
C         from the two vertexes of the edge are the same;
C         Find the smallest ring containing an edge
C         (2).ithere is one edge  that the topological distances  of the two 
C          vertexes of the edge to  the twoo edges of that edge are the same;
C         (3).record the smallest distance(rmin) and the corresponding vertex or edge
C         then find the smallest ring containing this  edge  and the recorded edge or 
C         the recorded vertex: two vertexes are inside this ring if they satisfy:
C         (1).the distance from the two points to the recorded point or edge are the same (r1);
C         (2).the  distance fromthe teo points to the original edge are the same (r2).
C         (3).r1+r2=rmin
      nfound=0
 600  continue
C     From edge i-----j
      do 700 i=1,natom-1
      if(kill(i).eq.1) goto 700
      do 800 j=i+1,natom
      if(kill(j).eq.1) goto 800
      if(killw(i).eq.1.and.killw(j).eq.1) goto 800
      if(ic(i,j).ne.1) goto 800
      nmin=natom
      nf=2
C     Identification vertex  k: {k:min(rik=rjk)}
      do 900 k=1,natom
      if(k.eq.i.or.k.eq.j) goto 900
      if(kill(k).eq.1) goto 900
      if(is(i,k).eq.0.or.is(j,k).eq.0) goto 900
      if(is(i,k).eq.is(j,k)) then
           if(is(i,k).lt.nmin) then
               nf=0
               nmin=is(i,k)
               keep=k
            end if
            goto 900
      end if
C     identification edge k----m: {k,m:min (rki=rmj) subject to rki<rkj and rmj<rmi}
      do 910 m=1,natom
      if(m.eq.i.or.m.eq.j) goto 910
      if(kill(m).eq.1) goto 910
      if(ic(k,m).ne.1) goto 910
      if(is(i,m).eq.0.or.is(j,m).eq.0) goto 910
      if(is(i,k).eq.is(j,m)) then
           if(is(i,k).lt.is(j,k).and.is(m,j).lt.is(m,i)) then    
              if(is(i,k).lt.nmin) then
                nf=1
                nmin=is(i,k)
                keepi=k
                keepj=m
               end if
           end if
      end if
 910  continue
 900  continue
C
      if(nf.eq.2) goto 800
      nfound=nfound+1
      if(nf.eq.0) nringsize(nfound)=2*nmin+1
      if(nf.eq.1) nringsize(nfound)=2*nmin+2
C
      if(nf.eq.0) then
          nringatm(nfound,1)=i
          nringatm(nfound,2)=j
          nringatm(nfound,3)=keep
      else
          nringatm(nfound,1)=i
          nringatm(nfound,2)=j
          nringatm(nfound,3)=keepi
          nringatm(nfound,4)=keepj
      end if
C
      if(nf.eq.0) then
          killw(i)=1
          killw(j)=1
          killw(keep)=1
      else            
          killw(i)=1
          killw(j)=1
          killw(keepi)=1
          killw(keepj)=1
      end if
C
C (i,j,keep) are on the same ring if nv=1
C (i,j,keepi,keepj) are on the same ring if ne=1
C
C then,find the the other atoms within the same ring: ni----nj
C if (nv.ge.1): r(i,ni)=r(j,nj)=r1 and r(keep,ni)=r(keep,nj),=r2: r1+r2=rmin 
C if (ne.ge.1): r(i,ni)=r(j,nj)=r1 and r(keepi,ni)=r(keepj,nj),=r2: r1+r2=rmin 
C
      if(nf.eq.0) then
        nn=3
      else
        nn=4
      end if
C
      do 930 ni=1,natom
      if(ni.eq.i) goto 930
      if(kill(ni).eq.1) goto 930
      if(is(ni,i).ge.nmin) goto 930
      do 940 nj=1,natom
      if(nj.eq.j) goto 940
      if(kill(nj).eq.1) goto 940
      if(is(nj,j).ge.nmin) goto 940
C
      if (nf.eq.0) then
          if(is(ni,i).eq.is(nj,j).and.is(ni,keep).eq.is(nj,keep)) then
             i1=is(ni,i)
             i2=is(ni,keep)
             ii=i1+i2
             if(ii.eq.nmin) then
              killw(ni)=1
              killw(nj)=1
              nringatm(nfound,nn+1)=ni
              nringatm(nfound,nn+2)=nj
              nn=nn+2
             end if
          end if
      else
          if(is(ni,i).eq.is(nj,j).and.is(ni,keepi).eq.is(nj,keepj)) then
             i1=is(ni,i)
             i2=is(ni,keepi)
             ii=i1+i2
             if(ii.eq.nmin) then
               killw(ni)=1
               killw(nj)=1
               nringatm(nfound,nn+1)=ni
               nringatm(nfound,nn+2)=nj
               nn=nn+2
             end if
          end if
      end if
C
 940  continue
 930  continue
 800  continue
 700  continue
C     write(6,*) 'Number of found rings',nfound
      nring=nfound
      return
      end
C***************************************
C     judge if an ring is aromatic,
C     only for 5 and 6 memebers rings
C***************************************
      subroutine ringarom(natom,number,
     & nring,nringatm,nringsize,nringarom)
      implicit  double precision(a-h,o-z)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension nringarom(natom),number(natom)
      character hyd*3,hybrid(natom)*3
      do 100 i=1,nring
      nringarom(i)=0
 100  continue
C
      do 200 i=1,nring
      nsi=nringsize(i)
      if(nsi.lt.5.or.nsi.gt.6) goto 200
C     for ring of size 5 and  6,check if each atom on this ring is SP2 atom
        ip=1
        do 300 j=1,nsi
        ia=nringatm(i,j)
        hyd=hybrid(ia)
        if(hyd(2:3).eq.'P2') then
           ip2=1
        else
           ip2=0
        end if
        ip=ip*ip2
 300    continue
        if(ip.eq.1) nringarom(i)=1
C
C     for ring of size 5,check if it has 4 SP2 C atom and one atom of lone pair
C     like pyrrole,furan,thiophene
C
      if(nsi.eq.5) then
        icp2=0
        ix=0
        do 400 j=1,nsi
        ia=nringatm(i,j)
        hyd=hybrid(ia)
        nmb=number(ia)
        if(nmb.eq.6.and.hyd(2:3).eq.'P2') icp2=icp2+1
        if(nmb.eq.7.or.nmb.eq.8.or.nmb.eq.16) ix=ix+1
 400    continue
        if(ix.eq.1.and.icp2.eq.4) nringarom(i)=1
      end if
 200  continue
      iprint=0
      if(iprint.eq.0) goto 600
      do  500 i=1,nring
      if(nringarom(i).eq.1) then
        write(6,*) 'Ring ',i,'is aromatici.'
       else
        write(6,*) 'Ring ',i,'is not aromatic.'
       end if
 500  continue
 600  continue
      return
      end 
C********************************************    
C     Fingerprint: according to function groups
C     from h-deleted model
C
C    Assuning that no free radicals exist, excess valence or 
C    less calence for an atom is because of formal charge
C
C********************************************
      subroutine groups(natom,number,nring,nbond
     & ,nringatm,nringsize,nringarom,ic,numh,iarom)
      implicit  double precision(a-h,o-z)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension nringarom(natom),number(natom)
      dimension nbond(natom,natom),ic(natom,natom)
      dimension numh(natom),iarom(natom)
      dimension nc(natom),linc(natom,8)
      dimension nhal(natom),ncc(natom),nod(natom),noh(natom),non(natom)
      dimension nn(natom),ns(natom),no(natom),nor(natom)
      dimension ncoo(natom)
      dimension nt(natom),nbd(natom)
c
C     connect information for  heavy atom 
      do 10  i=1,natom
        nu=0
        do  20 j=1,natom
         if(ic(i,j).eq.1) then
          nu=nu+1
          linc(i,nu)=j
         end if
         nc(i)=nu
20      continue
C     Total number of atoms bonded to atom i
      nt(i)=nc(i)+numh(i)
10    continue
C
      do 30 i=1,natom
C     total number of bonds connected to atom(i)
      nbd(i)=numh(i)
C     Number of halogen atoms bonded to atom i (-X)
       nhal(i)=0
C     Number of C  atoms bonded to atom i (-C)
       ncc(i)=0
C     Number of O atoms bonded to atom i
       no(i)=0
C     Number of O  atoms doube-bonded to atom i,(=O)
       nod(i)=0
C     Number of OH groups bonded to atom i,(-OH)
       noh(i)=0
C     Number of O(-) bonded to atom i,( -O(-) )
       non(i)=0
C     Number of -OR bonded to atom i (R ne H)
      nor(i)=0
C     Number of N  atoms bonded to atom i (-C)
       nn(i)=0
C     Number of S  atoms bonded to atom i (-C)
       ns(i)=0
C
      do 40 j=1,nc(i)
      ija=linc(i,j)
      nb=number(ija)
      nbd(i)=nbd(i)+nbond(i,ija)
      iprd=(nb-9)*(nb-17)*(nb-35)*(nb-53)
      if(iprd.eq.0) nhal(i)=nhal(i)+1
      if(nb.eq.6)ncc(i)=ncc(i)+1
      if(nb.eq.8) then
         no(i)=no(i)+1
         if(nbond(i,ija).eq.2) nod(i)=nod(i)+1
         if(nbond(i,ija).eq.1.and.numh(ija).eq.1) noh(i)=noh(i)+1
         if(nbond(i,ija).eq.1.and.nt(ija).eq.1) non(i)=non(i)+1
         if(nbond(i,ija).eq.1.and.nc(ija).eq.2) nor(i)=nor(i)+1
      end if
      if(nb.eq.7) nn(i)=nn(i)+1
      if(nb.eq.16) ns(i)=ns(i)+1
40    continue
C
30    continue
C
C     number of -COOH  or COO(-)  bonded to atom i
      do 50 i=1,natom
      ncoo(i)=0
      if(number(i).eq.8) goto 50
      do 60 j=1,nc(i)
      ija=linc(i,j)
      nb=number(ija)
      if(nb.eq.6.and.nod(ija).eq.1.and.no(ija).eq.2) ncoo(i)=ncoo(i)+1
60    continue
50    continue
C------Carbon cation [1-1]-------
C(C1)  RCH2(+)(Primary carbocation)
C(C2)  R2CH(+)(secondary carbocation)
C(C3)  R3C(+) (Tertiary carbocation)
      icn1=0
      icn2=0
      icn3=0
      do 100 i=1,natom
      if(number(i).ne.6.or.nbd(i).ne.3) goto 100
      if(nc(i).eq.1.and.numh(i).eq.2) icn1=1
      if(nc(i).eq.2.and.numh(i).eq.1) icn2=1
      if(nc(i).eq.3.and.numh(i).eq.0) icn3=1
100   continue
      write(16,1010) icn1
      write(16,1020) icn2
      write(16,1030) icn3
1010  format('fingerprint for primary carbocation:',i2)
1020  format('fingerprint for secondary carbocation:',i2)
1030  format('fingerprint for tertiary carbocation:',i2)
C
C-----X=F,Cl,Br,I [1-2]------- (organohalide)
C(X1)   R-F
C(X2)   R-Cl
C(X2)   R-Br
       ix=0
      do 200 i=1,natom
       if(number(i).ne.6.or.nhal(i).lt.1) goto 200
       ix=1
200   continue
      write(16,2010) ix       
2010  format('fingerprint for organohalide:',i2)
C
C-----N [1-3]----
C(N1)   RNH3(+),R2(NH2)(+),R3NH(+),R4N(+)    (Ammonium ions)
C
C(N2)   RNH2                                 (Primary aimine)
C(N3)   R2NH                                 (secondary amine)
C(N4)   R3N                                  (Tertiary amine) 
      iami=0
      ipa=0
      isa=0
      ita=0   
C
      do 300 i=1,natom
        if(number(i).ne.7) goto 300
       nch=nc(i)+numh(i)
       if(nch.eq.4) iami=1
       if(nbd(i).ne.3) goto 300
       if(ncc(i).eq.1.and.numh(i).eq.2) ipa=1
       if(ncc(i).eq.2.and.numh(i).eq.1) isa=1
       if(ncc(i).eq.3.and.numh(i).eq.0) ita=1
300   continue
      write(16,3010) iami
      write(16,3020) ipa
      write(16,3030) isa
      write(16,3040) ita
3010  format('fingerprint for amonium ion:',i2)
3020  format('fingerprint for primary amonium:',i2)
3030  format('fingerprint for secondary amonium: ',i2)
3040  format('fingerprint for tertiary amonium: ',i2)
C(N5) R-NO2                                (Nitro)
      ino2=0
      do 400 i=1,natom
      if(number(i).ne.7) goto 400
      if(ncc(i).eq.1.and.no(i).eq.2) then
         ino2=1
         goto 420
      end if
400   continue    
420   write(16,3050)ino2
3050  format('fingerprint for nitro:',i2)
C(N6) R-CN  (Nitrile)
      icn=0
      do 500 i=1,natom
      if(number(i).ne.6) goto 500
       if(nc(i).ne.2) goto 520
       if(ncc(i).ne.1) goto 520
       if(nn(i).ne.1) goto 520
       icn=1
500   continue
520   write(16,3060) icn
3060  format('fingerprint for nitrile:',i2)
C(N7) R-N=N-R                              (Diazo)
      inn=0
      do 600 i=1,natom-1
      nbi=number(i)
      if(nbi.ne.7.or.nc(i).ne.2) goto 600
      do 610 j=i+1,natom
      nbj=number(j)
      if(nbj.ne.7.or.nc(j).ne.2) goto 610
      if(nbond(i,j).ne.2) goto 610
      inn=1
      goto 620
610   continue
600   continue
620   write(16,3070) inn
3070  format('fingerprint for diazo:',i2)
C
C-------O [1.4]--------
C(O1)     RCH2-OH                  (Primary alcohol)
C(O2)     R2CH-OH                  (Secondary alcohol)
C(O3)     R3C-OH                  (Tertiary alcohol)
C(O4)     Ph-OH ( Phenol)
      ipa=0
      isa=0
      ita=0
      ipo=0
C
      do 700 i=1,natom
      ncci=ncc(i)
      if(number(i).ne.6.or.noh(i).eq.0) goto 700
      if(iarom(i).eq.1) then
         ipo=1
         goto 700
      end if
      if(ncci.eq.1.and.numh(i).eq.2) then
         ipa=1
         goto 700
      end if
      if(ncci.eq.2.and.numh(i).eq.1) then
             isa=1
             goto 700
      end if   
      if(ncci.eq.3.and.numh(i).eq.0) then
            ita=1
            goto 700
      end if
700   continue
      write(16,4010) ipo
      write(16,4020) ipa
      write(16,4030) isa
      write(16,4040) ita
4010  format('fingerprint for phenol(Ph-OH):',i2)
4020  format('fingerprint for primary alcohol:',i2)
4030  format('fingerprint for second alcohol:',i2)
4040  format('fingerprint for tertiary alcohol:',i2)
C
C(o5)     R-O-R                    (ether)
C(o6)     Ph-O-Ph
      iror=0
      ipop=0
      do 800 i=1,natom
       if(number(i).ne.8.or.nc(i).ne.2) goto 800
       j1=linc(i,1)
       j2=linc(i,2)
       if(number(j1).ne.6.or.number(j2).ne.6) goto 800
       if(iarom(j1).eq.1.and.iarom(j2).eq.1) then
        ipop=1
        goto 800
       end if
       nch1=nc(j1)+numh(j1)
       nch2=nc(j2)+numh(j2)
       if(nch1.eq.4.and.nch2.eq.4) then
         iror=1
       end if
 800  continue
      write(16,4050) ipop
      write(16,4060) iror
4050  format('fingerprint for Ph-O-Ph:',i2)
4060  format('fingerprint for ether(R-O-R):',i2)
C(O7)     R-CHO                    (Aldehyde)
C(O8)     R-CO-R                   (Ketone)
C(09)     RCOOH                    (Carboxylic acid)
C(O10)    RCOO(-)                  (Carboxylate ion or salt)
C(O11)    RCO(+)                   (Acyl cation)
       ircho=0
       ircor=0
       ircooh=0
       ircoo=0
       irco=0
       Do 900 i=1,natom
        if(number(i).ne.8.or.nc(i).ne.1) goto 900
        jo=linc(i,1)
        nccjo=ncc(jo)
        if(number(jo).ne.6.or.nbond(i,jo).ne.2) goto 900
        kj=linc(jo,1)
        if(numh(jo).eq.1.and.ncc(jo).eq.1) ircho=1
        if(nccjo.eq.2) ircor=1
        if(nccjo.eq.1.and.noh(jo).eq.1) ircooh=1
        if(nccjo.eq.1.and.non(jo).eq.1) ircoo=1
        if(nccjo.eq.1.and.nt(jo).eq.2) irco=1
900    continue
       write(16,4070)ircho
       write(16,4080)ircor
       write(16,4090)ircooh
       write(16,4100)ircoo
       write(16,4110)irco
4070   format('fingerprint for aldehyde(R-CHO):',i2)
4080   format('fingerprint for ketone(R-CO-R):',i2)
4090   format('fingerprint for carboxylic acid(R-COOH):',i2)
4100   format('fingerprint for carboxylate ion(R-COO(-)):',i2)
4110   format('fingerprint for acyl cation(R-CO(+)):',i2)
C(O12)    RCOOR                    (Ester)
       ircoor=0
       do 1000 i=1,natom
       ncci=ncc(i)
       if(number(i).ne.6) goto 1000
       if(ncci.eq.1.and.nod(i).eq.1.and.nor(i).eq.1) then
       ircoor=1
       goto 1100
       end if
1000   continue
1100   continue
       write(16,4120) ircoor
4120   format('fingerprint for ester(R-COOR):',i2)
C(O13)    RCO-O-COR                (Acid Anhydride)
       ircoocor=0
       do 1200 i=1,natom
       if(number(i).ne.8.or.ncc(i).ne.2) goto 1200
       j1=linc(i,1)
       j2=linc(i,2)
       if(no(j1).ne.2.or.no(j2).ne.2) goto 1200
       if(nod(j1).ne.1.or.nod(j2).ne.1) goto 1200
       ircoocor=1
1200   continue
       write(16,4130) ircoocor
4130   format('fingerprint for Acid anhydride(R-CO-O-COR):',i2)
C(O14)    RO(-)                    (Alkoxide ion)
       iro=0
       do 1300 i=1,natom
       if(number(i).ne.8.or.nt(i).ne.1) goto 1300
       ja=linc(i,1)
       nch=ncc(ja)+numh(ja)
       if(nch.eq.3.and.nbond(i,ja).eq.1) then
        iro=1
        goto 1310
       end if
1300   continue
1310   continue
       write(16,4140) iro
4140   format('fingerprint for Alkoxide ion(R-O(-)):',i2)
C(O15)    R-O-OH,R-O-O-R                  (Peroxide)
       irooh=0
       do 1400 i=1,natom
       if(number(i).ne.8.or.nc(i).ne.2) goto 1400
       if(ncc(i).eq.1.and.no(i).eq.1) then
       irooh=1
       end if
1400   continue
       write(16,4150) irooh
4150   format('fingerprint for peroxide(R-O-O-R):',i2)
C(O16)    -(C-O-C)ring-            (Epoxide)
       icocr=0
       do 1500 i=1,nring
       if(nringsize(i).ne.3) goto 1500
       do 1510 j=1,3
       ija=nringatm(i,j)
       if(number(ija).eq.8.and.nc(ija).eq.2) then
       icocr=1
       goto 1520
       end if
1510   continue
1500   continue
1520   continue
       write(16,4160) icocr
4160   format('fingerprint for epoxide(c-O-c ring):',i2)
C(017)    -C(OH)-C(OH)-            (Diol)
       idiol=0
       do 1600 i=1,natom
       if(number(i).ne.6.or.noh(i).ne.1.or.nt(i).ne.4) goto 1600
       do 1610 j=1,nc(i)
       ja=linc(i,j)
       if(number(ja).eq.6.and.noh(ja).eq.1.and.nt(ja).eq.4) then
       idiol=1
       goto 1620
       end if
1610   continue
1600   continue
1620   continue
       write(16,4170) idiol
4170   format('fingerprint for diol(C(OH)-C(OH)-):',i2)

C
C------Si [5]-----
C      -RSiH3,-R2SiH2,R3SiH,R4Si       (Organosilicon)
       isi=0
       do 1700 i=1,natom
       if(number(i).ne.14)goto 1700
       nch=ncc(i)+numh(i)
       if(nch.eq.4) then
       isi=1
       goto 1710
       end if
1700   continue
1710   continue
       write(16,5001) isi
5001   format('fingerprint for organosilicon:',i2)
C
C ------As [6]-----
C    (Organoarsenical)
       ias=0
       do 1800 i=1,natom
       if(number(i).eq.33) then
       ias=1
       end if
1800   continue
       write(16,6001) ias
6001   format('fingerprint for organoarsenical:',i2)
C
C------S   [7]------
C(S1) R-SH                 (Thiol)
C(S2) Ph-SH                (Thiophenol)
C(s3) R-S-R
C(S4) Ph-S-Ph
      irsh=0
      iphsh=0
      irsr=0
      iphsph=0
      do 1900 i=1,natom
      if(number(i).ne.16.or.nt(i).ne.2) goto 1900
      if(nc(i).eq.1.and.numh(i).eq.1) then
       ija=linc(i,1)
       if(iarom(ija).eq.1) then
         iphsh=1
       else
         irsh=1
       end if
       goto 1900
      end if
      if(nc(i).ne.2) goto 1900
      j1=linc(i,1)
      j2=linc(i,2)
      if(iarom(j1).eq.1.and.iarom(j2).eq.1) then
        iphsph=1
      else
        irsr=1
      end if
1900  continue
       write(16,7001)irsh 
       write(16,7002)iphsh 
       write(16,7003)irsr 
       write(16,7004)iphsph
7001   format('fingerprint for thiol(R-SH):',i2)
7002   format('fingerprint for thiophenol(Ph-SH):',i2)
7003   format('fingerprint for R-S-R:',i2)
7004   format('fingerprint for Ph-S-Ph:',i2)
C
C(S4)  R-SO3H or R-SO3(-)   Sulfonic Acid)
C(S5)  -C=S                 (Thioketone)
       isa=0
       ics=0
       do 2100 i=1,natom
       if(number(i).ne.16) goto 2100
       if(ncc(i).eq.1.and.no(i).eq.3) isa=1
       if(nc(i).eq.1.and.nt(i).eq.1) then
         j1=linc(i,1)
         if(nbond(i,j1).eq.2) ics=1
       end if
2100   continue
2110   continue
       write(16,7005)isa
       write(16,7006)ics
7005   format('fingerprint for sulfonic acid:',i2)
7006   format('fingerprint for thioketone(R-C=S):',i2)
C
C------P and  O-----
C         C-P(=O)(-OH)(-OH)                                      (Phosphonic acid)
C         C-P(=O)(-OH)-C                                    (Phosphinic Acid)
      ipho=0
      iphi=0
      do 2200 i=1,natom
      if(number(i).ne.15) goto 2200
      if(ncc(i).eq.1.and.no(i).eq.3) then
        ipho=1
      end if
      if(ncc(i).eq.2.and.no(i).eq.2) then
        iphi=1
      end if
2200  continue
       write(16,7007)ipho
       write(16,7008)iphi
7007   format('fingerprint for phosphonic acid:',i2)
7008   format('fingerprint for phosphinic acid:',i2)
C
C         C-O-P(=O)(-OH)-OH,C-O-P(=O)(-OH)-O (-)    (Organophosphate ester)
       iope=0
       do 2300 i=1,natom
        if(number(i).ne.8.or.nt(i).ne.2.or.nc(i).ne.2) goto 2300
        if(ncc(i).ne.1) goto 2300
        do 2310 j=1,nc(i)
        jo=linc(i,j)
        if(number(jo).eq.15.and.no(jo).eq.4)then
          iope=1
          goto 2320
        end if
2310    continue
2300   continue
2320   continue
       write(16,7009)iope
7009   format('fingerprint for organophophosphate ester:',i2)
C
C-----O and S  [1.8]-----
C      -C-C(=O)-S-C                   (Carboxylic thioester)
       ict=0
       do 2400 i=1,natom
       if(number(i).ne.6.or.nt(i).ne.3) goto 2400
       if(nod(i).ne.1) goto 2400
       do 2410 j=1,nc(i)
       ja=linc(i,j)
       nccja=ncc(ja)
       if(number(ja).eq.16.and.nccja.eq.2.and.nt(ja).eq.2) then
         ict=1 
         goto 2420
       end if
2410   continue
2400   continue
2420   continue
       write(16,8001)ict
8001   format('fingerprint for carboxylic thioester:',i2)
C       -C-O-S(=O)(=O)-O (-)           (Sulfate ester)
        ise=0
      do 2500 i=1,natom
      ncc(i)=ncc(i)
      if(number(i).ne.8.or.nc(i).ne.2) goto 2500
      if(ncci.ne.1) goto 2500
      do 2510 j=1,nc(i)
      jo=linc(i,j)
      if(number(jo).eq.16.and.no(jo).eq.4) then
        ise=1
        goto 2520
      end if
2510  continue
2500  continue 
2520  continue
       write(16,8002)ise
8002   format('fingerprint for sulfate ester:',i2)
C
C  -----P and O and S [1.9] ----
C       C-O-P(=S)(-OH)-O(-)    (Thiophosphate ester)
        ite=0
        do 2600 i=1,natom
        if(number(i).ne.15.or.nc(i).ne.4) goto 2600
        if(no(i).ne.3) goto 2600
        do 2610 j=1,nc(i)
        jp=linc(i,j)
        if(number(jp).eq.16.and.nbond(i,jp).eq.2) then
             ite=1
             goto 2620
        end if
2610    continue
2600    continue
2620    continue
       write(16,9001)ite
9001   format('fingerprint for thiophosphate ester:',i2)
     
C
C-------O and N [1.10]---
C(NO1) RCO(NH2),RC(=O)-NH-R (Amide)
       iam=0
       do 2700 i=1,natom
       if(number(i).ne.6.or.nt(i).ne.3) goto 2700
       if(nod(i).eq.1.and.nn(i).eq.1.and.ncc(i).eq.1) then
         iam=1
       end if
2700   continue
       write(16,1101)iam
1101   format('fingerprint for amide:',i2)
C(NO2) R-CH(NH2)(COOH) (alpha-Amino acida)
       iaac=0
       do 2800 i=1,natom
       if(number(i).ne.6.or.ncoo(i).ne.1.) goto 2800
       do 2810 j=1,nc(i)
       jia=linc(i,j)
       if(number(jia).eq.7.and.numh(jia).ge.2) then
        iaac=1
        goto 2820
       end if
2810   continue
2800   continue
2820   continue
       write(16,1102)iaac
1102   format('fingerprint for alpha-amino acid:',i2)
C(NO3) R2C(CN)(OH)     (Hydroxynitrile)
       ihn=0
       do 2900 i=1,natom
       if(number(i).ne.6.or.nt(i).ne.4) goto 2900
       if(noh(i).ne.1) goto 2900
       do 2910 j=1,nc(i)
       ij=linc(i,j)
       if(number(ij).eq.6.and.nn(ij).eq.1.and.nt(ij).eq.2) then
        ihn=1
        goto 2920
       end if
2910   continue
2900   continue
2920   continue
       write(16,1103)ihn
1103   format('fingerprint for hydroxynitrile:',i2)
C
C(NO4) C=N-OH  (Oxime)
       iox=0
       do 3100 i=1,natom
       if(number(i).ne.7.or.nt(i).ne.2) goto 3100
       if(noh(i).ne.1) goto 3100
       do 3110 j=1,nc(i)
       ij=linc(i,j)
       if(number(ij).eq.6.and.nbond(i,ij).eq.2) then
          iox=1
          goto 3120
       end if
3110   continue
3100   continue
3120   continue
       write(16,1104)iox
1104   format('fingerprint for oxime:',i2)
C(NO5) R-C-O-NO2 (Nitrate ester)
       ine=0
       do 3200 i=1,natom
       if(number(i).ne.8.or.nt(i).ne.2) goto 3200
       if(ncc(i).ne.1.or.nn(i).ne.1) goto 3200
       do 3210 j=1,nc(i)
       ji=linc(i,j)
       if(number(ji).eq.7.and.no(ji).eq.3) then
        ine=1
        goto 3220
       end if
3210   continue
3200   continue
3220   continue
       write(16,1105)ine
1105   format('fingerprint for nitrate ester:',i2)
C
C------O and  CL [1.11]-------
C(OCl1) RCOCL,RCOBR,RCOI  ( Acid halide)
      iah=0
      do 3300 i=1,natom
      if(number(i).ne.6) goto 3300
      if(nod(i).eq.1.and.nhal(i).eq.1) then
         iah=1
         goto 3310
      end if
3300  continue
3310  continue
       write(16,1111)iah
1111   format('fingerprint for acid halide(RCOX):',i2)
       write(16,*)
      return
      end
C***************************************************************
C     Fingerprint derived from ring
C
C     is(i,j): topological distactnce between atom i and j
C     isr(i,j): the shortest topological distance between ring i and j
C     ncr(i,j): number of common atoms between two rings
C     nclus(i): number of rings in cluster i (cluster: ring cluster according to fused rings)
C     nclusr(i,j) the jth ring in cluster i is ring nclusr(i,j)
C     indf(i,j): indicator for two ring i and j are fused diectly or indirectly (in one fused cluster)
C
C*************************************************************
      subroutine fingerring(natom,number,nbond,nring
     & ,nringatm,nringsize,nringarom,ic,numh,iarom)
      implicit  double precision(a-h,o-z)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension nringarom(natom),number(natom)
      dimension nbond(natom,natom),ic(natom,natom)
      dimension numh(natom),iarom(natom)
      dimension nc(natom),linc(natom,8)
      dimension is(natom,natom),isr(nring,nring),ncr(nring,nring)
      dimension kill(nring),nmaxfri(nring)
      dimension nclus(nring),indf(nring,nring)
      dimension nclusr(nring,nring)
C
      n3=0
      n4=0
      n5non=0
      n6non=0
      n7=0
C     number of n-member aromatic rings
      n5ar=0
      n6ar=0
C
      do 100 i=1,nring
      if(nsi.eq.3) n3=n3+1
      if(nsi.eq.4) n4=n4+1
      if(nsi.eq.7) n7=n7+1
      nsi=nringsize(i)
      if( nringarom(i).eq.0) then
        if(nsi.eq.5) n5non=n5non+1
        if(nsi.eq.6) n6non=n6non+1
      else
        if(nsi.eq.5) n5ar=n5ar+1
        if(nsi.eq.6) n6ar=n6ar+1
      end if
100   continue
C
      in3=0
      in4=0
      in7=0
      in5non=0
      in6non=0
      in5ar=0
      in6ar=0
C
      if(n3.gt.0)in3=1
      if(n4.gt.0)in4=1
      if(n7.gt.0)in5=1
      if(n5non.gt.0)in5non=1
      if(n6non.gt.0)in6non=1
      if(n5ar.gt.0)in6non=1
      if(n6ar.gt.0) in6ar=1
C
      write(6,1001) in3
      write(6,1002) in4
      write(6,1003) in7
      write(6,1004) in5non
      write(6,1005) in6non
      write(6,1006) in5ar
      write(6,1007) in6ar
C
C     heterocyclic rings
C     nhr:number of heterocyclic rings
C     nhrn:number of heterocyclic rings containing at least one N atom
C     nhro:number of heterocyclic rings containing at least one O atom
C     nhrs:number of heterocyclic rings containing at least one S atom
C
      nhr=0
      nhrn=0
      nhro=0
      nhrs=0
      do 200 i=1,nring
      do 300 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
C      write(6,*) 'i,j,ija,nb=',i,j,ija,nb
      if(nb.eq.7) then
       nhrn=nhrn+1
       goto 200
      end if
300   continue
200   continue
      do 400 i=1,nring
      do 500 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.eq.8) then
       nhro=nhro+1
       goto 400
      end if
500   continue
400   continue
      do 600 i=1,nring
      do 700 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.eq.16) then
       nhrs=nhrs+1
       goto 600
      end if
700   continue
600   continue
      do 800 i=1,nring
      do 900 j=1,nringsize(i)
      ija=nringatm(i,j)
      nb=number(ija)
      if(nb.ne.6) then
       nhr=nhr+1
       goto 800
      end if
900   continue
800   continue
      inhr=0
      inhrn=0
      inhro=0
      inhrs=0
      if(nhr.gt.0) inhr=1
      if(nhrn.gt.0) inhrn=1
      if(nhro.gt.0) inhro=1
      if(nhrs.gt.0) inhrs=1
      write(6,1101) nhr
      write(6,1102) nhrn
      write(6,1103) nhro
      write(6,1104) nhrs
1001  format('[1.12-1] Fingerprint for 3-member rings:',i2)
1002  format('[1.12-2] fingerprint for 4-member ings:',i2)
1003  format('[1.12-3] Fingerprint for 7-member rings:',i2)
1004  format('[1.12-4] Fingerprint for 5-member non-aromatic rings:',i2)
1005  format('[1.12-5] Fingerprint for 6-member non-aromatic rings:',i2)
1006  format('[1.12-6] Fingerprint for 5-member aromatic rings:',i2)
1007  format('[1.12-7] Fingerprint for 6-member aromatic rings:',i2)
1101  format('[1.12-8] Fingerprint for heterocyclic rings:',i2)
1102  format('[1-12-9] Fingerprint for N heterocyclic rings:',i2) 
1103  format('[1-12-10] Fingerprint for O heterocyclic rings:',i2) 
1104  format('[1-12-11] Fingerprint for S heterocyclic rings:',i2) 
C
C     fingerprint involving ring connection
C
      call topdis(natom,ic,is)
C
C     shortest topological distance between rings
C      
C     Initialization
C
      if(nring.le.1) goto 90 
      do i=1,nring
      do j=1,nring
         ncr(i,j)=0
      if(i.eq.j) then
         isr(i,j)=0 
      else
         isr(i,j)=natom
      end if
      end do
      end do
C
      do 1000 i=1,nring-1
      do 1100 j=i+1, nring
C
      ncom=0
      do 1010 ik=1,nringsize(i)
      ia=nringatm(i,ik)
      do 1110 jk=1,nringsize(j)
      ja=nringatm(j,jk)
      if(is(ja,ia).lt.isr(i,j)) then
        isr(i,j)=is(ja,ia)
        isr(j,i)=is(ja,ia)
      end if       
      if(ia.eq.ja) ncom=ncom+1
1110  continue
1010  continue
       ncr(i,j)=ncom
       ncr(j,i)=ncom
C
1100  continue
1000  continue
C
C      write(6,*) 'Shortest topological distance and #common atoms'
C      do i=1,nring
C      do j=1,nring
C      write(6,*) 'i,j=',i,j,'  isr=',isr(i,j),' ncr=',ncr(i,j)
C      end do
C      end do
C
C    Cluster of rings according to if rings are fused or not
C
      do 1210 k=1,nring
      kill(k)=0
1210  continue
C
      ncl=0
      do 1200 i=1,nring
      if(kill(i).eq.1) goto 1200
      ncl=ncl+1 
      kill(i)=1
      nfdt=1
      nclusr(ncl,1)=i
1220  newfd=0
      do  1400 k=1,nfdt
      nfk=nclusr(ncl,k)
      do 1300 j=1,nring
      if(kill(j).eq.1) goto 1300
      if(ncr(j,nfk).ge.2) then
        newfd=newfd+1
        nfdt=nfdt+1
        kill(j)=1
        nclusr(ncl,nfdt)=j
      end if
1300  continue
1400  continue
      if(newfd.ge.1) goto 1220
       nclus(ncl)=nfdt
1200  continue
C
C     indicator for two ring i and j are fused diectly or indirectly (in one fused cluster)
      do 1410 i=1,nring
      do 1420 j=1,nring
      indf(i,j)=0
1420  continue
1410  continue
C
       do 1430 i=1,ncl
       do 1440 j1=1,nclus(i)
       do 1450 j2=1,nclus(i)
       ja=nclusr(i,j1)
       jb=nclusr(i,j2)
       indf(ja,jb)=1
1450   continue
1440   continue
1430   continue
C      write(6,*) 'Number of cluster of fused rings',ncl
C      write(6,*) 'Cluster  -->    rings'
C      do i=1,ncl
C      write(6,*) i,'  -->  ',(nclusr(i,j),j=1,nclus(i))
C      end do
C
90    continue
C
C  fingerprint for ring linkage information
C
C     number of fused rings:6,5,4,3,2
      in6fr=0
      in5fr=0
      in4fr=0
      in3fr=0
      in2fr=0
      if(nring.le.1) goto 1590
      do  1500 i=1,ncl
      if(nclus(i).eq.6) in6fr=1
      if(nclus(i).eq.5) in5fr=1
      if(nclus(i).eq.4) in4fr=1
      if(nclus(i).eq.3) in3fr=1
      if(nclus(i).eq.2) in2fr=1
1500  continue
1590  continue
      write(6,1105) in6fr
      write(6,1106) in5fr
      write(6,1107) in4fr
      write(6,1108) in3fr
      write(6,1109) in2fr
1105  format('[1-12-12] Fingerprint for fused rings with 6 rings:',i2) 
1106  format('[1-12-13] Fingerprint for fused rings with 5 rings:',i2) 
1107  format('[1-12-14] Fingerprint for fused rings with 4 rings:',i2) 
1108  format('[1-12-15] Fingerprint for fused rings with 3 rings:',i2) 
1109  format('[1-12-16] Fingerprint for fused rings with 2 rings:',i2) 
C
C     indicator for two rings with only one common atom
      icr1=0
C     indicator for molecule containing rings seperated by n edges:n=1,2,3,4,5,6,7
      irs1=0
      irs2=0
      irs3=0
      irs4=0
      irs5=0
      irs6=0
      irs7=0
      if(nring.le.1) goto 1690
      do 1600 i=1,nring-1
      do 1610 j=i+1,nring
      if(ncr(i,j).eq.1.and.indf(i,j).eq.0) icr1=1
      if(isr(i,j).eq.1.and.indf(i,j).eq.0) irs1=1 
      if(isr(i,j).eq.2.and.indf(i,j).eq.0) irs2=1 
      if(isr(i,j).eq.3.and.indf(i,j).eq.0) irs3=1 
      if(isr(i,j).eq.4.and.indf(i,j).eq.0) irs4=1 
      if(isr(i,j).eq.5.and.indf(i,j).eq.0) irs5=1 
      if(isr(i,j).eq.6.and.indf(i,j).eq.0) irs6=1 
      if(isr(i,j).eq.7.and.indf(i,j).eq.0) irs7=1 
1610  continue
1600  continue
1690  continue
      write(6,1201) icr1
      write(6,1202) irs1
      write(6,1203) irs2
      write(6,1204) irs3
      write(6,1205) irs4
      write(6,1206) irs5
      write(6,1207) irs6
      write(6,1208) irs7
1201  format('[1-12-17] Fingerprint for two rings',
     & ' with only one common atom:',i2) 
1202  format('[1-12-18] Fingerprint for containing ',
     &'rings connected  by one non-ring edge:',i2)  
1203  format('[1-12-19] Fingerprint for containing ',
     &'rings connected  by 2 non-ring edges:',i2)  
1204  format('[1-12-20] Fingerprint for containing ',
     &'rings connected  by 3 non-ring edges:',i2)  
1205  format('[1-12-21] Fingerprint for containing ',
     &'rings connected  by 4 non-ring edges:',i2)  
1206  format('[1-12-22] Fingerprint for containing ',
     &'rings connected  by 5 non-ring edges:',i2)  
1207  format('[1-12-23] Fingerprint for containing ',
     &'rings connected  by 6 non-ring edges:',i2)  
1208  format('[1-12-24] Fingerprint for containing ',
     &'rings connected  by 7 non-ring edges:',i2)  
      return
      end
C****************************************************************************************
C     Program to dentify rotable bonds for conformation analysis  from unit-atom model
C--------------------------------------------------------------------
C     Definition of these rotable bonds:
C     (1). Single bonds
C     (2). A bond not on a ring
C     (3). not a terminal bond
C
C     (maybe) only heavy atom information input in cor(3,natm)
C-------------------------------------------------------------------
C     natm: only used to define arrays
C     natom: number of non-hydrogen atom
C     nx(i): vertex degree,i.e.,number of non-hydrogen atoms bonded to atom i
C     lx(i,j): the jth non-hydrogen atom bonded to atom i is lx(i,j)
C     cor(3,i): x,y,z coordinates for atom i
C     ic(i,j): adjacent matrix
C     is(i,j): topological distance between atoms i and j
C     irot(i): if an atom is or not on rotable bonds (=1:yes, =0: no)
C              and a bond is rotable only the two atoms connected the bond 
C              are rotable atoms.
C***************************************************************************
      subroutine numbrot(natom,number,nring,nbond
     & ,nringatm,nringsize,nringarom,ic,is,numh,iarom)
      implicit double precision(a-h,o-z)
      dimension nbond(natom,natom)
      dimension number(natom),nx(natom),lx(natom,8)
      dimension nc(natom),linc(natom,8),ic(natom,natom)
      dimension is(natom,natom),iringatm(natom),irot(natom)
      dimension nringatm(natom,natom),nringsize(natom)
C
      do 100 i=1,natom
      nt=0
      do 200 j=1,natom
      if(ic(i,j).eq.1) then
      nt=nt+1
      linc(i,nt)=j
      end if
  200 continue
      nc(i)=nt
  100 continue  
C
C
C     number of edges in this molecular graph
C
      ne=0
      do 500 i=1,natom-1
      do 600 j=i+1,natom
      if(ic(i,j).eq.1) ne=ne+1
600   continue
500   continue
C     
C     labels for atoms on rings or not
C     and atoms are on rotable bonds or not
C
      do 700 i=1,natom
      iringatm(i)=0
      irot(i)=1
700   continue
      if(nring.eq.0) goto 850
      do 750 i=1,nring
      do 760 j=1,nringsize(i)
      ia=nringatm(i,j)
      iringatm(ia)=1
760   continue
750   continue
C
      do 800 i=1,natom
      if(iringatm(i).eq.1) irot(i)=0
800   continue
850   continue
C
C     number of bonded heavy atoms
C
      do 900 i=1,natom
      nx(i)=0
      do 1000 j=1,natom
      if(number(j).eq.1) goto 1000
      if(ic(i,j).eq.1) then
        nx(i)=nx(i)+1
        lx(i,nx(i))=j
      end if
1000  continue
900   continue
C
C     An atom that has only one bonded non-hydrogen atom is considered
C     as not on rotable bond, because it is meaningless for conformation
C     analysis.
C     For nx=2 or 3,If the bond order is greater than 1,
C     then it is  considered as sp or sp2 hybrid and irot=0.
C     For nx=4 or more, it is considered as rotable (as C in X-CX3
C     or ligand bond in metal).
C
      do 1100 i=1,natom
      if(nx(i).lt.2) then
        irot(i)=0
      else
        if(nx(i).eq.2) then
          lc1=lx(i,1)
          lc2=lx(i,2)
          nb1=nbond(i,lc1)
          nb2=nbond(i,lc2)
          if(nb1.gt.1.or.nb2.gt.1) irot(i)=0
        end if
        if(nx(i).eq.3) then
          lc1=lx(i,1)
          lc2=lx(i,2)
          lc3=lx(i,3)
          nb1=nbond(i,lc1)
          nb2=nbond(i,lc2)
          nb3=nbond(i,lc3)
          if(nb1.gt.1.or.nb2.gt.1.or.nb3.gt.1) irot(i)=0
        end if
C       if(nx(i).ge.4) then
C          irot(i)=1
C        end if
         end if
1100  continue
C
C     now, judge rotable bonds
C
      nprint=0
C
      nrot=0
      if(nprint.eq.1) write(16,*)'Rotable bonds'
      do 1200 i=1,natom-1
      do 1300 j=i+1,natom
      if(ic(i,j).eq.0) goto 1300
      if(nx(i).eq.1.or.nx(j).eq.1) goto 1300
      if(irot(i).eq.0.and.irot(j).eq.0) goto 1300
      nrot=nrot+1
      if(nprint.eq.1) write(16,*) 'Atom',i,j
1300  continue
1200  continue
      write(16,1001) nrot
C
1400  continue
1001  format('Number of rotable bonds:',i5)
      return
      end
C*********************************************************************
C      Topological BCUT descriptors
C**********************************************************************
C     bd(i,j): bond order between atom i,j
C               =3.0: triple bond
C               =2.0: double bond
C               =1.5: aromatic bonds or resonance bonds
C               =1.0: single bond
C               =0.0: nonbonded
C     bm(i,j):  burden matrix
C     refer to: handbook of molecular descriptors,pp133
C-----------------------------------------------------------------------
      subroutine bcuttop(natom,number,ic,sand,bd,es,rvd,alpha,vdw)
      implicit double precision(a-h,o-z)
      dimension ic(natom,natom),bd(natom,natom)
      dimension nc(natom),bm(natom,natom),number(natom)
      dimension aw(100),eigen(natom)
      dimension w(natom,natom),sand(natom)
      dimension es(natom),rvd(natom)
      dimension alpha(natom),vdw(natom)
C
      write(16,*)
C     
C     Off-diagonal element of  Burden matrix
C     
      do 500 i=1,natom-1
      do 600 j=i+1,natom
      if(ic(i,j).eq.0) then
      bm(i,j)=0.001
      else
      bm(i,j)=0.1*bd(i,j)
      if(nc(i).eq.1.or.nc(j).eq.1) bm(i,j)=bm(i,j)+0.01
      end if
C
      bm(j,i)=bm(i,j)
 600  continue
 500  continue
C
C     BCUT Diagonal matrix elements as atomic mass
C
      call  atomw(aw)
      do  700 i=1,natom
      in=number(i)
      bm(i,i)=aw(in)
 700  continue
      do  i=1,natom
      do  j=1,natom
      w(i,j)=bm(i,j)
      end do
      end do
C
       call jacobi(natom,natom,w,eigen)
C
      do  800 i=1,5
      write(16,1002) i,eigen(i)
 800  continue
      do 900 i=1,5
      write(16,1003)i, eigen(natom-i)
 900  continue
C
C     BCUT by  Sanderson electronegativity 
      do 1100 i=1,natom
      bm(i,i)=sand(i)
 1100 continue
C
      do  i=1,natom
      do  j=1,natom
      w(i,j)=bm(i,j)
      end do
      end do
      call jacobi(natom,natom,w,eigen)
C
      do  1200 i=1,5
      write(16,1004) i,eigen(i)
 1200 continue
      do 1300 i=1,5
      write(16,1005) i, eigen(natom-i)
 1300 continue
C
C     Diagonal matrix elements  as van der waals radius
C
      do 1400  i=1,natom
      bm(i,i)=rvd(i)
 1400 continue
      do  i=1,natom
         do  j=1,natom
          w(i,j)=bm(i,j)
         end do
      end do
      call jacobi(natom,natom,w,eigen)
      do  1600 i=1,5
      write(16,1006) i,eigen(i)
 1600 continue
      do 1700 i=1,5
      write(16,1007)  i, eigen(natom-i)
 1700 continue
C
C     BCUT  Diagonal matrix elements  as  Estate values
C
      do 3100  i=1,natom
      bm(i,i)=es(i)
 3100 continue
      do  i=1,natom
      do  j=1,natom
      w(i,j)=bm(i,j)
      end do
      end do
      call jacobi(natom,natom,w,eigen)
C
      do  1800 i=1,5
      write(16,1008) i,eigen(i)
 1800 continue
      do 1900 i=1,5
      write(16,1009)  i, eigen(natom-i)
 1900 continue
C
C     BCUT Diagonal matrix elements  as  polarizability
C
      do 4100  i=1,natom
      bm(i,i)=alpha(i)
 4100 continue
      do  i=1,natom
         do  j=1,natom
          w(i,j)=bm(i,j)
         end do
      end do
      call jacobi(natom,natom,w,eigen)
C
      do  2100 i=1,5
      write(16,1011) i,eigen(i)
 2100 continue
      do 2200 i=1,5
      write(16,1012)  i, eigen(natom-i)
 2200 continue
C
C      Diagonal matrix elements  as  van dew Waals volume
C
       do 5100  i=1,natom
       bm(i,i)=vdw(i)
 5100  continue
       do  i=1,natom
       do  j=1,natom
       w(i,j)=bm(i,j)
       end do
       end do
       call jacobi(natom,natom,w,eigen)
C
      do  2300 i=1,5
      write(16,1013) i,eigen(i)
 2300 continue
      do 2400 i=1,5
      write(16,1014) i, eigen(natom-i)
 2400 continue
C
C
 1002 format('BCUT ',i2,' highest eigenvalues by atomic mass:',e12.4)
 1003 format('BCUT ',i2,' lowest  eigenvalues by atomic mass:',e12.4)
 1004 format('BCUT ',i2,' highest eigenvalues by EN :',e12.4)
 1005 format('BCUT ',i2,' lowest eigenvalues by  EN:',e12.4)
 1006 format('BCUT ',i2,' highest eigenvalues by VDW radius:',e12.4)
 1007 format('BCUT ',i2,' lowest eigenvalues by  VDW radius:',e12.4)
 1008 format('BCUT ',i2,' highest eigenvalues by Estate:',e12.4)
 1009 format('BCUT ',i2,' lowest eigenvalues by  Estate:',e12.4)
 1011 format('BCUT ',i2,' highest eigenvalues by polarizability:',e12.4)
 1012 format('BCUT ',i2,' lowest eigenvalues by polarizability:',e12.4)
 1013 format('BCUT ',i2,' highest eigenvalues by VDW volume:',e12.4)
 1014 format('BCUT ',i2,' lowest eigenvalues by VDW volume:',e12.4)
C
      return
      end
C*****************************************************************
C     Dagonalization of a real symmetric matrix
C     by the method of Jacobi rotations
C
C     Refer to: chapter 11.1 in <<NUmerical recipes in Fortran>>
C*****************************************************************
C     n    logical dimension of the matrix to be diagonalized
C     np   physical dimension of the matrix storage area
C     a    input with the matrix to be diagonalized; only
C             the upper triangle and diagonal are required
C     d    returned with the eigenvalues in descending order
C     v    returned with the eigenvectors of the matrix
C     b    temporary work vector
C     z    temporary work vector
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jacobi(n,np,a,d)
      implicit double precision(a-h,o-z)
      dimension a(np,np),d(np),v(np,np),b(np),z(np)
c
c
c     setup and initialization
c
C
      maxrot = 19500
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     print warning if not converged
c
   10 continue
      if (nrot .eq. maxrot) then
         write (6,20)
   20    format (' JACOBI  --  Matrix Diagonalization not Converged')
      return
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .gt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end
C**********************************************
C     Moreau-Broto topological descriptors   
C**********************************************
C     refer to:Handbook,p17-p20 
C     Definition:
C     ACd=SUMi SUMj delta(i,j)w(i) w(j)
C     where: delta(i,j)=1 if d(i,j)=d,otherwise=0
C     d(i,j) is the topological distance between atom i and j
C     w(i) is the atom properties for atom i
C---------------------------------------------------------------
      subroutine moreau(natom,is,number,sand,
     & es,rvd,alpha,vdw)
C---------------------------------------------
C     natom:      number of atoms
C     number(i):  atomic number of atom i
C     name(i):    atomic symbol of atom i
C     is(i,j):    topological distances between atom i ang j
C     sand(i):    Sanderson electronegativity scale
C     es(i):      estate value for atom i 
C     rvd(i):     radius of van der Waals for atom i
C     vdw(i):     van der Waals voulme for atom i
C     alpha(i):   polarizability for atom i
C---------------------------------------------------
      implicit  double precision(a-h,o-z)
      character name(natom)*2
      dimension number(natom)
      dimension is(natom,natom),rvd(natom)
      dimension aw(100),es(natom)
      dimension atmass(natom),sand(natom)
      dimension alpha(natom),vdw(natom)
C
C
C     atomic weight
C
      call  atomw(aw)
C
C     Atomic mass for atoms in the molecule
C
      do  100 i=1,natom
      inum=number(i)
      atmass(i)=aw(inum)
 100  continue
C
      ip=1
      call moreauw(natom,is,atmass,ip)
C
C
C
      ip=2
      call moreauw(natom,is,sand,ip)
C
      ip=3
      call moreauw(natom,is,rvd,ip)
C
      ip=4
      call moreauw(natom,is,es,ip)
C
      ip=5
      call moreauw(natom,is,alpha,ip)
C
      ip=6
      call moreauw(natom,is,vdw,ip)
C
      return
      end
C*************************************************************
C     Moran autocorrelation of a topological   structure  
C*************************************************************
C  refer to:Handbook,p17-p20 
C  note: when the lag id is too large, Sum (delta(i,j)) will be 0.0
C        or when all the atom are the same, 0.0 is given to ATS
C---------------------------------------------------------------
      subroutine moran(natom,is,number,sand,
     & es,rvd,alpha,vdw)
C------------------------------------------------------------
C     natom:      number of atoms
C     number(i):  atomic number of atom i
C     name(i):    atomic symbol of atom i
C     is(i,j):    topological distances between atom i ang j
C     sand(i):    Sanderson electronegativity scale
C     es(i):      estate value
C     rvd(i):     radius of ver der waals
C     vdw(i):     van der Waals volume for atom i
C     alpha(i):   polarizability of atom i
C---------------------------------------------
      implicit  double precision(a-h,o-z)
      character name(natom)*2
      dimension number(natom)
      dimension is(natom,natom)
      dimension aw(100),es(natom)
      dimension atmass(natom),sand(natom)
      dimension rvd(natom),vdw(natom),alpha(natom)
C
C     atomic  mass
C
      call  atomw(aw)
C
C     Atomic mass for atoms in the molecule and average atomic mass
C
      do  100 i=1,natom
      inum=number(i)
      atmass(i)=aw(inum)
 100  continue
C
      ip=1
      call moranw(natom,is,atmass,ip)
C
      ip=2
      call moranw(natom,is,sand,ip)
C
C
      ip=3
      call moranw(natom,is,rvd,ip)
C
      ip=4
      call moranw(natom,is,es,ip)
C
      ip=5
      call moranw(natom,is,alpha,ip)
C
      ip=6
      call moranw(natom,is,vdw,ip)
C
      return
      end
C*********************************************************
C     Geary autocorrelation of a topological   structure  
C*********************************************************
C     refer to:Handbook,p17-p20 
C     note: when the lag id is too large, Sum (delta(i,j)) will be 0.0
C     or when all the atom are the same, 0.0 is given to all ATS
C------------------------------------------------------------------
      subroutine  geary(natom,is,number,sand,es,rvd,alpha,vdw)
      implicit  double precision(a-h,o-z)
C------------------------------------------------------------------
C     natom:      number of atoms
C     number(i):  atomic number of atom i
C     name(i):    atomic symbol of atom i
C     is(i,j):    topological distances between atom i ang j
C     sand(i):    Sanderson electronegativity scale
C     es(i):      estate value for atom i
C     rvd(i):     vdw radius of atom i      
C     vdw(i):     van der Waals volume
C     alpha(i):   polarizability
C-------------------------------------------------------------------
      character name(natom)*2
      dimension number(natom)
      dimension is(natom,natom),es(natom)
      dimension aw(100)
      dimension atmass(natom),sand(natom)
      dimension rvd(natom),alpha(natom),vdw(natom)
C
C     atomic mass
C
      call  atomw(aw)
      do  100 i=1,natom
      inum=number(i)
      atmass(i)=aw(inum)
 100  continue
C
      ip=1
      call  gearyw(natom,is,atmass,ip)
C
      ip=2
      call  gearyw(natom,is,sand,ip)
C
      ip=3
      call  gearyw(natom,is,rvd,ip)
C
C
C     Estate  weighted Geary Topological descriptors
C
      ip=4
      call  gearyw(natom,is,es,ip)
C
      ip=5
      call gearyw(natom,is,alpha,ip)
C
      ip=6
      call gearyw(natom,is,vdw,ip)
C
C
      return
      end
C-------------------------------------------------------------
C     Moreau-Broto topological descriptors  weighted by W  
C-------------------------------------------------------
C     refer to:Handbook,p17-p20 
C     Definition:
C     ACd=SUMi SUMj delta(i,j)w(i) w(j)
C     where: delta(i,j)=1 if d(i,j)=d,otherwise=0
C     d(i,j) is the topological distance between atom i and j
C     w(i) is the atom properties for atom i
C---------------------------------------------------------------
      subroutine moreauw(natom,is,w,ip)
C------------------------------------------------------------
C     natom:      number of atoms
C     is(i,j):    topological distances between atom i ang j
C     w(i):       weighting factor
C------------------------------------------------------------
      implicit  double precision(a-h,o-z)
      dimension is(natom,natom)
      dimension wac(10),w(natom)
C
      do 200 id=1,10
      wac(id)=0.0d0
      wac0=0.0d0
      do 300 i=1,natom
      wac0=wac0+w(i)*w(i)
      do 400 j=1,natom
      if (is(i,j).ne.id) goto 400
      wac(id)=wac(id)+w(i)*w(j)
 400  continue
 300  continue
 200  continue
      write(16,*)
      if(ip.eq.1) write(16,1001)  0,wac0
      if(ip.eq.2) write(16,1002)  0,wac0
      if(ip.eq.3) write(16,1003)  0,wac0
      if(ip.eq.4) write(16,1004)  0,wac0
      if(ip.eq.5) write(16,1005)  0,wac0
      if(ip.eq.6) write(16,1006)  0,wac0
C
      do 500 id=1,10
      if(ip.eq.1) write(16,1001) id,wac(id)
      if(ip.eq.2) write(16,1002) id,wac(id)
      if(ip.eq.3) write(16,1003) id,wac(id)
      if(ip.eq.4) write(16,1004) id,wac(id)
      if(ip.eq.5) write(16,1005) id,wac(id)
      if(ip.eq.6) write(16,1006) id,wac(id)
 500  continue
 1001 format('Atomic mass  Moreau-Broto  Autocorrelation lagged ',
     &i5,':    ',e12.4)
 1002 format('Electronegativity  Moreau-Broto Autocorrelation lagged ',
     & i5,':   ',e12.4)
 1003 format('VDW radius  Moreau-Broto Autocorrelation lagged ',
     & i5,':   ',e12.4)
 1004 format('E-State  Moreau-Broto Autocorrelation lagged ',
     & i5,':   ',e12.4)
 1005 format('Polarizability  Moreau-Broto  Autocorrelation lagged ',
     & i5,':   ',e12.4)
 1006 format('VDW volume Moreau-Broto Autocorrelation lagged ',
     & i5,':   ',e12.4)
C
      return
      end
C---------------------------------------------------------
C     Moran autocorrelation  weighted by W 
C----------------------------------------
C  refer to:Handbook,p17-p20 
C  note: when the lag id is too large, Sum (delta(i,j)) will be 0.0
C        or when all the atom are the same, 0.0 is given to ATS
C---------------------------------------------------------------
      subroutine moranw(natom,is,w,ip)
C----------------------------------------
C     natom:      number of atoms
C     is(i,j):    topological distances between atom i ang j
C---------------------------------------------------------------
      implicit  double precision(a-h,o-z)
      dimension is(natom,natom),w(natom)
      dimension wac(10)
C
      avgw=0.0d0
      do  100 i=1,natom
      avgw=avgw+w(i)
 100  continue
      avgw=avgw/dble(natom)
C
      do 200 id=1,10
      a=0.0d0
      b=0.0d0
      is0=0
      do 300 i=1,natom
      b=b+(w(i)-avgw)*(w(i)-avgw)
      do 400 j=1,natom
      if (is(i,j).ne.id) goto 400
      is0=is0+1
      a=a+(w(i)-avgw)*(w(j)-avgw)
 400  continue
 300  continue
      if (is0.eq.0.or.dabs(b).lt.1.0e-8) then
         wac(id)=0.0d0
      else
         wac(id)=(a/dble(is0))/(b/dble(natom))
      end if
 200  continue
      do 500 id=1,10
      if (ip.eq.1) write(16,1001) id,wac(id)
      if (ip.eq.2) write(16,1002) id,wac(id)
      if (ip.eq.3) write(16,1003) id,wac(id)
      if (ip.eq.4) write(16,1004) id,wac(id)
      if (ip.eq.5) write(16,1005) id,wac(id)
      if (ip.eq.6) write(16,1006) id,wac(id)
 500  continue
 1001 format('Atomic mass Moran  Autocorrelation lagged '
     &,i5,':   ',e12.4)
 1002 format('EN  Moran  Autocorrelation lagged '
     &,i5,':   ',e12.4)
 1003 format('VDW radius Moran  Autocorrelation lagged '
     &,i5,':   ',e12.4)
 1004 format('E-State Moran autocorrelation lagged ',
     & i5,':   ',e12.4)
 1005 format('Polarizability Moran Autocorrelation lagged ',
     & i5,':   ',e12.4)
 1006 format('VDW volume Moran Autocorrelation lagged '
     & ,i5,':   ',e12.4)
C
      return
      end 
C-----------------------------------------------------------------
C     Geary autocorrelation   weithted by W 
C--------------------------------
C  refer to:Handbook,p17-p20 
C  note: when the lag id is too large, Sum (delta(i,j)) will be 0.0
C        or when all the atom are the same, 0.0 is given to all ATS
C------------------------------------------------------------------
      subroutine  gearyw(natom,is,w,ip)
C------------------------------------------------------------------
C     natom:      number of atoms
C      w(i):       weighting factor
C     is(i,j):    topological distances between atom i ang j
C-------------------------------------------------------------------
      implicit  double precision(a-h,o-z)
      dimension is(natom,natom),w(natom)
      dimension wac(10)
C
      avgw=0.0d0
      do  100 i=1,natom
      avgw=avgw+w(i)
 100  continue
      avgw=avgw/dble(natom)
C
C
      do 200 id=1,10
      a=0.0d0
      b=0.0d0
      is0=0
      do 300 i=1,natom
      b=b+(w(i)-avgw)*(w(i)-avgw)
      do 400 j=1,natom
      if (is(i,j).ne.id) goto 400
      is0=is0+1
      a=a+(w(i)-w(j))*(w(i)-w(j))
 400  continue
 300  continue
      if (is0.eq.0.or.dabs(b).lt.1.0e-8) then
         wac(id)=0.0d0
      else
         wac(id)=(a/dble(2*is0))/(b/dble(natom-1))
      end if
 200  continue
C
      do 500 id=1,10
      if(ip.eq.1) write(16,1001) id,wac(id)
      if(ip.eq.2) write(16,1002) id,wac(id)
      if(ip.eq.3) write(16,1003) id,wac(id)
      if(ip.eq.4) write(16,1004) id,wac(id)
      if(ip.eq.5) write(16,1005) id,wac(id)
      if(ip.eq.6) write(16,1006) id,wac(id)
 500  continue
C
 1001 format('Atomic mass Geary Autocorrelation lagged',i5,':   ',e12.4)
 1002 format('EN Geary Aotocorrelation lagged',i5,':   ',e12.4)
 1003 format('VDW radius  Geary Autocorrelation lagged',i5,':   ',e12.4)
 1004 format('Estate Geary Autocorrelation lagged',i5,':',e12.4)
 1005 format('Polarizability Geary Autocorrelation lagged'
     &,i5,':   ',e12.4)
 1006 format('VDW volume Geary Autocorrelation lagged',i5,':   ',e12.4)
C
      return
      end 
C********************************************************
C      transfer characters in name to nameout :
C      (1) characters in lower case to upercase
C      (2) remove character ' ' 
C*********************************************************
      subroutine trans(name,idm,nameout,ii)
      implicit double precision (a-h,o-z)
      character name*201, nameout*201
      do 100 i=1,idm
      im=ichar(name(i:i))
C     transfer lower case name to upercase name 
      if(im.ge.97.and.im.le.122) then
        im=im-32
        name(i:i)=char(im)
      end if
100   continue
C
C     Remove blank character 
C
      ii=0
      do 200 i=1,idm
       if(name(i:i).eq.' ') goto 200
         ii=ii+1
         nameout(ii:ii)=name(i:i)
200   continue
      return
      end 
C
C***************************************************
C      Topological polar surface area
C      Definition: Summation over solvent accessible surface area for N, O, S, P atoms 
C      refer to:
C         Peter Ertl, Bernhard Rohde, and Paul Selzer
C         Fast  Calculation of Molecular Polar Surface Area as a Sum 
C         of Fragment-Based Contributions and Its Application to the
C          Prediction of Drug Transport
C          J. Med. chem.,2000,43,3714
C***************************************************************
      subroutine tpsa(natom,number,ic,nbond,numh,iarom 
     & ,nring,nringatm,nringsize)
      implicit double precision(a-h,o-z)
      dimension ic(natom,natom),nbond(natom,natom)
      dimension nc(natom),number(natom),numh(natom)
      dimension linc(natom,8),iarom(natom)
      dimension nringatm(natom,natom),nringsize(natom)
      dimension i3mring(natom),iarbond(natom,natom)
C
C     nc(i): number of non-H atoms bonded to atoms
C     linc(i,j): the jth non-h atom bonded to atom i is atom linc(i,j)
C
      do 100 i=1,natom
        nu=0
        do 200 j=1,natom
         if(ic(i,j).eq.1) then
          nu=nu+1
          linc(i,nu)=j
         end if
         nc(i)=nu
200   continue
100   continue
C
C      write(16,*) 'natom=',natom
C      do 110 i=1,natom
C      write(16,*) i,number(i),'=',(linc(i,k),k=1,nc(i))
C110   continue
C      write(16,*) 'nring=',nring
C      if(nring.ge.1) then
C       do 120 i=1,nring
C        write(16,*) 'Atoms for ring ',i
C        write(16,*) (nringatm(i,j),j=1,nringsize(i))
C120    continue
C      end if
C
C     if atoms are in 3-member ring:i3mring(i)
C
      do 300 i=1,natom
       i3mring(i)=0
300   continue
      do 400 i=1,nring
      ns=nringsize(i)
      if(ns.ne.3) goto 400
      do 500 j=1,ns
      ja=nringatm(i,j)
      i3mring(ja)=1
500   continue
400   continue
C
C     if a bond is an aromatic bond? iarbond(i,j)
C
      do 600 i=1,natom
      do 700 j=1,natom
      iarbond(i,j)=0
      if(iarom(i).ne.1) goto 700
      if(iarom(j).ne.1) goto 700
      if(ic(i,j).ne.1) goto 700
      iarbond(i,j)=1
700   continue
600   continue
      xt=0.0
C
      do 1000 i=1,natom
C
C     ibd: number of aromatic bonds that are bonded to atom i
C     nbd1: number of nonaromatic single bonds that are bonded to atom i
C     nbd2: number of nonaromatic doublle bonds that are bonded to atom i
C     nbd2: number of nonaromatic triple  bonds that are bonded to atom i
C
      ibd=0
      nbd1=0
      nbd2=0
      nbd3=0
      do 900 j=1,nc(i)
      ja=linc(i,j)
      iabd=iarbond(i,ja)
      nbd=nbond(i,ja)
      if(iabd.eq.1) then
        ibd=ibd+1
      else
        if(nbd.eq.1) nbd1=nbd1+1
        if(nbd.eq.2) nbd2=nbd2+1
        if(nbd.eq.3) nbd3=nbd3+1
      end if
900   continue
C
C **************** N *****************************
C
      if (number(i).ne.7) goto 2000
C
C     for N atom that is not on an aromatic ring: type (1)-(18)
C
      if(iarom(i).eq.1) goto 1300
C
C    (1)[N](-*)(-*)-*   
C    (6) N1(-*)(-*)(-*)-1  
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1110
      if(ibd.eq.0.and.nbd1.eq.3.and.nbd2.eq.0.and.nbd3.eq.0) then
         if(i3mring(i).eq.1)then 
              xt=xt+3.01
C             write(16,*) '(6) [N]1-(-*)(-*)(-*)-1 ',i
             goto 1000
         else
             xt=xt+3.24
C             write(16,*) '(1) [N](-*)(-*)-* ',i
            goto 1000
         end if
      end if
C
1110  continue
C
C    (2) [N] (-*)=*
C         
      if(nc(i).ne.2.or.numh(i).ne.0) goto 1130
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+12.36
C         write(16,*) '(1) N(-*)=* ',i
         goto 1000
      end if
C
1130  continue
C
C    (3) [N]#*
C         
      if(nc(i).ne.1.or.numh(i).ne.0) goto 1140
      if(ibd.eq.0.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.1) then
         xt=xt+23.79
C         write(16,*) '(3) [N]#* ',i
         goto 1000
      end if
C
1140  continue
C
C    (4) [N] (-*)(=*)=*
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1150
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.2.and.nbd3.eq.0) then
         xt=xt+11.68
C         write(16,*) '(4) [N](-*)(=*)=* ',i
         goto 1000
      end if
C
1150  continue
C
C    (5) N(=*)#*  (middle atom nitrogen in azide group)
C         
      if(nc(i).ne.2.or.numh(i).ne.0) goto 1160
      if(ibd.eq.0.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.1) then
        xt=xt+13.60
C        write(16,*) '(5) [N](=*)#* ',i
        goto 1000
      end if
C
1160  continue
C    (7) [NH](-*)(-*)  
C    (8) [NH]1(-*)(-*)-1  
C         
      if(nc(i).ne.2.or.numh(i).ne.1) goto 1190
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.0.and.nbd3.eq.0) then
         if(i3mring(i).ne.1) then
           xt=xt+12.03
C           write(16,*) '(7) [NH](-*)(-*) ',i
           goto 1000
         else
           xt=xt+21.94
C           write(16,*) '(8) [NH]1-*-*-1 ',i
          goto 1000
        end if
      end if
1190  continue
C
C    (9) [NH]=*
C         
      if(nc(i).ne.1.or.numh(i).ne.1) goto 1200
      if(ibd.eq.0.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+23.85
C         write(16,*) '(9) [NH]=* ',i
         goto 1000
      end if
1200  continue
C
C    (10) [NH2]-*
C         
      if(nc(i).ne.1.or.numh(i).ne.2) goto 1210
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+26.02
C         write(16,*) '(10) [NH2]-* ',i
         goto 1000
      end if
1210  continue
C
C     (11) [N+](-*)(-*)(-*)-*  
C         
      if(nc(i).ne.4.or.numh(i).ne.0) goto 1220
      if(ibd.eq.0.and.nbd1.eq.4.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+0.00
C         write(16,*) '(11) [N+](-*)(-*)(-*)-* ',i
         goto 1000
      end if
1220  continue
C
C     (12) [N+](-*)(-*)=*  
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1230
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+3.01
C         write(16,*) '(12) [N+](-*)(-*)=* ',i
         goto 1000
      end if
1230  continue
C
C     (13) [N+](-*)#*  
C         
      if(nc(i).ne.2.or.numh(i).ne.0) goto 1240
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.1) then
       xt=xt+4.36
C       write(16,*) '(13) [N+](-*)#* ',i
       goto 1000
      end if
1240  continue
C
C     (14) [NH+](-*)(-*)-*  
C         
      if(nc(i).ne.3.or.numh(i).ne.1) goto 1250
      if(ibd.eq.0.and.nbd1.eq.3.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+4.44
C         write(16,*) '(14) [NH+](-*)(-*)-* ',i
         goto 1000
      end if
1250  continue
C
C     (15) [NH+](-*)=*  
C         
      if(nc(i).ne.2.or.numh(i).ne.1) goto 1260
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+13.97
C         write(16,*) '(15) [NH+](-*)=* ',i
         goto 1000
      end if
1260  continue
C
C     (16) [NH2+](-*)-*  
C         
      if(nc(i).ne.2.or.numh(i).ne.2) goto 1270
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+16.61
C         write(16,*) '(16) [NH2+](-*)-* ',i
         goto 1000
      end if
1270  continue
C
C     (17) [NH2+]=*  
C         
      if(nc(i).ne.1.or.numh(i).ne.2) goto 1280
      if(ibd.eq.0.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+25.59
C         write(16,*) '(17) [NH2+]=* ',i
         goto 1000
      end if
1280  continue
C
C     (18) [NH3+]-*  
C         
      if(nc(i).ne.1.or.numh(i).ne.3) goto 1290
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+27.64
C         write(16,*) '(18) [NH3+]-* ',i
         goto 1000
      end if
1290  continue
      write(16,*) 'Warning, no atom type for non aromatic N',i
      goto 1000
C
C     For N atom that is part of an aromatic ring
C
1300   continue
C
C
C     (19) [n](:*):* 
C         
      if(nc(i).ne.2.or.numh(i).ne.0) goto 1310
      if(ibd.eq.2.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+12.89
C         write(16,*) '(19) [n](:*):* ',i
         goto 1000
      end if
1310  continue
C
C     (20) [n](:*)(:*):* 
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1320
      if(ibd.eq.3.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+4.41
C         write(16,*) '(20) [n](:*)(:*):* ',i
         goto 1000
      end if
1320  continue
C
C     (21) [n](-*)(:*):* 
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1330
      if(ibd.eq.2.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+4.93
C         write(16,*) '(21) [n](-*)(:*):* ',i
         goto 1000
      end if
1330  continue
C
C     (22) [n](=*)(:*):* 
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1340
      if(ibd.eq.2.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+8.39
C         write(16,*) '(22) [n](=*)(:*):* ',i
         goto 1000
      end if
1340  continue
C
C     (23) [nH](:*):* 
C         
      if(nc(i).ne.2.or.numh(i).ne.1) goto 1350
      if(ibd.eq.2.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.0) then
        xt=xt+15.79
C        write(16,*) '(23) [nH](:*):* ',i
        goto 1000
      end if
1350  continue
C
C     (24) [n+](:*)(:*):* 
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1360
      if(ibd.eq.3.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.0) then
        xt=xt+4.10
C        write(16,*) '(24) [n+](:*)(:*):* ',i
        goto 1000
      end if
1360  continue
C
C     (25) [n+](-*)(:*):* 
C         
      if(nc(i).ne.3.or.numh(i).ne.0) goto 1370
      if(ibd.eq.2.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+3.88
C         write(16,*) '(25) [n+](-*)(:*):* ',i
         goto 1000
      end if
1370  continue
C
C     (26) [nH+](:*):* 
C         
      if(nc(i).ne.2.or.numh(i).ne.1) goto 1380
      if(ibd.eq.2.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.0) then
C         write(16,*) '(26) [nH+](:*):* ',i
         goto 1000
      end if
1380  continue
       write(16,*) 'Warning, no atom type for aromatic N atom ',i
       goto 1000
C
2000  continue
C
C    *********   O  atom *******************
C
      if(number(i).ne.8) goto 3000
      if(iarom(i).eq.1) goto 2100
C
C     for non-aromatic O atom
C
C     (27) [O](-*)-*
C     (28) [O}1-*-*-1
C
      if(nc(i).ne.2.or.numh(i).ne.0) goto 2110
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.0.and.nbd3.eq.0) then
         if(i3mring(i).eq.1) then
           xt=xt+12.53
C           write(16,*) '(28)[O]1-*-*-1 ',i
           goto 1000
        else
           xt=xt+9.23
C           write(16,*) '(27)[O](-*)-* ',i
          goto 1000
        end if
      end if
2110  continue
C
C     [29] [O]=*
C
      if(nc(i).ne.1.or.numh(i).ne.0) goto 2120
      if(ibd.eq.0.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.0) then
        xt=xt+17.07
C        write(16,*) '(29)[O]=* ',i
        goto 1000
      end if
2120  continue
C
C     [30] [OH]-*
C
      if(nc(i).ne.1.or.numh(i).ne.1) goto 2130
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
        xt=xt+20.23
C        write(16,*) '(30)[OH]-* ',i
        goto 1000
      end if
2130  continue
C
C     [31] [O-]-*
C
      if(nc(i).ne.1.or.numh(i).ne.0) goto 2140
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
        xt=xt+23.06
C        write(16,*) '(31)[O-]-* ',i
        goto 1000
      end if
2140  continue
C
2100  continue
C     for aromatic O atom
C
C     (32) [o](-*)-*
      xt=xt+13.14
C      write(16,*) ' (32)[o](-*)-* ',i
      goto 1000
C
3000  continue
      if(number(i).ne.16) goto 4000
C   
C***************  S   atom  ***********************
C
C
       if(iarom(i).eq.1) goto 3500
C
C      non-aromatic S atom
C     (33) [S](-*)-*
C
      if(nc(i).ne.2.or.numh(i).ne.0) goto 3110
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.0.and.nbd3.eq.0) then
          xt=xt+25.30
C         write(16,*) '(33) [S](-*)-* ',i
         goto 1000
      end if
3110  continue
C
C     (34) [S]=*
C
      if(nc(i).ne.1.or.numh(i).ne.0) goto 3120
      if(ibd.eq.0.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.0) then
         xt=xt+32.09
C         write(16,*) '(34) [S]=* ',i
         goto 1000
      end if
3120  continue
C
C     (35) [S](-*)(-*)=*
C
      if(nc(i).ne.3.or.numh(i).ne.0) goto 3130
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.1.and.nbd3.eq.0) then
        xt=xt+19.21
C        write(16,*) '(35) [S](-*)(-*)=* ',i
        goto 1000
      end if
3130  continue
C
C     (36) [S](-*)(-*)(=*)=*
C
      if(nc(i).ne.4.or.numh(i).ne.0) goto 3140
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.2.and.nbd3.eq.0) then
        xt=xt+8.38
C        write(16,*) '(36) [S](-*)(-*)(=*)=* ',i
        goto 1000
      end if
3140  continue
C
C     (37) [SH]-*
C
      if(nc(i).ne.1.or.numh(i).ne.1) goto 3150
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.0.and.nbd3.eq.0) then
        xt=xt+38.80
C        write(16,*) '(37) [SH]- ',i
        goto 1000
      end if
3150  continue
C      write(16,*) 'Warning, no atom type for noaromatic S atom',i
      goto 1000
C
3500  continue
C
C      aromatic S atom
C
C      (38) [s](:*):*
C
      if(nc(i).ne.2.or.numh(i).ne.0) goto 3160
      if(ibd.eq.2.and.nbd1.eq.0.and.nbd2.eq.0.and.nbd3.eq.0) then
         xt=xt+28.24
C         write(16,*) ' (38) [s](:*):* ',i
         goto 1000
      end if
3160  continue
C
C      (39) [s](=*)(:*):*
C
      if(nc(i).ne.3.or.numh(i).ne.0) goto 3170
      if(ibd.eq.2.and.nbd1.eq.0.and.nbd2.eq.1.and.nbd3.eq.0) then
        xt=xt+21.70
C        write(16,*) ' (39) [s](=*)(:*):* ',i
        goto 1000
      end if
3170  continue
C      write(16,*) 'warning, no atom type for aromatic S atom ',i
      goto 1000
C
C
4000  continue
C
C**************    P    atom    ********************
C
      if(number(i).ne.15) goto 1000
C
C     (40) [P](-*)(-*)-*
C
      if(nc(i).ne.3.or.numh(i).ne.0) goto 4110
      if(ibd.eq.0.and.nbd1.eq.3.and.nbd2.eq.0.and.nbd3.eq.0) then
        xt=xt+13.59
C        write(16,*) ' (40) [P](-*)(-*)-* ',i
        goto 1000
      end if
4110  continue
C
C     (41) [P](-*)=*
C
      if(nc(i).ne.2.or.numh(i).ne.0) goto 4120
      if(ibd.eq.0.and.nbd1.eq.1.and.nbd2.eq.1.and.nbd3.eq.0) then
        xt=xt+34.14
C        write(16,*) ' (41) [P](-*)=* ',i
        goto 1000
      end if
4120  continue
C
C     (42) [P](-*)(-*)(-*)=*
C
      if(nc(i).ne.4.or.numh(i).ne.0) goto 4130
      if(ibd.eq.0.and.nbd1.eq.3.and.nbd2.eq.1.and.nbd3.eq.0) then
        xt=xt+9.81
C        write(16,*) '(42) [P](-*)(-*)(-*)=* ',i
        goto 1000
      end if
4130  continue
C
C     (43) [PH](-*)(-*)=*
C
      if(nc(i).ne.3.or.numh(i).ne.1) goto 4140
      if(ibd.eq.0.and.nbd1.eq.2.and.nbd2.eq.1.and.nbd3.eq.0) then
        xt=xt+23.47
C        write(16,*) '(43) [PH](-*)(-*)=* ',i
        goto 1000
      end if
4140  continue
C
C
1000  continue
C
      write(16,1001) xt
1001  format('Topological polar surface area (TPSA):',e12.4)
      return 
      return
      end
