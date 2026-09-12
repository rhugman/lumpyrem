        program lumprem

        implicit none

        integer MAXDELAY
        parameter(MAXDELAY=500)
        integer, parameter :: MAX_LUMPREM_NAME=26
        integer, parameter :: MAXSSFSER=20

        integer i,numdrain,nummac,nstep,mxiter,numdays,noutdays,
     +  ierr,iline,iveg,nveg,idayold,ifin,imod,iday1,iday2,
     +  iday,jline,rechmod_call
        integer itemp
        integer ifile,ifail
        integer irrigcode_all
        integer narg
        integer numssfser,issfser
        integer ibeg,iend,nn,nb
        integer dds,mms,yys,dd,mm,yy,datespec,hhs,mins,sss
        integer nbucket
        integer lw(MAXDELAY),rw(MAXDELAY)
        integer lumpser(MAXSSFSER)
        double precision maxvol,rdelay,mdelay,ks,m,l,mflowmax,vol
        double precision tol,rtemp,
     +  cropfac1,gamma1,recharge,mrecharge,runoff,
     +  evapn,rainfall,rbal,deltav,oldvol,totd,totm,ototd,ototm,
     +  totrain,totevpn,totrech,totmrech,tottrech,totrun,totbal,
     +  potevapn,totepot,totirrigation,totgwithdrawal,totresidepot,
     +  totnetrecharge,totevpn_br,totepot_br
        double precision totdrainage_br,totoverflow_br
        double precision irrigation,gwithdrawal,irrigvolfrac
        double precision dfac,cdiff,gdiff,dday1,dday2,ddiff,iddiff
        double precision offset,factor1,factor2,power,dtemp,elevation
        double precision elevmin,elevmax
        double precision datum,dwt
        double precision cropfac_all,gamma_all
        double precision gwirrigfrac_all
        double precision maxvol_br,ks_br,l_br,m_br,vol_br,deltav_br
        double precision oldvol_br,drainage_br,overflow_br,extravol_br
        double precision allrecharge,epot_all_br,gamma_br
        double precision evapn_br,potevapn_br,gw_pot_evap
        double precision drainsub(maxdelay),macsub(maxdelay)
        double precision scaleser(MAXSSFSER),offsetser(MAXSSFSER),
     +                   minvalser(MAXSSFSER),maxvalser(MAXSSFSER)
        character*1 aa
        character*4 aext
        character*8 timeser
        character*10 aline
        character*15 lumprem_name(MAX_LUMPREM_NAME)
        character*15 ssfser(MAXSSFSER)
        character*15 adate,atime,anum
        character*20 atemp20,asite,asection
        character*25 atemp25
        character*100 atext
        character*200 atemp200a,atemp200b
        character*256 infle,outfle,cgfle,rainfle,epotfle,irrigfle,afile
        character*256 ssffle,csvfle,epotfle_br,recfle
        character*600 cline
        character*1000 errmsg

        integer, allocatable :: ioutday(:),vday(:)
        integer, allocatable :: irrigcode(:)
        double precision, allocatable :: cropfac(:),gamma(:),rain(:),
     +  epot(:),epot_br(:)
        double precision, allocatable :: cropfac_day(:),gamma_day(:)
        double precision, allocatable :: gwirrigfrac(:)
        double precision, allocatable :: valser(:,:)
        character*10, allocatable     :: dateser(:)

        data lumprem_name /
     +  'vol_upper      ','vol_lower      ','vol_drain      ',
     +  'vol_macro      ','del_vol_upper  ','del_vol_lower  ',
     +  'del_vol_drain  ',
     +  'del_vol_macro  ','rainfall       ','irrigation     ',
     +  'drain_upper    ','macro_upper    ','drain_lower    ',
     +  'overflow_lower ','total_rech     ',
     +  'gw_withdrawal  ','net_recharge   ','runoff         ',
     +  'pot_evap_upper ','evap_upper     ',
     +  'pot_evap_lower ','evap_lower     ',
     +  'gw_pot_evap    ',
     +  'balance        ','elevation      ','depth-to-water'/

C -- Obtain names of input files.

        narg=iargc()
        if(narg.eq.1)then
          write(6,9)
9         format(' LUMPREM is run using the command:')
          write(6,91)
91        format(/,'    LUMPREM infile outfile [csvoutfile]')
          write(6,911)
911       format(/,' where')
          write(6,912)
912       format(/,'    infile     is the name of its input file,')
          write(6,913)
913       format(  '    outfile    is the name of its tabular output ',
     +    'file, and')
          write(6,914)
914       format(  '    csvoutfile is the name of an optional CSV ',
     +    'output file.')
          write(6,915)
915       format(/,' Or provide no command line arguments ',
     +    'and let LUMPREM prompt for these.')
          stop
        else if(narg.ge.2)then
          call getarg(1,infle)
          call getarg(2,outfle)
          open(unit=10,file=infle,status='old',action='read',
     +    iostat=ierr)
          if(ierr.ne.0)then
            write(errmsg,15) trim(infle)
15          format('Cannot open file ',a,' to read LUMPREM ',
     +      'input dataset.')
            go to 9998
          end if
          if(narg.eq.2)then
            csvfle=' '
          else if(narg.eq.3)then
            call getarg(3,csvfle)
            nb=len_trim(csvfle)
            nn=nb-3
            if(nn.lt.1)nn=1
            aext=csvfle(nn:nb)
            call lowcase(aext)
            if(aext.ne.'.csv')then
              write(errmsg,16)
16            format('If third filename is supplied on LUMPREM command ',
     +        'line, it must have an extension of ".csv".')
              go to 9998
            end if
          else
            write(6,9)
            write(6,91)
            write(6,911)
            write(6,912)
            write(6,913)
            write(6,914)
            write(6,915)
            stop
          end if
        else if(narg.eq.0)then
40        write(6,50,advance='no')
50        format(' Enter name of LUMPREM input file: ')
          read(5,*) infle
          open(unit=10,file=infle,status='old',action='read',err=40)
69        write(6,70,advance='no')
70        format(' Enter name for LUMPREM tabular output file: ')
          read(5,*) outfle
701       write(6,702,advance='no')
702       format(' Enter name for CSV output file (<Enter> if none): ')
          read(5,'(a)') csvfle
          if(csvfle.ne.' ')then
            cline=csvfle
            ibeg=1
            iend=len_trim(cline)
            call getfile(ifail,cline,csvfle,ibeg,iend)
            if(ifail.ne.0) go to 701
            nb=len_trim(csvfle)
            nn=nb-3
            if(nn.lt.1)nn=1
            aext=csvfle(nn:nb)
            call lowcase(aext)
            if(aext.ne.'.csv')then
              write(6,703)
703           format(' *** Requires an extension ',
     +        'of ".csv" - try again ***')
              go to 701
            end if
          end if
          write(6,*)
        end if

C -- Initialize accumulators.

        totrain=0.0
        totevpn=0.0
        totevpn_br=0.0
        totrech=0.0
        totmrech=0.0
        tottrech=0.0
        totrun=0.0
        totbal=0.0
        totepot=0.0
        totepot_br=0.0
        totirrigation=0.0
        totgwithdrawal=0.0
        totresidepot=0.0
        totnetrecharge=0.0
        totdrainage_br=0.0
        totoverflow_br=0.0

C -- Initialize certain variables.

        do i=1,maxdelay
          drainsub(i)=0.0
          macsub(i)=0.0
        end do
        numssfser=0
        ssffle=' '
        epot_all_br=-1.0d301
        epotfle_br=' '
        iline=0
        recfle='lumprem_variables.rec'
        rechmod_call=0

C -- Read soil properties.

1       iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1
        call lowcase(cline)
        if(cline(1:7).ne.'* earth')then
          atext='* earth properties'
          go to 9550
        end if

2       iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 2
        call writint(aline,iline)
        call linesplit(ifail,2,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call drealread(ifail,cline(lw(1):rw(1)),maxvol)
        if(ifail.ne.0) go to 9400
        if(maxvol.le.0.0d0)then
          write(errmsg,99) trim(aline),trim(infle)
99        format('MAXVOL must be greater than zero at line ',a,
     +    ' of file ',a,'.')
          go to 9998
        end if
        call drealread(ifail,cline(lw(2):rw(2)),irrigvolfrac)
        if(ifail.ne.0) go to 9400
        if((irrigvolfrac.gt.1.0).or.(irrigvolfrac.lt.0.0))then
          write(errmsg,101) trim(aline),trim(infle)
101       format('IRRIGVOLFRAC must not be greater than 1.0 ',
     +    'or less than 0.0 at line ',a,' of file ',a,'.')
          go to 9998
        end if
        call linesplit(ifail,3,lw,rw,cline)
        if(ifail.ne.0)then
          maxvol_br=0.0d0
          extravol_br=0.0d0
          gamma_br=0.0d0
        else
          call drealread(ifail,cline(lw(3):rw(3)),maxvol_br)
          if(ifail.ne.0)then
            write(errmsg,5) trim(aline),trim(infle)
5           format('Cannot read value for MAXVOL_BR at line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          if(maxvol_br.lt.0.0d0) maxvol_br=0.0d0
          if(maxvol_br.gt.0.0d0)then
            call linesplit(ifail,4,lw,rw,cline)
            if(ifail.ne.0)then
              write(errmsg,511) trim(aline),trim(infle)
511           format('If MAXVOL_BR is greater than zero, then a value ',
     +        'must be supplied for EXTRAVOL_BR at line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
            call drealread(ifail,cline(lw(4):rw(4)),extravol_br)
            if(ifail.ne.0)then
              write(errmsg,512) trim(aline),trim(infle)
512           format('Cannot read value for EXTRAVOL_BR at line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
            if(extravol_br.lt.0.0d0)then
              write(errmsg,513) trim(aline),trim(infle)
513           format('EXTRAVOL_BR must be greater than or equal to ',
     +        'zero at line ',a,' of file ',a,'.')
              go to 9998
            end if
          else
            extravol_br=0.0d0
            gamma_br=0.0d0
          end if
          if(maxvol_br.gt.0.0d0)then
            call linesplit(ifail,5,lw,rw,cline)
            if(ifail.ne.0)then
              write(errmsg,5111) trim(aline),trim(infle)
5111          format('If MAXVOL_BR is non-zero then a value ',
     +        'must be supplied for GAMMA_BR at line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
            call drealread(ifail,cline(lw(5):rw(5)),gamma_br)
            if(ifail.ne.0)then
              write(errmsg,5121) trim(aline),trim(infle)
5121          format('Cannot read value for GAMMA_BR at line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
            if(gamma_br.gt.10.0) gamma_br=10.0d0
            if(gamma_br.lt.0.1) gamma_br=0.1
          end if
        end if

1011    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1011
        call linesplit(ifail,2,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call drealread(ifail,cline(lw(1):rw(1)),rdelay)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(2):rw(2)),mdelay)
        if(ifail.ne.0) go to 9400
        if(rdelay.gt.MAXDELAY-2)rdelay=MAXDELAY-2
        if(mdelay.gt.MAXDELAY-2)mdelay=MAXDELAY-2

1012    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1012
        call linesplit(ifail,4,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call drealread(ifail,cline(lw(1):rw(1)),ks)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(2):rw(2)),m)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(3):rw(3)),l)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(4):rw(4)),mflowmax)
        if(ifail.ne.0) go to 9400
        if(maxvol_br.le.0.0d0)then
          ks_br=0.0d0
          m_br=0.0d0
          l_br=0.0d0
        else
          call linesplit(ifail,7,lw,rw,cline)
          if(ifail.ne.0)then
            call writint(aline,iline)
            write(errmsg,71) trim(aline),trim(infle)
71          format('If MAXVOL_BR is supplied, and greater than zero, ',
     +      'then values must be provided for KS_BR, M_BR and L_BR ',
     +      'on line ',a,' of file ',a,'.')
            go to 9998
          end if
          call drealread(ifail,cline(lw(5):rw(5)),ks_br)
          if(ifail.ne.0) go to 9400
          call drealread(ifail,cline(lw(6):rw(6)),m_br)
          if(ifail.ne.0) go to 9400
          call drealread(ifail,cline(lw(7):rw(7)),l_br)
          if(ifail.ne.0) go to 9400
        end if

C -- Read volume to elevation conversion parameters

1013    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1013
        call lowcase(cline)
        if(cline(1:8).ne.'* volume') then
          atext='* volume to elevation'
          go to 9550
        end if

1014    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1014
        call linesplit(ifail,4,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call drealread(ifail,cline(lw(1):rw(1)),offset)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(2):rw(2)),factor1)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(3):rw(3)),factor2)
        if(ifail.ne.0) go to 9400
        call drealread(ifail,cline(lw(4):rw(4)),power)
        if(ifail.ne.0) go to 9400
        cline=cline(rw(4)+1:)
        cline=adjustl(cline)
        call lowcase(cline)
        nbucket=0
        nn=index(cline,'upper')
        if(nn.ne.0)then
          if(nbucket.eq.2)then
            call writint(aline,iline)
            go to 9998
          end if
          nbucket=1
          cline(nn:nn+4)=' '
        end if
        nn=index(cline,'lower')
        if(nn.ne.0)then
          if(nbucket.eq.1)then
            call writint(aline,iline)
            write(errmsg,721) trim(aline),trim(infle)
721         format('"lower" and "upper" cannot appear on the same ',
     +      'line at line ',a,' of file ',a,'.')
            go to 9998
            nbucket=2
            cline(nn:nn+4)=' '
          end if
          if(maxvol_br.le.0.0d0)then
            call writint(aline,iline)
            write(errmsg,722) trim(aline),trim(infle)
722         format('The string "lower" must not appear on line ',a,
     +      ' of file ',a,' unless MAXVOL_BR is greater than zero, ',
     +      'for this is what activates the lower LUMPREM bucket.')
            go to 9998
          end if
          nbucket=2
        end if
        if(nbucket.eq.0)nbucket=1
        call linesplit(ifail,1,lw,rw,cline)
        if(ifail.ne.0)then
          elevmin=-1.0d20
          elevmax=1.0d20
        else
          call drealread(ifail,cline(lw(1):rw(1)),elevmin)
          if(ifail.ne.0) then
            call writint(aline,iline)
            write(errmsg,7221) trim(aline),trim(infle)
7221        format('Cannot read value for ELEVMIN from line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          call linesplit(ifail,2,lw,rw,cline)
          if(ifail.ne.0)then
            elevmax=1.0d20
          else
            call drealread(ifail,cline(lw(2):rw(2)),elevmax)
            if(ifail.ne.0)then
              call writint(aline,iline)
              write(errmsg,7222) trim(aline),trim(infle)
7222          format('Cannot read value for ELEVMAX from line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
          end if
        end if

1015    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1015
        call lowcase(cline)
        if(cline(1:6).ne.'* topo') then
          atext='* topographic surface'
          go to 9550
        end if

1016    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1016
        call linesplit(ifail,1,lw,rw,cline)
        call drealread(ifail,cline(lw(1):rw(1)),datum)
        if(ifail.ne.0) go to 9400

C -- Read initial conditions.

1017    iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1017
        call lowcase(cline)
        if(cline(1:6).ne.'* init') then
          atext='* initial conditions'
          go to 9550
        end if

75      iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 75
        call linesplit(ifail,1,lw,rw,cline)
        call drealread(ifail,cline(lw(1):rw(1)),vol)
        if(ifail.ne.0) go to 9400
        if(vol.lt.0.0d0)then
          call writint(aline,iline)
          write(errmsg,751) trim(aline),trim(infle)
751       format('Initial volume must not be negative at line ',a,
     +    ' of file ',a,'.')
          go to 9998
        end if
        if(vol.gt.maxvol)then
          call writint(aline,iline)
          write(errmsg,752) trim(aline),trim(infle)
752       format('Initial volume of upper bucket must not exceed ',
     +    'MAXVOL at line ',a,' of file ',a,'.')
          go to 9998
        end if
        if(maxvol_br.le.0.0d0)then
          vol_br=0.0d0
        else
          call linesplit(ifail,2,lw,rw,cline)
          if(ifail.ne.0)then
            call writint(aline,iline)
            write(errmsg,77) trim(aline),trim(infle)
77          format('If MAXVOL_BR is supplied, and greater than zero, ',
     +      'then a value must be provided for VOL_BR ',
     +      'on line ',a,' of file ',a,'.')
            go to 9998
          end if
          call drealread(ifail,cline(lw(2):rw(2)),vol_br)
          if(ifail.ne.0) go to 9400
          if(vol_br.gt.maxvol_br+extravol_br)then
            call writint(aline,iline)
            write(errmsg,771) trim(aline),trim(infle)
771         format('The value supplied for VOL_BR at line ',a,
     +      ' of file ',a,' must not exceed MAXVOL_BR+EXTRAVOL_BR.')
            go to 9998
          end if
          if(vol_br.lt.0.0d0)then
            call writint(aline,iline)
            write(errmsg,751) trim(aline),trim(infle)
            go to 9998
          endif
        end if

772     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 772
        call linesplit(ifail,2,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call intread(ifail,cline(lw(1):rw(1)),numdrain)
        if(ifail.ne.0) go to 9400
        call intread(ifail,cline(lw(2):rw(2)),nummac)
        if(ifail.ne.0) go to 9400
        if((numdrain.gt.maxdelay).or.(nummac.gt.maxdelay)) go to 9000
        if((numdrain.le.0).or.(nummac.le.0))then
          call writint(aline,iline)
          write(errmsg,7721) trim(aline),trim(infle)
7721      format('NRBUF and NMBUF must both be greater than zero at',
     +    ' line ',a,' of file ',a,'.')
          go to 9998
        end if

773     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 773
        call linesplit(ifail,numdrain,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        do i=1,numdrain
          call drealread(ifail,cline(lw(i):rw(i)),drainsub(i))
          if(ifail.ne.0) go to 9420
          if(drainsub(i).lt.0.0d0) go to 9020
        end do

774     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 774
        call linesplit(ifail,nummac,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        do i=1,nummac
          call drealread(ifail,cline(lw(i):rw(i)),macsub(i))
          if(ifail.ne.0) go to 9440
          if(macsub(i).lt.0.0d0) go to 9020
        end do
        totd=0.0
        totm=0.0
        do 150 i=1,maxdelay
          totd=totd+drainsub(i)
          totm=totm+macsub(i)
150     continue

C -- Read solution parameters.

775     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 775
        call lowcase(cline)
        if(cline(1:7).ne.'* solut') then
          atext='* solution parameters'
          go to 9550
        end if

776     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 776
        call linesplit(ifail,3,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call intread(ifail,cline(lw(1):rw(1)),nstep)
        if(ifail.ne.0) go to 9460
        call intread(ifail,cline(lw(2):rw(2)),mxiter)
        if(ifail.ne.0) go to 9460
        call drealread(ifail,cline(lw(3):rw(3)),tol)
        if(ifail.ne.0) go to 9460
        if((nstep.le.0).or.(mxiter.le.0).or.(tol.le.0.0d0))then
          call writint(aline,iline)
          write(errmsg,7761) trim(aline),trim(infle)
7761      format('NSTEP, MXITER and TOL must all be greater than ',
     +    'zero at line ',a,' of file ',a,'.')
          go to 9998
        end if

C -- Read timing information

777     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 777
        call lowcase(cline)
        if(cline(1:8).ne.'* timing') then
          atext='* timing information'
          go to 9550
        end if

778     iline=iline+1
        read(10,'(a)',err=9400,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 778
        call linesplit(ifail,2,lw,rw,cline)
        if(ifail.ne.0) go to 9350
        call intread(ifail,cline(lw(1):rw(1)),numdays)
        if(ifail.ne.0) go to 9480
        call intread(ifail,cline(lw(2):rw(2)),noutdays)
        if(ifail.ne.0) go to 9480
        allocate(rain(numdays),epot(numdays),epot_br(numdays),
     +  ioutday(noutdays),stat=ierr)
        if(ierr.ne.0) go to 9050
        allocate(cropfac_day(numdays),gamma_day(numdays),stat=ierr)
        if(ierr.ne.0) go to 9050
        allocate(irrigcode(numdays),gwirrigfrac(numdays),stat=ierr)
        if(ierr.ne.0) go to 9050
        iline=iline+1
        read(10,*,err=9495,end=9495) (ioutday(i),i=1,noutdays)

C -- Read input data filenames. For two of these files parameter values may be provided instead of filenames.

        atext='* data filenames'
779     iline=iline+1
        read(10,'(a)',err=9550,end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 779
        call lowcase(cline)
        if(cline(1:6).ne.'* data') go to 9550
        jline=0
        asection='data filenames'

C -- Either a filename or crop factor and gamma values may be read from the next line.

        cgfle=' '
781     iline=iline+1
        jline=jline+1
        read(10,'(a)',end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 781
        cline=adjustl(cline)
        ifile=0
        if((cline(1:1).eq.'''').or.(cline(1:1).eq.'"'))ifile=1
        if(ifile.eq.0)then
          call linesplit(ifail,2,lw,rw,cline)
          if(ifail.ne.0)then
            ifile=1
            go to 155
          else
            call drealread(ifail,cline(lw(1):rw(1)),cropfac_all)
            if(ifail.ne.0) then
              ifile=1
              go to 155
            end if
            call drealread(ifail,cline(lw(2):rw(2)),gamma_all)
            if(ifail.ne.0)then
              call writint(aline,jline)
              write(errmsg,153) trim(aline),trim(infle)
153           format('Error: cannot read GAMMA value from line ',a,
     +        ' of "data filenames" section of file ',a,'.')
              go to 9998
            end if
          end if
        end if
155     continue
        if(ifile.eq.1)then
          ibeg=1
          iend=len_trim(cline)
          call getfile(ifail,cline,cgfle,ibeg,iend)
          if(ifail.ne.0)then
            call writint(aline,jline)
            write(errmsg,1531) trim(aline),trim(infle)
1531        format('Cannot read name of vegetation file on line ',a,
     +      ' of "data filenames" section of file ',a,'.')
            go to 9998
          end if
        end if

782     iline=iline+1
        jline=jline+1
        read(10,'(a)',err=9500,end=9500) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 782
        ibeg=1
        iend=len_trim(cline)
        call getfile(ifail,cline,rainfle,ibeg,iend)
        if(ifail.ne.0) go to 9500

783     iline=iline+1
        jline=jline+1
        read(10,'(a)',err=9520,end=9520) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 783
        ibeg=1
        iend=len_trim(cline)
        call getfile(ifail,cline,epotfle,ibeg,iend)
        if(ifail.ne.0) go to 9520

C -- Either a filename or irrigation code and groundwater extraction fraction may be read from the next line.

        irrigfle=' '
784     iline=iline+1
        jline=jline+1
        read(10,'(a)',end=9370) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 784
        cline=adjustl(cline)
        ifile=0
        if((cline(1:1).eq.'''').or.(cline(1:1).eq.'"'))ifile=1
        if(ifile.eq.0)then
          call linesplit(ifail,2,lw,rw,cline)
          if(ifail.ne.0)then
            ifile=1
            go to 158
          else
            call intread(ifail,cline(lw(1):rw(1)),irrigcode_all)
            if(ifail.ne.0) then
              ifile=1
              go to 158
            end if
            call drealread(ifail,cline(lw(2):rw(2)),gwirrigfrac_all)
            if(ifail.ne.0)then
              call writint(aline,jline)
              write(errmsg,156) trim(aline),trim(infle)
156           format('Cannot read value for GWIRRIGFRAC from line ',a,
     +        ' of "data filenames" section of file ',a,'.')
              go to 9998
            end if
          end if
        end if
158     continue
        if(ifile.eq.1)then
          ibeg=1
          iend=len_trim(cline)
          call getfile(ifail,cline,irrigfle,ibeg,iend)
          if(ifail.ne.0)then
            call writint(aline,jline)
            write(errmsg,1581) trim(aline),trim(infle)
1581        format('Cannot read name of irrigation file on line ',a,
     +      ' of "data filenames" section of file ',a,'.')
            go to 9998
          end if
        end if

! -- Now we look at extraction from the lower bucket.
!    This can be a single number, or it can be the name of a file.

        do
          iline=iline+1
          jline=jline+1
          read(10,'(a)',end=1251) cline
          call comstrip(cline)
          if(cline.ne.' ') exit
        end do
        if(cline(1:1).eq.'*') go to 1259
        call linesplit(ifail,1,lw,rw,cline)
        call drealread(ifail,cline(lw(1):rw(1)),epot_all_br)
        if(ifail.eq.0)then
          if(epot_all_br.lt.0.0d0)then
            call writint(aline,jline)
            write(errmsg,1510) trim(aline),trim(infle)
1510        format('EPOT_ALL_BR must be zero or greater at line ',a,
     +      ' of "data filenames" section of file ',a,'.')
            go to 9998
          end if
        else
          ibeg=1
          iend=len_trim(cline)
          call getfile(ifail,cline,epotfle_br,ibeg,iend)
          if(ifail.ne.0) then
            call writint(aline,jline)
            write(errmsg,1520) trim(aline),trim(infle)
1520        format('Cannot read name of lower bucket evaporation ',
     +      'file from line ',a,
     +      ' of "data filenames" section of file ',a,'.')
            go to 9998
          end if
          atemp200a=epotfle
          atemp200b=epotfle_br
          call lowcase(atemp200a)
          call lowcase(atemp200b)
          if(atemp200a.eq.atemp200b)then
            write(errmsg,1521) trim(atemp200a)
1521        format('EPOTFILE and EPOTFILE_BR must be different. They ',
     +      'are both "',a,'".')
            go to 9998
          end if
        end if

! -- We see if we need to read a "ssf file" section

        do
          iline=iline+1
          read(10,'(a)',end=1251) cline
          call comstrip(cline)
          if(cline.ne.' ') exit
        end do
        call lowcase(cline)
        cline=adjustl(cline)
        if(cline(1:1).ne.'*') go to 9600
1259    continue
        cline=adjustl(cline(2:))
        if(cline(1:4).ne.'ssf ') go to 9600
        jline=0
        asection='ssf file'

1260    iline=iline+1
        jline=jline+1
        read(10,'(a)',err=9620,end=9640) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1260
        ibeg=1
        iend=len_trim(cline)
        call getfile(ifail,cline,ssffle,ibeg,iend)
        if(ifail.ne.0)then
          write(errmsg,1265) trim(infle)
1265      format('Error reading name of site sample file ',
     +    'from "ssf file" section of file ',a,'.')
          go to 9998
        end if

1270    iline=iline+1
        jline=jline+1
        read(10,'(a)',err=9620,end=9640) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1270
        call linesplit(ifail,3,lw,rw,cline)
        if(ifail.ne.0) go to 9660
        adate=cline(lw(1):rw(1))
        call lowcase(adate)
        if(adate.eq.'dd/mm/yyyy')then
          datespec=1
        else if(adate.eq.'mm/dd/yyyy')then
          datespec=2
        else
          call writint(aline,jline)
          write(errmsg,1280) trim(aline),trim(asection),trim(infle)
1280      format('Date specification must be "dd/mm/yyyy" or ',
     +    '"mm/dd/yyyy" at line ',a,' of "',a,
     +    '" section of file ',a,'.')
          go to 9998
        end if
        adate=cline(lw(2):rw(2))
        call getdate2(ifail,adate,dds,mms,yys,datespec)
        if(ifail.ne.0)then
          call writint(aline,jline)
          write(errmsg,1290) trim(aline),trim(asection),trim(infle)
1290      format('Cannot read simulation starting date from line ',
     +    a,' of "',a,'" section of file ',a,'.')
          go to 9998
        end if
        atime=cline(lw(3):rw(3))
        call gettime(ifail,atime,hhs,mins,sss)
        if(ifail.ne.0)then
          call writint(aline,jline)
          write(errmsg,1301) trim(aline),trim(asection),trim(infle)
1301      format('Cannot read simulation starting time from line ',
     +    a,' of "',a,'" section of file ',a,'.')
          go to 9998
        end if

1311    iline=iline+1
        jline=jline+1
        read(10,'(a)',err=9620,end=9640) cline
        call comstrip(cline)
        if(cline.eq.' ') go to 1311
        call linesplit(ifail,1,lw,rw,cline)
        call intread(ifail,cline(lw(1):rw(1)),numssfser)
        if(ifail.ne.0) go to 9620
        if(numssfser.le.0) then
          numssfser=0
          go to 1251
        end if
        if(numssfser.gt.MAXSSFSER)then
          call writint(anum,MAXSSFSER)
          call writint(aline,jline)
          write(errmsg,1320) trim(anum),trim(aline),trim(asection),
     +    trim(infle)
1320      format('NUMSERIES can be a maximum of ',a,' at line ',a,
     +    ' of "',a,'" section of file ',a,'.')
          go to 9998
        end if
        do issfser=1,numssfser
1330      continue
          iline=iline+1
          jline=jline+1
          read(10,'(a)',err=9620,end=9640) cline
          call comstrip(cline)
          if(cline.eq.' ') go to 1330
          cline=adjustl(cline)
          call lowcase(cline)
          ibeg=1
          iend=len_trim(cline)
          call getfile(ifail,cline,atemp20,ibeg,iend)
          if(ifail.ne.0)then
            call writint(aline,jline)
            write(errmsg,1340) trim(aline),trim(asection),trim(infle)
1340        format('Cannot read LUMPREM series name at line ',a,
     +      ' of "',a,'" section of file ',a,'.')
            go to 9998
          end if
          atemp20=adjustl(atemp20)
          do i=1,MAX_LUMPREM_NAME
            if(atemp20.eq.lumprem_name(i)) go to 1345
          end do
          call writint(aline,jline)
          write(errmsg,1342) trim(atemp20),trim(aline),trim(asection),
     +    trim(infle)
1342      format('Series name "',a,'" at line ',a,' of "',a,
     +    '" section of file ',a,
     +    ' is not a LUMPREM output series name.')
          go to 9998
1345      lumpser(issfser)=i
          write(6,*) trim(lumprem_name(i))    !debug
          cline=adjustl(cline(iend+1:))
          ibeg=1
          iend=len_trim(cline)
          call getfile(ifail,cline,atemp25,ibeg,iend)
          if(ifail.ne.0)then
            call writint(aline,jline)
            write(errmsg,1350) trim(aline),trim(asection),trim(infle)
1350        format('Cannot read SSF series name at line ',a,
     +      ' of "',a,'" section of file ',a,'.')
            go to 9998
          end if
          atemp25=adjustl(atemp25)
          nn=len_trim(atemp25)
          if(nn.gt.15)then
            call writint(aline,jline)
            write(errmsg,1355) trim(atemp25),trim(aline),trim(asection),
     +      trim(infle)
1355        format('SSF series name "',a,'" exceeds 15 characters ',
     +      'at line ',a,' of "',a,'" section of file ',a,'.')
            go to 9998
          end if
          do i=1,nn
            if(atemp25(i:i).eq.' ')then
              call writint(aline,jline)
              write(errmsg,1396) trim(atemp25),trim(aline),
     +        trim(asection),trim(infle)
1396          format('SSF series name "',a,'" contains a space at ',
     +        'line ',a,' of "',a,'" section of file ',a,'.')
              go to 9998
            end if
          end do
          atemp25=adjustl(atemp25)
          ssfser(issfser)=atemp25
          if(issfser.gt.1)then
            do i=1,issfser-1
              if(atemp25.eq.ssfser(i))then
                call writint(aline,jline)
                write(errmsg,1397) trim(atemp25),trim(aline),
     +          trim(asection),trim(infle)
1397            format('Series name "',a,
     +          '" for SSF file is duplicated at line ',a,
     +          ' of "',a,'" section of file ',a,'.')
                go to 9998
              end if
            end do
          end if
          cline=adjustl(cline(iend+1:))
          call linesplit(ifail,4,lw,rw,cline)
          if(ifail.ne.0) go to 9660
          do i=1,4
            call drealread(ifail,cline(lw(i):rw(i)),dtemp)
            if(ifail.ne.0) then
              call writint(aline,jline)
              write(errmsg,1360) trim(aline),trim(asection),trim(infle)
1360          format('Error reading either SCALE, OFFSET, MINVAL ',
     +        'or MAXVAL from line ',a,' of "',a,
     +        '" section of file ',a,'.')
              go to 9998
            end if
            if(i.eq.1)then
              scaleser(issfser)=dtemp
            else if(i.eq.2)then
              offsetser(issfser)=dtemp
            else if(i.eq.3)then
              minvalser(issfser)=dtemp
            else if(i.eq.4)then
              maxvalser(issfser)=dtemp
            end if
          end do
          if(maxvalser(issfser).le.minvalser(issfser))then
            call writint(aline,jline)
            write(errmsg,1361) trim(aline),trim(asection),trim(infle)
1361        format('Maximum value does not exceed minimum value at ',
     +      'line ',a,' of "',a,'" section of file ',a,'.')
            go to 9998
          end if
        end do

1251    continue
        if((epot_all_br.lt.-1.0d300).and.(epotfle_br.eq.' '))then
          if(maxvol_br.gt.0.0d0)then
            write(errmsg,1252) trim(infle)
1252        format('If MAXVOL_BR is greater than zero, then a value ',
     +      'must be supplied for EPOT_ALL_BR (normally zero) or ',
     +      'EPOTFILE_BR in the "data filenames" section ',
     +      'of file ',a,'.')
            go to 9998
          end if
        end if
        close(unit=10)
        write(6,159) trim(infle)
159     format(' - file ',a,' read ok.')

! -- The LUMPREM variables file is now recorded.

        open(unit=30,file=recfle,action='write',iostat=ierr)
        if(ierr.ne.0)then
          write(errmsg,2510) trim(recfle)
2510      format('Cannot record LUMPREM variables file ',a,'.')
          go to 9998
        end if
        write(30,2520)
2520    format(' Variables used by LUMPREM')
        write(30,*)
        write(30,2530)
2530    format('* earth properties',/)
        write(30,2540) maxvol
2540    format('  MAXVOL't30,1pg14.7)
        write(30,2550) irrigvolfrac
2550    format('  IRRIGVOLFRAC't30,1pg14.7)
        if(maxvol_br.gt.0.0d0)then
          write(30,2560) maxvol_br
2560      format('  MAXVOL_BR't30,1pg14.7)
          write(30,2570) extravol_br
2570      format('  EXTRAVOL_BR't30,1pg14.7)
          write(30,2580) gamma_br
2580      format('  GAMMA_BR't30,1pg14.7)
        end if
        write(30,2590) rdelay
2590    format('  RDELAY't30,1pg14.7)
        write(30,2600) mdelay
2600    format('  MDELAY't30,1pg14.7)
        write(30,2610) ks
2610    format('  KS't30,1pg14.7)
        write(30,2620) m
2620    format('  M't30,1pg14.7)
        write(30,2621) l
2621    format('  L't30,1pg14.7)
        if(maxvol_br.gt.0.0d0)then
          write(30,2630) ks_br
2630      format('  KS_BR't30,1pg14.7)
          write(30,2640) m_br
2640      format('  M_BR't30,1pg14.7)
          write(30,2650) l_br
2650      format('  L_BR't30,1pg14.7)
        end if
        write(30,2660)
2660    format(/,'* volume to elevation',/)
        write(30,2670) offset
2670    format('  OFFSET't30,1pg14.7)
        write(30,2680) factor1
2680    format('  FACTOR1't30,1pg14.7)
        write(30,2690) factor2
2690    format('  FACTOR2't30,1pg14.7)
        write(30,2700) power
2700    format('  POWER't30,1pg14.7)
        write(30,2710) elevmin
2710    format('  ELEVMIN't30,1pg14.7)
        write(30,2720) elevmax
2720    format('  ELEVMAX't30,1pg14.7)
        if(maxvol_br.gt.0.0d0)then
          if(nbucket.eq.1)then
            atemp20='upper'
          else
            atemp20='lower'
          end if
          atemp20=adjustl(atemp20)
          write(30,2730) trim(atemp20)
2730      format('  BUCKET',t31,a)
        end if
        write(30,2740)
2740    format(/,'* topographic surface',/)
        write(30,2750) datum
2750    format('  SURFACE't30,1pg14.7)
        write(30,2760)
2760    format(/,' * initial conditions',/)
        write(30,2770) vol
2770    format('  VOL't30,1pg14.7)
        if(maxvol_br.gt.0.0d0)then
          write(30,2780) vol_br
2780      format('  VOL_BR't30,1pg14.7)
        end if
        write(30,2790) numdrain
2790    format('  NRBUF't30,i5)
        write(30,2800) nummac
2800    format('  NMBUF't30,i5)
        write(30,2810)
2810    format('  RBUF --->')
        write(30,2820) (drainsub(i),i=1,numdrain)
2820    format(2x,20(1pg14.7))
        write(30,2830)
2830    format('  MBUF --->')
        write(30,2820) (macsub(i),i=1,nummac)
        write(30,2840)
2840    format(/,'* timing information',/)
        write(30,2850) numdays
2850    format('  NUMDAYS',t30,i6)
        write(30,2860) noutdays
2860    format('  NOUTDAYS',t30,i6)
        write(30,2870)
2870    format('  OUTDAY --->')
        write(30,2880) (ioutday(i),i=1,noutdays)
2880    format(30i7)
        write(30,2890)
2890    format(/,'* data filenames',/)
        if(cgfle.eq.' ')then
          write(30,2900) cropfac_all
2900      format('  CROPFAC_ALL',t30,1pg14.7)
          write(30,2912) gamma_all
2912      format('  GAMMA_ALL',t30,1pg14.7)
        else
          write(30,2921) trim(cgfle)
2921      format('  VEGFILE',t30,a)
        end if
        write(30,2930) trim(rainfle)
2930    format('  RAINFILE',t30,a)
        write(30,2940) trim(epotfle)
2940    format('  EPOTFILE',t30,a)
        if(irrigfle.eq.' ')then
          write(30,2950) irrigcode_all
2950      format('  IRRIGCODE_ALL',t30,i5)
          write(30,2960) gwirrigfrac_all
2960      format('  GWIRRIGFRAC_ALL',t30,1pg14.7)
        else
          write(30,2970) trim(irrigfle)
2970      format('  IRRIGFILE',t30,a)
        end if
        if(maxvol_br.gt.0.0d0)then
          if(epotfle_br.ne.' ')then
            write(30,2980) trim(epotfle_br)
2980        format('  EPOTFILE_BR',t30,a)
          else
            write(30,2990) epot_all_br
2990        format('  EPOT_ALL_BR',t30,1pg14.7)
          end if
        end if
        if(ssffle.ne.' ')then
          write(30,29911)
29911     format(/,'* ssf file',/)
          write(30,29921) trim(ssffle)
29921     format('  SSFFILE',t30,a)
          if(datespec.eq.1)then
            atemp20='dd/mm/yyyy'
          else
            atemp20='mm/dd/yyyy'
          end if
          write(30,2993) trim(atemp20)
2993      format('  DATE_FORMAT',t30,a)
          if(datespec.eq.1)then
            write(atemp20,2994) dds,mms,yys
          else
            write(atemp20,2994) mms,dds,yys
2994        format(i2.2,'/',i2.2,'/',i4.4)
          end if
          atemp20=adjustl(atemp20)
          write(30,2995) trim(atemp20)
2995      format('  START_DATE',t30,a)
          write(atemp20,2996) hhs,mins,sss
2996      format(i2.2,':',i2.2,':',i2.2)
          atemp20=adjustl(atemp20)
          write(30,2997) trim(atemp20)
2997      format('  START_TIME',t30,a)
          write(30,2998) numssfser
2998      format('  NUMSSFSER',t30,i2)
          if(numssfser.gt.0)then
            cline=' '
            cline(3:)='LUMPSERNAME'
            cline(23:)='SSFSERNAME'
            cline(43:)='SCALE'
            cline(63:)='OFFSET'
            cline(83:)='MINVAL'
            cline(103:)='MAXVAL'
            write(30,'(a)') trim(cline)
            do i=1,len_trim(cline)
              if(cline(i:i).ne.' ')cline(i:i)='-'
            end do
            write(30,'(a)') trim(cline)
            do i=1,numssfser
              atemp20=adjustl(lumprem_name(lumpser(i)))
              cline(3:)=trim(atemp20)
              atemp20=adjustl(ssfser(i))
              cline(23:)=trim(atemp20)
              write(atemp20,2999) scaleser(i)
2999          format(1pg14.7)
              atemp20=adjustl(atemp20)
              cline(43:)=trim(atemp20)
              write(atemp20,2999) offsetser(i)
              atemp20=adjustl(atemp20)
              cline(63:)=trim(atemp20)
              write(atemp20,2999) minvalser(i)
              atemp20=adjustl(atemp20)
              cline(83:)=trim(atemp20)
              write(atemp20,2999) maxvalser(i)
              atemp20=adjustl(atemp20)
              cline(103:)=trim(atemp20)
              write(30,'(a)') trim(cline)
            end do
          end if
        end if

        close(unit=30)
        write(6,2050) trim(recfle)

C -- Open the vegetation and climate files.

        if(cgfle.ne.' ')
     +  open(unit=10,file=cgfle,status='old',err=9150)
        open(unit=11,file=rainfle,status='old',err=9200)
        open(unit=12,file=epotfle,status='old',err=9250)
        if(irrigfle.ne.' ')
     +  open(unit=13,file=irrigfle,status='old',err=9300)
        if(maxvol_br.gt.0.0d0)then
          if(epotfle_br.ne.' ')
     +    open(unit=14,file=epotfle_br,status='old',err=9320)
        end if

C -- The vegetation file is read.

        if(cgfle.eq.' ')then
          do iday=1,numdays
            cropfac_day(iday)=cropfac_all
          end do
          do iday=1,numdays
            gamma_day(iday)=gamma_all
          end do
          go to 475
        end if

        iline=1
330     continue
          read(10,*,end=350)
          iline=iline+1
        go to 330
350     allocate(vday(iline),cropfac(iline),gamma(iline),stat=ierr)
        if(ierr.ne.0) then
          write(errmsg,360)
360       format('Cannot allocate sufficient memory to continue ',
     +    'exection.')
          go to 9998
        endif
        rewind(unit=10)
        iveg=0
        iline=0
        afile=cgfle
380     continue
          iline=iline+1
          read(10,'(a)',end=450) cline
          if(cline.eq.' ') go to 380
          iveg=iveg+1
          call linesplit(ifail,3,lw,rw,cline)
          if(ifail.ne.0) go to 9590
          call intread(ifail,cline(lw(1):rw(1)),vday(iveg))
          if(ifail.ne.0) go to 9570
          call drealread(ifail,cline(lw(2):rw(2)),cropfac(iveg))
          if(ifail.ne.0) go to 9570
          if(cropfac(iveg).lt.0.0d0)then
            call writint(aline,iline)
            write(errmsg,381) trim(aline),trim(cgfle)
381         format('Crop factor is negative at line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          call drealread(ifail,cline(lw(3):rw(3)),gamma(iveg))
          if(ifail.ne.0) go to 9570
          if(gamma(iveg).lt.0.0d0)then
            call writint(aline,iline)
            write(errmsg,382) trim(aline),trim(cgfle)
382         format('Gamma is negative at line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          if(iveg.eq.1)then
            if(vday(iveg).ne.1) then
              write(errmsg,400) trim(cgfle)
400           format('The first day cited in vegetation file ',
     +        a,' must be day one of simulation.')
              go to 9998
            end if
          else
            if(vday(iveg).le.vday(iveg-1))then
              call writint(aline,iline)
              write(errmsg,401) trim(aline),trim(cgfle)
401           format('Consecutive entries in vegetation file ',
     +        'pertain to non-increasing days at line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
          end if
        go to 380
450     close (unit=10)
        nveg=iveg
        if(nveg.eq.0)then
          write(errmsg,451) trim(cgfle)
451       format('Vegetation file ',a,' is empty.')
          go to 9998
        end if
        if(vday(nveg).lt.numdays)then
          write(errmsg,470) trim(cgfle)
470       format('The last day cited in vegetation file ',a,' must ',
     +    'be greater than or equal to model simulation time.')
          go to 9998
        end if
        write(6,159) trim(cgfle)

C -- We write vegetation data to daily arrays.

        if(nveg.eq.1)then
          do iday=1,numdays
           cropfac_day(iday)=cropfac(1)
          end do
          do iday=1,numdays
           gamma_day(iday)=gamma(1)
          end do
        else
          do iveg=1,nveg-1
            dday1=dble(vday(iveg))
            dday2=dble(vday(iveg+1))
            cropfac1=cropfac(iveg)
            gamma1=gamma(iveg)
            ddiff=dday2-dday1
            iddiff=1.0d0/ddiff
            cdiff=cropfac(iveg+1)-cropfac(iveg)
            gdiff=gamma(iveg+1)-gamma(iveg)
            do iday=vday(iveg),vday(iveg+1)
              if(iday.lt.1) cycle
              if(iday.gt.numdays) go to 475
              dfac=(dble(iday)-dday1)*iddiff
              cropfac_day(iday)=cropfac1+cdiff*dfac
              gamma_day(iday)=gamma1+gdiff*dfac
            end do
          end do
        end if
475     continue

C -- The rainfall file is read.

        afile=rainfle
        iline=0
        idayold=0
550     continue
          iline=iline+1
          read(11,'(a)',end=700) cline
          if(cline.eq.' ') go to 550
          call linesplit(ifail,2,lw,rw,cline)
          if(ifail.ne.0) go to 9590
          call intread(ifail,cline(lw(1):rw(1)),iday)
          if(ifail.ne.0) go to 9570
          call drealread(ifail,cline(lw(2):rw(2)),rtemp)
          if(ifail.ne.0) go to 9570
          if(rtemp.lt.0.0d0)then
            call writint(aline,iline)
            write(errmsg,551) trim(aline),trim(rainfle)
551         format('Rainfall is less than zero at line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          if(iday.lt.1) then
            call writint(aline,iline)
            write(errmsg,570) trim(aline),trim(rainfle)
570         format('Zero or negative day cited at line ',a,
     +      ' of rainfall file ',a,'.')
            go to 9998
          end if
          if(iday.le.idayold)then
            call writint(aline,iline)
            write(errmsg,590) trim(aline),trim(infle)
590         format('Days are not in correct order at line ',a,
     +      ' of rainfall file ',a,'.')
            go to 9998
          end if
          if(iday.le.numdays) rain(iday)=rtemp
          if(idayold.lt.numdays)then
            if(idayold.lt.iday-1)then
              do 610 i=idayold+1,min(iday-1,numdays)
                rain(i)=0.0
610           continue
            end if
          end if
          idayold=iday
        go to 550
700     close(unit=11)
        if(idayold.lt.numdays)then
          do 750 i=idayold+1,numdays
            rain(i)=0.0
750       continue
        end if
        write(6,159) trim(rainfle)

C -- The evaporation file is read.

        iline=0
        idayold=0
        afile=epotfle
850     continue
          iline=iline+1
          read(12,'(a)',end=1000) cline
          if(cline.eq.' ') go to 850
          call linesplit(ifail,2,lw,rw,cline)
          if(ifail.ne.0) go to 9590
          call intread(ifail,cline(lw(1):rw(1)),iday)
          if(ifail.ne.0) go to 9570
          call drealread(ifail,cline(lw(2):rw(2)),rtemp)
          if(ifail.ne.0) go to 9570
          if(iday.lt.1) then
            call writint(aline,iline)
            write(errmsg,870) trim(aline),trim(afile)
            go to 9998
          end if
          if(idayold.eq.0)then
            if(iday.ne.1) then
              write(errmsg,880) trim(afile)
              go to 9998
            end if
          end if
          if(iday.le.idayold)then
            call writint(aline,iline)
            write(errmsg,890) trim(aline),trim(afile)
            go to 9998
          end if
          if(iday.le.numdays) epot(iday)=rtemp
          if(iday.eq.1) go to 920
          if(idayold.le.numdays)then
            if(idayold.lt.iday-1) then
              rtemp=epot(idayold)
              do 910 i=idayold+1,min(iday-1,numdays)
                epot(i)=rtemp
910           continue
            end if
          end if
920       idayold=iday
        go to 850
1000    close(unit=12)
        if(idayold.lt.numdays)then
          rtemp=epot(idayold)
          do 1050 i=idayold+1,numdays
            epot(i)=rtemp
1050      continue
        end if
        write(6,159) trim(epotfle)

C -- The irrigation file is read.

        if(irrigfle.eq.' ')then
          do iday=1,numdays
            irrigcode(iday)=irrigcode_all
          end do
          do iday=1,numdays
            gwirrigfrac(iday)=gwirrigfrac_all
          end do
          go to 2310
        end if

        iline=0
        idayold=0
        afile=irrigfle
2200    continue
          iline=iline+1
          read(13,'(a)',end=2300) cline
          if(cline.eq.' ') go to 2200
          call linesplit(ifail,3,lw,rw,cline)
          if(ifail.ne.0) go to 9590
          call intread(ifail,cline(lw(1):rw(1)),iday)
          if(ifail.ne.0) go to 9570
          call intread(ifail,cline(lw(2):rw(2)),itemp)
          if(ifail.ne.0) go to 9570
          call drealread(ifail,cline(lw(3):rw(3)),dtemp)
          if(ifail.ne.0) go to 9570
          if(iday.lt.1) then
            call writint(aline,iline)
            write(errmsg,2210) trim(aline),trim(irrigfle)
2210        format('Zero or negative day at line ',a,
     +      ' of irrigation file ',a,'.')
            go to 9998
          end if
          if(idayold.eq.0)then
            if(iday.ne.1) then
              write(errmsg,2220) trim(irrigfle)
2220          format('First line of irrigation file ',a,
     +        ' must pertain to day one.')
              go to 9998
            end if
          end if
          if(iday.le.idayold)then
            call writint(aline,iline)
            write(errmsg,2230) trim(aline),trim(irrigfle)
2230        format('Days are not in increasing order at line ',a,
     +      ' of irrigation file ',a,'.')
            go to 9998
          end if
          if(iday.le.numdays) then
            irrigcode(iday)=itemp
            gwirrigfrac(iday)=dtemp
          end if
          if((itemp.ne.0).and.(itemp.ne.1))then
            call writint(aline,iline)
            write(errmsg,2240) trim(aline),trim(irrigfle)
2240        format('Irrigation code must be zero or one at line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          if((dtemp.lt.0.0).or.(dtemp.gt.1.0))then
            call writint(aline,iline)
            write(errmsg,2250) trim(aline),trim(irrigfle)
2250        format('Groundwater irrigation source fraction must ',
     +      'not be less than zero or greater than one at line ',a,
     +      ' of file ',a,'.')
            go to 9998
          end if
          if(iday.eq.1) go to 2290
          if(idayold.le.numdays)then
            if(idayold.lt.iday-1) then
              itemp=irrigcode(idayold)
              dtemp=gwirrigfrac(idayold)
              do i=idayold+1,min(iday-1,numdays)
                irrigcode(i)=itemp
                gwirrigfrac(i)=dtemp
              end do
            end if
          end if
2290      continue
          idayold=iday
        go to 2200
2300    close(unit=13)
        if(idayold.lt.numdays)then
          itemp=irrigcode(idayold)
          dtemp=gwirrigfrac(idayold)
          do i=idayold+1,numdays
            irrigcode(i)=itemp
            gwirrigfrac(i)=dtemp
          end do
        end if
        write(6,159) trim(irrigfle)
2310    continue

C -- The below-root-zone-evaporation file is read.

        if((maxvol_br.gt.0.0d0).and.(epotfle_br.ne.' '))then
          iline=0
          idayold=0
          afile=epotfle_br
1850      continue
            iline=iline+1
            read(14,'(a)',end=2001) cline
            if(cline.eq.' ') go to 1850
            call linesplit(ifail,2,lw,rw,cline)
            if(ifail.ne.0) go to 9590
            call intread(ifail,cline(lw(1):rw(1)),iday)
            if(ifail.ne.0) go to 9570
            call drealread(ifail,cline(lw(2):rw(2)),rtemp)
            if(ifail.ne.0) go to 9570
            if(iday.lt.1) then
              call writint(aline,iline)
              write(errmsg,870) trim(aline),trim(afile)
870           format('Error: zero or negative day at line ',a,
     +        ' of file ',a,'.')
              go to 9998
            end if
            if(idayold.eq.0)then
              if(iday.ne.1) then
                write(errmsg,880) trim(afile)
880            format('Error: first line of file ',a,' must pertain',
     +         ' to day one.')
                go to 9998
              end if
            end if
            if(iday.le.idayold)then
              call writint(aline,iline)
              write(errmsg,890) trim(aline),trim(afile)
890           format('Error: days are not in correct order at ',
     +        'line ',a,' of file ',a,'.')
              go to 9998
            end if
            if(iday.le.numdays) epot_br(iday)=rtemp
            if(iday.eq.1) go to 1920
            if(idayold.le.numdays)then
              if(idayold.lt.iday-1) then
                rtemp=epot_br(idayold)
                do 1910 i=idayold+1,min(iday-1,numdays)
                  epot_br(i)=rtemp
1910            continue
              end if
            end if
1920        idayold=iday
          go to 1850
2001      close(unit=14)
          if(idayold.lt.numdays)then
            rtemp=epot_br(idayold)
            do 2051 i=idayold+1,numdays
              epot_br(i)=rtemp
2051        continue
          end if
          write(6,159) trim(epotfle_br)
        else
          epot_br=epot_all_br    ! An array
        end if

C -- The header is written to the LUMPREM2 output file and initial conditions
C    are recorded.

        if(nbucket.eq.1)then
          dtemp=min(vol,maxvol)
        else
          dtemp=min(vol_br,maxvol_br)
        end if
        if(dtemp.lt.1.0e-10) dtemp=1.0e-10                       ! arbitrary - needed for a negative power
        elevation=offset+factor1*dtemp+factor2*(dtemp**power)
        if(elevation.lt.elevmin) elevation=elevmin
        if(elevation.gt.elevmax) elevation=elevmax
        dwt=datum-elevation
        open(unit=10,file=outfle,action='write',iostat=ierr)
        if(ierr.ne.0)then
          write(errmsg,1249) trim(outfle)
1249      format('Cannot write to file ',a,'.')
          go to 9998
        end if
        write(10,1250) (lumprem_name(i),i=1,MAX_LUMPREM_NAME)
1250    format('  days  ',1000a15)
        write(10,1261) 0,vol,vol_br,totd,totm,0.0,0.0,0.0,0.0,0.0,0.0,
     +  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     +  elevation,dwt
1261    format(i6,1000(1x,1pg14.7))
        if(csvfle.ne.' ')then
          open(unit=20,file=csvfle,action='write',iostat=ierr)
          if(ierr.ne.0)then
            write(errmsg,1249) trim(csvfle)
            go to 9998
          end if
          write(20,1263)'days',(trim(lumprem_name(i)),
     +    i=1,MAX_LUMPREM_NAME)
1263      format(a,1000(',',a))
          write(cline,1262) 0,vol,vol_br,totd,totm,
     +    0.0,0.0,0.0,0.0,0.0,0.0,
     +    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     +    elevation,dwt
1262      format(i6,26(',',1pg14.7))
          call remspace(cline)
          write(20,'(a)') trim(cline)
        end if

! -- Memory is allocated if required for SSF file output.

        if(numssfser.ne.0)then
          allocate(dateser(0:noutdays),valser(0:noutdays,numssfser),
     +    stat=ierr)
          if(ierr.ne.0) go to 9050
          do issfser=1,numssfser
            select case(lumpser(issfser))
              case(1)
                valser(0,issfser)=vol
              case(2)
                valser(0,issfser)=vol_br
              case(3)
                valser(0,issfser)=totd
              case(4)
                valser(0,issfser)=totm
              case(5)
                valser(0,issfser)=0.0
              case(6)
                valser(0,issfser)=0.0
              case(7)
                valser(0,issfser)=0.0
              case(8)
                valser(0,issfser)=0.0
              case(9)
                valser(0,issfser)=0.0
              case(10)
                valser(0,issfser)=0.0
              case(11)
                valser(0,issfser)=0.0
              case(12)
                valser(0,issfser)=0.0
              case(13)
                valser(0,issfser)=0.0
              case(14)
                valser(0,issfser)=0.0
              case(15)
                valser(0,issfser)=0.0
              case(16)
                valser(0,issfser)=0.0
              case(17)
                valser(0,issfser)=0.0
              case(18)
                valser(0,issfser)=0.0
              case(19)
                valser(0,issfser)=0.0
              case(20)
                valser(0,issfser)=0.0
              case(21)
                valser(0,issfser)=0.0
              case(22)
                valser(0,issfser)=0.0
              case(23)
                valser(0,issfser)=0.0
              case(24)
                valser(0,issfser)=0.0
              case(25)
                valser(0,issfser)=elevation
              case(26)
                valser(0,issfser)=dwt
            end select
          end do

! -- Evaluate the date strings and the single time string.

          write(timeser,1400) hhs,mins,sss
1400      format(i2.2,':',i2.2,':'i2.2)
          timeser=adjustl(timeser)
          if(datespec.eq.1)then
            write(dateser(0),1410) dds,mms,yys
          else
            write(dateser(0),1410) mms,dds,yys
          end if
          dateser(0)=adjustl(dateser(0))
          do i=1,noutdays
            call newdate(ioutday(i),dds,mms,yys,dd,mm,yy)
            if(datespec.eq.1)then
              write(dateser(i),1410) dd,mm,yy
            else
              write(dateser(i),1410) mm,dd,yy
            end if
            dateser(i)=adjustl(dateser(i))
1410        format(i2.2,'/',i2.2,'/',i4.4)
          end do
        end if

C -- Subroutine RECHMOD is called once for every time at which output is
C    required.

        ifin=0
        imod=0
        oldvol=vol
        oldvol_br=vol_br
        ototd=totd
        ototm=totm
1300    continue
          imod=imod+1
          if(imod.gt.noutdays) go to 2000
          if(imod.eq.1)then
            iday1=0
          else
            iday1=idayold
          end if
          iday2=ioutday(imod)
          if(iday2.ge.numdays)then
            iday2=numdays
            ifin=1
          end if
          rechmod_call=rechmod_call+1
          call rechmod(rechmod_call,
     +    mxiter,tol,iday2-iday1,maxdelay,nstep,maxvol,
     +    irrigvolfrac,cropfac_day(iday1+1),gamma_day(iday1+1),ks,l,m,
     +    mflowmax,rdelay,
     +    mdelay,vol,drainsub,macsub,rain(iday1+1),epot(iday1+1),
     +    recharge,mrecharge,runoff,evapn,rainfall,potevapn,
     +    irrigation,gwithdrawal,
     +    irrigcode(iday1+1),gwirrigfrac(iday1+1),
     +    maxvol_br,extravol_br,ks_br,l_br,m_br,vol_br,drainage_br,
     +    overflow_br,epot_br(iday1+1),gamma_br,evapn_br,potevapn_br)

C -- Some aspects of the water balance are calculated.

          deltav=vol-oldvol
          deltav_br=vol_br-oldvol_br
          totd=0.0
          totm=0.0
          do 1700 i=1,maxdelay
            totd=totd+drainsub(i)
            totm=totm+macsub(i)
1700      continue
          if(maxvol_br.gt.0.0d0)then
            allrecharge=mrecharge+drainage_br+overflow_br
          else
            allrecharge=mrecharge+recharge
          end if
          rbal=rainfall+irrigation-deltav-deltav_br-allrecharge-
     +         runoff-evapn-(totd-ototd)-(totm-ototm)-evapn_br

C -- Elevation and depth to water table are calculated.

          if(nbucket.eq.1)then
            dtemp=min(vol,maxvol)
          else
            dtemp=min(vol_br,maxvol_br)
          end if
          if(dtemp.lt.1.0e-10) dtemp=1.0e-10                       ! arbitrary - needed for a negative power
          elevation=offset+factor1*dtemp+factor2*(dtemp**power)
          if(elevation.lt.elevmin) elevation=elevmin
          if(elevation.gt.elevmax) elevation=elevmax
          dwt=datum-elevation

C -- Data is written to the output file.

!          gw_pot_evap=potevapn-evapn-evapn_br
          gw_pot_evap=potevapn-evapn                     ! We don't subtract evap from the lower bucket
          write(10,1261) iday2,vol,vol_br,totd,totm,deltav,deltav_br,
     +    (totd-ototd),
     +    (totm-ototm),rainfall,irrigation,recharge,
     +    mrecharge,drainage_br,overflow_br,allrecharge,gwithdrawal,
     +    allrecharge-gwithdrawal,runoff,
     +    potevapn,evapn,potevapn_br,evapn_br,
     +    gw_pot_evap,rbal,elevation,dwt
          idayold=iday2
          if(csvfle.ne.' ')then
          write(cline,1262) iday2,vol,vol_br,totd,totm,deltav,deltav_br,
     +    (totd-ototd),
     +    (totm-ototm),rainfall,irrigation,recharge,
     +    mrecharge,drainage_br,overflow_br,allrecharge,gwithdrawal,
     +    allrecharge-gwithdrawal,runoff,
     +    potevapn,evapn,potevapn_br,evapn_br,
     +    gw_pot_evap,rbal,elevation,dwt
          call remspace(cline)
          write(20,'(a)') trim(cline)
          end if

! -- If necessary data is registered for SSF file output.

          if(numssfser.ne.0)then
            do issfser=1,numssfser
              select case(lumpser(issfser))
                case(1)
                  valser(imod,issfser)=vol
                case(2)
                  valser(imod,issfser)=vol_br
                case(3)
                  valser(imod,issfser)=totd
                case(4)
                  valser(imod,issfser)=totm
                case(5)
                  valser(imod,issfser)=deltav
                case(6)
                  valser(imod,issfser)=deltav_br
                case(7)
                  valser(imod,issfser)=(totd-ototd)
                case(8)
                  valser(imod,issfser)=(totm-ototm)
                case(9)
                  valser(imod,issfser)=rainfall
                case(10)
                  valser(imod,issfser)=irrigation
                case(11)
                  valser(imod,issfser)=recharge
                case(12)
                  valser(imod,issfser)=mrecharge
                case(13)
                  valser(imod,issfser)=drainage_br
                case(14)
                  valser(imod,issfser)=overflow_br
                case(15)
                  valser(imod,issfser)=allrecharge
                case(16)
                  valser(imod,issfser)=gwithdrawal
                case(17)
                  valser(imod,issfser)=allrecharge-gwithdrawal
                case(18)
                  valser(imod,issfser)=runoff
                case(19)
                  valser(imod,issfser)=potevapn
                case(20)
                  valser(imod,issfser)=evapn
                case(21)
                  valser(imod,issfser)=potevapn_br
                case(22)
                  valser(imod,issfser)=evapn_br
                case(23)
                  valser(imod,issfser)=gw_pot_evap
                case(24)
                  valser(imod,issfser)=rbal
                case(25)
                  valser(imod,issfser)=elevation
                case(26)
                  valser(imod,issfser)=dwt
              end select
            end do
          end if

C -- Accumulators are updated.

          totrain=totrain+rainfall
          totrech=totrech+recharge
          totmrech=totmrech+mrecharge
          tottrech=tottrech+allrecharge
          totdrainage_br=totdrainage_br+drainage_br
          totoverflow_br=totoverflow_br+overflow_br
          totrun=totrun+runoff
          totevpn=totevpn+evapn
          totevpn_br=totevpn_br+evapn_br
          totbal=totbal+rbal
          totepot=totepot+potevapn
          totepot_br=totepot_br+potevapn_br
          totresidepot=totresidepot+gw_pot_evap
          totirrigation=totirrigation+irrigation
          totgwithdrawal=totgwithdrawal+gwithdrawal
          totnetrecharge=tottrech-totgwithdrawal

          oldvol=vol
          oldvol_br=vol_br
          ototd=totd
          ototm=totm
          if(ifin.eq.1) go to 2000
        go to 1300

C -- Grand totals are recorded.

2000    write(10,*)
        write(10,1800) totrain,totirrigation,totrech,totmrech,
     +  totdrainage_br,totoverflow_br,
     +  tottrech,totgwithdrawal,totnetrecharge,totrun,
     +  totepot,totevpn,totepot_br,totevpn_br,totresidepot,totbal
1800    format(1x,'total',8(1x,14x),16(1x,1pg14.7))
        close(unit=10)
        write(6,2050) trim(outfle)
2050    format(' - file ',a,' written ok.')
        if(csvfle.ne.' ')then
!          write(20,*)
!          write(cline,1801) totrain,totirrigation,totrech,totmrech,   ! We don't do this now as
!     +    totdrainage_br,totoverflow_br,                              ! people forget its there and
!     +    tottrech,totgwithdrawal,totnetrecharge,totrun,              ! it messes up the plotting when
!     +    totepot,totevpn,totepot_br,totevpn_br,totresidepot,totbal   ! imported into EXCEL
!1801      format('total',8(','),16(',',1pg14.7),',,')
!          call remspace(cline)
!          write(20,'(a)') trim(cline)
          close(unit=20)
          write(6,2050) trim(csvfle)
        end if

! -- If necessary the site sample file is written.

        if(numssfser.gt.0)then
          open(unit=20,file=ssffle,action='write',iostat=ierr)
          if(ierr.ne.0)then
            write(errmsg,1249) trim(ssffle)
            go to 9998
          end if
          do issfser=1,numssfser
            asite=ssfser(issfser)
            do i=0,noutdays
              dtemp=
     +        valser(i,issfser)*scaleser(issfser)+offsetser(issfser)
              dtemp=min(dtemp,maxvalser(issfser))
              dtemp=max(dtemp,minvalser(issfser))
              write(atemp20,'(1pg17.10)') dtemp
              atemp20=adjustl(atemp20)
              write(20,1450) trim(asite),trim(dateser(i)),
     +        trim(timeser),trim(atemp20)
1450          format(1x,a,t25,a,t39,a,t55,a)
            end do
          end do
          write(6,2050) trim(ssffle)
          close(unit=20,iostat=ierr)
        end if

        go to 9999

C -- Some error messages.

9000    call writint(aline,iline)
        write(errmsg,9010) trim(aline),trim(infle)
9010    format('Line ',a,' of file ',a,
     +  ': too many entries proposed for subsurface delay buffers: ',
     +  'redimension MAXDELAY and recompile LUMPREM.')
        go to 9998
9020    write(errmsg,9030)
9030    format('All elements of the RBUF and MBUF arrays must be ',
     +  '0.0 or greater.')
        go to 9998
9050    write(errmsg,9060)
9060    format('Cannot allocate sufficient memory for LUMPREM ',
     +  'execution.')
        go to 9998
9150    write(errmsg,9160) trim(cgfle)
9160    format('Cannot open vegetation file ',a,'.')
        go to 9998
9200    write(errmsg,9210) trim(rainfle)
9210    format('Cannot open rainfall file ',a,'.')
        go to 9998
9250    write(errmsg,9260)trim(epotfle)
9260    format('Cannot open evaporation file ',a,'.')
        go to 9998
9300    write(errmsg,9310) trim(irrigfle)
9310    format('Cannot open irrigation file ',a,'.')
        go to 9998
9320    write(errmsg,9330) trim(epotfle_br)
9330    format('Cannot open below-root-zone evaporation file ',a,'.')
        go to 9998
9350    continue
        call writint(aline,iline)
        write(errmsg,9360) trim(aline),trim(infle)
9360    format('Insufficient entries on line ',a,' of file ',a,'.')
        go to 9998
9370    continue
        write(errmsg,9380) trim(infle)
9380    format('Premature end encountered to file ',a,'.')
        go to 9998
9400    continue
        call writint(aline,iline)
        write(errmsg,9410) trim(aline),trim(infle)
9410    format('Error reading data on line ',a,' of file ',a,'.')
        go to 9998
9420    call writint(aline,iline)
        write(errmsg,9430) trim(aline),trim(infle)
9430    format('Error reading initial contents of recharge ',
     +  'delay buffer from line ',a,' of file ',a,'.')
        go to 9998
9440    call writint(aline,iline)
        write(errmsg,9450) trim(aline),trim(infle)
9450    format('Error reading initial contents of macropore ',
     +  'recharge delay buffer from line ',a,' of file ',a,'.')
        go to 9998
9460    call writint(aline,iline)
        write(errmsg,9470) trim(aline),trim(infle)
9470    format('Error reading solution parameters from line ',a,
     +  ' of file ',a,'.')
        go to 9998
9480    call writint(aline,iline)
        write(errmsg,9490) trim(aline),trim(infle)
9490    format('Error reading NUMDAYS or NOUTDAYS from line ',a,
     +  ' of file ',a,'.')
        go to 9998
9495    call writint(aline,iline)
        write(errmsg,9496) trim(aline),trim(infle)
9496    format('Cannot read OUTDAY array. Entries for this array ',
     +  'should begin on line ',a,' of file ',a,'. No comments ',
     +  'are allowed on this line or any subsequent lines ',
     +  'intermixed with array elements. This allows the array ',
     +  'to be recorded on multiple lines.')
        go to 9998
9500    call writint(aline,jline)
        write(errmsg,9510) trim(aline),trim(infle)
9510    format('Error reading name of rainfall file ',
     +  'from line ',a,' of "data filenames" section of ',
     +  'file ',a,'.')
        go to 9998
9520    call writint(aline,jline)
        write(errmsg,9530) trim(aline),trim(infle)
9530    format('Error reading name of potential evaporation file ',
     +  'from line ',a,' of "data filenames" section of ',
     +  'file ',a,'.')
        go to 9998
9550    call writint(aline,iline)
        write(errmsg,9560) trim(atext),trim(aline),trim(infle)
9560    format('Text string "',a,'" expected on or about line ',a,
     +  ' of file ',a,'.')
        go to 9998
9570    call writint(aline,iline)
        write(errmsg,9580) trim(aline),trim(afile)
9580    format('Error reading data on line ',a,' of file ',a,'.')
        go to 9998
9590    call writint(aline,iline)
        write(errmsg,9595) trim(aline),trim(afile)
9595    format('Insufficient entries on line ',a,' of file ',a,'.')
        go to 9998
9600    write(errmsg,9610) trim(infle)
9610    format('"* ssf file" or end of file ',a,' expected.')
        go to 9998
9620    call writint(aline,jline)
        write(errmsg,9630) trim(aline),trim(asection),trim(infle)
9630    format('Error reading data from line ',a,' of "',a,
     +  '" section of file ',a,'.')
        go to 9998
9640    write(errmsg,9650) trim(infle)
9650    format('Unexpected end encountered to file ',a,' while ',
     +  'reading "ssf file" section of this file.')
        go to 9998
9660    call writint(aline,jline)
        write(errmsg,9670) trim(aline),trim(asection),trim(infle)
9670    format('Insufficient entries on line ',a,' of "',a,
     +  '" section of file ',a,'.')
        go to 9998

9998    call writerr(errmsg)

9999    continue
        deallocate(ioutday,vday,cropfac,gamma,rain,epot,stat=ierr)
        deallocate(cropfac_day,gamma_day,stat=ierr)
        deallocate(irrigcode,gwirrigfrac,stat=ierr)
        deallocate(valser,dateser,stat=ierr)
        deallocate(epot_br,stat=ierr)

        end



        subroutine comstrip(cline)

! -- This subroutine removes text that includes and follows the comment symbol.

        implicit none
        character*(*) cline
        integer nn,ll,icount,jcount,i

        if(cline.eq.' ') return
        cline=adjustl(cline)
        nn=index(cline,'#')
        if(nn.eq.0) return
        if(nn.eq.1)then
          cline=' '
          return
        end if

        ll=len_trim(cline)
        icount=0
        jcount=0
        do i=1,ll
          if(cline(i:i).eq.'''') then
            icount=1-icount
          else if(cline(i:i).eq.'"') then
            jcount=1-jcount
          else if(cline(i:i).eq.'#') then
            if((cline(i-1:i-1).eq.' ').or.                 ! Notice that we require the space before #
     +         (cline(i-1:i-1).eq.char(9))) then           ! This allows a filename to have a # in it (mostly)
               if((icount.eq.0).and.(jcount.eq.0)) then
                 cline(i:)=' '
                 return
               end if
            end if
          end if
        end do

        return
        end



        subroutine remspace(cline)

! -- REMSPACE removes spaces from a string (rather inefficiently).

        implicit none
        character*(*), intent(inout) :: cline

        integer :: i,nb,nn
        character*3 :: aa

        aa='   '
        do i=3,1,-1
          do
            nb=len_trim(cline)
            nn=index(cline(1:nb),aa(1:i))
            if(nn.eq.0) exit
            cline(nn:)=cline(nn+i:)
          end do
        end do

        return
        end


        SUBROUTINE WRITERR(ERRMSG)

C -- Subroutine WRITERR formats and writes an error message.

        INTEGER J,JEND,I,NBLC,ITAKE
        CHARACTER*(*) ERRMSG

        write(6,*)
        ITAKE=0
        J=0
        NBLC=len_trim(ERRMSG)
5       JEND=J+78-ITAKE
        IF(JEND.GE.NBLC) GO TO 100
        DO 10 I=JEND,J+1,-1
        IF(ERRMSG(I:I).EQ.' ') THEN
          IF(ITAKE.EQ.0) THEN
            WRITE(6,'(A)',ERR=200) ' '//ERRMSG(J+1:I)
            ITAKE=4
          ELSE
            WRITE(6,'(A)',ERR=200) '     '//ERRMSG(J+1:I)
          END IF
          J=I
          GO TO 5
        END IF
10      CONTINUE
        IF(ITAKE.EQ.0)THEN
          WRITE(6,'(A)',ERR=200) ' '//ERRMSG(J+1:JEND)
          ITAKE=4
        ELSE
          WRITE(6,'(A)',ERR=200) '     '//ERRMSG(J+1:JEND)
        END IF
        J=JEND
        GO TO 5
100     JEND=NBLC
        IF(ITAKE.EQ.0)THEN
          WRITE(6,'(A)',ERR=200) ' '//ERRMSG(J+1:JEND)
        ELSE
          WRITE(6,'(A)',ERR=200) '     '//ERRMSG(J+1:JEND)
        END IF
        RETURN

200     STOP
        END



        subroutine writint(atemp,ival)

c       Subroutine WRITINT writes an integer to a character variable.

        integer*4 ival
        character*6 afmt
        character*(*) atemp

        afmt='(i   )'
        write(afmt(3:5),'(i3)') len(atemp)
        write(atemp,afmt)ival
        atemp=adjustl(atemp)
        return
        end



        SUBROUTINE DREALREAD(IFAIL,CLINE,RTEMP)

! -- Subroutine REALREAD reads a real number from a string.

        INTEGER IFAIL
        double precision RTEMP
        CHARACTER*8 AFMT
        CHARACTER*(*) CLINE

        IFAIL=0
        AFMT='(F   .0)'
        WRITE(AFMT(3:5),'(I3)') len_trim(CLINE)
        READ(CLINE,AFMT,ERR=100) RTEMP
        RETURN

100     IFAIL=1
        RETURN
        END



        SUBROUTINE INTREAD(IFAIL,CLINE,iTEMP)

C -- Subroutine INTREAD reads a real number from a string.

        INTEGER IFAIL
        integer iTEMP
        CHARACTER*6 AFMT
        CHARACTER*(*) CLINE

        IFAIL=0
        AFMT='(i   )'
        WRITE(AFMT(3:5),'(I3)') LEN(CLINE)
        READ(CLINE,AFMT,ERR=100) iTEMP
        RETURN

100     IFAIL=1
        RETURN
        END




         SUBROUTINE LINESPLIT(IFAIL,NUM,LW,RW,CLINE)

C -- Subroutine LINESPLIT splits a string into blank-delimited fragments.

        INTEGER IFAIL,NW,NBLC,J,I
        INTEGER NUM,NBLNK
        INTEGER LW(NUM),RW(NUM)
        CHARACTER*(*) CLINE

        IFAIL=0
        NW=0
        NBLC=LEN_TRIM(CLINE)
        IF((NBLC.NE.0).AND.(INDEX(CLINE,CHAR(9)).NE.0)) THEN
          CALL TABREM(CLINE)
          NBLC=LEN_TRIM(CLINE)
        ENDIF
        IF(NBLC.EQ.0) THEN
          IFAIL=-1
          RETURN
        END IF
        J=0
5       IF(NW.EQ.NUM) RETURN
        DO 10 I=J+1,NBLC
          IF((CLINE(I:I).NE.' ').AND.(CLINE(I:I).NE.',').AND.
     +    (ICHAR(CLINE(I:I)).NE.9)) GO TO 20
10      CONTINUE
        IFAIL=1
        RETURN
20      NW=NW+1
        LW(NW)=I
        DO 30 I=LW(NW)+1,NBLC
          IF((CLINE(I:I).EQ.' ').OR.(CLINE(I:I).EQ.',').OR.
     +    (ICHAR(CLINE(I:I)).EQ.9)) GO TO 40
30      CONTINUE
        RW(NW)=NBLC
        IF(NW.LT.NUM) IFAIL=1
        RETURN
40      RW(NW)=I-1
        J=RW(NW)
        GO TO 5

        END


        SUBROUTINE TABREM(CLINE)

C -- Subroutine TABREM removes tabs from a string.

        INTEGER I
        CHARACTER*(*) CLINE

        DO 10 I=1,LEN(CLINE)
10      IF(ICHAR(CLINE(I:I)).EQ.9) CLINE(I:I)=' '

        RETURN
        END


        SUBROUTINE LOWCASE(ASTRNG)

! -- Subroutine LOWCASE converts a string to lower case.

        INTEGER I,J
        CHARACTER*(*) ASTRNG

        DO 10 I=1,len_trim(ASTRNG)
        J=ICHAR(ASTRNG(I:I))
        IF((J.GE.65).AND.(J.LE.90)) ASTRNG(I:I)=CHAR(J+32)
10      CONTINUE
        RETURN
        END



        subroutine getdate2(ifail,adate,dd,mm,yy,datespec)

! -- Subroutine getdate2 reads a date from a string.

        logical leap
        integer ifail,dd,mm,yy,n,datespec,itemp
        character*(*) adate
        character*(15) tdate

        ifail=0
        adate=adjustl(adate)
        tdate=adate(1:min(15,len(adate)))
        tdate=adjustl(tdate)
        n=index(tdate,'/')
        if(n.eq.0) go to 9999
        if(n.eq.1) go to 9999
        if(n.eq.2)then
          read(tdate(1:1),'(i1)',err=9999) itemp
        else if(n.eq.3)then
          read(tdate(1:2),'(i2)',err=9999) itemp
        else
          go to 9999
        end if
        if(datespec.eq.1)then
          dd=itemp
        else
          mm=itemp
        end if
        tdate=tdate(n+1:)
        n=index(tdate,'/')
        if(n.eq.0) go to 9999
        if(n.eq.1) go to 9999
        if(n.eq.2) then
          read(tdate(1:1),'(i1)',err=9999) itemp
        else if(n.eq.3) then
          read(tdate(1:2),'(i2)',err=9999) itemp
        else
          go to 9999
        end if
        if(datespec.eq.1)then
          mm=itemp
        else
          dd=itemp
        end if
        tdate=tdate(n+1:)
        n=len_trim(tdate)
        if(n.ne.4) go to 9999
        read(tdate(1:4),'(i4)',err=9999) yy

        if((mm.lt.1).or.(mm.gt.12)) go to 9999
        if((mm.eq.1).or.(mm.eq.3).or.(mm.eq.5).or.
     +  (mm.eq.7).or.(mm.eq.8).or.(mm.eq.10).or.(mm.eq.12))then
          if(dd.gt.31) go to 9999
        else if((mm.eq.4).or.(mm.eq.6).or.(mm.eq.9).or.
     +  (mm.eq.11))then
          if(dd.gt.30) go to 9999
        else
          if(leap(yy)) then
            if(dd.gt.29) go to 9999
          else
            if(dd.gt.28) go to 9999
          end if
        end if

        return

9999    ifail=1
        return
        end


        logical function leap(year)

! -- Function LEAP returns .true. if a year is a leap year.

        integer year

        leap = ( mod(year,4).eq.0 .and. mod(year,100).ne.0 ) .or.
     +  ( mod(year,400).eq.0 .and. year.ne.0 )

        return
        end



	subroutine gettime(ifail,atime,hh,mm,ss)

! -- Subroutine gettime reads a time; note that it does not destroy the time string.

	integer ifail,hh,mm,ss,n,i
	character*15 tatime
	character*(*) atime

	if(atime.eq.' ') go to 9999
	do i=1,len(atime)
	  if(atime(i:i).ne.' ') go to 10
	end do
10	tatime=atime(i:)

	ifail=0
	tatime=adjustl(tatime)
	n=index(tatime,':')
	if(n.eq.0) then
	  n=index(tatime,'.')
	  if(n.eq.0) go to 9999
	end if
	if(n.eq.1) go to 9999
	if(n.eq.2)then
	  read(tatime(1:1),'(i1)',err=9999) hh
	else if(n.eq.3)then
	  read(tatime(1:2),'(i2)',err=9999) hh
	else
	  go to 9999
	end if
	tatime=tatime(n+1:)
	n=index(tatime,':')
	if(n.eq.0)then
	  n=index(tatime,'.')
	  if(n.eq.0) go to 9999
	end if
	if(n.eq.0) go to 9999
	if(n.eq.1) go to 9999
	if(n.eq.2) then
	  read(tatime(1:1),'(i1)',err=9999) mm
	else if(n.eq.3) then
	  read(tatime(1:2),'(i2)',err=9999) mm
	else
	  go to 9999
	end if
	tatime=tatime(n+1:)
	n=len_trim(tatime)
	if(n.gt.2) go to 9999
	if(n.eq.1)then
	  read(tatime(1:n),'(i1)',err=9999) ss
	else
	  read(tatime(1:n),'(i2)',err=9999) ss
	end if

	if((mm.lt.0).or.(mm.gt.59)) go to 9999
	if((hh.lt.0).or.(hh.gt.23)) go to 9999
	if((ss.lt.0).or.(ss.gt.59)) go to 9999

	return

9999	ifail=1
	return
	end



        subroutine getfile(ifail,cline,filename,ibeg,iend)

!-- Subroutine getfile extracts a filename from a string.

! -- Arguments are as follows:-
!       ifail: returned as zero if filename successfully read
!       cline: a character string containing the file name
!       filename: the name of the file read from the string
!       ibeg: character position at which to begin search for filename
!       iend: on input  - character position at which to end search for filename
!             on output - character postion at which filename ends


        integer, intent(out)               :: ifail
        integer, intent(in)                :: ibeg
        integer, intent(inout)             :: iend
        character (len=*), intent(in)      :: cline
        character (len=*), intent(out)     :: filename

        integer                            :: i,j,k
        character (len=1)                  :: aa

        ifail=0
        do i=ibeg,iend
          aa=cline(i:i)
          if((aa.ne.' ').and.(aa.ne.',').and.(aa.ne.char(9)))go to 50
        end do
        ifail=1
        return

50      if((aa.eq.'"').or.(aa.eq.''''))then
          do j=i+1,iend
            if(cline(j:j).eq.aa) go to 60
          end do
          ifail=1
          return
60        iend=j
          if(i+1.gt.j-1)then
            ifail=1
            return
          else
            filename=cline(i+1:j-1)
          end if
        else
          do j=i+1,iend
            if((cline(j:j).eq.' ').or.(cline(j:j).eq.',')
     +     .or.(cline(j:j).eq.char(9)))then
              k=j-1
              go to 100
            end if
          end do
          k=iend
100       filename=cline(i:k)
          if(cline(k:k).eq.'"')then
            ifail=1
            return
          else if(cline(k:k).eq.'''')then
            ifail=1
            return
          end if

          iend=k
        end if
        filename=adjustl(filename)
        return

        end



        subroutine newdate(ndays,day1,mon1,year1,day2,mon2,year2)

! -- Subroutine NEWDATE evaluates the date after NDAYS days have elapsed from
!    a provided date. NDAYS may be negative.

! -- Arguments are as follows:-
!       ndays:            elapsed number of days
!       day1,mon1,year1:  days, month and year of first date
!       day2,mon2,year2:  days, month and year of second date


        implicit none

        integer, intent(in)     :: ndays,day1,mon1,year1
        integer, intent(out)    :: day2,mon2,year2

        integer  :: yearref,newdays,idays,iyear,jdays,i
        integer, dimension(12) :: monthdays

        integer numdays
        logical leap

        data monthdays /31,28,31,30,31,30,31,31,30,31,30,31/

! -- First a reference date is chosen. This is the beginning of the first
! -- year. Alternatively the reference date is the beginning of a year prior
! -- to the likely calculated date if NDAYS is negative.

        if(ndays.ge.0) then
          yearref=year1
        else
          yearref=year1-abs(ndays)/365-1
        end if
        newdays=numdays(31,12,yearref-1,day1,mon1,year1)
        newdays=ndays+newdays
        if(newdays.lt.0) call sub_error('NEWDATE')

! -- Next days are counted, starting at the new reference date.

        idays=0
        iyear=yearref
        do
          jdays=idays+365
          if(leap(iyear)) jdays=jdays+1
          if(jdays.ge.newdays) go to 20
          iyear=iyear+1
          idays=jdays
        end do
        call sub_error('NEWDATE')
20      year2=iyear

        do i=1,12
          jdays=idays+monthdays(i)
          if((i.eq.2).and.(leap(year2))) jdays=jdays+1
          if(jdays.ge.newdays) go to 40
          idays=jdays
        end do
        call sub_error('NEWDATE')
40      mon2=i
        day2=newdays-idays
!       if((day2.le.0).or.(mon2.le.0).or.(year2.le.0)) call sub_error('NEWDATE')
        if((day2.le.0).or.(mon2.le.0)) call sub_error('NEWDATE')
        return

        end




        integer function numdays(DR,MR,YR,D,M,Y)

! -- Function numdays calculates the number of days between dates
!    D-M-Y and DR-MR-YR. If the former preceeds the latter the answer is
!    negative.

! -- Arguments are as follows:-
!       dr,mr,yr:     days, months and years of first date
!       d,m,y:        days, months and years of second date
!       numdays returns the number of elapsed days


        integer, intent(in)     :: dr,mr,yr,d,m,y

        INTEGER FLAG,I,J,DA(12),YE,ME,DE,YL,ML,DL
        logical leap

        DATA DA /31,28,31,30,31,30,31,31,30,31,30,31/

! --    THE SMALLER OF THE TWO DATES IS NOW CHOSEN TO DO THE COUNTING FROM.

        IF(Y.LT.YR)GO TO 10
        IF((Y.EQ.YR).AND.(M.LT.MR)) GO TO 10
        IF((Y.EQ.YR).AND.(M.EQ.MR).AND.(D.LT.DR)) GO TO 10
        FLAG=0
        YE=YR
        ME=MR
        DE=DR
        YL=Y
        ML=M
        DL=D
        GO TO 20
10      FLAG=1
        YE=Y
        ME=M
        DE=D
        YL=YR
        ML=MR
        DL=DR

! --    IN THE ABOVE THE POSTSCRIPT "E" STANDS FOR EARLIER DATE, WHILE
!       "L" STANDS FOR THE LATER DATE.

20      numdays=0
        IF((ME.EQ.ML).AND.(YL.EQ.YE))THEN
        numdays=DL-DE
        IF(FLAG.EQ.1) numdays=-numdays
        RETURN
        END IF

        DO 30 J=ME,12
        IF((ML.EQ.J).AND.(YE.EQ.YL))GOTO 40
        numdays=numdays+DA(J)
        IF((J.EQ.2).AND.(leap(ye)))numdays=numdays+1
30      CONTINUE
        GO TO 50
40      numdays=numdays+DL-DE
        IF(FLAG.EQ.1)numdays=-numdays
        RETURN

50      DO 60 I=YE+1,YL
        DO 70 J=1,12
        IF((YL.EQ.I).AND.(ML.EQ.J))GO TO 80
        numdays=numdays+DA(J)
        IF((J.EQ.2).AND.(leap(i))) numdays=numdays+1
70      CONTINUE
60      CONTINUE
        call sub_error('NUMDAYS')
        RETURN

80      numdays=numdays+DL-DE
        IF(FLAG.EQ.1) numdays=-numdays

        RETURN
        end


        subroutine sub_error(subname)

! -- Subroutine sub_error names the subroutine causing a run-time error.

! -- Arguments are as follows:-
!       subname:  name of offending subroutine


        character (len=*)               ::subname

        write(6,10) trim(subname)
10      format(/,' *** PROGRAMMING ERROR CALLING SUBROUTINE ',a,' ***')
        stop

       end

