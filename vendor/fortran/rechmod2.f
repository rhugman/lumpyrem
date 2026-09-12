        subroutine rechmod(icall,
     +  mxiter,tol,ndays,subdim,nstep,maxvol,
     +  irrigvolfrac,cropfac_day,gamma_day,ks,l,m,mflowmax,rdelay,
     +  mdelay,vol,drainsub,macsub,rain,epot,recharge,mrecharge,
     +  runoff,evapn,rainfall,potevapn,irrigation,gwithdrawal,
     +  irrigcode,gwirrigfrac,
     +  maxvol_br,extravol_br,ks_br,l_br,m_br,vol_br,drainage_br,
     +  overflow_br,epot_br,gamma_br,evapn_br,potevapn_br)

! -- The recharge to groundwater is RECHARGE if MAXVOL_BR is zero.
!    Otherwise it is RECHARGE_BR.
!    RECHARGE is inherited from the older version of this subroutine
!    where there was no lower bucket.

        implicit none

C -- Subroutine arguments are first declared.

        integer icall,ndays,subdim,nstep,mxiter
        integer irrigcode(ndays)
        double precision maxvol,irrigvolfrac,ks,l,m,
     +  mflowmax,vol,recharge,mrecharge,runoff,rdelay,mdelay,
     +  evapn,tol,drainage,evap,rainfall,potevapn
        double precision maxvol_br,ks_br,l_br,m_br,vol_br,drainage_br,
     +  overflow_br,extravol_br
        double precision daydrainage_br,dayoverflow_br,recharge_for_br,
     +  tvol_br,otvol_br,totvol_br
        double precision irrigation,gwithdrawal,vvol
        double precision rain(ndays),epot(ndays),drainsub(subdim)
        double precision macsub(subdim)
        double precision cropfac_day(ndays),gamma_day(ndays)
        double precision gwirrigfrac(ndays)
        double precision epot_br(ndays),gamma_br,evapn_br,potevapn_br

C -- Other variables are declared.

        integer irdelay,imdelay,i,iday,istep,iter
        integer ircode
        integer iflag_br
        double precision train,tstep,tepot,cropfac,dayrech,daymrech,
     +  dayevap,dayevap_br,
     +  dayrunoff,tvol,vd,otvol,rtemp1,rtemp2,tempvol,frdelay,
     +  fmdelay,gamma
        double precision gwfrac,dayirrigation,daygwithdrawal,dirrig,
     +  odirrig
        double precision cropfac_br,tepot_br,rtemp3

C -- Initialization

        vvol=irrigvolfrac*maxvol

C -- Variables accumulated over a single RECHMOD call are zeroed.

        rainfall=0.0
        recharge=0.0        ! drainage from upper bucket
        drainage_br=0.0     ! drainage from lower bucket
        overflow_br=0.0
        mrecharge=0.0
        runoff=0.0
        evapn=0.0
        potevapn=0.0
        evapn_br=0.0
        potevapn_br=0.0
        irrigation=0.0
        gwithdrawal=0.0
        tstep=1.0/dble(nstep)
        cropfac_br=1.0      ! hardwired
        recharge_for_br=0.0d0
        iflag_br=1

C -- Any water in the recharge or macropore flow recharge buffers that
C    should already have been assigned to recharge is now claimed.

        irdelay=int(rdelay)
        frdelay=rdelay-irdelay
        irdelay=irdelay+1
        imdelay=int(mdelay)
        fmdelay=mdelay-imdelay
        imdelay=imdelay+1
        if(irdelay+1.le.subdim) then
          do 100 i=irdelay+1,subdim
            recharge=recharge+drainsub(i)
            drainsub(i)=0.0
100       continue
        end if
        if(maxvol_br.gt.0.0d0) recharge_for_br=recharge
        if(imdelay+1.le.subdim) then
          do 150 i=imdelay+1,subdim
            mrecharge=mrecharge+macsub(i)
            macsub(i)=0.0
150       continue
        end if

C -- Cycle through the days.

        do 1000 iday=1,ndays

C -- First the amount of rainfall and potential evaporation per timestep
C    is established.

          train=rain(iday)*tstep
          tepot=epot(iday)
          tepot_br=epot_br(iday)
          cropfac=cropfac_day(iday)
          gamma=gamma_day(iday)
          ircode=irrigcode(iday)
          gwfrac=gwirrigfrac(iday)

C -- Next daily accumulators are zeroed.

          dayrech=0.0
          daydrainage_br=0.0
          dayoverflow_br=0.0
          daymrech=0.0
          dayevap=0.0
          dayevap_br=0.0
          dayrunoff=0.0
          dayirrigation=0.0
          daygwithdrawal=0.0

          do 700 istep=1,nstep
            tvol=vol
            otvol=vol
            odirrig=0.0d0
            do 400 iter=1,mxiter
c            vd=tvol/maxvol                        !fully implicit
              vd=(tvol+vol)*0.5/maxvol             !Crank-Nicholson
              if(vd.gt.1.0) then
                vd=1.0
              else if(vd.lt.0.0)then
                vd=0.0
              end if
              rtemp1=drainage(vd,ks,m,l)*tstep
              rtemp2=evap(vd,tepot,cropfac,gamma)*tstep
              tvol=vol-rtemp1-rtemp2+train
              if(tvol.lt.0.0) tvol=0.0
              dirrig=0.0
              if(ircode.ne.0)then
                if(tvol.lt.vvol)then
                  dirrig=vvol-tvol
                  tvol=vvol
                end if
              end if
              tempvol=tvol
              if(tvol.gt.maxvol) tvol=maxvol
              if(iter.eq.1) go to 390
              if(ircode.eq.0)then
                if(abs(tvol-otvol).le.tol*maxvol) go to 600
              else
                if(abs(tvol-otvol).le.tol*maxvol) then
                if(abs(dirrig-odirrig).le.tol*(dirrig+odirrig)*0.5)
     +            go to 600
                end if
              end if
390           otvol=tvol
              odirrig=dirrig
400         continue
            write(6,510) icall,iday,istep
510         format(' RECHMOD upper bucket: call',i6,' day',i6,
     +      ' step',i6,': itn limit exceeded')
C           stop
600         continue
!            write(6,610) icall,iday,istep,iter                      !debug
!610         format(' Main  bucket call',i4,' day',i4,' step',i4,    !debug
!     +      ':  iterations = ',i4)                                  !debug

            if(tempvol.eq.0.0d0)then
              rtemp3=rtemp1+rtemp2
              if(vol.lt.rtemp3)then
                if(vol.eq.0.0d0)then
                  rtemp1=0.0d0
                  rtemp2=0.0d0
                else
                  if(rtemp3.gt.0.0d0)then
                    rtemp2=rtemp2*vol/rtemp3
                    rtemp1=rtemp1*vol/rtemp3
                  end if
                end if
              end if
            end if

            dayrech=dayrech+rtemp1
            dayevap=dayevap+rtemp2
            dayirrigation=dayirrigation+dirrig
            daygwithdrawal=daygwithdrawal+gwfrac*dirrig
            if(tempvol.gt.maxvol) then
              rtemp1=tempvol-maxvol
              if(rtemp1.lt.mflowmax*tstep) then
                daymrech=daymrech+rtemp1
              else
                daymrech=daymrech+mflowmax*tstep
                dayrunoff=dayrunoff+rtemp1-mflowmax*tstep
              end if
              tempvol=maxvol
            end if
            vol=tempvol
700       continue

          runoff=runoff+dayrunoff
          evapn=evapn+dayevap
          irrigation=irrigation+dayirrigation
          gwithdrawal=gwithdrawal+daygwithdrawal

C -- Next the recharge and macropore recharge delay buffers are shifted,
C    emptied from the bottom, and filled from the top.

          recharge=recharge+drainsub(irdelay)
          mrecharge=mrecharge+macsub(imdelay)
          if(iflag_br.eq.1)then
            recharge_for_br=drainsub(irdelay)
          else
            recharge_for_br=recharge_for_br+drainsub(irdelay)
            iflag_br=0
          end if
          if(irdelay.gt.1) then
            do 800 i=irdelay,2,-1
              drainsub(i)=drainsub(i-1)
800         continue
          end if
          if(imdelay.gt.1) then
            do 850 i=imdelay,2,-1
              macsub(i)=macsub(i-1)
850         continue
          end if
          drainsub(1)=dayrech
          macsub(1)=daymrech
          recharge=recharge+(1.0-frdelay)*drainsub(irdelay)
          recharge_for_br=
     +    recharge_for_br+(1.0-frdelay)*drainsub(irdelay)
          drainsub(irdelay)=frdelay*drainsub(irdelay)
          mrecharge=mrecharge+(1.0-fmdelay)*macsub(imdelay)
          macsub(imdelay)=fmdelay*macsub(imdelay)
          rainfall=rainfall+rain(iday)
          potevapn=potevapn+epot(iday)

! -- We now, optionally, route daily recharge through the below-root-zone bucket.

          if(maxvol_br.gt.0.0d0)then
            totvol_br=maxvol_br+extravol_br
            do 1700 istep=1,nstep
              tvol_br=vol_br
              otvol_br=vol_br
              do 1400 iter=1,mxiter
c               vd=tvol_br/maxvol_br                        !fully implicit
                vd=(tvol_br+vol_br)*0.5/maxvol_br           !Crank-Nicholson
                if(vd.gt.1.0) then
                  vd=1.0
                else if(vd.lt.0.0)then
                  vd=0.0
                end if
                rtemp1=drainage(vd,ks_br,m_br,l_br)*tstep
                rtemp2=evap(vd,tepot_br,cropfac_br,gamma_br)*tstep
                tvol_br=vol_br-rtemp1-rtemp2+recharge_for_br*tstep
                if(tvol_br.lt.0.0) tvol_br=0.0
                tempvol=tvol_br
                if(tvol_br.gt.totvol_br) tvol_br=totvol_br     ! Is this necessary?
                if(iter.eq.1) go to 1390
                if(abs(tvol_br-otvol_br).le.tol*maxvol_br) go to 1600
1390            otvol_br=tvol_br
1400          continue
              write(6,1510) icall,iday,istep
1510          format(' RECHMOD lower bucket: call',i6,' day',i6,
     +        ' step',i6,': itn limit exceeded')
C             stop
1600          continue
!              write(6,1610) icall,iday,istep,iter                     !debug
!1610          format(' Lower bucket call',i4,' day',i4,' step',i4,    !debug
!     +        ':  iterations = ',i4)                                  !debug

! -- We see if there is a mass balance error because of too much extraction.

              if(tempvol.eq.0.0d0)then
                rtemp3=rtemp1+rtemp2
                if(vol_br.lt.rtemp3)then
                  if(vol_br.eq.0.0d0)then
                    rtemp1=0.0d0
                    rtemp2=0.0d0
                  else
                    if(rtemp3.gt.0.0d0)then
                      rtemp2=rtemp2*vol_br/rtemp3
                      rtemp1=rtemp1*vol_br/rtemp3
                    end if
                  end if
                end if
              end if
              daydrainage_br=daydrainage_br+rtemp1
              dayevap_br=dayevap_br+rtemp2
              if(tempvol.gt.totvol_br) then
                rtemp1=tempvol-totvol_br
                dayoverflow_br=dayoverflow_br+rtemp1
                tempvol=totvol_br
              end if
              vol_br=tempvol
1700        continue
            drainage_br=drainage_br+daydrainage_br
            overflow_br=overflow_br+dayoverflow_br
            evapn_br=evapn_br+dayevap_br
            potevapn_br=potevapn_br+epot_br(iday)
          end if

1000    continue

        return
        end



      double precision function drainage(vd,ks,m,l)

C -- Function DRAINAGE calculates drainage as a function of the amount
C    of water in the upper moisture store.

      implicit none

      double precision vd,ks,m,l,rtemp

      if(vd.le.0.0)then
        drainage=0.0
      else if(vd.ge.1.0)then
        drainage=ks
      else
        rtemp=1.0-vd**(1.0/m)
        rtemp=rtemp**m
        rtemp=(1.0-rtemp)
        drainage=ks*(vd**l)*rtemp*rtemp
      end if

      return
      end



      double precision function evap(vd,epot,cropfac,gamma)

C -- Subroutine evap calculates evaporation rate as a function amount of
C    water in the upper moisture store.

      implicit none

      double precision vd,epot,cropfac,gamma,rtemp

      if(vd.le.0.0)then
        evap=0.0
      else if(vd.ge.1.0)then
        evap=cropfac*epot
      else
        rtemp=exp(-gamma*vd)
        evap=cropfac*epot*(1.0-rtemp)/(1.0-2.0*exp(-gamma)+rtemp)
      end if
      return
      end

