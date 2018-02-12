       subroutine NMAmodel(y1, sd1, npt1, x1, ids1, iarm1, igroup1,
     1    nrep, nthin, nwarm, iseed, conf,
     2    DICinfo, alogcpo, betapo, gampo, tau2po, eRho, seRho,
     3    Rholow, Rhoupp)
c      Bayesian NMA model
c      Hao Li, UConn
c      March 28,2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       parameter (tol=2.220446049250313E-14)
       parameter (epi= 3.141592741012573)
       double precision y1(ns),sd1(ns),x1(ns,nx), conf
       integer ids1(ns),iarm1(ns), npt1(ns), igroup1(nT)
       integer nrep, nthin, nwarm, iseed
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision pOmega(nT,nT), pRho(nT,nT)
       double precision alam(K),Rgam(ns),tau2(nT), Rho(nT,nT)
       double precision xmean(nx),xsd(nx),x2mean(nx)
       double precision alow(2),aupp(2)
       double precision betapo(nx,4), gampo(nT,4)
       double precision tau2po(nT,4)
c
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns), igroup(nT)
       integer icount1, idim
       double precision seqbeta(nx,nrep),seqgam(nT,nrep)
       double precision seqOmega(nrep,nT,nT), seqsig2(ns,nrep)
       double precision seqRho(nrep,nT,nT),seqRgam(ns,nrep)
       double precision ahpd(nrep)
c......Beta
       double precision sbeta(nx),s2beta(nx)
       double precision ebeta(nx),sebeta(nx)
c......Gamma Fixed
       double precision sgam(nT),s2gam(nT)
       double precision egam(nT),segam(nT)
c......sig2
       double precision ssig2(ns),esig2(ns)
c       double precision s2sig2(ns),sesig2(ns),sig2low(ns),sig2upp(ns)
c......Omega
       double precision sOmega(nT,nT),s2Omega(nT,nT)
       double precision eOmega(nT,nT),seOmega(nT,nT),Omegalow(nT,nT)
       double precision Omegaupp(nT,nT)
c......Rho
       double precision sRho(nT,nT),s2Rho(nT,nT),Rholow(nT,nT)
       double precision eRho(nT,nT),seRho(nT,nT),Rhoupp(nT,nT)
c......barDIC and DICbar
       double precision yk(4), xx(4), ywk(4), EI(nT,4), B2(4,nT)
       double precision B3(4,4), SigI(4,4)
       double precision sum1, barDIC, DICbar, DICinfo(5)
       double precision ein, detf
       double precision CI2(2,2), CIinv2(2,2),CIchol2(2,2)
       double precision CI3(3,3), CIinv3(3,3),CIchol3(3,3)
       double precision CI4(4,4), CIinv4(4,4),CIchol4(4,4)
c......CPO and LPML
       double precision g(K,nrep), alogcpo(K), alpml, gmax(K), sumrep
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecRgam/Rgam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecigroup/igroup       
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE
c......common variable
       y=y1
       sd=sd1
       npt=npt1
       ids=ids1
       iarm=iarm1
       igroup=igroup1  
       do i1=1,ns
        do j1=1,nx
         x(j1,i1)=x1(i1,j1)
        enddo
       enddo    
c......calculate # of arms for each trial
         do kk=1,K
          narms(kk)=0
         enddo
         do i=1,ns
          iarm(i)=iarm(i)+1
          narms(ids(i))=narms(ids(i))+1
         enddo
c......Standardize X
         do j=1,nx
          xmean(j)=0.0d0
          x2mean(j)=0.0d0
          do i=1,ns
           xmean(j)=xmean(j)+x(j,i)
           x2mean(j)=x2mean(j)+x(j,i)**2
          enddo
          xmean(j)=xmean(j)/real(ns)
          xsd(j)=dsqrt( (x2mean(j)-
     1     real(ns)*xmean(j)**2)/real(ns-1) )
         enddo
         do j=1,nx
           do i=1,ns
             x(j,i)=(x(j,i)-xmean(j))/xsd(j)
           enddo
         enddo
c
         do kk=1,K
           do j1=1,nT
            do j2=1,nT
             randE(kk,j1,j2)=0.0d0
            enddo
           enddo
           do i=1,ns
             if (ids(i) .eq. kk) then
               do j=1,nT
               if (iarm(i) .eq. j) randE(kk,j,j)=1.0d0
               enddo
             endif
           enddo
         enddo
c......set initial values
        do j1=1,nT
         do j2=1,nT
          pOmega(j1,j2)=0.0d0
         enddo
         pOmega(j1,j1)=2.67d0
        enddo
       do j=1,nx
        beta(j)=0.0d0
        sbeta(j)=0.0d0
        s2beta(j)=0.0d0
        ebeta(j)=0.0d0
        sebeta(j)=0.0d0
        do j1=1,4
         betapo(j,j1)=0.0d0
        enddo
       enddo
       do jj=1,ns
         Rgam(jj)=0.0d0
       enddo
       do i=1,ns
        sig2(i)=sd(i)**2
        ssig2(i)=0.0d0
c        s2sig2(i)=0.0d0
        esig2(i)=0.0d0
c        sesig2(i)=0.0d0
c        sig2low(i)=0.0d0
c        sig2upp(i)=0.0d0
       enddo
       do j1=1,nT
        gam(j1)=0.0d0
        sgam(j1)=0.0d0
        s2gam(j1)=0.0d0
        egam(j1)=0.0d0
        segam(j1)=0.0d0
        do j=1,4
         gampo(j1,j)=0.0d0
        enddo
        do j2=1,nT
         Omega(j1,j2)=pOmega(j1,j2)
         pRho(j1,j2)=0.0d0
         sOmega(j1,j2)=0.0d0
         s2Omega(j1,j2)=0.0d0
         eOmega(j1,j2)=0.0d0
         seOmega(j1,j2)=0.0d0
         Omegalow(j1,j2)=0.0d0
         Omegaupp(j1,j2)=0.0d0
         sRho(j1,j2)=0.0d0
         s2Rho(j1,j2)=0.0d0
         eRho(j1,j2)=0.0d0
         seRho(j1,j2)=0.0d0   
        enddo
       enddo
       beta(1)=0.05113
       beta(2)=-1.38866
       beta(3)=1.09817
       beta(4)=-0.85855
       beta(5)=-1.12056
       beta(6)=-1.14133
       beta(7)=-0.22435
       beta(8)=3.63453
       beta(9)=-2.09322
       beta(10)=1.07858  
       gam(1)=0.80566
       gam(2)=-40.76753
       gam(3)=-45.07127
       gam(4)=-28.27232
       gam(5)=-44.14054
       gam(6)=-28.13203
       gam(7)=-19.19989
       gam(8)=-47.21824
       gam(9)=-51.31234
       gam(10)=-48.46266
       gam(11)=-47.71443
c
       do kk=1,K
        alam(kk)=1.0d0
       enddo
c
       do i1=1,nT
        tau2(i1)=Omega(i1,i1)
        do i2=1,4
        tau2po(i1,i2)=0.0d0
        enddo
       enddo
       do j1=1,nT
        do j2=1,nT
         Rho(j1,j2)= Omega(j1,j2)/dsqrt(Omega(j1,j1)*Omega(j2,j2))
        enddo
       enddo
c......set initial values for partial Rho
        do j=1,nT
         pRho(j,j)=1.0d0
        enddo
c......warm up Gibbs
c       call rnset( iseed )
       do i1=1,nwarm
        call gibbs( iseed )
       enddo
c       call rnget( iseed )
c......initial values
       do 213 i=1,nrep
        do j1=1,nx
         seqbeta(j1,i)=0.0d0
        enddo
        do j2=1,nT
         seqgam(j2,i)=0.0d0
        enddo
        do j3=1,ns
         seqRgam(j3,i)=0.0d0
        enddo
        do j4=1,nT
         do j5=1,nT
          seqOmega(i,j4,j5)=0.0d0
          seqRho(i,j4,j5)=0.0d0
         enddo
        enddo
        do j6=1,ns
         seqsig2(j6,i)=0.0d0
        enddo
 213   continue
c
       barDIC=0.0d0
       DICbar=0.0d0
       do kk=1,K
        do jj=1,nrep
         g(kk,jj)=0.0d0
        enddo
        alogcpo(kk)=0.0d0
       enddo
       alpml=0.0d0
c......write gibbs samplings
c       call rnset( iseed )
       do 220 i=1,nrep
        do ithin = 1,nthin
         call gibbs( iseed )
        enddo
        do 230 j1=1,nx
         seqbeta(j1,i)=beta(j1)
         sbeta(j1)=sbeta(j1)+beta(j1)
         s2beta(j1)=s2beta(j1)+beta(j1)**2
 230    continue
        do 231 j2=1,nT
         seqgam(j2,i)=gam(j2)
         sgam(j2)=sgam(j2)+gam(j2)
         s2gam(j2)=s2gam(j2)+gam(j2)**2
 231    continue
c
        do j3=1,ns
         seqRgam(j3,i)=Rgam(j3)
c          sRgam(kk,j3)=sRgam(kk,j3)+Rgam(kk,j3)
c          s2Rgam(kk,j3)=s2Rgam(kk,j3)+Rgam(kk,j3)**2
        enddo
c
        do 232 j6=1,ns
c         seqsig2(j6,i)=sig2(j6)
         ssig2(j6)=ssig2(j6)+sig2(j6)
c         s2sig2(j6)=s2sig2(j6)+sig2(j6)**2
 232    continue
        do j4=1,nT
         do j5=1,nT
          seqOmega(i,j4,j5)=Omega(j4,j5)
          sOmega(j4,j5)=sOmega(j4,j5)+Omega(j4,j5)
          s2Omega(j4,j5)=s2Omega(j4,j5)+Omega(j4,j5)**2
          seqRho(i,j4,j5)=Rho(j4,j5)
          sRho(j4,j5)=sRho(j4,j5)+Rho(j4,j5)
          s2Rho(j4,j5)=s2Rho(j4,j5)+Rho(j4,j5)**2
         enddo
        enddo
c.......compute barDIC
        sum1=0.0d0
        do kk=1,K
         ein=0.0d0
         detf=0.0d0
         idim=narms(kk)
c............Set initial values
         do i1=1,4
          do i2=1,4
           B3(i1,i2)=0.0d0
           SigI(i1,i2)=0.0d0
          enddo
          do j1=1,nT
           EI(j1,i1)=0.0d0
           B2(i1,j1)=0.0d0
          enddo
          yk(i1)=0.0d0
          xx(i1)=0.0d0
          ywk(i1)=0.0d0
         enddo
c..........................................
         sum1=sum1 + real(narms(kk))*dlog(2.0d0*epi)
         g(kk,i)=0.50d0*real(narms(kk))*dlog(2.0d0*epi)
c
         icount1=0
         do j=1,ns
          if (ids(j) .eq. kk) then
           icount1=icount1+1
           SigI(icount1,icount1)=sig2(j)/real(npt(j))
           yk(icount1)=y(j)
          endif
         enddo
c
         icount1=0
         do j3=1,ns
          if (ids(j3) .eq. kk) then
           do j4=1,nT
            if (iarm(j3) .eq. j4) then
             icount1=icount1+1
             do j5=1,nT
               EI(j5,icount1)=randE(kk,j5,j4)
             enddo
            endif
           enddo
           endif
          enddo
c
        do j1=1,4
         do j2=1,nT
          do j3=1,nT
           B2(j1,j2)=B2(j1,j2)+EI(j3,j1)*Omega(j3,j2)
          enddo
         enddo
        enddo
c        
        do j1=1,4
         do j2=1,4
          do j3=1,nT
           B3(j1,j2)=B3(j1,j2)+B2(j1,j3)*EI(j3,j2)
          enddo
         enddo
        enddo
c
         icount1=0
         do j6=1,ns
          if (ids(j6) .eq. kk) then
           icount1=icount1+1
           do i5=1,nx
            xx(icount1)=xx(icount1)+x(i5,j6)*beta(i5)
           enddo
          endif
         enddo
         do i6=1,idim
          do i7=1,nT
           xx(i6)=xx(i6)+EI(i7,i6)*gam(i7)
          enddo
          ywk(i6)=yk(i6)-xx(i6)
         enddo
c.........narms(kk) = 2
        if (idim .eq. 2) then
         do i1=1,idim
          do i2=1,idim
           CI2(i1,i2)=0.0d0
           CIinv2(i1,i2)=0.0d0
           CIchol2(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            CIinv2(i3,i4)=CIinv2(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)
           enddo
          enddo
          call inverse(CIinv2, CI2, idim)
c          call DLFTDS(idim,CIinv2,idim,CIchol2,idim)
c          call DLFDDS(idim,CIchol2,idim,detf1,detf2)
         detf=DMGT(tol,idim,CIinv2) 
         sum1=sum1 + dlog(detf)
         g(kk,i)=g(kk,i) + 0.50d0*dlog(detf)
         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*CI2(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c.........narms(kk) = 3
        if (idim .eq. 3) then
         do i1=1,idim
          do i2=1,idim
           CI3(i1,i2)=0.0d0
           CIinv3(i1,i2)=0.0d0
           CIchol3(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            CIinv3(i3,i4)=CIinv3(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)
           enddo
          enddo
          call inverse(CIinv3, CI3, idim)
c          call DLFTDS(idim,CIinv3,idim,CIchol3,idim)
c          call DLFDDS(idim,CIchol3,idim,detf1,detf2)
          detf=DMGT(tol,idim,CIinv3) 
         sum1=sum1 + dlog(detf)
         g(kk,i)=g(kk,i) + 0.50d0*dlog(detf)
         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*CI3(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c.........narms(kk) = 4
        if (idim .eq. 4) then
         do i1=1,idim
          do i2=1,idim
           CI4(i1,i2)=0.0d0
           CIinv4(i1,i2)=0.0d0
           CIchol4(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            CIinv4(i3,i4)=CIinv4(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)
           enddo
          enddo
          call inverse(CIinv4, CI4, idim)
c          call DLFTDS(idim,CIinv4,idim,CIchol4,idim)
c          call DLFDDS(idim,CIchol4,idim,detf1,detf2)
          detf=DMGT(tol,idim,CIinv4) 
         sum1=sum1 + dlog(detf)
         g(kk,i)=g(kk,i) + 0.50d0*dlog(detf)
         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*CI4(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c          
         sum1=sum1 + ein
         g(kk,i)=g(kk,i) + 0.50d0*ein
c........kk=1,K end
        enddo
       barDIC=barDIC+sum1
c......end gibbs
 220   continue
c
       barDIC=barDIC/real(nrep)
c.......compute LPML
       do kk=1,K
        gmax(kk)=g(kk,1)
        do j1=2,nrep
         if (gmax(kk) .lt. g(kk,j1)) gmax(kk) = g(kk,j1)
        enddo
        sumrep=0.0d0
        do j1=1,nrep
         sumrep=sumrep+dexp(g(kk,j1)-gmax(kk))
        enddo
        alogcpo(kk)=-gmax(kk)-dlog( sumrep / real(nrep) )
        alpml=alpml+alogcpo(kk)
       enddo
       alpml = alpml - alogcpo(2)-alogcpo(10)
     1        - alogcpo(11)

c......compute posterior mean and sd for parameters
       do j=1,nx
        ebeta(j)=sbeta(j)/real(nrep)
        sebeta(j)=dsqrt( (s2beta(j)-
     1      real(nrep)*ebeta(j)**2)/real(nrep-1) )
       enddo
       do j=1,nT
        egam(j)=sgam(j)/real(nrep)
        segam(j)=dsqrt( (s2gam(j)-
     1      real(nrep)*egam(j)**2)/real(nrep-1) )
       enddo
       do i=1,ns
        esig2(i)=ssig2(i)/real(nrep)
c        sesig2(i)=dsqrt( (s2sig2(i)-
c     1      real(nrep)*esig2(i)**2)/real(nrep-1) )
       enddo
       do j1=1,nT
        do j2=j1,nT
         eOmega(j1,j2)=sOmega(j1,j2)/real(nrep)
         seOmega(j1,j2)=dsqrt( (s2Omega(j1,j2)-
     1      real(nrep)*eOmega(j1,j2)**2)/real(nrep-1) )
         eOmega(j2,j1)=eOmega(j1,j2)
         seOmega(j2,j1)=seOmega(j1,j2)
         eRho(j1,j2)=sRho(j1,j2)/real(nrep)
         seRho(j1,j2)=dsqrt( (s2Rho(j1,j2)-
     1      real(nrep)*eRho(j1,j2)**2)/real(nrep-1) )
         eRho(j2,j1)=eRho(j1,j2)
         seRho(j2,j1)=seRho(j1,j2)
        enddo
       enddo
c
c.......compute DICbar
        sum1=0.0d0
        do kk=1,K
         ein=0.0d0
         detf=0.0d0
         idim=narms(kk)
c............Set initial values
         do i1=1,4
          do i2=1,4
           B3(i1,i2)=0.0d0
           SigI(i1,i2)=0.0d0
          enddo
          do j1=1,nT
           EI(j1,i1)=0.0d0
           B2(i1,j1)=0.0d0
          enddo
          yk(i1)=0.0d0
          xx(i1)=0.0d0
          ywk(i1)=0.0d0
         enddo
c..........................................
         sum1=sum1 + real(narms(kk))*dlog(2.0d0*epi)
         icount1=0
         do j=1,ns
          if (ids(j) .eq. kk) then
           icount1=icount1+1
           SigI(icount1,icount1)=esig2(j)/real(npt(j))
           yk(icount1)=y(j)
          endif
         enddo
c
         icount1=0
         do j3=1,ns
          if (ids(j3) .eq. kk) then
           do j4=1,nT
            if (iarm(j3) .eq. j4) then
             icount1=icount1+1
             do j5=1,nT
               EI(j5,icount1)=randE(kk,j5,j4)
             enddo
            endif
           enddo
           endif
          enddo
c
        do j1=1,4
         do j2=1,nT
          do j3=1,nT
           B2(j1,j2)=B2(j1,j2)+EI(j3,j1)*eOmega(j3,j2)
          enddo
         enddo
        enddo
c        
        do j1=1,4
         do j2=1,4
          do j3=1,nT
           B3(j1,j2)=B3(j1,j2)+B2(j1,j3)*EI(j3,j2)
          enddo
         enddo
        enddo
c
         icount1=0
         do j6=1,ns
          if (ids(j6) .eq. kk) then
           icount1=icount1+1
           do i5=1,nx
            xx(icount1)=xx(icount1)+x(i5,j6)*ebeta(i5)
           enddo
          endif
         enddo
         do i6=1,idim
          do i7=1,nT
           xx(i6)=xx(i6)+EI(i7,i6)*egam(i7)
          enddo
          ywk(i6)=yk(i6)-xx(i6)
         enddo
c
c.........narms(kk) = 2
        if (idim .eq. 2) then
         do i1=1,idim
          do i2=1,idim
           CI2(i1,i2)=0.0d0
           CIinv2(i1,i2)=0.0d0
           CIchol2(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            CIinv2(i3,i4)=CIinv2(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)
           enddo
          enddo
          call inverse(CIinv2, CI2, idim)
c          call DLFTDS(idim,CIinv2,idim,CIchol2,idim)
c          call DLFDDS(idim,CIchol2,idim,detf1,detf2)
          detf=DMGT(tol,idim,CIinv2) 
         sum1=sum1 + dlog(detf)
         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*CI2(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c.........narms(kk) = 3
        if (idim .eq. 3) then
         do i1=1,idim
          do i2=1,idim
           CI3(i1,i2)=0.0d0
           CIinv3(i1,i2)=0.0d0
           CIchol3(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            CIinv3(i3,i4)=CIinv3(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)
           enddo
          enddo
          call inverse(CIinv3, CI3, idim)
c          call DLFTDS(idim,CIinv3,idim,CIchol3,idim)
c          call DLFDDS(idim,CIchol3,idim,detf1,detf2)
          detf=DMGT(tol,idim,CIinv3) 
         sum1=sum1 + dlog(detf)
         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*CI3(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c.........narms(kk) = 4
        if (idim .eq. 4) then
         do i1=1,idim
          do i2=1,idim
           CI4(i1,i2)=0.0d0
           CIinv4(i1,i2)=0.0d0
           CIchol4(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            CIinv4(i3,i4)=CIinv4(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)
           enddo
          enddo
          call inverse(CIinv4, CI4, idim)
c          call DLFTDS(idim,CIinv4,idim,CIchol4,idim)
c          call DLFDDS(idim,CIchol4,idim,detf1,detf2)
          detf=DMGT(tol,idim,CIinv4) 
         sum1=sum1 + dlog(detf)
         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*CI4(j1,j2)*ywk(j2)
          enddo
         enddo
        endif         
c        
         sum1=sum1  + ein
        enddo
       DICbar=sum1
       pD=barDIC-DICbar
       DIC=DICbar+2.0d0*pD
       DICinfo(1)=DIC
       DICinfo(2)=pD
       DICinfo(3)=barDIC
       DICinfo(4)=DICbar
       DICinfo(5)=alpml
c......compute 95% HPD for beta     
       do 333 j=1,nx
c       write(*,*) '95% HPD for beta(j)=',j
        do ii=1,2
          alow(ii)=0.0d0
          aupp(ii)=0.0d0
        enddo
        do i=1,nrep
         ahpd(i)=seqbeta(j,i)
        enddo
        call hpd(nrep,conf,ahpd,alow,aupp)
        betapo(j,1)=ebeta(j)
        betapo(j,2)=sebeta(j)
        betapo(j,3)=alow(1)
        betapo(j,4)=aupp(1)
 333   continue
c......compute 95% HPD for gam
       do 334 j=1,nT
        do ii=1,2
          alow(ii)=0.0d0
          aupp(ii)=0.0d0
        enddo
        do i=1,nrep
         ahpd(i)=seqgam(j,i)
        enddo
        call hpd(nrep,conf,ahpd,alow,aupp)
        gampo(j,1)=egam(j)
        gampo(j,2)=segam(j)
        gampo(j,3)=alow(1)
        gampo(j,4)=aupp(1)
 334   continue
c......compute 95% HPD for Omega
       do i=1,nrep
        ahpd(i)=0.0d0
       enddo
       do 335 j1=1,nT
        do 336 j2=j1,nT
         do ii=1,2
          alow(ii)=0.0d0
          aupp(ii)=0.0d0
         enddo
         do i=1,nrep
          ahpd(i)=seqOmega(i,j1,j2)
         enddo
         call hpd(nrep,conf,ahpd,alow,aupp)
         Omegalow(j1,j2)=alow(1)
         Omegaupp(j1,j2)=aupp(1)
 336    continue
 335   continue
c......compute 95% HPD for Rho
       do i=1,nrep
        ahpd(i)=0.0d0
       enddo
       do 337 j1=1,nT
        do 338 j2=j1,nT
         do ii=1,2
          alow(ii)=0.0d0
          aupp(ii)=0.0d0
         enddo
         do i=1,nrep
          ahpd(i)=seqRho(i,j1,j2)
         enddo
         call hpd(nrep,conf,ahpd,alow,aupp)
         Rholow(j1,j2)=alow(1)           
         Rhoupp(j1,j2)=aupp(1)
         Rholow(j2,j1)=Rholow(j1,j2)
         Rhoupp(j2,j1)=Rhoupp(j1,j2)
 338    continue
 337   continue
c......compute 95% HPD for tau2
       do i1=1,nT
       tau2po(i1,1)=eOmega(i1,i1)
       tau2po(i1,2)=seOmega(i1,i1)
       tau2po(i1,3)=Omegalow(i1,i1)
       tau2po(i1,4)=Omegaupp(i1,i1)
       enddo
c......compute 95% HPD for sig2
c       do 347 j=1,ns
c        do i=1,nrep
c         ahpd(i)=seqsig2(j,i)
c        enddo
c        do ii=1,2
c         alow(ii)=0.0d0
c         aupp(ii)=0.0d0
c        enddo
c        call hpd(nrep,conf,ahpd,alow,aupp)
c        sig2low(j)=alow(1)
c        sig2upp(j)=aupp(1)
c 347   continue
       return
       end
       
c       include 'gibbsNMAGroup1norm.f'
c       include 'hpd.f'
c       include 'optim1.f'
c..........................................
c used to generate standard gamma random variate
c       include 'gengam.f'
c       include 'rexpo.f'
c       include 'rnorm.f'
c       include 'runif.f'     
c........................................... 
c       include 'inverse.f'
c       include 'determinant.f'
c............................................      
c......cholesky decomposition used for multivariate normal
c       include 'cholesky.f' 

       subroutine gibbs( iseed )
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision alam(K),Rgam(ns),tau2(nT)
       double precision Rho(nT,nT), pRho(nT,nT)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns), igroup(nT),iseed
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecRgam/Rgam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecigroup/igroup
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE
c......generate variances: sig2
       call Gsig2( iseed )

c......generate beta
       call Gbeta( iseed )

c......generate Omega
       call GOmega( iseed )
c       write(*,*) 'after Omega'

c......generate Rgam
       call GRgam( iseed )
c       write(*,*) 'after Rgam'

       return
       end


       subroutine Gsig2( iseed )
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision alam(K),Rgam(ns)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns), iseed
       double precision shape,scale,einb,ein,r
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecRgam/Rgam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vecrandE/randE
c
       do i=1,ns
        shape=real(npt(i))/2.0d0+0.00001d0
        scale=0.0d0
        r=0.0d0
        einb=0.0d0
        ein=0.0d0
        do j1=1,nx
         einb=einb+beta(j1)*x(j1,i)
        enddo
        do i2=1,nT
         if (iarm(i) .eq. i2) then
          einb=einb+gam(i2)  
         endif
        enddo
        einb=einb+Rgam(i)
        ein=real(npt(i))*(y(i)-einb)**2
        scale=scale+0.5d0*(ein+((sd(i)**2)*(real(npt(i))-1.0d0)))
     1        +0.00001d0
c        call rnset( iseed )        
c        call DRNGAM(1,shape,r)
        r = gengam(shape, iseed)       
c        call rnget( iseed )
        sig2(i)=scale/r
       enddo
       return
       end



       subroutine Gbeta( iseed )
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision alam(K)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns),iseed
c
       double precision SigBeta(21,21), SigBetaInv(21,21), B4(21,21)
       double precision Rsig(21,21),RV(1,21), r1(21)
       double precision bhat(21),bstar(21), B6(21), yk(4)
       double precision EI(nT,4), B2(4,nT), B3(4,4), W(21,4), B5(21,4)
       double precision AI2(2,2), AIinv2(2,2)
       double precision AI3(3,3), AIinv3(3,3)
       double precision AI4(4,4), AIinv4(4,4)
       double precision tmp, SigBetavec(231), bn, upper(231)
       integer idim, icount, nn
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vecrandE/randE
c       external DCHFAC
c......Sigma_beta initial value
       do jj=1,21
        do jk=1,21
          SigBeta(jj,jk)=0.0d0
          SigBetaInv(jj,jk)=0.0d0
          Rsig(jj,jk)=0.0d0
        enddo
        SigBetaInv(jj,jj)=0.00001d0
        bstar(jj)=0.0d0
        bhat(jj)=0.0d0
        RV(1,jj)=0.0d0
        r1(jj)=0.0d0
       enddo

c......Sigma_beta amd Mu_beta calculate
c...........Set initial value
       do 113 kk=1,K
        idim=narms(kk)
        do i1=1,4
          yk(i1)=0.0d0
          do j1=1,nT
           EI(j1,i1)=0.0d0
           B2(i1,j1)=0.0d0
          enddo
          do j2=1,21
           W(j2,i1)=0.0d0
           B5(j2,i1)=0.0d0
          enddo
          do i2=1,4
           B3(i1,i2)=0.0d0
          enddo
        enddo
        do jj=1,21
         B6(jj)=0.0d0
         do jk=1,21
          B4(jj,jk)=0.0d0
         enddo
        enddo

c...........Calculate Inverse of sum of covariance matrix
c.......E'_k Omega E_k
        icount=0
        do j3=1,ns
          if (ids(j3) .eq. kk) then
           do j4=1,nT
            if (iarm(j3) .eq. j4) then
             icount=icount+1
             do j5=1,nT
               EI(j5,icount)=randE(kk,j5,j4)
             enddo
            endif
           enddo
          endif
        enddo
c
        do j1=1,4
         do j2=1,nT
          do j3=1,nT
           B2(j1,j2)=B2(j1,j2)+EI(j3,j1)*Omega(j3,j2)
          enddo
         enddo
        enddo
c        
        do j1=1,4
         do j2=1,4
          do j3=1,nT
           B3(j1,j2)=B3(j1,j2)+B2(j1,j3)*EI(j3,j2)
          enddo
         enddo
        enddo
c.......W transpose
        icount=0
        do j6=1,ns
         if (ids(j6) .eq. kk) then
          icount=icount+1
          do i5=1,nx
            W(i5,icount)=x(i5,j6)
          enddo
         endif
        enddo
        do j7=1,idim
         do j8=1,nT
          W(j8+nx,j7)=EI(j8,j7)
         enddo
        enddo
c......narms(kk) = 2
        if (idim .eq. 2) then
        do i1=1,idim
         do i2=1,idim
          AI2(i1,i2)=0.0d0
          AIinv2(i1,i2)=0.0d0
         enddo
        enddo
        icount=1
        do j=1,ns
         if (ids(j) .eq. kk) then
          AI2(icount,icount)=sig2(j)/real(npt(j))
          yk(icount)=y(j)
          icount=icount+1
         endif
        enddo
        do i1=1,idim
         do j1=1,idim
          AI2(i1,j1)=AI2(i1,j1) + B3(i1,j1)/alam(kk)
         enddo
        enddo
        call inverse(AI2, AIinv2, idim)
        do i1=1,21
         do j1=1,idim
          do i2=1,idim
           B5(i1,j1)=B5(i1,j1)+W(i1,i2)*AIinv2(i2,j1)
          enddo
         enddo
        enddo
        endif
c......narms(kk) = 3
        if (idim .eq. 3) then
        do i1=1,idim
         do i2=1,idim
          AI3(i1,i2)=0.0d0
          AIinv3(i1,i2)=0.0d0
         enddo
        enddo
        icount=1
        do j=1,ns
         if (ids(j) .eq. kk) then
          AI3(icount,icount)=sig2(j)/real(npt(j))
          yk(icount)=y(j)
          icount=icount+1
         endif
        enddo
        do i1=1,idim
         do j1=1,idim
          AI3(i1,j1)=AI3(i1,j1) + B3(i1,j1)/alam(kk)
         enddo
        enddo
        call inverse(AI3, AIinv3, idim)
        do i1=1,21
         do j1=1,idim
          do i2=1,idim
           B5(i1,j1)=B5(i1,j1)+W(i1,i2)*AIinv3(i2,j1)
          enddo
         enddo
        enddo
       endif
c......narms(kk) = 4
        if (idim .eq. 4) then
        do i1=1,idim
         do i2=1,idim
          AI4(i1,i2)=0.0d0
          AIinv4(i1,i2)=0.0d0
         enddo
        enddo
        icount=1
        do j=1,ns
         if (ids(j) .eq. kk) then
          AI4(icount,icount)=sig2(j)/real(npt(j))
          yk(icount)=y(j)
          icount=icount+1
         endif
        enddo
        do i1=1,idim
         do j1=1,idim
          AI4(i1,j1)=AI4(i1,j1) + B3(i1,j1)/alam(kk)
         enddo
        enddo
        call inverse(AI4, AIinv4, idim)
        do i1=1,21
         do j1=1,idim
          do i2=1,idim
           B5(i1,j1)=B5(i1,j1)+W(i1,i2)*AIinv4(i2,j1)
          enddo
         enddo
        enddo
       endif
c
        do i1=1,21
         do j1=1,21
          do i2=1,4
           B4(i1,j1)=B4(i1,j1)+ B5(i1,i2)*W(j1,i2)
          enddo
         enddo
        enddo
c...............Sigma_beta_inverse
        do j9=1,21
         do j10=1,21
           SigBetaInv(j9,j10)=SigBetaInv(j9,j10)+B4(j9,j10)
         enddo
        enddo
c...............bstar

        do j11=1,21
          do j12=1,4
           B6(j11)=B6(j11)+B5(j11,j12)*yk(j12)
          enddo
          bstar(j11)=bstar(j11) + B6(j11)
        enddo
 113   continue
c......Sigma_beta
       call inverse(SigBetaInv, SigBeta, 21)
c......Mu_beta
       do i=1,21
        do j=1,21
         bhat(i)=bhat(i)+SigBeta(i,j)*bstar(j)
        enddo
       enddo
c.......generate multivariate normal
        lda=21
        ldr=21
c        call DCHFAC(21,SigBeta,lda,tol,irank,Rsig,ldr)
c       
       tmp=0.0d0
       do i1=1,lda
        do j1=1,i1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        SigBetavec(ij)=SigBeta(i1, j1)       
        enddo
       enddo
       bn=real(lda)*(real(lda)+1.0d0)*0.5d0
       nn=int(bn)
       call cholesky(SigBetavec, lda, nn, upper, nullty, ifault) 
c       
       tmp=0.0d0
       do j1=1,lda
        do i1=1,j1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        Rsig(i1, j1)=upper(ij)       
        enddo
       enddo                 
c    
        do i1=1,21
         r1(i1)=r8_normal_01_sample(iseed)
        enddo                
        do i1=1,21
         do j1=1,21
         RV(1,i1) = RV(1,i1) + Rsig(j1,i1) * r1(j1) 
         enddo
        enddo              
        do j=1,nx
         beta(j)=RV(1,j)+ bhat(j)
c         write(*,*) j, beta(j)
        enddo
        do jj=1,nT
         gam(jj)=RV(1,jj+nx)+ bhat(jj+nx)
c         write(*,*) jj, gam(jj)
        enddo
       return
       end


       subroutine GOmega( iseed )
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns), tau2(nT)
       double precision Rho(nT,nT), pRho(nT,nT)
       double precision beta(nx),gam(nT), Omega(nT,nT), sig2(ns)
       double precision alam(K)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns),igroup(nT),iseed
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecigroup/igroup       
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE

       call Gtau2( iseed )
c
       call GRho( iseed )
c
       do j1=1,nT
        Omega(j1,j1)=tau2(j1)
       enddo
c
       do j1=1,nT
        do j2=1,nT
         Omega(j1,j2)=Rho(j1,j2)*dsqrt(Omega(j1,j1)*Omega(j2,j2))
        enddo
       enddo
       return
       end


       subroutine Gtau2( iseed )
c      Hao Li, Uconn
c      March 28, 2016
c       use linear_operators
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT), Omega(nT,nT), sig2(ns)
       double precision alam(K),tau2(nT), Rho(nT,nT), pRho(nT,nT)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns), igroup(nT),iseed
c
       double precision start(1), step(1), ximax(1), plmax
       double precision e1, fl(3, 1), cl(5, 3)
       double precision dl(5, 1), al(3, 3), el(3, 1), alinv(3, 3)
       double precision sigmaa, xistar, r1(1), r2(1), rat1
       integer konvge, kcount, ifault1, numres, icount, igrp, iflag
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecigroup/igroup
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE
       common  /vecgrp/igrp
       external bNloglike
c
       do 474 ii=1,nT
        ximax(1)=0.0d0
        plmax=0.0d0
        rat1=0.0d0
        r1(1)=0.0d0
        r2(1)=0.0d0
        nopt = 1
        reqmin = 1.0d-20
        konvge = 5
        kcount = 1000
        step(1)=0.2d0
        if (ii .eq. 1) then 
c........simulate
        igrp=ii
        start(1)=dlog(tau2(ii))
        call nelmin(bNloglike, nopt, start, ximax, plmax, reqmin,
     1         step, konvge, kcount, icount, numres, ifault1)
        step1=0.50d0
 55     do i=1, 5
          e1 = real(i-3)
          cl(i, 1)=(ximax(1)+ e1*step1)**2
          cl(i, 2)= ximax(1)+ e1*step1
          cl(i, 3)=1.0d0
          dl(i, 1)=bNloglike(ximax(1) + e1*step1)
        end do
c
        do ni=1, 5
         if (ni .ne. 3) then
          if(dl(ni, 1) .le. plmax) then
           step1=step1*1.50d0
           GOTO 55
          endif
         endif
        enddo
c
        al = matmul(transpose(cl),cl)
        el = matmul(transpose(cl),dl)
c        call lin_sol_self(al, el, fl)
        ial=3
        call inverse(al, alinv, ial)
        fl = matmul(alinv, el)        
c
        sigmaa = dsqrt(1.0d0/(fl(1, 1)*2.0d0))
        r1(1)=r8_normal_01_sample(iseed)
        xistar = sigmaa * r1(1)  + ximax(1)
c
         rat1 = -bNloglike(xistar) + bNloglike(start(1))
     1            -0.50d0*(start(1) - ximax(1))**2/sigmaa**2
     2           + 0.50d0*(xistar - ximax(1))**2/sigmaa**2
c
         if (rat1 .ge. 0.0d0) then
           do j1=1,nT
            if (igroup(j1) .eq. igroup(ii)) tau2(j1) = dexp(xistar)
           enddo
         else
         r2(1)=r8_uniform_01_sample(iseed)
          if ( dlog(r2(1)) .le. rat1 ) then
           do j1=1,nT
            if (igroup(j1) .eq. igroup(ii)) tau2(j1) = dexp(xistar)
           enddo
          endif
         endif
c         
c        
        else  
        iflag=0
        do jj=1,ii-1
         if (igroup(jj) .eq. igroup(ii)) iflag=1
        enddo        
        if (iflag .eq. 0) then   
c........simulate
        igrp=ii
        start(1)=dlog(tau2(ii))
        call nelmin(bNloglike, nopt, start, ximax, plmax, reqmin,
     1         step, konvge, kcount, icount, numres, ifault1)
        step1=0.50d0
 65     do i=1, 5
          e1 = real(i-3)
          cl(i, 1)=(ximax(1)+ e1*step1)**2
          cl(i, 2)= ximax(1)+ e1*step1
          cl(i, 3)=1.0d0
          dl(i, 1)=bNloglike(ximax(1) + e1*step1)
        end do
c
        do ni=1, 5
         if (ni .ne. 3) then
          if(dl(ni, 1) .le. plmax) then
           step1=step1*1.50d0
           GOTO 65
          endif
         endif
        enddo
c
        al = matmul(transpose(cl),cl)
        el = matmul(transpose(cl),dl)
c        call lin_sol_self(al, el, fl)
        ial=3
        call inverse(al, alinv, ial)
        fl = matmul(alinv, el)        
c
        sigmaa = dsqrt(1.0d0/(fl(1, 1)*2.0d0))
        r1(1)=r8_normal_01_sample(iseed)
        xistar = sigmaa * r1(1)  + ximax(1)
c
         rat1 = -bNloglike(xistar) + bNloglike(start(1))
     1            -0.50d0*(start(1) - ximax(1))**2/sigmaa**2
     2           + 0.50d0*(xistar - ximax(1))**2/sigmaa**2
c
         if (rat1 .ge. 0.0d0) then
           do j1=1,nT
            if (igroup(j1) .eq. igroup(ii)) tau2(j1) = dexp(xistar)
           enddo
         else
          r2(1)=r8_uniform_01_sample(iseed)
          if ( dlog(r2(1)) .le. rat1 ) then
           do j1=1,nT
            if (igroup(j1) .eq. igroup(ii)) tau2(j1) = dexp(xistar)
           enddo
          endif             
         endif
c         
        endif
       endif
 474   continue
        return
        end

c..... - log likelihood
       double precision function bNloglike(xi)
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       parameter (tol=2.220446049250313E-14)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision alam(K), tau2(nT), Rho(nT,nT), pRho(nT,nT)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns),igroup(nT)

       double precision xi, ximat(nT,nT), V(nT)
       double precision yk(4), xx(4), ywk(4), EI(nT,4), B2(4,nT)
       double precision B3(4,4), SigI(4,4)
       double precision ein, sum1, detf
c       double precision detf1, detf2
       integer icount, idim, igrp
       double precision AI2(2,2), AIinv2(2,2), AIchol2(2,2)
       double precision AI3(3,3), AIinv3(3,3), AIchol3(3,3)
       double precision AI4(4,4), AIinv4(4,4), AIchol4(4,4)
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecigroup/igroup
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE
       common  /vecgrp/igrp
c       external DLFTDS, DLFDDS
c
c       V=tau2
       do i1=1,nT
       V(i1)=tau2(i1)
       if (igroup(i1) .eq. igroup(igrp)) V(i1) = dexp(xi)
       enddo
        do j1=1,nT
         do j2=1,nT
          ximat(j1,j2)=dsqrt(V(j1))*Rho(j1,j2)*dsqrt(V(j2))
         enddo
        enddo
       sum1=0.0d0
       do 313 kk=1,K
         idim=narms(kk)
         ein=0.0d0
         detf=0.0d0
c............Set initial values
         do i1=1,4
          do i2=1,4
           B3(i1,i2)=0.0d0
           SigI(i1,i2)=0.0d0
          enddo
          do j1=1,nT
           EI(j1,i1)=0.0d0
           B2(i1,j1)=0.0d0
          enddo
          yk(i1)=0.0d0
          xx(i1)=0.0d0
          ywk(i1)=0.0d0
         enddo
c..........................................
         icount=0
         do j=1,ns
          if (ids(j) .eq. kk) then
           icount=icount+1
           SigI(icount,icount)=sig2(j)/real(npt(j))
           yk(icount)=y(j)
          endif
         enddo
c
         icount=0
         do j3=1,ns
          if (ids(j3) .eq. kk) then
           do j4=1,nT
            if (iarm(j3) .eq. j4) then
             icount=icount+1
             do j5=1,nT
               EI(j5,icount)=randE(kk,j5,j4)
             enddo
            endif
           enddo
          endif
         enddo
c
        do j1=1,4
         do j2=1,nT
          do j3=1,nT
           B2(j1,j2)=B2(j1,j2)+EI(j3,j1)*ximat(j3,j2)
          enddo
         enddo
        enddo
c        
        do j1=1,4
         do j2=1,4
          do j3=1,nT
           B3(j1,j2)=B3(j1,j2)+B2(j1,j3)*EI(j3,j2)
          enddo
         enddo
        enddo
c
         icount=0
         do j6=1,ns
          if (ids(j6) .eq. kk) then
           icount=icount+1
           do i5=1,nx
            xx(icount)=xx(icount)+x(i5,j6)*beta(i5)
           enddo
          endif
         enddo
c
         do i6=1,idim
          do i7=1,nT
           xx(i6)=xx(i6)+EI(i7,i6)*gam(i7)
          enddo
          ywk(i6)=yk(i6)-xx(i6)
         enddo
c
c.........narms(kk) = 2
        if (idim .eq. 2) then
         do i1=1,idim
          do i2=1,idim
           AI2(i1,i2)=0.0d0
           AIinv2(i1,i2)=0.0d0
           AIchol2(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            AIinv2(i3,i4)=AIinv2(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)/alam(kk)
           enddo
          enddo
          call inverse(AIinv2, AI2, idim)
c          call DLFTDS(idim,AIinv2,idim,AIchol2,idim)
c          call DLFDDS(idim,AIchol2,idim,detf1,detf2)
          detf=DMGT(tol,idim,AIinv2) 
          sum1=sum1 - 0.50d0*dlog(detf)

         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*AI2(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c.........narms(kk) = 3
        if (idim .eq. 3) then
         do i1=1,idim
          do i2=1,idim
           AI3(i1,i2)=0.0d0
           AIinv3(i1,i2)=0.0d0
           AIchol3(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            AIinv3(i3,i4)=AIinv3(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)/alam(kk)
           enddo
          enddo
          call inverse(AIinv3, AI3, idim)
c          call DLFTDS(idim,AIinv3,idim,AIchol3,idim)
c          call DLFDDS(idim,AIchol3,idim,detf1,detf2)
          detf=DMGT(tol,idim,AIinv3) 
          sum1=sum1 - 0.50d0*dlog(detf)

         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*AI3(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c.........narms(kk) = 4
        if (idim .eq. 4) then
         do i1=1,idim
          do i2=1,idim
           AI4(i1,i2)=0.0d0
           AIinv4(i1,i2)=0.0d0
           AIchol4(i1,i2)=0.0d0
          enddo
         enddo
          do i3=1,idim
           do i4=1,idim
            AIinv4(i3,i4)=AIinv4(i3,i4)+SigI(i3,i4)
     1                   +B3(i3,i4)/alam(kk)
           enddo
          enddo
          call inverse(AIinv4, AI4, idim)
c          call DLFTDS(idim,AIinv4,idim,AIchol4,idim)
c          call DLFDDS(idim,AIchol4,idim,detf1,detf2)
          detf=DMGT(tol,idim,AIinv4) 
          sum1=sum1 - 0.50d0*dlog(detf)

         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*AI4(j1,j2)*ywk(j2)
          enddo
         enddo
        endif
c
        sum1 = sum1 - 0.50d0*ein      
 313   continue
c
       sum1=sum1 - ( (0.5d0+1.0d0)*xi + 0.5d0/dexp(xi) ) + xi
       bNloglike = -sum1
       return
       end


       subroutine GRho( iseed )
c      Hao Li, Uconn
c      March 28, 2016
c       use linear_operators
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT), Omega(nT,nT), sig2(ns)
       double precision alam(K),tau2(nT), Rho(nT,nT), pRho(nT,nT)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns),iseed
c
       double precision start(1), step(1), zprhomax(1), plmax
       double precision e1, fl(3, 1), cl(5, 3)
       double precision dl(5, 1), al(3, 3), el(3, 1), alinv(3, 3)
       double precision sigmaa, zprhostar, r1(1), r2(1), rat1, step1
       integer konvge, kcount, ifault1, numres, icount
       integer index1, index2
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE
       common  /vecindex1/index1
       common  /vecindex2/index2
       external cNloglike
c
       do 405 iR=1, nT
        do 410 iC=iR+1, nT
         index1=iR
         index2=iC
         zprhomax(1)=0.0d0
         plmax=0.0d0
         e1=0.0d0
         rat1=0.0d0
         r1(1)=0.0d0
         r2(1)=0.0d0
         nopt = 1
         reqmin = 1.0d-20
         konvge = 5
         kcount = 1000
         step(1)=0.20d0
         start(1)=0.50d0*dlog((1.0d0+pRho(iR,iC))/(1.0d0-pRho(iR,iC)))
         call nelmin(cNloglike, nopt, start, zprhomax, plmax, reqmin,
     1         step, konvge, kcount, icount, numres, ifault1)
         step1=0.50d0
 56      do i=1, 5
          e1 = real(i-3)
          cl(i, 1)=(zprhomax(1)+ e1*step1)**2
          cl(i, 2)= zprhomax(1)+ e1*step1
          cl(i, 3)=1.0d0
          dl(i, 1)=cNloglike(zprhomax(1) + e1*step1)
         end do
c
        do ni=1, 5
         if (ni .ne. 3) then
          if(dl(ni, 1) .le. plmax) then
           step1=step1*1.20d0
           GOTO 56
          endif
         endif
        enddo
c
         al = matmul(transpose(cl),cl)
         el = matmul(transpose(cl),dl)
c         call lin_sol_self(al, el, fl)
        ial=3
        call inverse(al, alinv, ial)
        fl = matmul(alinv, el) 
c
         sigmaa = dsqrt(1.0d0/(fl(1, 1)*2.0d0))
c         call rnset( iseed )
         r1(1)=r8_normal_01_sample(iseed)
c         call rnget( iseed )
         zprhostar = sigmaa * r1(1)  + zprhomax(1)
c
          rat1 = -cNloglike(zprhostar) + cNloglike(start(1))
     1            -0.50d0*(start(1) - zprhomax(1))**2/sigmaa**2
     2           + 0.50d0*(zprhostar - zprhomax(1))**2/sigmaa**2

          if (rat1 .ge. 0.0d0) then
            pRho(iR,iC) = (dexp(2.0d0*zprhostar)-1.0d0)/
     1                    (dexp(2.0d0*zprhostar)+1.0d0)
            call prhoTorho(nT, pRho, Rho)
          else
c           call rnset(iseed)
           r2(1)=r8_uniform_01_sample( iseed )
c           call rnget(iseed)
           if (dlog(r2(1)) .le. rat1) then
             pRho(iR,iC) = (dexp(2.0d0*zprhostar)-1.0d0)/
     1                     (dexp(2.0d0*zprhostar)+1.0d0)
             call prhoTorho(nT, pRho, Rho)
           endif
          endif
 410    continue
 405   continue
       return
       end


       subroutine prhoTorho(nT, pRho, Rho)
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       double precision  Rho(nT,nT), pRho(nT,nT)
       double precision rvec1(nT-2), rvec3(nT-2), rr11, rr13, rr33
       double precision rmat1(1,1), rmatinv1(1,1)
       double precision rmat2(2,2), rmatinv2(2,2)
       double precision rmat3(3,3), rmatinv3(3,3)
       double precision rmat4(4,4), rmatinv4(4,4)
       double precision rmat5(5,5), rmatinv5(5,5)
       double precision rmat6(6,6), rmatinv6(6,6)
       double precision rmat7(7,7), rmatinv7(7,7)
       double precision rmat8(8,8), rmatinv8(8,8)
       double precision rmat9(9,9), rmatinv9(9,9)
c
        do j=1,nT
         Rho(j,j)=pRho(j,j)
        enddo
c
        do j=1,nT-1
         Rho(j,j+1)=pRho(j,j+1)
         Rho(j+1,j)=Rho(j,j+1)
        enddo
c
        do 555 iL=2,nT-1
         do 556 iR=1,nT-iL
          do i1=1,nT-2
           rvec1(i1)=0.0d0
           rvec3(i1)=0.0d0
          enddo
          do i1=1,iL-1
           rvec1(i1)=Rho(iR,iR+i1)
           rvec3(i1)=Rho(iR+i1,iR+iL)
          enddo
          rr11=0.0d0
          rr13=0.0d0
          rr33=0.0d0
c.......iL=2         
         if (iL .eq. 2) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat1(i1,j1)=0.0d0
           rmatinv1(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat1(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat1, rmatinv1, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv1(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv1(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv1(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=3
         if (iL .eq. 3) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat2(i1,j1)=0.0d0
           rmatinv2(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat2(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat2, rmatinv2, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv2(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv2(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv2(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=4
         if (iL .eq. 4) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat3(i1,j1)=0.0d0
           rmatinv3(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat3(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat3, rmatinv3, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv3(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv3(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv3(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=5
         if (iL .eq. 5) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat4(i1,j1)=0.0d0
           rmatinv4(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat4(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat4, rmatinv4, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv4(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv4(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv4(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=6
         if (iL .eq. 6) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat5(i1,j1)=0.0d0
           rmatinv5(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat5(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat5, rmatinv5, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv5(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv5(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv5(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=7
         if (iL .eq. 7) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat6(i1,j1)=0.0d0
           rmatinv6(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat6(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat6, rmatinv6, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv6(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv6(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv6(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=8
         if (iL .eq. 8) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat7(i1,j1)=0.0d0
           rmatinv7(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat7(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat7, rmatinv7, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv7(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv7(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv7(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif         
c.......iL=9
         if (iL .eq. 9) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat8(i1,j1)=0.0d0
           rmatinv8(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat8(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat8, rmatinv8, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv8(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv8(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv8(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=10
         if (iL .eq. 10) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat9(i1,j1)=0.0d0
           rmatinv9(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat9(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat9, rmatinv9, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv9(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv9(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv9(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c
        Rho(iR,iR+iL) = rr13 +
     1                pRho(iR,iR+iL)*dsqrt((1.0d0-rr11)*(1.0d0-rr33))
        Rho(iR+iL,iR) = Rho(iR,iR+iL)

 556    continue
 555   continue

        return
        end



       double precision function cNloglike(zprho)
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       parameter (tol=2.220446049250313E-14)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision alam(K), tau2(nT), Rho(nT,nT), pRho(nT,nT)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns)
c
       double precision zprho, temppRho(nT,nT), tempRho(nT,nT)
       double precision yk(4), xx(4), ywk(4), EI(nT,4), B2(4,nT)
       double precision B3(4,4), SigI(4,4), tempOmega(nT,nT)
       double precision ein, sum1, detf
c       double precision detf1, detf2
       integer icount, idim
       double precision AI2(2,2), AIinv2(2,2), AIchol2(2,2)
       double precision AI3(3,3), AIinv3(3,3), AIchol3(3,3)
       double precision AI4(4,4), AIinv4(4,4), AIchol4(4,4)
       integer iL
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vectau2/tau2
       common  /vecRho/Rho
       common  /vecpRho/pRho
       common  /vecrandE/randE
c
       common  /vecindex1/index1
       common  /vecindex2/index2
c       external DLFTDS, DLFDDS
c
       do j1=1,nT
        do j2=1,nT
         temppRho(j1,j2)=pRho(j1,j2)
        enddo
       enddo
       temppRho(index1, index2)= (dexp(2.0d0*zprho)-1.0d0)/
     1                        (dexp(2.0d0*zprho)+1.0d0)

       temppRho(index2, index1)= temppRho(index1, index2)
c
       call prhoTorho(nT, temppRho, tempRho)
c
       sum1=0.0d0
       do 400 kk=1,K
         idim=narms(kk)
         ein=0.0d0
         detf=0.0d0
c         detf1=0.0d0
c         detf2=0.0d0
c............Set initial values
         do i1=1,4
          do i2=1,4
           B3(i1,i2)=0.0d0
           SigI(i1,i2)=0.0d0
          enddo
          do j1=1,nT
           EI(j1,i1)=0.0d0
           B2(i1,j1)=0.0d0
          enddo
          yk(i1)=0.0d0
          xx(i1)=0.0d0
          ywk(i1)=0.0d0
         enddo
c..........................................
         icount=0
         do j=1,ns
          if (ids(j) .eq. kk) then
           icount=icount+1
           SigI(icount,icount)=sig2(j)/real(npt(j))
           yk(icount)=y(j)
          endif
         enddo
c
         icount=0
         do j3=1,ns
          if (ids(j3) .eq. kk) then
           do j4=1,nT
            if (iarm(j3) .eq. j4) then
             icount=icount+1
             do j5=1,nT
               EI(j5,icount)=randE(kk,j5,j4)
             enddo
            endif
           enddo
          endif
         enddo
c
         do j1=1,nT
          do j2=1,nT
           tempOmega(j1,j2)=tempRho(j1,j2)*dsqrt(tau2(j1)*tau2(j2))
          enddo
         enddo
c
        do j1=1,4
         do j2=1,nT
          do j3=1,nT
           B2(j1,j2)=B2(j1,j2)+EI(j3,j1)*tempOmega(j3,j2)
          enddo
         enddo
        enddo
c        
        do j1=1,4
         do j2=1,4
          do j3=1,nT
           B3(j1,j2)=B3(j1,j2)+B2(j1,j3)*EI(j3,j2)
          enddo
         enddo
        enddo
c
         icount=0
         do j6=1,ns
          if (ids(j6) .eq. kk) then
           icount=icount+1
           do i5=1,nx
            xx(icount)=xx(icount)+x(i5,j6)*beta(i5)
           enddo
          endif
         enddo
         do i6=1,idim
          do i7=1,nT
           xx(i6)=xx(i6)+EI(i7,i6)*gam(i7)
          enddo
          ywk(i6)=yk(i6)-xx(i6)
         enddo
c.........narms(kk) = 2
       if (idim .eq. 2) then
         do i1=1,idim
          do i2=1,idim
           AI2(i1,i2)=0.0d0
           AIinv2(i1,i2)=0.0d0
           AIchol2(i1,i2)=0.0d0
          enddo
         enddo
         do i3=1,idim
          do i4=1,idim
           AIinv2(i3,i4)=AIinv2(i3,i4)+SigI(i3,i4)
     1                   + B3(i3,i4)/alam(kk)
          enddo
         enddo
         call inverse(AIinv2, AI2, idim)
c         call DLFTDS(idim,AIinv2,idim,AIchol2,idim)
c         call DLFDDS(idim,AIchol2,idim,detf1,detf2)
          detf=DMGT(tol,idim,AIinv2) 
         sum1=sum1 - 0.50d0*dlog(detf)

         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*AI2(j1,j2)*ywk(j2)
          enddo
         enddo
        sum1 = sum1 - 0.50d0*ein
        endif
c.........narms(kk) = 3
       if (idim .eq. 3) then
         do i1=1,idim
          do i2=1,idim
           AI3(i1,i2)=0.0d0
           AIinv3(i1,i2)=0.0d0
           AIchol3(i1,i2)=0.0d0
          enddo
         enddo
         do i3=1,idim
          do i4=1,idim
           AIinv3(i3,i4)=AIinv3(i3,i4)+SigI(i3,i4)
     1                   + B3(i3,i4)/alam(kk)
          enddo
         enddo
         call inverse(AIinv3, AI3, idim)
c         call DLFTDS(idim,AIinv3,idim,AIchol3,idim)
c         call DLFDDS(idim,AIchol3,idim,detf1,detf2)
          detf=DMGT(tol,idim,AIinv3) 
         sum1=sum1 - 0.50d0*dlog(detf)

         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*AI3(j1,j2)*ywk(j2)
          enddo
         enddo
        sum1 = sum1 - 0.50d0*ein
        endif
c.........narms(kk) = 4
       if (idim .eq. 4) then
         do i1=1,idim
          do i2=1,idim
           AI4(i1,i2)=0.0d0
           AIinv4(i1,i2)=0.0d0
           AIchol4(i1,i2)=0.0d0
          enddo
         enddo
         do i3=1,idim
          do i4=1,idim
           AIinv4(i3,i4)=AIinv4(i3,i4)+SigI(i3,i4)
     1                   + B3(i3,i4)/alam(kk)
          enddo
         enddo
         call inverse(AIinv4, AI4, idim)
c         call DLFTDS(idim,AIinv4,idim,AIchol4,idim)
c         call DLFDDS(idim,AIchol4,idim,detf1,detf2)
          detf=DMGT(tol,idim,AIinv4) 
         sum1=sum1 - 0.50d0*dlog(detf)

         do j1=1,idim
          do j2=1,idim
           ein=ein+ywk(j1)*AI4(j1,j2)*ywk(j2)
          enddo
         enddo
        sum1 = sum1 - 0.50d0*ein
        endif
        
 400   continue
c
       iL = index2 - index1
       sum1=sum1 + 0.50d0*(10.0d0-real(iL))*dlog(1.0d0-
     1    ((dexp(2.0d0*zprho)-1.0d0)/(dexp(2.0d0*zprho)+1.0d0))**2)+
     2     2.0d0*zprho - 2.0d0*dlog(dexp(2.0d0*zprho)+1.0d0)
c
       cNloglike = - sum1
        return
        end


       subroutine GRgam( iseed )
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       parameter (ns=73, K=29, nx=10, nT=11)
       double precision y(ns),sd(ns),x(nx,ns)
       double precision beta(nx),gam(nT),Omega(nT,nT),sig2(ns)
       double precision alam(K),Rgam(ns)
       double precision randE(1:K,1:nT,1:nT)
       integer ids(ns),iarm(ns),narms(K),npt(ns),iseed
c
       double precision yk(4), xx(4), ywk1(4)
       double precision ywk2(4),EI(nT,4), B2(4,nT), r1(4)
       double precision B3(4,4), SigIinv(4,4), Rgamhat(4),RV(1,4)
       double precision AI2(2,2), AIinv2(2,2),AIchol2(2,2)
       double precision AI3(3,3), AIinv3(3,3),AIchol3(3,3)
       double precision AI4(4,4), AIinv4(4,4),AIchol4(4,4)
       double precision AI2vec(3), AI3vec(6), AI4vec(10)
       double precision upper2(3), upper3(6), upper4(10) 
       double precision tmp, bn
       integer idim, icount, isum, nn
c
       common  /vecy/y
       common  /vecsd/sd
       common  /vecnpt/npt
       common  /vecx/x
       common  /vecids/ids
       common  /veciarm/iarm
       common  /vecnarms/narms
       common  /vecsig2/sig2
       common  /vecbeta/beta
       common  /vecgam/gam
       common  /vecRgam/Rgam
       common  /vecalam/alam
       common  /vecOmega/Omega
       common  /vecrandE/randE
c       external DCHFAC
c
       do 557 kk=1,K
        isum=0
        idim=narms(kk)
c............Set initial values
        do i1=1,4
          do j1=1,nT
           EI(j1,i1)=0.0d0
           B2(i1,j1)=0.0d0
          enddo
          yk(i1)=0.0d0
          xx(i1)=0.0d0
          ywk1(i1)=0.0d0
          ywk2(i1)=0.0d0
          Rgamhat(i1)=0.0d0
          RV(1,i1)=0.0d0
          r1(i1)=0.0d0
          do j2=1,4
           B3(i1,j2)=0.0d0
           SigIinv(i1,j2)=0.0d0
          enddo
        enddo
c.......Sigma_Rgamma:  AI
         icount=0
         do j=1,ns
          if (ids(j) .eq. kk) then
           icount=icount+1
           SigIinv(icount,icount)=real(npt(j))/sig2(j)
           yk(icount)=y(j)
          endif
         enddo
c
        icount=0
        do j3=1,ns
          if (ids(j3) .eq. kk) then
           do j4=1,nT
            if (iarm(j3) .eq. j4) then
             icount=icount+1
             do j5=1,nT
               EI(j5,icount)=randE(kk,j5,j4)
             enddo
            endif
           enddo
          endif
        enddo
c
        do j1=1,4
         do j2=1,nT
          do j3=1,nT
           B2(j1,j2)=B2(j1,j2)+EI(j3,j1)*Omega(j3,j2)
          enddo
         enddo
        enddo
c        
        do j1=1,4
         do j2=1,4
          do j3=1,nT
           B3(j1,j2)=B3(j1,j2)+B2(j1,j3)*EI(j3,j2)
          enddo
         enddo
        enddo
c
        icount=0
        do j6=1,ns
         if (ids(j6) .eq. kk) then
          icount=icount+1
          do i5=1,nx
            xx(icount)=xx(icount)+x(i5,j6)*beta(i5)
          enddo
         endif
        enddo
c
        do i6=1,idim
         do i7=1,nT
          xx(i6)=xx(i6)+EI(i7,i6)*gam(i7)
         enddo
         ywk1(i6)=yk(i6)-xx(i6)
        enddo
         do i8=1,idim
          do i9=1,idim
           ywk2(i8)=ywk2(i8)+SigIinv(i8,i9)*ywk1(i9)
          enddo
         enddo
c.......narms(kk) = 2
       if (idim .eq. 2) then
        do i1=1,idim
          do i2=1,idim
           AI2(i1,i2)=0.0d0
           AIinv2(i1,i2)=0.0d0
           AIchol2(i1,i2)=0.0d0
          enddo
        enddo
        do i3=1,idim
          do i4=1,idim
            AI2(i3,i4)=B3(i3,i4)
          enddo
        enddo
        call inverse(AI2, AIinv2, idim)
        do i3=1,idim
          do i4=1,idim
            AIinv2(i3,i4)=alam(kk)*AIinv2(i3,i4)+SigIinv(i3,i4)
          enddo
        enddo
        call inverse(AIinv2, AI2, idim)
c.........Mu_Rgamma
        do j8=1,idim
         do j9=1,idim
          Rgamhat(j8)=Rgamhat(j8)+AI2(j8,j9)*ywk2(j9)
         enddo
        enddo
c.......generate multivariate normal
        lda=idim
        ldr=idim
c        call DCHFAC(lda,AI2,lda,tol,irank,AIchol2,ldr)    
c   
       tmp=0.0d0
       do i1=1,lda
        do j1=1,i1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        AI2vec(ij)=AI2(i1, j1)       
        enddo
       enddo
       bn=real(lda)*(real(lda)+1.0d0)*0.5d0
       nn=int(bn)
       call cholesky(AI2vec, lda, nn, upper2, nullty, ifault) 
c       
       tmp=0.0d0
       do j1=1,lda
        do i1=1,j1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        AIchol2(i1, j1)=upper2(ij)       
        enddo
       enddo     
c
        do i1=1,idim
         r1(i1)=r8_normal_01_sample(iseed)
        enddo
        do i1=1,idim
         do j1=1,idim
         RV(1,i1)=RV(1,i1) + AIchol2(j1,i1)*r1(j1)
         enddo
        enddo             
       endif
c.......narms(kk) = 3
       if (idim .eq. 3) then
        do i1=1,idim
          do i2=1,idim
           AI3(i1,i2)=0.0d0
           AIinv3(i1,i2)=0.0d0
           AIchol3(i1,i2)=0.0d0
          enddo
        enddo
        do i3=1,idim
          do i4=1,idim
            AI3(i3,i4)=B3(i3,i4)
          enddo
        enddo
        call inverse(AI3, AIinv3, idim)
        do i3=1,idim
          do i4=1,idim
            AIinv3(i3,i4)=alam(kk)*AIinv3(i3,i4)+SigIinv(i3,i4)
          enddo
        enddo
        call inverse(AIinv3, AI3, idim)
c.........Mu_Rgamma
        do j8=1,idim
         do j9=1,idim
          Rgamhat(j8)=Rgamhat(j8)+AI3(j8,j9)*ywk2(j9)
         enddo
        enddo
c.......generate multivariate normal
        lda=idim
        ldr=idim
c        call DCHFAC(lda,AI3,lda,tol,irank,AIchol3,ldr)
c   
       tmp=0.0d0
       do i1=1,lda
        do j1=1,i1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        AI3vec(ij)=AI3(i1, j1)       
        enddo
       enddo
       bn=real(lda)*(real(lda)+1.0d0)*0.5d0
       nn=int(bn)
       call cholesky(AI3vec, lda, nn, upper3, nullty, ifault) 
c       
       tmp=0.0d0
       do j1=1,lda
        do i1=1,j1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        AIchol3(i1, j1)=upper3(ij)       
        enddo
       enddo     
c
        do i1=1,idim
         r1(i1)=r8_normal_01_sample(iseed)
        enddo
        do i1=1,idim
         do j1=1,idim
          RV(1,i1)=RV(1,i1) + AIchol3(j1,i1)*r1(j1)
         enddo
        enddo 
       endif
c.......narms(kk) = 4
       if (idim .eq. 4) then
        do i1=1,idim
          do i2=1,idim
           AI4(i1,i2)=0.0d0
           AIinv4(i1,i2)=0.0d0
           AIchol4(i1,i2)=0.0d0
          enddo
        enddo
        do i3=1,idim
          do i4=1,idim
            AI4(i3,i4)=B3(i3,i4)
          enddo
        enddo
        call inverse(AI4, AIinv4, idim)
        do i3=1,idim
          do i4=1,idim
            AIinv4(i3,i4)=alam(kk)*AIinv4(i3,i4)+SigIinv(i3,i4)
          enddo
        enddo
        call inverse(AIinv4, AI4, idim)
c.........Mu_Rgamma
        do j8=1,idim
         do j9=1,idim
          Rgamhat(j8)=Rgamhat(j8)+AI4(j8,j9)*ywk2(j9)
         enddo
        enddo
c.......generate multivariate normal
        lda=idim
        ldr=idim
c        call DCHFAC(lda,AI4,lda,tol,irank,AIchol4,ldr)
c   
       tmp=0.0d0
       do i1=1,lda
        do j1=1,i1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        AI4vec(ij)=AI4(i1, j1)       
        enddo
       enddo
       bn=real(lda)*(real(lda)+1.0d0)*0.5d0
       nn=int(bn)
       call cholesky(AI4vec, lda, nn, upper4, nullty, ifault) 
c       
       tmp=0.0d0
       do j1=1,lda
        do i1=1,j1
        tmp=tmp + 1.0d0
        ij=int(tmp)
        AIchol4(i1, j1)=upper4(ij)       
        enddo
       enddo     
c
        do i1=1,idim
         r1(i1)=r8_normal_01_sample(iseed)
        enddo
        do i1=1,idim
         do j1=1,idim
          RV(1,i1)=RV(1,i1) + AIchol4(j1,i1)*r1(j1)
         enddo
        enddo 
       endif
c
        if (kk .eq. 1) then
          isum=0
        else
          do jj=2,kk
           isum=isum + narms(jj-1)
          enddo
        endif
        do j=1,idim
          Rgam(isum +j)=Rgamhat(j)+RV(1,j)
        enddo
 557   continue
       return
       end


       subroutine hpd(n,alpha,x,alow,aupp)
c      computing 100(1-alpha)% HPD and credible intervals for x 
c      use Chen-Shao HPD Estimation Algorithm 
c      (see page 219 of Monte Carlo Methods in Bayesian Computation, Springer-Verlag, 2000) 
c      ming-hui chen
c      july 23, 2001 at wpi
c      input:
c        alpha: confidence level,  0 < alpha < 1
c        n = mcmc sample size
c        x(n): a univariate vector of mcmc sample 
c      output: 
c        (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c        (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible interval
c
       implicit double precision (a-h,o-z)
       double precision x(n)
       double precision aupp(2),alow(2)

       q1=(alpha/2.0d0)*real(n)
       q2=(1.0d0-alpha/2.0d0)*real(n)
       nq1=int(q1+0.5d0)
       nq2=int(q2+0.5d0)
       nq=nq2-nq1
       do 145 i=1,n-1
        do 178 j=i+1,n
         if (x(i) .gt. x(j)) then
          temp=x(i)
          x(i)=x(j)
          x(j)=temp
         endif
 178    continue
 145   continue
       do 120 j=1,n-nq
        pdiff1=x(j)
        pdiff2=x(j+nq)
        wb=pdiff2-pdiff1
        if (j .eq. 1) then
         whpd=wb
         aupp1=pdiff2
         alow1=pdiff1
        else
         if (whpd .gt. wb) then
          whpd=wb
          aupp1=pdiff2
          alow1=pdiff1
         endif
        endif
 120   continue
       alow(1)=alow1       
       aupp(1)=aupp1       
       alow(2)=x(nq1)
       aupp(2)=x(nq2)

       return
       end
      
c This file contains two versions of the Nelder & Mead simplex algorithm
c for function minimization.   The first is that published in the journal
c of Applied Statistics.   This does not include the fitting of a quadratic
c surface, which provides the only satisfactory method of testing whether
c a minimum has been found.   The search for a minimum is liable to
c premature termination.
c The second version is one which has been developed jointly by CSIRO and
c Rothamsted, and does include the quadratic-surface fitting.
c
c
      subroutine nelmin(fn, n, start, xmin, ynewlo, reqmin, step,
     #     konvge, kcount, icount, numres, ifault)
      implicit double precision (a-h,o-z)
c
c     Simplex function minimisation procedure due to Nelder+Mead(1965),
c     as implemented by O'Neil
c
c     The Nelder-Mead Simplex Minimisation Procedure
c
c
c        Purpose :: To find the minimum value of a user-specified 
c                   function
c
c        Formal parameters ::
c         
c            fn :        : The name of the function to be minimized.
c             n :  input : The number of variables over which we are
c                        : minimising
c         start :  input : Array; Contains the coordinates of the
c                          starting point.
c                 output : The values may be overwritten.
c          xmin : output : Array; Contains the coordinates of the
c                        : minimum.
c        ynewlo : output : The minimum value of the function.
c        reqmin :  input : The terminating limit for the variance of
c                        : function values.
c          step :  input : Array; Determines the size and shape of the
c                        : initial simplex.  The relative magnitudes of
c                        : its n elements should reflect the units of
c                        : the n variables.
c        konvge :  input : The convergence check is carried out every
c                        : konvge iterations.
c        kcount :  input : Maximum number of function evaluations.
c        icount : output : Function evaluations performed
c        numres : output : Number of restarts.
c        ifault : output : 1 if reqmin, n, or konvge has illegal value;
c                        : 2 if terminated because kcount was exceeded
c                        :   without convergence;
c                        : 0 otherwise.
c
c        All variables and arrays are to be declared in the calling
c        program as double precision.
c
c
c        Auxiliary algorithm :: The double precision function
c        subprogram fn(a) calculates the function value at point a.
c        a is double precision with n elements.
c
c
c        Reference :: Nelder,J.A. and Mead,R.(1965).  A simplex method
c        for function minimization.  Computer J., Vol.7,308-313.
c
c************************************************************************
c
      double precision start(n), xmin(n), ynewlo, reqmin, step(n),
     1   p(20,21), pstar(20), p2star(20), pbar(20), y(21),
     2   dn, dnn, z, ylo, rcoeff, ystar, ecoeff, y2star, ccoeff,
     3   rq, x, del, fn, one, half, zero, eps
      external fn
c
      data rcoeff/1.0d0/, ecoeff/2.0d0/, ccoeff/5.0d-1/
      data one/1.0d0/, half/0.5d0/, zero/0.0d0/, eps/0.001d0/
c        reflection, extension and contraction coefficients.
c
c       validity checks on input parameters.
c
      ifault=1
      if(reqmin .le. zero .or. n .lt. 1 .or. n .gt. 20
     #   .or. konvge .lt. 1) return
      ifault=2
      icount=0
      numres=0
c
      jcount=konvge                                                         
      dn=float(n)                                                          
      nn=n+1                                                                
      dnn=float(nn)                                                        
      del=one
      rq=reqmin*dn
c
c        construction of initial simplex.                                   
c
   10 do i=1,n                                                            
       p(i,nn)=start(i)  
      enddo                                                    
      y(nn)=fn(start)
      do 40 j=1,n                                                            
        x=start(j)
        start(j)=start(j)+step(j)*del                                         
        do i=1,n                                                            
         p(i,j)=start(i)   
        enddo                                                    
        y(j)=fn(start)
        start(j)=x
   40 continue
      icount=icount+nn
c                                                                           
c       simplex construction complete                                       
c                                                                           
c       find highest and lowest y values.  ynewlo (=y(ihi) ) indicates       
c       the vertex of the simplex to be replaced.                           
c                                                                           
   43 ylo=y(1)
      ilo=1   
      do 47 i=2,nn
        if(y(i).ge.ylo) goto 47
        ylo=y(i)                                                              
        ilo=i                                                                 
   47 continue
   50 ynewlo=y(1)
      ihi=1
      do 70 i=2,nn
        if(y(i) .le. ynewlo) goto 70
        ynewlo=y(i)
        ihi=i                                                                 
   70 continue
c
c      calculate pbar,the centroid of the simplex vertices                  
c          excepting that with y value ynewlo.                              
c
      do 90 i=1,n                                                            
        z=zero
        do j=1,nn                                                           
         z=z+p(i,j)
        enddo
        z=z-p(i,ihi)                                                          
        pbar(i)=z/dn                                                          
   90 continue
c
c      reflection through the centroid                                      
c
      do i=1,n
       pstar(i)=pbar(i) + rcoeff * (pbar(i) - p(i,ihi))
      enddo
      ystar=fn(pstar)
      icount=icount+1                                                       
      if (ystar.ge.ylo) goto 140
c
c      successful reflection,so extension                                   
c
       do i=1,n
        p2star(i)=pbar(i) + ecoeff * (pstar(i)-pbar(i))
       enddo
      y2star=fn(p2star)
      icount=icount+1                                                       
c
c       check extension
c
      if(y2star .ge. ystar) goto 133
c
c       retain extension or contraction.                                    
c
      do i=1,n
       p(i,ihi)=p2star(i) 
      enddo                                                   
      y(ihi)=y2star                                                         
      goto 230
c
c     retain reflection
c
  133 do i=1,n
       p(i,ihi)=pstar(i)
      enddo
      y(ihi)=ystar
      goto 230
c
c     no extension
c
  140 l=0
      do 150 i=1,nn
        if (y( i).gt.ystar) l=l+1                                             
  150 continue                                                              
      if (l.gt.1) goto 133
      if (l.eq.0) goto 170
c
c     contraction on the reflection side of the centroid.                   
c
      do i=1,n
       p2star(i) = pbar(i) + ccoeff * (pstar(i) - pbar(i))
      enddo
      y2star = fn(p2star)
      icount=icount+1
      if(y2star .le. ystar) goto 182
c
c        retain reflection
c
      do i=1,n
       p(i,ihi)=pstar(i)
      enddo
      y(ihi)=ystar                                                          
      goto 230
c
c      contraction on the  y(ihi) side of the centriod.                     
c
  170 do i=1,n
       p2star(i)=pbar(i) + ccoeff * (p(i,ihi) - pbar(i))
      enddo
      y2star=fn(p2star)
      icount=icount+1                                                       
      if (y2star .gt. y(ihi)) goto 188
c
c        retain contraction
c
  182 do i=1,n
       p(i,ihi) = p2star(i)
      enddo
      y(ihi)=y2star
      goto 230
c
c       contract whole simplex.                                             
c
  188 do 200 j=1,nn
        do 190 i=1,n
          p(i,j)=(p(i,j)+p(i,ilo))*half
          xmin(i)=p(i,j) 
  190   continue
        y(j)=fn(xmin)
  200 continue
      icount=icount+nn
      if (icount .gt. kcount) go to 260
      goto 43
c
c        Check if ylo improved
c
  230 if (y(ihi) .ge. ylo) goto 235
      ylo=y(ihi)
      ilo=ihi
  235 jcount=jcount-1
      if(jcount .ne. 0) goto 50
c
c     check to see if minimum reached.                                      
c
      if (icount.gt.kcount) goto 260
      jcount=konvge
      z=zero
      do i=1, nn
       z = z+y(i)
      enddo
      x=z / dnn
      z=zero
      do i=1,nn
       z = z+(y(i)-x) ** 2
      enddo
      if (z .gt. rq) goto 50
c
c       factorial tests to check that ynewlo is a local minimum.             
c
  260 do i=1,n
       xmin(i)=p(i,ilo)
      enddo
      ynewlo=y(ilo)
      if (icount.gt.kcount) return
      do 280 i=1,n
        del=step(i)*eps
        xmin(i)=xmin(i)+del
        z=fn(xmin)
        icount=icount+1
        if (z.lt.ynewlo) goto 290
        xmin(i)=xmin(i)-del-del
        z=fn(xmin)
        icount=icount+1
        if (z.lt.ynewlo) goto 290
        xmin(i)=xmin(i)+del
  280 continue
      ifault = 0
      return
c
c     restart procedure
c
  290 do i=1,n
       start(i) = xmin(i)
      enddo
      del=eps
      numres = numres + 1
      goto 10
      end
c
      function gengam ( a , iseed)
c*********************************************************************72
c
cc gengam samples the standard Gamma distribution.
c
c  Discussion:
c
c    This procedure corresponds to algorithm GD in the reference.
c
c    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 April 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Generating Gamma Variates by a Modified Rejection Technique,
c    Communications of the ACM,
c    Volume 25, Number 1, January 1982, pages 47-54.
c
c  Parameters:
c
c    Input, double precision A, the shape parameter of the standard gamma
c    distribution. 0 .lt. A.
c
c    Output, double precision gengam, a random deviate
c    from the distribution.
c
      implicit double precision (a-h,o-z)
      integer iseed
      double precision a
      double precision a1
      parameter ( a1 =  0.3333333D+00 )
      double precision a2
      parameter ( a2 = -0.2500030D+00 )
      double precision a3
      parameter ( a3 =  0.2000062D+00 )
      double precision a4
      parameter ( a4 = -0.1662921D+00 )
      double precision a5
      parameter ( a5 =  0.1423657D+00 )
      double precision a6
      parameter ( a6 = -0.1367177D+00 )
      double precision a7
      parameter ( a7 =  0.1233795D+00 )
      double precision b
      double precision c
      double precision d
      double precision e
      double precision, parameter :: e1 = 1.0D+00
      double precision, parameter :: e2 = 0.4999897D+00
      double precision, parameter :: e3 = 0.1668290D+00
      double precision, parameter :: e4 = 0.0407753D+00
      double precision, parameter :: e5 = 0.0102930D+00
      double precision p
      double precision q
      double precision q0
      double precision, parameter :: q1 =  0.04166669D+00
      double precision, parameter :: q2 =  0.02083148D+00
      double precision, parameter :: q3 =  0.00801191D+00
      double precision, parameter :: q4 =  0.00144121D+00
      double precision, parameter :: q5 = -0.00007388D+00
      double precision, parameter :: q6 =  0.00024511D+00
      double precision, parameter :: q7 =  0.00024240D+00
      double precision r
      double precision r8_exponential_01_sample
      double precision gengam
      double precision r8_normal_01_sample
      double precision r8_uniform_01_sample
      double precision s
      double precision s2
      double precision si
      double precision, parameter :: sqrt32 = 5.656854D+00
      double precision t
      double precision u
      double precision v
      double precision w
      double precision x

      if ( 1.0D+00 .le. a ) then

        s2 = a - 0.5D+00
        s = sqrt ( s2 )
        d = sqrt32 - 12.0D+00 * s
c
c  Immediate acceptance.
c
        t = r8_normal_01_sample ( iseed )
        x = s + 0.5D+00 * t
        gengam = x * x

        if ( 0.0D+00 .le. t ) then
          return
        end if
c
c  Squeeze acceptance.
c
        u = r8_uniform_01_sample ( iseed )
        if ( d * u .le. t * t * t ) then
          return
        end if

        r = 1.0D+00 / a
        q0 = (((((( q7 
     &    * r + q6 ) 
     &    * r + q5 ) 
     &    * r + q4 ) 
     &    * r + q3 ) 
     &    * r + q2 ) 
     &    * r + q1 ) 
     &    * r
c
c  Approximation depending on size of parameter A.
c
        if ( 13.022D+00 .lt. a ) then
          b = 1.77D+00
          si = 0.75D+00
          c = 0.1515D+00 / s
        else if ( 3.686D+00 .lt. a ) then
          b = 1.654D+00 + 0.0076D+00 * s2
          si = 1.68D+00 / s + 0.275D+00
          c = 0.062D+00 / s + 0.024D+00
        else
          b = 0.463D+00 + s + 0.178D+00 * s2
          si = 1.235D+00
          c = 0.195D+00 / s - 0.079D+00 + 0.16D+00 * s
        end if
c
c  Quotient test.
c
        if ( 0.0D+00 .lt. x ) then

          v = 0.5D+00 * t / s

          if ( 0.25D+00 .lt. abs ( v ) ) then
            q = q0 - s * t + 0.25D+00 * t * t 
     &        + 2.0D+00 * s2 * log ( 1.0D+00 + v )
          else
            q = q0 + 0.5D+00 * t * t * (((((( a7 
     &        * v + a6 ) 
     &        * v + a5 ) 
     &        * v + a4 ) 
     &        * v + a3 ) 
     &        * v + a2 ) 
     &        * v + a1 ) 
     &        * v
          end if

          if ( log ( 1.0D+00 - u ) .le. q ) then
            return
          end if

        end if

 13    continue

          e = r8_exponential_01_sample ( iseed )
          u = 2.0D+00 * r8_uniform_01_sample ( iseed ) - 1.0D+00

          if ( 0.0D+00 .le. u ) then
            t = b + abs ( si * e )
          else
            t = b - abs ( si * e )
          end if
c
c  Possible rejection.
c
          if ( t .lt. -0.7187449D+00 ) then
            go to 13
          end if
c
c  Calculate V and quotient Q.
c
          v = 0.5D+00 * t / s

          if ( 0.25D+00 .lt. abs ( v ) ) then
            q = q0 - s * t + 0.25D+00 * t * t 
     &        + 2.0D+00 * s2 * log ( 1.0D+00 + v )
          else
            q = q0 + 0.5D+00 * t * t * (((((( a7 
     &        * v + a6 ) 
     &        * v + a5 ) 
     &        * v + a4 ) 
     &        * v + a3 ) 
     &        * v + a2 ) 
     &        * v + a1 ) 
     &        * v
          end if
c
c  Hat acceptance.
c
          if ( q .le. 0.0D+00 ) then
            go to 13
          end if

          if ( 0.5D+00 .lt. q ) then
            w = exp ( q ) - 1.0D+00
          else
            w = (((( e5 * q + e4 ) * q + e3 ) * q + e2 ) * q + e1 ) * q
          end if
c
c  May have to sample again.
c
          if ( c * abs ( u ) .le. w * exp ( e - 0.5D+00 * t * t ) ) then
            go to 23
          end if

        go to 13

 23     continue

        x = s + 0.5D+00 * t
        gengam = x * x

        return
c
c  Method for A .lt. 1.
c
      else

        b = 1.0D+00 + 0.3678794D+00 * a

 33     continue

          p = b * r8_uniform_01_sample ( iseed )

          if ( p .lt. 1.0D+00 ) then

            gengam = exp ( log ( p ) / a )

            if ( gengam .le. 
     &        r8_exponential_01_sample ( iseed ) ) then
              return
            end if

            go to 33

          end if

          gengam = - log ( ( b - p ) / a )

          if ( ( 1.0D+00 - a ) * log ( gengam ) .le. 
     &      r8_exponential_01_sample ( iseed ) ) then
            go to 40
          end if

        go to 33

40      continue

      end if

      return
      end



      function r8_exponential_01_sample ( iseed )

c*********************************************************************72
c
cc R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 April 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
c
      implicit double precision (a-h,o-z)
      integer iseed
      double precision r
      double precision r8_exponential_01_sample
      double precision r8_uniform_01_sample

      r = r8_uniform_01_sample ( iseed )

      r8_exponential_01_sample = - log ( r )

      return
      end
      
      function r8_normal_01_sample ( iseed )

c*********************************************************************72
c
cc R8_NORMAL_01_SAMPLE returns a unit pseudonormal R8.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    The Box-Muller method is used, which is efficient, but
c    generates two values at a time.
c
c    Typically, we would use one value and save the other for the next call.
c    However, the fact that this function has saved memory makes it difficult
c    to correctly handle cases where we want to re-initialize the code,
c    or to run in parallel.  Therefore, we will instead use the first value
c    and DISCARD the second.
c
c    EFFICIENCY must defer to SIMPLICITY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_NORMAL_01_SAMPLE, a sample of the standard
c    normal PDF.
c
      implicit double precision (a-h,o-z)
      integer iseed
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01_sample
      double precision r8_uniform_01_sample
      double precision x

      r1 = r8_uniform_01_sample ( iseed )
      r2 = r8_uniform_01_sample ( iseed )

      x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )

      r8_normal_01_sample = x

      return
      end
      
      function r8_uniform_01_sample ( iseed )
c*********************************************************************72
c
cc R8_UNIFORM_01_SAMPLE generates a uniform random deviate from [0,1].
c
c  Discussion:
c
c    This function should be the only way that the package accesses random
c    numbers.
c
c    Setting OPTION to 0 accesses the R8_UNI_01() function in RNGLIB,
c    for which there are versions in various languages, which should result
c    in the same values being returned.  This should be the only function
c    that calls a function in RNGLIB.
c
c    Setting OPTION to 1 in the FORTRAN90 version calls the system
c    RNG "random_number()".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 August 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Output, double precision R8_UNIFORM_01_SAMPLE, a random deviate
c    from the distribution.
      implicit double precision (a-h,o-z)
      integer na,np,iseed,nb15,nb16,nxhi,nxalo,leftlo,nfhi,nk      
      data na/16807/,nb15/32768/,nb16/65536/,np/2147483647/
      nxhi=int(real(iseed)/real(nb16))
      nxalo=(iseed-nxhi*nb16)*na
      leftlo=int(real(nxalo)/real(nb16))
      nfhi=nxhi*na+leftlo
      nk=int(real(nfhi)/real(nb15))
      iseed=(((nxalo-leftlo*nb16)-np)+(nfhi-nk*nb15)*nb16)+nk  
      if(iseed.lt.0) iseed=iseed+np
      r8_uniform_01_sample=real(iseed)*4.656612875d-10
      return
      end

          subroutine inverse(w,c,n)
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed 
        ! during the calculation
        !===========================================================
        implicit double precision (a-h,o-z)
        integer n
        double precision a(n,n), c(n,n), w(n,n)
        double precision L(n,n), U(n,n), b(n), d(n), x(n)
        double precision coeff
        integer i, j, k
        
        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0
        do i1=1,n
         do j1=1,n
          a(i1,j1)=w(i1,j1)
         enddo
        enddo
        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do
        
        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do
        
        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            c(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        end subroutine inverse
        
        !*****************************************************
        !* Calculate the determinant of a real square matrix *
        !* A(n,n) by Gauss method with full pivoting.        *
        !* ------------------------------------------------- *
        !* Ref.: "Algèbre - Algorithmes et programmes en     *
        !*        Pascal By Jean-Louis Jardrin, Dunod -      *
        !*        Bordas Paris, 1988 p. 76-79".              *
        !* ------------------------------------------------- *
        !* SAMPLE RUN:                                       *
        !* (Calculate the determinant of matrix:             *
        !*            10 18  1 14 22                         *
        !*             4 12 25  8 16                         *
        !*            23  6 19  2 15                         *
        !*            17  5 13 21  9                         *
        !*            11 24  7 20  3  )                      *
        !*                                                   *
        !* Input size of square real matrix: 5               *
        !*                                                   *
        !* Line 1                                            *
        !* Element 1: 10                                     *
        !* Element 2: 18                                     *
        !* Element 3: 1                                      *
        !* Element 4: 14                                     *
        !* Element 5: 22                                     *
        !*                                                   *
        !* Line 2                                            *
        !* Element 1: 4                                      *
        !* Element 2: 12                                     *
        !* Element 3: 25                                     *
        !* Element 4:  8                                     *
        !* Element 5: 16                                     *
        !*                                                   *
        !* Line 3                                            *
        !* Element 1: 23                                     *
        !* Element 2:  6                                     *
        !* Element 3: 19                                     *
        !* Element 4:  2                                     *
        !* Element 5: 15                                     *
        !*                                                   *
        !* Line 4                                            *
        !* Element 1: 17                                     *
        !* Element 2:  5                                     *
        !* Element 3: 13                                     *
        !* Element 4: 21                                     *
        !* Element 5:  9                                     *
        !*                                                   *
        !* Line 5                                            *
        !* Element 1: 11                                     *
        !* Element 2: 24                                     *
        !* Element 3:  7                                     *
        !* Element 4: 20                                     *
        !* Element 5:  3                                     *
        !*                                                   *
        !* Determinant = -0.468000E+07                       *
        !*                                                   *
        !*             F90 Version with dynamic allocations  *
        !*                 By Jean-Pierre Moreau, Paris.     *
        !*                      (www.jpmoreau.fr)            *
        !*****************************************************

c          eps=2.220446049250313E-14          
c          det=DMGT(eps,n,A)       

        
        
        !The subroutine TSRGT applies to input real square matrix A(n,n) the upper
        !triangularization algorithm of Gauss method with full pivoting and keeps
        !trace of successive transformations done in integer vectors KP and LP.
        !-----------------------------------------------------------------------------
        !  Input parameters:
        !  eps        precision (double precision)
        !  n          size of A matrix (integer)
        !  A          pointer to input real square matrix (double precision)
        !  Output parameters:
        !  it         flag=1 if A matrix ok, =0 if A matrix is singular (integer)
        !  C          pointer to table storing main diagonal elements and supra-
        !             diagonal elements of upper triangular matrix and the multi-
        !             plying coefficients used during triangularization process
        !  KP         table storing informations concerning the column exchanges
        !             during process (integer)
        !  LP         table storing informations concerning the line exchanges
        !             during process (integer)
        !-----------------------------------------------------------------------------
        !The table C is first initialized to A matrix, then receives at each step k
        !of the triangularization process, usefull elements of A matrix at step k for
        !k=1,2,...n.
        !The variables po(double precision), lo and ko(integer) store respectively pivot at step k,
        !its line number and its column number.
        !------------------------------------------------------------------------------
        Subroutine TSRGT(eps, n, A, it, C, Kp, Lp)
        implicit double precision (a-h,o-z)        
          double precision eps
          integer n,it
          double precision A(n,n), C(n,n)
          integer Kp(n),Lp(n) 
          double precision  po,t0
          C=A; it=1; k=1
          do while (it==1.and.k<n)
            po=C(k,k); lo=k; ko=k
            do i=k, n
              do j=k, n
                if (dabs(C(i,j))>dabs(po)) then
                  po=C(i,j); lo=i; ko=j
                end if
              end do
            end do
            Lp(k)=lo; Kp(k)=ko
            if (dabs(po)<eps) then
              it=0
            else
              if (lo.ne.k) then
                do j=k, n
                  t0=C(k,j); C(k,j)=C(lo,j); C(lo,j)=t0
                end do
              end if
              if (ko.ne.k) then
                do i=1, n
                  t0=C(i,k); C(i,k)=C(i,ko); C(i,ko)=t0
                end do
              end if 
              do i=k+1, n
                C(i,k)=C(i,k)/po
                do j=k+1, n
                  C(i,j)=C(i,j)-C(i,k)*C(k,j)
                end do 
              end do
              k=k+1
            end if
          end do
          if (it==1.and.dabs(C(n,n))<eps)  it=0
          return
        End !TSRGT
        
        !The function DMGT returns the determinant of a real square matrix
        !A(n,n) by Gauss method with full pivoting.
        !----------------------------------------------------------------------------
        !  Input parameters:
        !  eps        precision (double precision)
        !  n          size of A matrix (integer)
        !  A          pointer to input real square matrix
        !  Output parameters:
        !  None
        !-----------------------------------------------------------------------------
        !The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
        !Output variables are it(integer), C(n,n), Kp(n) and Lp(n).
        !If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
        !at location i,j (j>=i) the corresponding element of the upper triangular matrix.
        !Tables Lp and Kp contain informations relative to exchanges of line or column
        !that occured during the process. For instance, the element number k of Lp is
        !an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
        !The number of exchanges of lines and columns is stored in integer L. the
        !determinant of A matrix is stored in d0 (double precision).
        !-----------------------------------------------------------------------------
        Function DMGT(eps, n, A)
        implicit double precision (a-h,o-z)        
        integer n
        double precision eps, A(n,n)
        double precision d0
c        
        double precision C(n,n)
        integer Kp(n), Lp(n)
        call TSRGT(eps,n,A,it,C,Kp,Lp)  !call triangularization subroutine
          if (it==0) then
            d0=0.d0  !matrix singular, det=0
          else       !matrix regular, det<>0
            d0=1.d0
            do k=1, n
              d0=d0*C(k,k)
            end do
            l=0
            do k=1, n-1
              if (Lp(k).ne.k)  l=l+1
              if (Kp(k).ne.k)  l=l+1
            end do
            if (MOD(l,2).ne.0) d0=-d0  !l is odd
          end if
          DMGT=d0   !return determinant
          return
        End
        
        !End of file deter.f90
        
      subroutine cholesky ( a, n, nn, u, nullty, ifault )

c*********************************************************************72
c
cc CHOLESKY computes the Cholesky factorization of a PDS matrix.
c
c  Discussion:
c
c    For a positive definite symmetric matrix A, the Cholesky factor U
c    is an upper triangular matrix such that A = U' * U.
c
c    This routine was originally named "CHOL", but that conflicted with
c    a built in MATLAB routine name.
c
c    The missing initialization "II = 0" has been added to the code.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    Michael Healy
c    Modifications by AJ Miller.
c    Modifications by John Burkardt
c
c  Reference:
c
c    Michael Healy,
c    Algorithm AS 6:
c    Triangular decomposition of a symmetric matrix,
c    Applied Statistics,
c    Volume 17, Number 2, 1968, pages 195-197.
c
c  Parameters:
c
c    Input, double precision A((N*(N+1))/2), a positive definite matrix 
c    stored by rows in lower triangular form as a one dimensional array, 
c    in the sequence
c    A(1,1),
c    A(2,1), A(2,2),
c    A(3,1), A(3,2), A(3,3), and so on.
c
c    Input, integer N, the order of A.
c
c    Input, integer NN, the dimension of the array used to store A, 
c    which should be at least (N*(N+1))/2.
c
c    Output, double precision U((N*(N+1))/2), an upper triangular matrix,
c    stored by columns, which is the Cholesky factor of A.  The program is
c    written in such a way that A and U can share storage.
c
c    Output, integer NULLTY, the rank deficiency of A.  If NULLTY is zero,
c    the matrix is judged to have full rank.
c
c    Output, integer IFAULT, an error indicator.
c    0, no error was detected;
c    1, if N < 1;
c    2, if A is not positive semi-definite.
c    3, if NN < (N*(N+1))/2.
c
c  Local Parameters:
c
c    Local, double precision ETA, should be set equal to the smallest positive 
c    value such that 1.0 + ETA is calculated as being greater than 1.0 in the 
c    accuracy being used.
c
       implicit double precision (a-h,o-z)
      double precision a(nn)
      double precision eta
      parameter ( eta = 1.0D-09 )
      integer n
      integer ifault
      integer nn
      integer nullty
      double precision u(nn)
      double precision w
      double precision x
c
      ifault = 0
      nullty = 0
c
      if ( n .le. 0 ) then
        ifault = 1
        return
      end if

      if ( nn .lt. ( n * ( n + 1 ) ) / 2 ) then
        ifault = 3
        return
      end if

      j = 1
      k = 0
      ii = 0
c
c  Factorize column by column, ICOL = column number.
c
      do icol = 1, n

        ii = ii + icol
        x = eta * eta * a(ii)
        l = 0
        kk = 0
c
c  IROW = row number within column ICOL.
c
        do irow = 1, icol

          kk = kk + irow
          k = k + 1
          w = a(k)
          m = j

          do i = 1, irow - 1
            l = l + 1
            w = w - u(l) * u(m)
            m = m + 1
          end do

          l = l + 1

          if ( irow .eq. icol ) then
            go to 27
          end if

          if ( u(l) .ne. 0.0D+00 ) then

            u(k) = w / u(l)

          else

            u(k) = 0.0D+00

            if ( abs ( x * a(k) ) .lt. w * w ) then
              ifault = 2
              return
            end if

          end if

        end do
c
c  End of row, estimate relative accuracy of diagonal element.
c
  27   continue

        if ( abs ( w ) .le. abs ( eta * a(k) ) ) then

          u(k) = 0.0D+00
          nullty = nullty + 1
       
        else

          if ( w .lt. 0.0D+00 ) then
            ifault = 2
            return
          end if

          u(k) = sqrt ( w )

        end if

        j = j + icol

      end do

      return
      end                           





       
