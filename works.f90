    program works


    implicit none
    integer N,step,ij, sc,i,j,tt,k,scc,ph
    integer K_max,K_mean,Length
    parameter (N=10000,Length=100,step=100000,sc=50)
    integer L(Length+2,Length+2),nei(N,4),ndeg(N)

    integer ns(N),nei2(N,8),s1(N),s2(N),nc(N)

    real(kind=8) z,x,y,fc,fcs,fcss,C,D,com,Q1,Q2,gamma,alpha,beta,epsilon,r,mad,lambda,b,etal,etag,et,es,ess
    real(kind=8) s(N),reward(N),pay(N),Q(N,4,2),mem(N,9),env(N,9),e(N),fcl(N)
    character(len=100)::filename,filename1,filename2,filename3

    alpha=0.1d0
    gamma=0.9d0
    epsilon=0.1d0
    lambda=0.8d0
    beta=0.5d0

    write(filename,'(''beta'',f4.2)')beta
    open(11,file="fc\"//trim(adjustl(filename))//".txt",status="unknown")
 open(12,file="e\"//trim(adjustl(filename))//".txt",status="unknown")
 open(13,file="scfc\"//trim(adjustl(filename))//".txt",status="unknown")
 open(14,file="sce\"//trim(adjustl(filename))//".txt",status="unknown")
    do k=0,26
        r=3.2d0+k/20d0
        !b=1.0+l/100d0
        call regular4(L)
        do i=2,Length+1
            do j=2,Length+1
                nei(L(i,j),1)=L(i,j-1)
                nei(L(i,j),2)=L(i,j+1)
                nei(L(i,j),3)=L(i-1,j)
                nei(L(i,j),4)=L(i+1,j)
                nei2(L(i,j),1)=L(i-1,j-1)
                nei2(L(i,j),2)=L(i-1,j)
                nei2(L(i,j),3)=L(i-1,j+1)
                nei2(L(i,j),4)=L(i,j-1)
                nei2(L(i,j),5)=L(i,j+1)
                nei2(L(i,j),6)=L(i+1,j-1)
                nei2(L(i,j),7)=L(i+1,j)
                nei2(L(i,j),8)=L(i+1,j+1)
            enddo
        enddo
        fcss=0d0
        ess=0d0
        do scc=1,sc
            write(filename1,*)k
            write(filename2,*)scc

            call random_seed()
            fcs=0d0
            es=0d0
            !===================set initial strategy
            Q=0.0d0
            e=epsilon
            do i=1,N
                call random_number(z)
                if(z<0.50d0)then
                    ns(i)=1
                else
                    ns(i)=0
                endif
            enddo
            do i=1,N
                mem(i,9)=ns(i)
                do j=1,8
                    mem(i,j)=ns(nei2(i,j))
                enddo
            enddo
            s2=1

            do tt=1,step
                reward=0d0
                !================================策略更新
                fcl=0d0
                s1=s2
                do i=1,n
                    C=Q(i,s1(i),2)
                    D=Q(i,s1(i),1)
                    call random_number(z)
                    if (z<e(i).or.c==d) then
                        call random_number(y)
                        if (y<0.5d0) then
                            ns(i)=1
                        else
                            ns(i)=0
                        end if
                    else
                        if (c>d) then
                            ns(i)=1
                        else
                            ns(i)=0
                        endif
                    endif

                enddo
                fc=real(sum(ns))/N
                
                !======================PDG
                !            pay=0d0
                !do  i=1,n
                !                        do ij=1,4
                !                            j=nei(i,ij)
                !                            if (ns(i)==1.and.ns(j)==1)then
                !                                pay(i)=pay(i)+1d0
                !                            elseif (ns(i)==1.and.ns(j)==0)then
                !                                pay(i)=pay(i)+1d0-b
                !                            elseif (ns(i)==0.and.ns(j)==1)then
                !                                pay(i)=pay(i)+b
                !                            elseif (ns(i)==0.and.ns(j)==0)then
                !                                pay(i)=pay(i)
                !                            endif
                !                        enddo
                !                enddo
                !!=================================计算PGG收益
                do  i=1,N
                    nc(i)=ns(i)
                    do j=1,4
                        nc(i)=nc(i)+ns(nei(i,j))
                    enddo
                    reward(i)=nc(i)*r/5d0
                enddo
                
                do i=1,N
                    pay(i)=reward(i)
                    fcl(i)=real(nc(i))/5d0
                    do j=1,4
                        pay(i)=pay(i)+reward(nei(i,j))
                        fcl(i)=fcl(i)+real(nc(nei(i,j)))/5d0
                    enddo
                    pay(i)=pay(i)/5d0-ns(i)
                enddo
                
                do i=1,N
                    env(i,9)=ns(i)
                    mad=abs(env(i,9)-mem(i,9))/9d0
                    do j=1,8
                        env(i,j)=ns(nei2(i,j))
                        mad=mad+abs(env(i,j)-mem(i,j))/9d0
                    enddo
                    etal=tanh((1d0-mad*2)*5d0)
                    etag=tanh((fcl(i)-fc)*5d0)
                    !etag=tanh((sum(env(i,:))/9d0-fc)*5d0)
                    if (etal>0.and.etag>0) then
                        s2(i)=1
                    elseif (etal>0.and.etag<0) then
                        s2(i)=2
                    elseif (etal<0.and.etag>0) then
                        s2(i)=3
                    else
                        s2(i)=4
                    endif

                e(i)=beta*(epsilon**(1d0+etal))+(1d0-beta)*(epsilon**(1d0+etag))
                enddo
                mem=lambda*mem+(1d0-lambda)*env
                et=sum(e)/N                
                !=============q表更新
                !s2=1
                do i=1,n
                    Q1=Q(i,s1(i),ns(i)+1)
                    Q2=maxval(Q(i,s2(i),:))
                    Q(i,s1(i),ns(i)+1)=Q1+alpha*(pay(i)+gamma*Q2-Q1)
                enddo



                if (mod(tt,1000)==0) then
                    write(*,*) tt,fc,et
                endif
                if (tt>90000)then
                    fcs=fcs+fc/10000d0
                    es=es+et/10000d0
                endif
            ENDDO
        write (13,*)k,scc,fcs
        write (14,*)k,scc,es

            fcss=fcss+fcs/real(sc)
            ess=ess+es/real(sc)
            
        enddo

        write (11,*)fcss
        write (12,*)ess
    enddo




    contains

    !=========生成规则格子
    subroutine regular4(L)
    implicit none
    integer Length,i,j
    parameter (Length=100)
    integer L(Length+2,Length+2)

    do i=2,Length+1
        do j=2,Length+1
            L(i,j)=(i-2)*Length+(j-1)
        enddo
    enddo

    do j=2,Length+1
        L(1,j)=L(Length+1,j)
        L(Length+2,j)=L(2,j)
    enddo

    do i=2,Length+1
        L(i,1)=L(i,Length+1)
        L(i,Length+2)=L(i,2)
    enddo
    L(1,1)=Length*Length
    L(1,Length+2)=Length*(Length-1)+1
    L(Length+2,1)=Length
    L(Length+2,Length+2)=1
    return
    end






    end program works

