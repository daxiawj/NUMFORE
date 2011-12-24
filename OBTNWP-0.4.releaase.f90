!*********************************************************!
!-------------This Program is BSD Licensed!---------------!
!*********************************************************!
!       Copyleft (c) 2005, Wang Jun @ Atmos Sci NJU       !
!               All rights reserved.                      !
!      For more infomation about the CopyLeft, Please     ! 
!     refer to the files LICENSE.txt and COPYLEFT.txt.    ! 
!*********************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------��ѹԭʼԤ������-------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
implicit none
integer i,j,irec,step,IntHour,IntTime
integer,parameter::n=16,m=20,k=2
integer,parameter::hr_Euler=1,hr_NinSmth=11,hr_TimSmth=6
integer,parameter::hr_Prd=12,hr_IntTime=2
!����һЩ����������ʱҪ�ã������޸ģ���hr_PrdΪһ������ʱ�䣬�˴�Ϊ12Сʱ
!hr_IntTime=2��ʾ�������ڣ��˴���Ϊ24Сʱ
real dt,dd    !���ֲ���������࣬���Ǵ�init.ini�ļ��õ���ֵ
real, parameter::S=0.5
real mij(m,n),fij(m,n),z0(m,n)
real u(m,n,3),v(m,n,3),h(m,n,3)
real du(m,n),dv(m,n),dh(m,n)
open(1,file='rz.txt')
open(2,file='rm.txt')
open(3,file='rf.txt')
do j=1,n
  read(1,*)(z0(i,j),i=1,20)
  read(2,*)(mij(i,j),i=1,20)
  read(3,*)(fij(i,j),i=1,20)
enddo
close(1)
close(2)
close(3)
do j=1,n
  do i=1,m
    z0(i,j)=z0(i,j)*1.0E1
  enddo
enddo
write(6,*)"Staring initializing......"
call Init(dd,dt,m,n,u(1:m,1:n,k),v(1:m,1:n,k),fij,mij,z0,h(1:m,1:n,k))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!�������ԭʼ��
open(1,file='inituvz.dat',form='unformatted',access='direct',recl=n*m)
  write(1,rec=1)((u(i,j,k),i=1,m),j=1,n)
  write(1,rec=2)((v(i,j,k),i=1,m),j=1,n)
  write(1,rec=3)((h(i,j,k),i=1,m),j=1,n)
close(1)
!!!!!!!��ʼ����
u(1:m,1:n,k-1)=u(1:m,1:n,k) !��ʼ������ʱ���
v(1:m,1:n,k-1)=v(1:m,1:n,k)
h(1:m,1:n,k-1)=h(1:m,1:n,k)
u(1:m,1:n,k+1)=u(1:m,1:n,k)
v(1:m,1:n,k+1)=v(1:m,1:n,k)
h(1:m,1:n,k+1)=h(1:m,1:n,k)
step=0
write(6,*)"Starting Integrating......"
open(1,file='date.dat',form='unformatted',access='direct',recl=n*m)
do IntTime=1,hr_IntTime  !һ��Ҫ���м������ڵĻ���
  IntHour=1
  step=0
  do while(int(step*dt)<=hr_Prd*60*60) !һ����������
    step=step+1
    if(int(step*dt)<=60*60)then        !�ж��Ƿ����ŷ����
      call Euler(m,n,u,v,h,du,dv,dh,fij,mij,dd,dt)
    else
      call TriJmp(m,n,u,v,h,du,dv,dh,fij,mij,dd,dt)
    endif
    if(int((step-1)*dt)<hr_TimSmth*60*60.and.int((step+1)*dt)>=hr_TimSmth*60*60)then
      !�ж��Ƿ����ʱ��ƽ��
      call TimSmth(m, n, U, s)
      call TimSmth(m, n, V, s)
      call TimSmth(m, n, H, s)
    elseif(int((step-1)*dt)<2*hr_TimSmth*60*60.and.int((step+1)*dt)>=2*hr_TimSmth*60*60)then
      !�ж��Ƿ����ʱ��ƽ��
      call TimSmth(m, n, U, s)
      call TimSmth(m, n, V, s)
      call TimSmth(m, n, H, s)
    endif
!    if(int((step-1)*dt)<(IntHour-1)*60*60.and.int(step*dt)>=(IntHour-1)*60*60)then
    if(real(step)*real(dt)/(60.0*60.0)>real(IntHour))then
      !�ж��Ƿ��ֻ��ִﵽ��һ��Сʱ
      irec=(IntTime-1)*hr_Prd+IntHour
      write(6,*)"Hour=",IntHour
      write(6,*)"Start the",(IntTime-1)*hr_Prd+IntHour,"Hour(s) Integration...."
      if(IntHour<=hr_NinSmth)then  !�ж��Ƿ�þŵ�ƽ��
        call NinSmth(m, n, U(:,:,k), s)
        call NinSmth(m, n, U(:,:,k-1), s)
        call NinSmth(m, n, U(:,:,k+1), s)
        call NinSmth(m, n, V(:,:,k), s)
        call NinSmth(m, n, V(:,:,k-1), s)
        call NinSmth(m, n, V(:,:,k+1), s)
        call NinSmth(m, n, H(:,:,k), s)
        call NinSmth(m, n, H(:,:,k-1), s)
        call NinSmth(m, n, H(:,:,k+1), s)
      else
        !�ж��Ƿ�Ӧ�ý���һ�����ƽ��
        call FivSmth(m, n, U(:,:,k), s)
        call FivSmth(m, n, U(:,:,k-1), s)
        call FivSmth(m, n, U(:,:,k+1), s)
        call FivSmth(m, n, V(:,:,k), s)
        call FivSmth(m, n, V(:,:,k-1), s)
        call FivSmth(m, n, V(:,:,k+1), s)
        call FivSmth(m, n, H(:,:,k), s)
        call FivSmth(m, n, H(:,:,k-1), s)
        call FivSmth(m, n, H(:,:,k+1), s)
      endif
      !���ÿ��Сʱ�Ļ��ֽ��
      write(1,rec=(irec*3-2))((u(i,j,k),i=1,m),j=1,n)
      write(1,rec=(irec*3-1))((v(i,j,k),i=1,m),j=1,n)
      write(1,rec=(irec*3))((h(i,j,k),i=1,m),j=1,n)
      IntHour=IntHour+1
    endif
  enddo
enddo
close(1)
write(6,*)"The Integration is done!"
end program



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------��ʼ��------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init(dd,dt,m,n,u,v,fij,mij,z0,h)
integer m,n
real dd,dt,u(m,n),v(m,n),fij(m,n),mij(m,n),z0(m,n),h(m,n)
real,parameter::g=9.8
open(1,file='init.ini')
read(1,*)dd,dt
close(1)
!���ʼ��U,V,Hֵ���м��������߽�ֱ���ǰ����
do j=1,n
  do i=1,m
    fij(i,j)=fij(i,j)*1.0E-4
  enddo
enddo
do i=1, m
  do j=2, n-1
    u(i, j)=(-1.0*mij(i, j)/fij(i, j))*g*(z0(i, j+1)-z0(i, j-1))/(2.0*dd)
  enddo
  u(i, 1)=-1.0*mij(i, 1)/fij(i,1)*g*(z0(i, 2)-z0(i, 1))/dd
  u(i, n)=-1.0*mij(i, n)/fij(i,n)*g*(z0(i, n)-z0(i, n-1))/dd
enddo
do j=1, n
  do i=2,m-1
    v(i, j)=mij(i, j)/fij(i, j)*g*(z0(i+1, j)-z0(i-1, j))/(2.0*dd)
  end do
  v(1, j)=mij(1, j)/fij(1, j)*g*(z0(2, j)-z0(1, j))/dd
  v(m, j)=mij(n, j)/fij(m, j)*g*(z0(m, j)-z0(m-1, j))/dd
end do
do j=1,n
  do i=1,m
    h(i, j)=z0(i, j)
  enddo
enddo
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------��������ֵ---------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TrendEGH(m,n,u,v,h,du,dv,dh,fij,mij,dd,dt)
integer m, n
real mij(m,n),fij(m,n),dd,dt
real u(m,n),v(m,n),h(m,n)
real du(m,n),dv(m,n),dh(m,n)
real uu(m,n),uv(m,n),fstar(m,n)
real vv(m,n),vu(m,n)
real hu(m,n),hv(m,n),huv(m,n)
real,parameter::z0=2500.0
real,parameter::g=9.8
k=2
du=0.0; dv=0.0; dh=0.0; fstar=0.0
uu=0.0; uv=0.0; hu=0.0; hv=0.0; huv=0.0
do j=2,n-1
  do i=2,m-1
    !Ϊ���⹫ʽ�����������˼�����ʱ����
    uu(i,j)=(u(i+1,j)**2-u(i-1,j)**2)/(4.0*dd)
    uv(i,j)=((v(i,j+1)+v(i,j))*(u(i,j+1)-u(i,j))+(v(i,j)+v(i,j-1))*(u(i,j)-u(i,j-1)))/(4.0*dd)
    vv(i,j)=(v(i,j+1)**2-v(i,j-1)**2)/(4.0*dd)
    vu(i,j)=((u(i+1,j)+u(i,j))*(v(i+1,j)-v(i,j))+(u(i,j)+u(i-1,j))*(v(i,j)-v(i-1,j)))/(4.0*dd)
    hu(i,j)=((u(i+1,j)+u(i,j))*(h(i+1,j)/mij(i+1,j)-h(i,j)/mij(i,j))+(u(i,j)&
    &+u(i-1,j))*(h(i,j)/mij(i,j)-h(i-1,j)/mij(i-1,j)))/(4.0*dd)
    hv(i,j)=((v(i,j+1)+v(i,j))*(h(i,j+1)/mij(i,j+1)-h(i,j)/mij(i,j))+(v(i,j)&
    &+v(i,j -1))*(h(i,j)/mij(i,j)-h(i,j-1)/mij(i,j-1)))/(4.0*dd)
    huv(i,j)=(h(i,j)-z0)/mij(i,j)*((u(i+1,j)-u(i-1,j))/(2.0*dd)+(v(i,j+1)-v(i,j-1))/(2.0*dd))
    fstar(i,j)=(fij(i,j)+u(i,j)*(mij(i,j+1)-mij(i,j-1))/(2.0*dd)-v(i,j)*(mij(i+1,j)-mij(i-1,j))/(2.0*dd))
  enddo
enddo
do j=2,n-1
  do i=2,m-1
    du(i,j)=-1.0*mij(i,j)*(uu(i,j)+uv(i,j)+g*(h(i+1,j)-h(i-1,j))/(2.0*dd))+v(i,j)*fstar(i,j)
    dv(i,j)=-1.0*mij(i,j)*(vv(i,j)+vu(i,j)+g*(h(i,j+1)-h(i,j-1))/(2.0*dd))-u(i,j)*fstar(i,j)
    dh(i,j)=-1.0*mij(i,j)**2*(hu(i,j)+hv(i,j)+huv(i,j))
  enddo
enddo
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------ŷ�����-----------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Euler(m,n,u,v,h,du,dv,dh,fij,mij,dd,dt)
integer m, n, k
real mij(m,n),fij(m,n)
real u(m,n,3),v(m,n,3),h(m,n,3)
real du(m,n),dv(m,n),dh(m,n),dt,dd
k=2
!ע�⼸��ʱ������ֵ�Ľ���
u(1:m,1:n,k)=u(1:m,1:n,k+1)
v(1:m,1:n,k)=v(1:m,1:n,k+1)
h(1:m,1:n,k)=h(1:m,1:n,k+1)
call TrendEGH(m,n,u(1:m,1:n,k),v(1:m,1:n,k),h(1:m,1:n,k),du,dv,dh,fij,mij,dd,dt)
do j=2, n-1
  do i=2, m-1
    u(i,j,k+1)=u(i,j,k)+du(i,j)*dt
    v(i,j,k+1)=v(i,j,k)+dv(i,j)*dt
    h(i,j,k+1)=h(i,j,k)+dh(i,j)*dt
  enddo
enddo
call TrendEGH(m,n,u(1:m,1:n,k+1),v(1:m,1:n,k+1),h(1:m,1:n,k+1),du,dv,dh,fij,mij,dd,dt)
do j=2, n-1
  do i=2, m-1
    u(i,j,k+1)=u(i,j,k)+du(i,j)*dt
    v(i,j,k+1)=v(i,j,k)+dv(i,j)*dt
    h(i,j,k+1)=h(i,j,k)+dh(i,j)*dt
  enddo
enddo
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------�����𲽷�----------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TriJmp(m,n,u,v,h,du,dv,dh,fij,mij,dd,dt)
integer m, n
integer k
real mij(m,n),fij(m,n)
real u(m,n,3),v(m,n,3),h(m,n,3)
real du(m,n),dv(m,n),dh(m,n),dt,dd
k=2
u(1:m,1:n,k-1)=u(1:m,1:n,k)
v(1:m,1:n,k-1)=v(1:m,1:n,k)
h(1:m,1:n,k-1)=h(1:m,1:n,k)
call TrendEGH(m,n,u(1:m,1:n,k),v(1:m,1:n,k),h(1:m,1:n,k),du,dv,dh,fij,mij,dd,dt)
do j=2, n-1
  do i=2, m-1
    u(i,j,k)=u(i,j,k-1)+du(i,j)*dt*0.5
    v(i,j,k)=v(i,j,k-1)+dv(i,j)*dt*0.5
    h(i,j,k)=h(i,j,k-1)+dh(i,j)*dt*0.5
  enddo
enddo
call TrendEGH(m,n,u(1:m,1:n,k),v(1:m,1:n,k),h(1:m,1:n,k),du,dv,dh,fij,mij,dd,dt)
do j=2, n-1
  do i=2, m-1
    u(i,j,k)=u(i,j,k-1)+du(i,j)*dt
    v(i,j,k)=v(i,j,k-1)+dv(i,j)*dt
    h(i,j,k)=h(i,j,k-1)+dh(i,j)*dt
  enddo
enddo
call TrendEGH(m,n,u(1:m,1:n,k),v(1:m,1:n,k),h(1:m,1:n,k),du,dv,dh,fij,mij,dd,dt)
do j=2, n-1
  do i=2, m-1
    u(i,j,k+1)=u(i,j,k-1)+du(i,j)*dt*2.0
    v(i,j,k+1)=v(i,j,k-1)+dv(i,j)*dt*2.0
    h(i,j,k+1)=h(i,j,k-1)+dh(i,j)*dt*2.0
  enddo
enddo
end subroutine 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------�ŵ�ƽ��-----------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NinSmth(m, n, Var, s)
integer m, n
real S
real Temp(m,n), Var(m,n)
Temp=Var
do i=1, m
  if(i==2.or.i==m-1)then
    do j=2,n-1
     Var(i,j)=Temp(i,j)+S*(1.0-S)/2.0*(Temp(i+1,j)+Temp(i,j+1)+Temp(i-1,j)+Temp(i,j-1)-4.0*Temp(i,j))
     Var(i,j)=Var(i,j)+S*S/4.0*(Temp(i+1,j+1)+Temp(i-1,j+1)+Temp(i-1,j-1)+Temp( i+1,j-1)-4.0*Temp(i,j))
    end do
  endif
enddo
do j=1, n
  if(j==2.or.j==n-1)then
    do i=2,m-2
     Var(i,j)=Temp(i,j)+S*(1.0-S)/2.0*(Temp(i+1,j)+Temp(i,j+1)+Temp(i-1,j)+Temp(i,j-1)-4.0*Temp(i,j))
     Var(i,j)=Var(i,j)+S*S/4.0*(Temp(i+1,j+1)+Temp(i-1,j+1)+Temp(i-1,j-1)+Temp( i+1,j-1)-4.0*Temp(i,j))
    end do
  endif
enddo
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------���ƽ��------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FivSmth(m, n, Var, s)
integer m, n
real s
real Temp(m,n), Var(m,n)
Temp=Var
do i=2,m-1
  do j=2,n-1
     Var(i, j)=Temp(i, j)+s/4.0*(Temp(i+1,j)+Temp(i,j+1)+Temp(i-1,j)+Temp(i,j-1)-4.0*Temp(i,j))
  end do
end do
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------ʱ��ƽ��------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TimSmth(m, n, Var, s)
integer m, n
integer,parameter::k=2
real Var(m, n, 3), S
do i=2,m-1
  do j=2,n-1
    Var(i,j,k)=(1.0-S)*Var(i,j,k)+S/2.0*(Var(i,j,k+1)+Var(i,j,k-1))
  end do
end do
end subroutine 
