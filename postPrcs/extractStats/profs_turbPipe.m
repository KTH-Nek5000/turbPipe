%*******************************************************************************
%  The script: - reads in 'data_turbPipe_stat.mat'
%              - transfers quantities from (x,y,z) to cylindrical coordinates (r,\theta,z).
%              - does averaging over azimuthal direction \theta.
%              - writes out 1D profiles of the flow statistics.  
%  Note: set whichEnv at Ln.17
%*******************************************************************************
%-------------------------------------------------------------------------------
close all; 
clear; 
clc
double precision;
format long;

%Choose the environment your are using
%'matlab' or 'octave'
whichEnv='matlab';   

%% Load data
load './tempData/data_turbPipe_stat.mat'
load './tempData/angles.mat'
load './tempData/geo.mat'
pipe=struct;

%%re-order x,y
npoints=length(x);
k=0;
for i=1:nR
    for j=1:nTh
        k=k+1;
        X(i,j)=x(k);
        Y(i,j)=y(k);
    end
end
k=0;
for j=1:nTh
    for i=1:nR
        k=k+1;
        x1(k,1)=X(i,j);
        y1(k,1)=Y(i,j);
    end
end

%transpose the imported data tensors
U=U'; V=V'; W=W';
uu=uu'; uv=uv'; uw=uw'; vv=vv'; vw=vw'; ww=ww';
dUdx=dUdx'; dUdy=dUdy'; dUdz=dUdz';
dVdx=dVdx'; dVdy=dVdy'; dVdz=dVdz';
dWdx=dWdx'; dWdy=dWdy'; dWdz=dWdz';
Pxx=Pxx'; Pxy=Pxy'; Pxz=Pxz'; Pyy=Pyy'; Pyz=Pyz'; Pzz=Pzz'; 
Dxx=Dxx'; Dxy=Dxy'; Dxz=Dxz'; Dyy=Dyy'; Dyz=Dyz'; Dzz=Dzz'; 
Cxx=Cxx'; Cxy=Cxy'; Cxz=Cxz'; Cyy=Cyy'; Cyz=Cyz'; Czz=Czz'; 
Txx=Txx'; Txy=Txy'; Txz=Txz'; Tyy=Tyy'; Tyz=Tyz'; Tzz=Tzz'; 
VDxx=VDxx'; VDxy=VDxy'; VDxz=VDxz'; VDyy=VDyy'; VDyz=VDyz'; VDzz=VDzz'; 
PTxx=PTxx'; PTxy=PTxy'; PTxz=PTxz'; PTyy=PTyy'; PTyz=PTyz'; PTzz=PTzz'; 
PSxx=PSxx'; PSxy=PSxy'; PSxz=PSxz'; PSyy=PSyy'; PSyz=PSyz'; PSzz=PSzz'; 
uuu=uuu'; uvv=uvv'; uuw=uuw';
uuv=uuv';           uvw=uvw';
uuw=uuw';           uww=uww';
          vvv=vvv'; vvw=vvw';
                    vww=vww';
                    www=www'; 

%% Rotation to cylindrical coord.
n2=0;
n3=0;
for i=1:nTh
        
    n1=n2+1;
    n2=nR*i+n3;
        
    pipe.r1(i,:)  = sqrt(x1(n1:n2).^2+y1(n1:n2).^2);   
             
    R = [cos(angle(i)) sin(angle(i)) 0;
        -sin(angle(i)) cos(angle(i)) 0;
                 0             0            1];

    for jj=n1:n2
        j=jj-n1+1;

        %Mean velocities. Tensors of Rank 1.
        prod = R*[U(jj); V(jj); W(jj)]; 
            
        pipe.Ur(i,j) = prod(1);
        pipe.Ut(i,j) = prod(2);
        pipe.Uz(i,j) = prod(3); 
        
        %Reynolds stress tensor. Tensors of Rank 2.
        S = [uu(jj) uv(jj) uw(jj);   
             uv(jj) vv(jj) vw(jj); 
             uw(jj) vw(jj) ww(jj)];
        prod = R*S*R';

        pipe.urur(i,j) = prod(1,1);
        pipe.urut(i,j) = prod(1,2);
        pipe.uruz(i,j) = prod(1,3);
        pipe.utut(i,j) = prod(2,2);
        pipe.uzuz(i,j) = prod(3,3);
        pipe.utuz(i,j) = prod(2,3);
            
        %Velocity gradient tensor. Tensor of Rank 2.     
        prod = R*[dUdx(jj) dUdy(jj) dUdz(jj);
                  dVdx(jj) dVdy(jj) dVdz(jj);
                  dWdx(jj) dWdy(jj) dWdz(jj)]*R';
        
        pipe.dUrdr(i,j)  = prod(1,1);
        pipe.dUrdt(i,j)  = prod(1,2);
        pipe.dUrdz(i,j)  = prod(1,3);
        pipe.dUtdr(i,j)  = prod(2,1);
        pipe.dUtdt(i,j)  = prod(2,2);
        pipe.dUtdz(i,j)  = prod(2,3);
        pipe.dUzdr(i,j)  = prod(3,1);
        pipe.dUzdt(i,j)  = prod(3,2);   
        pipe.dUzdz(i,j)  = prod(3,3);
 
        %Production tensor. Tensor of Rank 2.
        S = [Pxx(jj) Pxy(jj) Pxz(jj); 
             Pxy(jj) Pyy(jj) Pyz(jj); 
             Pxz(jj) Pyz(jj) Pzz(jj)];

        prod=R*S*R';        
        pipe.Prr(i,j) = prod(1,1);
        pipe.Ptt(i,j) = prod(2,2);
        pipe.Pzz(i,j) = prod(3,3);
        pipe.Prt(i,j) = prod(1,2);
        pipe.Prz(i,j) = prod(1,3);
        pipe.Ptz(i,j) = prod(2,3);
 
        %Dissipation tensor. Tensor of Rank 2.
        S = [Dxx(jj) Dxy(jj) Dxz(jj); 
             Dxy(jj) Dyy(jj) Dyz(jj);  
             Dxz(jj) Dyz(jj) Dzz(jj)];

        prod=R*S*R';
        pipe.Drr(i,j) = prod(1,1);
        pipe.Dtt(i,j) = prod(2,2);
        pipe.Dzz(i,j) = prod(3,3);
        pipe.Drt(i,j) = prod(1,2);
        pipe.Drz(i,j) = prod(1,3);
        pipe.Dtz(i,j) = prod(2,3);
 
        %Mean convection tensor. Tensor of Rank 2.
        S = [Cxx(jj) Cxy(jj) Cxz(jj); 
             Cxy(jj) Cyy(jj) Cyz(jj);  
             Cxz(jj) Cyz(jj) Czz(jj)];

        prod=R*S*R';   
        pipe.Crr(i,j) = prod(1,1);
        pipe.Ctt(i,j) = prod(2,2);
        pipe.Czz(i,j) = prod(3,3);
        pipe.Crt(i,j) = prod(1,2);
        pipe.Crz(i,j) = prod(1,3);
        pipe.Ctz(i,j) = prod(2,3);
    
        %Turbulent transport tensor. Tensor of Rank 2.
        S = [Txx(jj) Txy(jj) Txz(jj); 
             Txy(jj) Tyy(jj) Tyz(jj);  
             Txz(jj) Tyz(jj) Tzz(jj)];

        prod=R*S*R'; 
        pipe.Trr(i,j) = prod(1,1);
        pipe.Ttt(i,j) = prod(2,2);
        pipe.Tzz(i,j) = prod(3,3);
        pipe.Trt(i,j) = prod(1,2);
        pipe.Trz(i,j) = prod(1,3);
        pipe.Ttz(i,j) = prod(2,3);

        %Viscous diffusion tensor. Tensor of Rank 2.
        S = [VDxx(jj) VDxy(jj) VDxz(jj); 
             VDxy(jj) VDyy(jj) VDyz(jj);  
             VDxz(jj) VDyz(jj) VDzz(jj)];

        prod=R*S*R';    
        pipe.VDrr(i,j) = prod(1,1);
        pipe.VDtt(i,j) = prod(2,2);
        pipe.VDzz(i,j) = prod(3,3);
        pipe.VDrt(i,j) = prod(1,2);
        pipe.VDrz(i,j) = prod(1,3);
        pipe.VDtz(i,j) = prod(2,3);

        %Pressure transport tensor. Tensor of Rank 2.   
        S = [PTxx(jj) PTxy(jj) PTxz(jj); 
             PTxy(jj) PTyy(jj) PTyz(jj);  
             PTxz(jj) PTyz(jj) PTzz(jj)];

        prod=R*S*R';    
        pipe.PTrr(i,j) = prod(1,1);
        pipe.PTtt(i,j) = prod(2,2);
        pipe.PTzz(i,j) = prod(3,3);
        pipe.PTrt(i,j) = prod(1,2);
        pipe.PTrz(i,j) = prod(1,3);
        pipe.PTtz(i,j) = prod(2,3); 

        %Pressure strain tensor. Tensor of Rank 2.
        S = [PSxx(jj) PSxy(jj) PSxz(jj); 
             PSxy(jj) PSyy(jj) PSyz(jj);  
             PSxz(jj) PSyz(jj) PSzz(jj)];

        prod=R*S*R';    
        pipe.PSrr(i,j) = prod(1,1);
        pipe.PStt(i,j) = prod(2,2);
        pipe.PSzz(i,j) = prod(3,3);
        pipe.PSrt(i,j) = prod(1,2);
        pipe.PSrz(i,j) = prod(1,3);
        pipe.PStz(i,j) = prod(2,3); 
    
            %Budget for each component of the Reynolds stress tensor 
%            %Without mean convection
%             pipe.Srr(i,j,k) = pipe.Prr(i,j,k)+pipe.Drr(i,j,k)+pipe.Trr(i,j,k)+pipe.VDrr(i,j,k)+pipe.Pirr(i,j,k); 
%             pipe.Stt(i,j,k) = pipe.Ptt(i,j,k)+pipe.Dtt(i,j,k)+pipe.Ttt(i,j,k)+pipe.VDtt(i,j,k)+pipe.Pitt(i,j,k);
%             pipe.Szz(i,j,k) = pipe.Pzz(i,j,k)+pipe.Dzz(i,j,k)+pipe.Tzz(i,j,k)+pipe.VDzz(i,j,k)+pipe.Pizz(i,j,k);
%             pipe.Srt(i,j,k) = pipe.Prt(i,j,k)+pipe.Drt(i,j,k)+pipe.Trt(i,j,k)+pipe.VDrt(i,j,k)+pipe.Pirt(i,j,k); 
%             pipe.Srz(i,j,k) = pipe.Prz(i,j,k)+pipe.Drz(i,j,k)+pipe.Trz(i,j,k)+pipe.VDrz(i,j,k)+pipe.Pirz(i,j,k);
%             pipe.Stz(i,j,k) = pipe.Ptz(i,j,k)+pipe.Dtz(i,j,k)+pipe.Ttz(i,j,k)+pipe.VDtz(i,j,k)+pipe.Pitz(i,j,k);
%            %With mean convection
%             pipe.Scrr(i,j,k) = pipe.Prr(i,j,k)+pipe.Drr(i,j,k)+pipe.Trr(i,j,k)+pipe.VDrr(i,j,k)+pipe.Pirr(i,j,k)-pipe.Crr(i,j,k); 
%             pipe.Sctt(i,j,k) = pipe.Ptt(i,j,k)+pipe.Dtt(i,j,k)+pipe.Ttt(i,j,k)+pipe.VDtt(i,j,k)+pipe.Pitt(i,j,k)-pipe.Ctt(i,j,k);
%             pipe.Sczz(i,j,k) = pipe.Pzz(i,j,k)+pipe.Dzz(i,j,k)+pipe.Tzz(i,j,k)+pipe.VDzz(i,j,k)+pipe.Pizz(i,j,k)-pipe.Czz(i,j,k);
%             pipe.Scrt(i,j,k) = pipe.Prt(i,j,k)+pipe.Drt(i,j,k)+pipe.Trt(i,j,k)+pipe.VDrt(i,j,k)+pipe.Pirt(i,j,k)-pipe.Crt(i,j,k); 
%             pipe.Scrz(i,j,k) = pipe.Prz(i,j,k)+pipe.Drz(i,j,k)+pipe.Trz(i,j,k)+pipe.VDrz(i,j,k)+pipe.Pirz(i,j,k)-pipe.Crz(i,j,k);
%             pipe.Sctz(i,j,k) = pipe.Ptz(i,j,k)+pipe.Dtz(i,j,k)+pipe.Ttz(i,j,k)+pipe.VDtz(i,j,k)+pipe.Pitz(i,j,k)-pipe.Ctz(i,j,k);

        %Skewness tensor. Tensor of Rank 3.    
        R3_tensor(:,:,1) = [uuu(jj) uvv(jj) uuw(jj);...
                            uuv(jj) uvv(jj) uvw(jj);...
                            uuw(jj) uvw(jj) uww(jj)];
        R3_tensor(:,:,2) = [uuv(jj) uvv(jj) uvw(jj);...
                            uvv(jj) vvv(jj) vvw(jj);...
                            uvw(jj) vvw(jj) vww(jj)];
        R3_tensor(:,:,3) = [uuw(jj) uvw(jj) uww(jj);...
                            uvw(jj) vvw(jj) vww(jj);...
                            uww(jj) vww(jj) www(jj)];

        aabc(1:3,1:3,1:3)=0;
        adef=R3_tensor(:,:,:);
        for aa=1:3
        for bb=1:3    
        for cc=1:3    
        for dd=1:3
        for ee=1:3    
        for ff=1:3
            aabc(aa,bb,cc)=aabc(aa,bb,cc)+R(aa,dd)*R(bb,ee) ...
               *R(cc,ff)*adef(dd,ee,ff);
        end    
        end
        end
        end    
        end
        end

        pipe.ururur(i,j) = aabc(1,1,1);
        pipe.ututut(i,j) = aabc(2,2,2);
        pipe.uzuzuz(i,j) = aabc(3,3,3);
        pipe.ururut(i,j) = aabc(1,2,1);
        pipe.ururuz(i,j) = aabc(1,3,1);
        pipe.urutut(i,j) = aabc(2,2,1);
        pipe.ututuz(i,j) = aabc(2,3,2);
        pipe.uruzuz(i,j) = aabc(3,3,1);
        pipe.utuzuz(i,j) = aabc(3,3,2);
        pipe.urutuz(i,j) = aabc(2,3,1);    
          
    end       
end    
n3=n2;

%% Average over the azimuthal direction
if whichEnv=='matlab'
   clearvars -except pipe nR nu yplus X Y    %matlab
elseif whichEnv=='octave'
   clear -x pipe nR nu yplus X Y   %octave
end
pipe2=struct;

for j=1:nR     
    %Mean velocities. Tensors of Rank 1.        
    pipe2.Ur(j) = mean(pipe.Ur(:,j));
    pipe2.Ut(j) = mean(pipe.Ut(:,j));       
    pipe2.Uz(j) = mean(pipe.Uz(:,j));

    %Reynolds stress tensor. Tensors of Rank 2.
    pipe2.urur(j) = mean(pipe.urur(:,j));
    pipe2.urut(j) = mean(pipe.urut(:,j));
    pipe2.uruz(j) = mean(pipe.uruz(:,j));
    pipe2.utut(j) = mean(pipe.utut(:,j));
    pipe2.uzuz(j) = mean(pipe.uzuz(:,j));
    pipe2.utuz(j) = mean(pipe.utuz(:,j));
        
    %Velocity gradient tensor. Tensor of Rank 2.
    pipe2.dUrdr(j) = mean(pipe.dUrdr(:,j));
    pipe2.dUrdt(j) = mean(pipe.dUrdt(:,j));
    pipe2.dUrdz(j) = mean(pipe.dUrdz(:,j));
    pipe2.dUtdr(j) = mean(pipe.dUtdr(:,j));
    pipe2.dUtdt(j) = mean(pipe.dUtdt(:,j));
    pipe2.dUtdz(j) = mean(pipe.dUtdz(:,j));
    pipe2.dUzdr(j) = mean(pipe.dUzdr(:,j));
    pipe2.dUzdt(j) = mean(pipe.dUzdt(:,j)); 
    pipe2.dUzdz(j) = mean(pipe.dUzdz(:,j));    
         
    %Production tensor. Tensor of Rank 2.     
    pipe2.Prr(j) = mean(pipe.Prr(:,j));
    pipe2.Ptt(j) = mean(pipe.Ptt(:,j));
    pipe2.Pzz(j) = mean(pipe.Pzz(:,j));
    pipe2.Prt(j) = mean(pipe.Prt(:,j));
    pipe2.Prz(j) = mean(pipe.Prz(:,j));
    pipe2.Ptz(j) = mean(pipe.Ptz(:,j));

    %Dissipation tensor. Tensor of Rank 2.
    pipe2.Drr(j) = mean(pipe.Drr(:,j));
    pipe2.Dtt(j) = mean(pipe.Dtt(:,j));
    pipe2.Dzz(j) = mean(pipe.Dzz(:,j));
    pipe2.Drt(j) = mean(pipe.Drt(:,j));
    pipe2.Drz(j) = mean(pipe.Drz(:,j));
    pipe2.Dtz(j) = mean(pipe.Dtz(:,j));
  
    %Mean convection tensor. Tensor of Rank 2.
    pipe2.Crr(j) = mean(pipe.Crr(:,j));
    pipe2.Ctt(j) = mean(pipe.Ctt(:,j));
    pipe2.Czz(j) = mean(pipe.Czz(:,j));
    pipe2.Crt(j) = mean(pipe.Crt(:,j));
    pipe2.Crz(j) = mean(pipe.Crz(:,j));
    pipe2.Ctz(j) = mean(pipe.Ctz(:,j));
  
    %Turbulent transport tensor. Tensor of Rank 2.
    pipe2.Trr(j) = mean(pipe.Trr(:,j));
    pipe2.Ttt(j) = mean(pipe.Ttt(:,j));
    pipe2.Tzz(j) = mean(pipe.Tzz(:,j));
    pipe2.Trt(j) = mean(pipe.Trt(:,j));
    pipe2.Trz(j) = mean(pipe.Trz(:,j));
    pipe2.Ttz(j) = mean(pipe.Ttz(:,j));

    %Viscous diffusion tensor. Tensor of Rank 2.
    pipe2.VDrr(j) = mean(pipe.VDrr(:,j));
    pipe2.VDtt(j) = mean(pipe.VDtt(:,j));
    pipe2.VDzz(j) = mean(pipe.VDzz(:,j));
    pipe2.VDrt(j) = mean(pipe.VDrt(:,j));
    pipe2.VDrz(j) = mean(pipe.VDrz(:,j));
    pipe2.VDtz(j) = mean(pipe.VDtz(:,j));
 
    %Pressure transport tensor. Tensor of Rank 2.   
    pipe2.PTrr(j) = mean(pipe.PTrr(:,j));
    pipe2.PTtt(j) = mean(pipe.PTtt(:,j));
    pipe2.PTzz(j) = mean(pipe.PTzz(:,j));
    pipe2.PTrt(j) = mean(pipe.PTrt(:,j));
    pipe2.PTrz(j) = mean(pipe.PTrz(:,j));
    pipe2.PTtz(j) = mean(pipe.PTtz(:,j));

    %Pressure strain tensor. Tensor of Rank 2.
    pipe2.PSrr(j) = mean(pipe.PSrr(:,j));
    pipe2.PStt(j) = mean(pipe.PStt(:,j));
    pipe2.PSzz(j) = mean(pipe.PSzz(:,j));
    pipe2.PSrt(j) = mean(pipe.PSrt(:,j));
    pipe2.PSrz(j) = mean(pipe.PSrz(:,j));
    pipe2.PStz(j) = mean(pipe.PStz(:,j));

%         
%         %Budget for each component of the Reynolds stress tensor 
%         %Without mean convection
%         pipe2.Srr(j,k) = mean(pipe.Srr(:,j,k));
%         pipe2.Stt(j,k) = mean(pipe.Stt(:,j,k));
%         pipe2.Szz(j,k) = mean(pipe.Szz(:,j,k));
%         pipe2.Srt(j,k) = mean(pipe.Srt(:,j,k));
%         pipe2.Srz(j,k) = mean(pipe.Srz(:,j,k));
%         pipe2.Stz(j,k) = mean(pipe.Stz(:,j,k));
%         %With mean convection
%         pipe2.Scrr(j,k) = mean(pipe.Scrr(:,j,k));
%         pipe2.Sctt(j,k) = mean(pipe.Sctt(:,j,k));
%         pipe2.Sczz(j,k) = mean(pipe.Sczz(:,j,k));
%         pipe2.Scrt(j,k) = mean(pipe.Scrt(:,j,k));
%         pipe2.Scrz(j,k) = mean(pipe.Scrz(:,j,k));
%         pipe2.Sctz(j,k) = mean(pipe.Sctz(:,j,k));
%               
   %Skewness tensor. Tensor of Rank 3. 
   pipe2.ururur(j) = mean(pipe.ururur(:,j));
   pipe2.ututut(j) = mean(pipe.ututut(:,j));
   pipe2.uzuzuz(j) = mean(pipe.uzuzuz(:,j));
   pipe2.ururut(j) = mean(pipe.ururut(:,j));
   pipe2.ururuz(j) = mean(pipe.ururuz(:,j));
   pipe2.urutut(j) = mean(pipe.urutut(:,j));
   pipe2.ututuz(j) = mean(pipe.ututuz(:,j));
   pipe2.uruzuz(j) = mean(pipe.uruzuz(:,j));
   pipe2.utuzuz(j) = mean(pipe.utuzuz(:,j));
   pipe2.urutuz(j) = mean(pipe.urutuz(:,j));
end    

%Compute TKE budget terms
P_k  = 0.5*(pipe2.Prr  + pipe2.Ptt  + pipe2.Pzz);   %production
T_k  = 0.5*(pipe2.Trr  + pipe2.Ttt  + pipe2.Tzz);   %turbulent transport
PS_k = 0.5*(pipe2.PSrr + pipe2.PStt + pipe2.PSzz);  %pressure-strain
PT_k = 0.5*(pipe2.PTrr + pipe2.PTtt + pipe2.PTzz);  %pressure-transport
Pi_k =     (PT_k - PS_k);
VD_k = 0.5*(pipe2.VDrr + pipe2.VDtt + pipe2.VDzz);  %viscous diffusion
D_k  = 0.5*(pipe2.Drr  + pipe2.Dtt  + pipe2.Dzz);   %dissipation
C_k  = 0.5*(pipe2.Crr  + pipe2.Ctt  + pipe2.Czz);   %mean convection

%Write postprocessed results in a file
% make an array of radii
r_=zeros(nR,1);
for i=1:nR
    r_(i,1)=sqrt(X(i,1)^2+Y(i,1)^2);
end
rMax=max(r_);
r_=rMax-r_;  %Changing the discrete radius from wall toward center

%%Compute uTau
uTau=sqrt(-pipe2.dUzdr(end)*nu);
%method 2
tauW=(pipe2.Uz(end-1)-pipe2.Uz(end))/(r_(end-1)-r_(end));
uTau2=sqrt(nu*tauW);
fprintf('uTau (method1)=%g\n',uTau)
fprintf('uTau (method2)=%g\n',uTau2)
ReTau=uTau*rMax/nu;

% velocity profiles + derivatives
fOut=fopen('./statsResults/turbPipe_meanVelocity.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# Mean velocity and their derivatives \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r\t Ur\t Ut\t Uz\t dUrdr\t dUrdt\t dUrdz\t dUtdr\t dUtdt\t dUtdz\t dUzdr\t dUzdt\t dUzdz \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t %g\t%g\t%g\t %g\t%g\t%g\t   \n',...
            r_(i,1),pipe2.Ur(i),pipe2.Ut(i),pipe2.Uz(i),...
            pipe2.dUrdr(i),pipe2.dUrdt(i),pipe2.dUrdz(i),... 
            pipe2.dUtdr(i),pipe2.dUtdt(i),pipe2.dUtdz(i),... 
            pipe2.dUzdr(i),pipe2.dUzdt(i),pipe2.dUzdz(i))
end
fclose(fOut);

%Reynolds stress components
fOut=fopen('./statsResults/turbPipe_ReyStress.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# Reynolds stress components \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r\t urp\t urut\t uruz\t utp\t uzp\t utuz \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t \n',...
            r_(i,1),sqrt(pipe2.urur(i)),pipe2.urut(i),pipe2.uruz(i),...
               sqrt(pipe2.utut(i)),sqrt(pipe2.uzuz(i)),pipe2.utuz(i))
end
fclose(fOut);

%TKE budget terms
fOut=fopen('./statsResults/turbPipe_kBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# TKE budget terms: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t P_k+ \t T_k+ \t PS_k+ \t PT_k+ \t VD_k+ \t D_k+ \t C_k+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu/uTau^4.;
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g\n',...
            r_(i,1)*ReTau,P_k(i)*fac,T_k(i)*fac,PS_k(i)*fac,PT_k(i)*fac,VD_k(i)*fac,D_k(i)*fac,C_k(i)*fac)
end
fclose(fOut);

%budget terms of <urur>
fOut=fopen('./statsResults/turbPipe_rrBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <urur>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t P_rr+ \t T_rr+ \t PS_rr+ \t PT_rr+ \t VD_rr+ \t D_rr+ \t C_rr+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu/uTau^4.;
P_rr  = pipe2.Prr;   %production
T_rr  = pipe2.Trr;   %turbulent transport
PS_rr = pipe2.PSrr;  %pressure-strain
PT_rr = pipe2.PTrr;  %pressure-transport
VD_rr = pipe2.VDrr;  %viscous diffusion
D_rr  = pipe2.Drr;   %dissipation
C_rr  = pipe2.Crr;   %mean convection
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n',...
            r_(i,1)*ReTau,P_rr(i)*fac,T_rr(i)*fac,PS_rr(i)*fac,PT_rr(i)*fac,VD_rr(i)*fac,D_rr(i)*fac,C_rr(i)*fac)
end
fclose(fOut);

%budget terms of <utut>
fOut=fopen('./statsResults/turbPipe_ttBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <utut>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t P_tt+ \t T_tt+ \t PS_tt+ \t PT_tt+ \t VD_tt+ \t D_tt+ \t C_tt+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu/uTau^4.;
P_tt  = pipe2.Ptt;   %production
T_tt  = pipe2.Ttt;   %turbulent transport
PS_tt = pipe2.PStt;  %pressure-strain
PT_tt = pipe2.PTtt;  %pressure-transport
VD_tt = pipe2.VDtt;  %viscous diffusion
D_tt  = pipe2.Dtt;   %dissipation
C_tt  = pipe2.Ctt;   %mean convection
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g\n',...
            r_(i,1)*ReTau,P_tt(i)*fac,T_tt(i)*fac,PS_tt(i)*fac,PT_tt(i)*fac,VD_tt(i)*fac,D_tt(i)*fac,C_tt(i)*fac)
end
fclose(fOut);

%budget terms of <uzuz>
fOut=fopen('./statsResults/turbPipe_zzBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <uzuz>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t P_zz+ \t T_zz+ \t PS_zz+ \t PT_zz+ \t VD_zz+ \t D_zz+ \t C_zz+  \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu/uTau^4.;
P_zz  = pipe2.Pzz;   %production
T_zz  = pipe2.Tzz;   %turbulent transport
PS_zz = pipe2.PSzz;  %pressure-strain
PT_zz = pipe2.PTzz;  %pressure-transport
VD_zz = pipe2.VDzz;  %viscous diffusion
D_zz  = pipe2.Dzz;   %dissipation
C_zz  = pipe2.Czz;   %mean convection
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n',...
            r_(i,1)*ReTau,P_zz(i)*fac,T_zz(i)*fac,PS_zz(i)*fac,PT_zz(i)*fac,VD_zz(i)*fac,D_zz(i)*fac,C_zz(i)*fac)
end
fclose(fOut);

%budget terms of <uruz>
fOut=fopen('./statsResults/turbPipe_rzBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <uruz>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# uTau = %g \n',uTau)
fprintf(fOut,'# nu = %g \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t P_rz+ \t T_rz+ \t PS_rz+ \t PT_rz+ \t VD_rz+ \t D_rz+ \t C_rz+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu/uTau^4.;
P_rz  = pipe2.Prz;   %production
T_rz  = pipe2.Trz;   %turbulent transport
PS_rz = pipe2.PSrz;  %pressure-strain
PT_rz = pipe2.PTrz;  %pressure-transport
VD_rz = pipe2.VDrz;  %viscous diffusion
D_rz  = pipe2.Drz;   %dissipation
C_rz  = pipe2.Crz;   %mean convection
for i=1:nR
    fprintf(fOut,'%g\t%g\t%g\t%g\t %g\t%g\t%g\t%g \n',...
            r_(i,1)*ReTau,P_rz(i)*fac,T_rz(i)*fac,PS_rz(i)*fac,PT_rz(i)*fac,VD_rz(i)*fac,D_rz(i)*fac,C_rz(i)*fac)
end
fclose(fOut);
