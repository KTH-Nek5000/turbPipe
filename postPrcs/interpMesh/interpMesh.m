%***************************************************************************
% Generating interpolation mesh (2d) over the circular surface of the pipe
%***************************************************************************
% Note: Set the SETTINGS
%***************************************************************************
close all; 
clear all; 
clc
double precision;
format long;

%------------------------------------------
%------ SETTINGS --------------------------
%% Simulation parameters
ReTau = 360;   %Nominal ReTau (is used only for dr+_1)
Rmax=1.0;      %Pipe radius
rho=1.0;       %Fluid density 
nu=1.0/5850.0; %Kinematic viscosity 

%% Specifications of the mesh to be used for interpolation
nR=140      %Number of points in radial direction
nTh=81      %Number of points in azimuthal direction
drp1=0.5;   %First r+ off the wall 
compressedMesh='True';   %1: Mesh compressed radially toward the wall
%------------------------------------------

%Points in the radial direction
if compressedMesh     
   %Non-uniform mesh in the wall-normal direction
   dr1=(drp1/ReTau)*Rmax;  %Distance from the wall of the first off-wall node
   gam=3.0;  %Grid compression controller >1
   xi=linspace(0.0,1.0,nR-1);
   R_=zeros(1,nR);
   for i=1:nR-1          
       R_(1,i+1)=dr1+(Rmax-dr1)*(1.0+(tanh(gam*((xi(1,i))/(1)-1.0)))/tanh(gam));
   end
   R_(1,1)=0.0;   %center of the pipe      
   %Reverse the order of points: to be from center toward the wall
   R_=Rmax-R_;
   for i=1:nR
       R(i)=R_(nR-i+1);
   end
else
   R = linspace(0.0,Rmax,nR);   %Uniform mesh in radial direction
end

%Discrete Azimuth angles
th = linspace(0,2.*pi,nTh);

%Number of mesh points
np = length(th)*length(R);

%Cartesian coordinates
k=0;
for i=1:nR
    for j=1:nTh
        xx(j,i) = R(i)*cos(th(j));
        yy(j,i) = R(i)*sin(th(j));
        theta(j,i) = th(j);
    end
end
x = xx(:);
y = yy(:);
angle = theta(:);

%Plot the mesh
figure
plot(x,y,'o-b'); 

%****** save data
save(['../extractStats/tempData/angles.mat'],'angle');
save(['../extractStats/tempData/geo.mat'],'R','nR','nTh','ReTau','Rmax','nu','rho','x','y');
%%
%Save x data points in Fortran binary format
fid=fopen('../ppNekOutData/ZSTAT/x.fort','w','ieee-le'); %saleh

%First write 4 bytes integer
data=np;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write real-8
data=x;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%Save y data points in Fortran binary format
fid=fopen('../ppNekOutData/ZSTAT/y.fort','w','ieee-le');

%First write 4 bytes integer
data=np;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write real-8
data=y;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

