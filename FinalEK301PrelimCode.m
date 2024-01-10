%Load input file---------------------------------------------
filename = input('Please type the name of the file you would like to load without its extension:\n','s');
 filename=strcat(filename,'.mat');
 load(filename)
tstart=tic;
x = X;
y = Y;
[J,G]= size(Sx);
M=(2*J)-3;
changex=zeros(M,1);
changey=zeros(M,1);
ux=zeros(M,1);
uy=zeros(M,1);
Ax=zeros(J,M);
Ay=zeros(J,M);
%Calculates coefficients for matrix A-----------------------------------
for m=1:M
 vec=find(C(:,m));
 changex(m)=abs(x(vec(2))-x(vec(1)));
 changey(m)=abs(y(vec(2))-y(vec(1)));
 r(m)=magnitude(changex(m),changey(m));
 ux(m)=changex(m)/r(m);
 uy(m)=changey(m)/r(m);
 if x(vec(1))<x(vec(2))
 Ax(vec(1),m)=ux(m);
 Ax(vec(2),m)=-ux(m);
 else
 Ax(vec(1),m)=-ux(m);
 Ax(vec(2),m)=ux(m);
 end
 if y(vec(1))<y(vec(2))
 Ay(vec(1),m)=uy(m);
 Ay(vec(2),m)=-uy(m);
 else
 Ay(vec(1),m)=-uy(m);
 Ay(vec(2),m)=uy(m);
 end
end
A=[Ax Sx;Ay Sy]; %Creates matrix A
T=A\L; %Solves equilibrium equations

%--------------result analysis-----------------------------
%tension/compression
Tm=T(1:length(T)-3);
tension=find(Tm>0);
compression=find(Tm<0);

zeroforce=find(Tm==0);
[Fb,U]=buckle(r(compression)); 
comprat=T(compression)./Fb'; 


maxcomprat=abs(min(comprat)); 
n=find(min(comprat)==comprat); 
critmember=compression(n); 

mtl=(Fb(n(1))'./abs(T(critmember(1)))).*Fb; 

mtlactual=(r(critmember(1))/(r(critmember(1))-.5))*mtl; 
cost=(10*J)+sum(r);
loadcostrat=mtl/cost;
telapsed=toc(tstart);

shortest=min(r);
longest=max(r);
shortestm=find(r==min(r));
longestm=find(r==max(r));

%----print results---------------------------------------------------
clc
fprintf('\n\n\t\t RESULTS\n\n')

fprintf('Truss Dimensions:\n')

fprintf('The length of the truss is %.1f in.\n',max(x))
fprintf('The height of the truss is %.1f in.\n',max(y)-min(y))

fprintf('Member Number\tJoint-Joint Length (cm)\n')
for m=1:M
 fprintf('%d\t\t%.1f\n',m,r(m))
end
fprintf('Shortest Members:\nMembers %s are each %.1fin.\n',num2str(shortestm),shortest)
fprintf('Longest Members:\nMembers %s are each %.1fin.\n\n',num2str(longestm),longest)

fprintf('/n Member & Support Forces:\n')
z=length(T);
k=1;
for m=1:M
 if abs(T(m))<=(10^-5)
 str=sprintf('Member %d is a zero force member.\n',m);
 elseif T(m)<0

 str=sprintf('Member %d: %.3f oz (C), Buckling Strength is %.3foz.\n',m,abs(T(m)),Fb(k));
 k=k+1;
 elseif T(m)>0
 str=sprintf('Member %d: %.3f oz (T).\n',m,abs(T(m)));
 end
 fprintf(str)
end
fprintf('Sx1 =  %.3foz \n',T(z-2))
fprintf('Sy1 = %.3foz \n',T(z-1))
fprintf('Sy2 = %.3foz \n\n',T(z))
fprintf('Failure Data:\n')
if length(critmember)>1
 fprintf('The critical members are:\n %s.\n',num2str(critmember'))
else
fprintf('The critical member is member %d.\n',critmember)
end
fprintf('Cost Data:\n')
fprintf('The total cost of the truss is $%.2f.\n',cost)
fprintf('The theoretical maximum load to cost ratio of the truss is %f(oz/$).\n\n',loadcostrat(end))
choice5=menu('Would you like to save your truss for next time?','Save Now','Quit');
switch choice5
 case 1
 filename=input('Please enter the file name you would like for the truss variables file, without the extension:\n','s');
 filename=strcat(filename,'.mat');
 save(filename,'J','M','x','y','C')
 fprintf('File Saved.\n')
end
%---functions-----------------------------------------------------------
function [F,U] = buckle(l)
%calculates the buckling strength (in oz) of a length of acrylic strip
%Uses data from the class testing labs and returns the error of the buckling strength in oz
F=3654.533.*(l.^-2.119);
U=1.685;
end
function r = magnitude(x,y)
%finds the magnitude of a vector based on its x and y components
r=((abs(x).^2)+(abs(y).^2)).^(1/2);
end
