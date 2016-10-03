%%%% copyright photrack ag, Oct 2016
%%%% By: Beat Lüthi
%%%%
%%%% this script transforms iMoMo Dischagre App measurements into 4 files:
%%%%
%%%% geometry.txt
%%%% m_corrds.txt
%%%%
%%%% discharge_geometry.txt
%%%% discharge_m_corrds.txt
%%%%
%%%% Instructions:
%%%% 
%%%% all user input goes betwen the lines 19 and 36
%%%% explanations are right at these lines

global l12 l14 l23 l34 w1 w2 w3 w4 dist2water water_level render

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Samu Furrow Oct 3 2016
l12=0.6; % intermarker dist between marker 1->2 etc
l34=0.53; % 
l14=1.13; %
l23=1.03; %
w1=0.00; % leave at 0!
w2=0.00; % 
w3=0.0; %
w4=0.0; %
load_path='C:\Users\Beat\Dropbox\discharge_app_sites\Samu_Forrow'; % path where the measBeam.txt file is
save_path='C:\Users\Beat\Dropbox\discharge_app_sites\Samu_Forrow/';% output path for the 4 files

al_z=0; % leave at 0
al_y=0; % leave at 0
T=[0.41 0.5 0]; % you can arbitrarily offset the coordinate system
water_level=0.085; % sets the water level
he_tol=1; % for the search algorithm, unit is metric. Distance allowed for markers above water level.
%%%% end of Samu Furrow Oct 3 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(load_path)
if exist('beamMeas.txt')
    data=load('beamMeas.txt');
    ind_beam=find(data(:,1)==-1);
    for i=1:length(ind_beam)
        be=ind_beam(i)+1;
        if i<length(ind_beam)
            en=ind_beam(i+1)-1;
        else
            en=length(data);
        end
        %load dist2water meas
        count=1;
        dist2water{i}(count,1)=data(ind_beam(i),2);
        dist2water{i}(count,2)=data(ind_beam(i),3);
        for j=be:en
            if isnan(data(j,2))==0
                count=count+1;
                dist2water{i}(count,1)=data(j,1);
                dist2water{i}(count,2)=data(j,2);
            end
        end
        aa0=1;
    end
end


he(1)=water_level;
he(2)=he(1);
he(3)=he(2);
he(4)=he(1);

ms=0.5*(l14+l23);%marker_scale
min_y12=min(w1,w2);
max_y34=-min(w3,w4);


p1=[-0.5*l12  0.5*ms he(1)];
p2=[ 0.5*l12  0.5*ms he(2)];
p3=[-0.5*l34 -0.5*ms he(3)];
p4=[ 0.5*l34 -0.5*ms he(4)];

vec=[p1 p2 p3 p4];
render=0;
initial_quality=eval_points_beam(vec);

LB=[-ms  min_y12 he(1)-he_tol  -0.5*ms  min_y12  he(2)-he_tol -ms -ms      he(3)-he_tol  0 -ms he(4)-he_tol];
UB=[  0.5*ms  ms      he(1)+he_tol  ms ms       he(2)+he_tol   0  max_y34 he(3)+he_tol ms  max_y34 he(4)+he_tol];
initial_quality = eval_points_beam(vec);
vec             = fminsearchbnd(@eval_points_beam, vec,LB,UB);%,optimset('maxiter',100000),'maxfunevals',100,
render=1;
current_quality = eval_points_beam(vec);

p1=vec(1:3);
p2=vec(4:6);
p3=vec(7:9);
p4=vec(10:12);
[p1;p2;p3;p4];

[LB;vec;UB]

p=[p1;p2;p3;p4];
Ry=[cosd(al_y) 0 sind(al_y);0 1 0;-sind(al_y) 0 cosd(al_y)];
p=p*Ry;
Rz=[cosd(al_z) sind(al_z) 0;-sind(al_z) cosd(al_z) 0;0 0 1];
p=p*Rz;

% figure(1);
% scatter3(p(:,1),p(:,2),p(:,3))
% axis equal

le=length(p);
p=p+repmat(T,le,1);

figure(1);hold off;
scatter3(p(1,1),p(1,2),p(1,3),'r','filled')
hold on;
scatter3(p(2,1),p(2,2),p(2,3),'g','filled')
scatter3(p(3,1),p(3,2),p(3,3),'b','filled')
scatter3(p(4,1),p(4,2),p(4,3),'c','filled')


%%%%%and now do geometry 
if exist('beamMeas.txt')
    geo=[];
    all_geo=[];
    interZ=[];
    col='bck';
    data=load('beamMeas.txt');
    ind_beam=find(data(:,1)==-1);
    for i=1:length(ind_beam)
        be=ind_beam(i)+1;
        if i<length(ind_beam)
            en=ind_beam(i+1)-1;
        else
            en=length(data);
        end
        vec=[p(data(ind_beam(i),3),:)-p(data(ind_beam(i),2),:)];
        vec=vec/norm(vec);
        %now build profile points
        for j=be:en
            geo{i}(j-be+1,1:3)=p(data(ind_beam(i),2),:)+data(j,1)*vec+data(j,3)*[0 0 -1];
        end
        %scatter3(geo{i}(:,1),geo{i}(:,2),geo{i}(:,3),col(i),'filled')
        all_geo=[all_geo;geo{i}];
    end
    %now build for each i an interp at all positions, so then we can take
    %the mean
    [tmp,k]=sort(all_geo(:,2),'ascend');%sort by y
    yz=[all_geo(k,2) all_geo(k,3)];
    [unique_y,ic,ia]=unique(yz(:,1));
    for i=1:length(ind_beam)
        interZ(1:length(unique_y),i)=interp1(geo{i}(:,2),geo{i}(:,3),unique_y);    
    end
    for i=1:length(unique_y)
       interGeo(i,1)=unique_y(i);
       interGeo(i,2)=nanmean_photrack(interZ(i,:));
    end
   
    meanX=mean(p(:,1));
    scatter3(repmat(meanX,length(unique_y),1),interGeo(:,1),interGeo(:,2))
    geometry=[repmat(meanX,length(unique_y),1) interGeo(:,1) interGeo(:,2)]
    rg=reshape(geometry(:,2:3),length(geometry)*2,1);
else
    rg=reshape(geometry,length(geometry)*2,1);
end

xlabel('x')
ylabel('y')
zlabel('z')
axis equal
title(['quality: ',num2str(current_quality)]);
box on;
grid off;


coords=p
cd(save_path)
save('m_coords.txt','coords','-ascii');
save('geometry.txt','geometry','-ascii');

rg=reshape(geometry(:,2:3),length(geometry)*2,1);
fid=fopen('discharge_geometry.txt','w');
for i=1:length(rg)
    fprintf(fid,'%f',rg(i));
    if i<length(rg)
       fprintf(fid,',',rg(i)); 
    end
end
fclose(fid);

rm=reshape(coords,length(coords)*3,1);
fid=fopen('discharge_m_coords.txt','w');
for i=1:length(rm)
    fprintf(fid,'%f',rm(i));
    if i<length(rm)
       fprintf(fid,',',rm(i)); 
    end
end
fclose(fid);

