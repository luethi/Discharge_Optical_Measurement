%%%% copyright photrack ag, Oct 2016
%%%% By: Beat Lüthi
%%%%
%%%% this script transforms Disto measurements into 4 files:
%%%%
%%%% geometry.txt
%%%% m_corrds.txt
%%%%
%%%% discharge_geometry.txt
%%%% discharge_m_corrds.txt
%%%%
%%%% Instructions:
%%%% 
%%%% all user input goes betwen the lines 17 and 31
%%%% explanations are right at these lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Wildbach 2016-09-15
disto_path='C:\Users\Beat\Dropbox\river_run_off\photrack_example_cases\Wildbach\2016-09-15 10-43-32_disto'
al_z=0; % leave at 0
al_y=0; % you can rotate the frame of reference if the slope is very steep. The e_x line should be parallel to the bottom of the channel/river/furrow
T=[0.3945 1.726+0.0274 1.265-0.4202+0.524]; % you can arbitrarily offset the coordinate system
i_ref=[1 2]; %list of reference points, actually always 1,2
iw=[1 2 3]; %list of waterline points (within the disto.txt file)
im=[4 5 6 7]; %list of marker points (within the disto.txt file)
ic=[]; %index of camera point, does not apply for app
ir=[]; %list of other points that you want to render anyway, for whatever reason

load_path='C:\Users\Beat\Dropbox\river_run_off\photrack_example_cases\Wildbach\2016-09-15 10-43-32_disto'; % path where the disto.txt & measBeam.txt file is
save_path='C:\Users\Beat\Dropbox\river_run_off\photrack_example_cases\Wildbach\2016-09-15 10-43-32_disto/';% output path for the 4 files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(load_path)
p=load('disto.txt');

Ry=[cosd(al_y) 0 sind(al_y);0 1 0;-sind(al_y) 0 cosd(al_y)];
p=p*Ry;
Rz=[cosd(al_z) sind(al_z) 0;-sind(al_z) cosd(al_z) 0;0 0 1];
p=p*Rz;


le=length(p);
tmp=[[1:le]' p]


p=p+repmat(T,le,1);
if exist('T2')>0
    p(ir',:)=p(ir',:)+repmat(T2,length(ir),1);
end

figure(1)
hold on;box on;
scatter3(p(ir,1),p(ir,2),p(ir,3),'k','filled')
scatter3(p(iw,1),p(iw,2),p(iw,3),'m','filled')
scatter3(p(im,1),p(im,2),p(im,3),'g','filled')
scatter3(p(ic,1),p(ic,2),p(ic,3),'r','filled')
%axis equal

mean_waterline=mean(p(iw,3));


if exist('beamMeas.txt')
    geo=[];
    all_geo=[];
    interZ=[];
    col='bck'
    data=load('beamMeas.txt')
    ind_beam=find(data(:,1)==-1);
    for i=1:length(ind_beam)
        be=ind_beam(i)+1;
        if i<length(ind_beam)
            en=ind_beam(i+1)-1;
        else
            en=length(data);
        end
        vec=[p(im(data(ind_beam(i),3)),:)-p(im(data(ind_beam(i),2)),:)];
        vec=vec/norm(vec);
        %now build profile points
        for j=be:en
            geo{i}(j-be+1,1:3)=p(im(data(ind_beam(i),2)),:)+data(j,1)*vec+data(j,3)*[0 0 -1];
        end
        %scatter3(geo{i}(:,1),geo{i}(:,2),geo{i}(:,3),col(i),'filled')
        all_geo=[all_geo;geo{i}];
    end
    %now build for each i an interp at all positions, so then we can take
    %the mean
    [tmp,k]=sort(all_geo(:,2),'ascend');%sort by y
    yz=[all_geo(k,2) all_geo(k,3)];
    for i=1:length(ind_beam)
        interZ(1:length(yz),i)=interp1(geo{i}(:,2),geo{i}(:,3),yz(:,1));    
    end
    for i=1:length(yz)
       interGeo(i,1)=yz(i,1);
       interGeo(i,2)=nanmean_photrack(interZ(i,:));
    end
    %%%%special
    if exist('G2')>0
       ind=find(interGeo(:,1)>=2.362 & interGeo(:,1)<=19.96);
       interGeo(ind,2)=interGeo(ind,2)+G2;
    end
    %%%%end special
    meanX=mean(p(im,1));
    figure(1)
    scatter3(repmat(meanX,length(yz),1),interGeo(:,1),interGeo(:,2))
    geometry=[repmat(meanX,length(yz),1) interGeo(:,1) interGeo(:,2)]
    
    cd(save_path);
    save_geo=geometry(:,2:3);
    save('geometry.txt','save_geo','-ascii');
    
    figure(2);hold on;box on
    plot(save_geo(:,1),save_geo(:,2))
    scatter(p(im,2), p(im,3),'g','filled')
    scatter(p(iw,2), p(iw,3),'m','filled')
    scatter(p(ir,2), p(ir,3),'k','filled')
    xlabel('y')
    ylabel('z')
    title(['mean water line: ',num2str(round(mean_waterline*1000)/1000),', edit by changing T(3), ~ln 22']);
    
    rg=reshape(geometry(:,2:3),length(geometry)*2,1);
    fid=fopen('discharge_geometry.txt','w');
    for i=1:length(rg)
        fprintf(fid,'%f',rg(i));
        if i<length(rg)
            fprintf(fid,',',rg(i));
        end
end
fclose(fid);
end

figure(1)
xlabel('x')
ylabel('y')
zlabel('z')
title(['mean water line: ',num2str(round(mean_waterline*1000)/1000),', edit by changing T(3), ~ln 22']);

m_coords=[p(im,1) p(im,2) p(im,3)]
cam_pos=[p(ic,1) p(ic,2) p(ic,3)]

cd(save_path);
rm=reshape(m_coords,length(m_coords)*3,1);
fid=fopen('discharge_m_coords.txt','w');
for i=1:length(rm)
    fprintf(fid,'%f',rm(i));
    if i<length(rm)
       fprintf(fid,',',rm(i)); 
    end
end
fclose(fid);



