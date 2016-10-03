function quality=eval_points_beam(vec)

global l12 l14 l23 l34 w1 w2 w3 w4 dist2water water_level render

p1=vec(1:3);
p2=vec(4:6);
p3=vec(7:9);
p4=vec(10:12);
p=[p1;p2;p3;p4];

% % % here follows a set of conditions
con=[];
if isnan(l12)==0
    con=[con abs(norm(p2-p1)-l12)/abs(l12)];%];%
end
if isnan(l34)==0
    con=[con abs(norm(p4-p3)-l34)/abs(l34)];%];%
end
if isnan(l14)==0
    con=[con abs(norm(p4-p1)-l14)/abs(l14)];%];%
end
if isnan(l23)==0
    con=[con abs(norm(p3-p2)-l23)/abs(l23)];%];%
end

%here wwe check the w1..w4 constraints
error_12=abs((w1-w2)-(p1(2)-p2(2)))/max((abs(w1)+abs(w2)),0.5);%;%
error_34=abs((w3-w4)-(-p3(2)+p4(2)))/max((abs(w3)+abs(w4)),0.5);%;%
con=[con error_12 error_34];

%now comes the dist to water part
for i=1:length(dist2water)
    ind1=dist2water{i}(1,1);
    ind2=dist2water{i}(1,2);
    vec=p(ind2,:)-p(ind1,:);
    vec=vec/norm(vec);
    error_water=0;
    for j=2:length(dist2water{i})
        tmp=p(ind1,:)+dist2water{i}(j,1)*vec-dist2water{i}(j,2)*[0 0 1];
        error_water=error_water+abs(tmp(3)-water_level);
    end
    error_water=error_water/(length(dist2water{i})-1);
    con=[con error_water/l34];
    aa=1;
end

quality=sum(con);%0.5*(mean(con)+max(con));


measured=[l12 l34 l14 l23 w1-w2 w3-w4];
modeled =[norm(p2-p1) norm(p4-p3) norm(p4-p1) norm(p3-p2) p1(2)-p2(2) -p3(2)+p4(2)];
compar=[measured;modeled];

if render==1
    aa=1;
end