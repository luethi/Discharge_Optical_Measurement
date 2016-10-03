function result=nanmean_photrack(vec)

result=NaN;

try
    ind=find(isnan(vec)==0 & vec~=Inf);
    result=mean(vec(ind));
catch
end