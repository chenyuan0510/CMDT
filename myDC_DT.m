function rlt=myDC_DT(trainy,train,test,if_continuous,randnum)
% for binary classfication
% rlt.predy=predy;
% rlt.myscore=myscore;
% subset.introduce_site=x_site(subset.introduce_site);
% rlt.subset=subset;
% rlt.site_info=site_info;
% rlt.randnum=randnum;
% rlt.del_position=del_position;

x_site=1:size(train,2);
if nargin<5
    [site_info,test_x,x,~,randnum,del_position]=backtracking_site_component(trainy,train,test,if_continuous);
else
    [site_info,test_x,x,~,randnum,del_position]=backtracking_site_component(trainy,train,test,if_continuous,randnum);
end
if any(del_position)
    x(:,del_position==1)=[];site_info(del_position==1')=[];if_continuous(del_position==1)=[];test_x(:,del_position==1)=[];
    x_site(del_position==1)=[];
end
subset=introduce_site(trainy,x,site_info,if_continuous);
[predy,myscore]=DR_pred_backtrack(trainy,subset(1).table,subset(1).site_info,subset(1).introduce_site,test_x,1);
rlt.predy=predy;
rlt.myscore=myscore;
subset(1).introduce_site=x_site(subset(1).introduce_site);
rlt.subset=subset;
rlt.site_info=site_info;
rlt.randnum=randnum;
rlt.del_position=del_position;
