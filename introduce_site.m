function subset=introduce_site(y,x,site_info,if_continuous)
% 调用方法：subset=introduce_site(trainy,train,site_info);
% the site_info is the output of "backtracking_site"
% see also "DR_pred_backtrack"
tabl=tabulate(y);
H0=-sum((tabl(:,3)/100).*mylog((tabl(:,3)/100)));
site_num=length(site_info);
mytable=cell(site_num,1);
myxx=cell(site_num,1);
% pval=nan(site_num,1);
myH=nan(site_num,1);
Hx=nan(site_num,1);
for myi=1:site_num
    table=zeros(2,length(unique(site_info{myi})));
    xx=cell(2,length(unique(site_info{myi})));
    for j=1:length(unique(y))
        for i=1:length(unique(site_info{myi}))
            tem_bases=find(site_info{myi}==i);
            for k=1:length(tem_bases)
                %             table(j,i)=table(j,i)+sum(x(y==j,firstsite)==tem_bases(k));
                temindex=find(y==j);
                xx{j,i}=[xx{j,i};temindex(x(y==j,myi)==tem_bases(k))];
            end
            table(j,i)=length(xx{j,i});
        end
    end
    myH(myi)=sum((-sum((table./repmat(sum(table,1),size(table,1),1)).*...
        mylog(table./repmat(sum(table,1),size(table,1),1)),1)).*(sum(table,1)./sum(sum(table))));
    Hx(myi)=-sum((sum(table,1)./sum(sum(table))).*mylog(sum(table,1)./sum(sum(table))));
    %     [pval(myi),~]=chi2test(table);
    mytable{myi}=table;
    myxx{myi}=xx;
end
Gain_H=H0-myH;
index1=find(Gain_H>=mean(Gain_H));
[~,indp]=max(Gain_H(index1)./Hx(index1));
% [~,indp]=sort(pval,'descend');
ini_site=1:site_num;
int_count=1;
while ~isempty(ini_site)
    [temtable3,last_info,introduce_site]=introduce_backtracking_site(x,...
        site_info,index1(indp(1)),mytable{ini_site(index1(indp(1)))},myxx{ini_site(index1(indp(1)))},if_continuous);
    subset(int_count).table=temtable3;
    subset(int_count).site_info=last_info;
    subset(int_count).introduce_site=ini_site(introduce_site);
    ini_site(introduce_site)=[];
    Gain_H(introduce_site)=[];
    Hx(introduce_site)=[];
    if length(ini_site)>1
        index1=find(Gain_H>=mean(Gain_H));
    elseif length(ini_site)==1
        index1=1;
    else
        index1=[];
    end
    [~,indp]=max(Gain_H(index1)./Hx(index1));
    %     pval(introduce_site)=[];
    %     [~,indp]=sort(pval,'descend');
    x(:,introduce_site)=[];
    site_info(introduce_site)=[];
    int_count=int_count+1;
end
subset=subset';