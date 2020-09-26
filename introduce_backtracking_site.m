function [table,last_info,introduce_site]=introduce_backtracking_site(x,site_info,firstsite,table,xx,if_continuous)
% y label must 1 ,2 ,3 .....
% 按小格（2*2）是否显著确定是否要回溯(可以加上每列样本数大于5的限制)
% see also "introduce_site"
% if_continuous(1:end)=0;
introduce_site=firstsite;
last_info{1}=site_info{firstsite};
site=1:size(x,2);
site(firstsite)=[];
is_go=1;
% row=sum(sum(table));
% chi2pva=(-sum((sum(table,2)/row).*mylog(sum(table,2)/row))-sum((sum(table)/row).*...
% mylog(sum(table)/row))+sum(sum((table/row).*mylog(table/row))));
%####
% [chi2pva,~]=chi2test(table);
%####
% H0=sum((-sum((table./repmat(sum(table,1),size(table,1),1)).*...
%         mylog(table./repmat(sum(table,1),size(table,1),1)),1)).*(sum(table,1)./sum(sum(table))));
intr_count=1;
while is_go
    %####
    H0=sum((-sum((table./repmat(sum(table,1),size(table,1),1)).*...
        mylog(table./repmat(sum(table,1),size(table,1),1)),1)).*(sum(table,1)./sum(sum(table))));
    intr_site=0;
    myH=nan(length(site),1);
    Hx=nan(length(site),1);
    temxx3=cell(length(site),1);
    temtable3=cell(length(site),1);
    tem_info3=cell(length(site),1);
    for i=1:length(site)
        for j=1:size(table,2)
            if any(table(:,j)==0)
                if j==1
                    temxx2=xx(:,j);
                    temtable2=table(:,j);
                    tem_info{j}=zeros(1,size(site_info{1},2));%%%%%%%%%%%%
                else
                    temxx2=[temxx2,xx(:,j)];
                    temtable2=[temtable2,table(:,j)];
                    tem_info{j}=zeros(1,size(site_info{1},2));%%%%%%%%%%%%%%
                end
            else
                %                 if i==3;
                %                     a=0;
                %                 end
                [temxx,tem_info{j},useful,temtable]=backtracking_sub(xx(:,j),site_info{site(i)},x(:,site(i)),if_continuous(site(i)));
                if j==1
                    if useful==1
                        temxx2=temxx;
                        temtable2=temtable;
                    else
                        temxx2=xx(:,j);
                        temtable2=table(:,j);
                    end
                else
                    if useful==1
                        temxx2=[temxx2,temxx];
                        temtable2=[temtable2,temtable];
                    else
                        temxx2=[temxx2,xx(:,j)];
                        temtable2=[temtable2,table(:,j)];
                    end
                end
            end
        end
        % using mutual information value
        %         row=sum(sum(temtable2));
        %         chi2pva2=(-sum((sum(temtable2,2)/row).*mylog(sum(temtable2,2)/row))-sum((sum(temtable2)/row).*...
        %         mylog(sum(temtable2)/row))+sum(sum((temtable2/row).*mylog(temtable2/row))));
        % using chi2test value
        %         [chi2pva2,~]=chi2test(temtable2);
        %         if chi2pva2>chi2pva
        %             intr_site=i;
        %             chi2pva=chi2pva2;
        %             temxx3=temxx2;
        %             temtable3=temtable2;
        %             tem_info3=tem_info;
        %         end
        % using gain of ch2testvalue
        %         [chi2pva2,~]=chi2test(temtable2);
        %         increment1=1 - chi2cdf((chi2pva2-chi2pva),(size(table,1)-1)*(size(temtable2,2)-size(table,2)));
        %         if increment1<0.01 && increment1<increment
        %             intr_site=i;
        %             increment=increment1;
        %             temxx3=temxx2;
        %             temtable3=temtable2;
        %             tem_info3=tem_info;
        %         end
        %#############
        % using gain ratio of information
        myH(i)=sum((-sum((temtable2./repmat(sum(temtable2,1),size(temtable2,1),1)).*...
            mylog(temtable2./repmat(sum(temtable2,1),size(temtable2,1),1)),1)).*(sum(temtable2,1)./sum(sum(temtable2))));
        Hx(i)=-sum((sum(temtable2,1)./sum(sum(temtable2))).*mylog(sum(temtable2,1)./sum(sum(temtable2))));
        temxx3{i}=temxx2;
        temtable3{i}=temtable2;
        tem_info3{i}=tem_info;
        if myH(i)<H0
            intr_site=i;
        end
    end
    if intr_site==0
        is_go=0;
    else
        Gain_H=H0-myH;
        index1=find(Gain_H>=mean(Gain_H(Gain_H>0)));
        if isempty(index1)
            [~,intr_site]=max(Gain_H);
        else
            [~,indp]=max(Gain_H(index1)./Hx(index1));
            intr_site=index1(indp);
        end
        intr_count=intr_count+1;
        introduce_site(intr_count)=site(intr_site);
        site(intr_site)=[];
        last_info{intr_count}=tem_info3{intr_site};
        xx=temxx3{intr_site};
        table=temtable3{intr_site};
    end
end
end

function [xx,tem_info,useful,table]=backtracking_sub(sub_xx,tem_info,temx,if_continuous)
table=nan(size(sub_xx,1),length(unique(tem_info)));
xx=cell(size(sub_xx,1),length(unique(tem_info)));
for j=1:size(sub_xx,1)
    for i=1:length(unique(tem_info))
        tem_bases=find(tem_info==i);
        for k=1:length(tem_bases)
            xx{j,i}=[xx{j,i};sub_xx{j}(temx(sub_xx{j})==tem_bases(k))];
        end
        table(j,i)=length(xx{j,i});
    end
end
tem_info0=tem_info;
col=size(table,2);
tem_info=1:col;
if col==2
    [~,pValue] = chi2test(table);
    if pValue>0.01% && all(sum(table,1)>=5)
        useful=0;
    else
        useful=1;
        infolabel=unique(tem_info0);
        for i=1:col
            info2{i}=find(tem_info0==infolabel(i));
        end
        tem_info=nan(size(tem_info0));
        for i=1:length(info2)
            tem_info(info2{i})=i;
        end
        count=0;
    end
else
    if if_continuous
        chi20=0;
        chi2_pValue0=1;
        for i=1:col-1
            temtable=nan(size(table,1),2);
            temtable(:,1)=sum(table(:,1:i),2);
            temtable(:,2)=sum(table(:,i+1:end),2);
            [chi2,chi2_pValue] = chi2test(temtable);
            if chi2>chi20
                chi20=chi2;
                chi2_pValue0=chi2_pValue;
                tem_comb=1:i;% 要用的
                tem_table=temtable;% 要用的
            end
        end
        if chi2_pValue0<=0.01
            useful=1;
            count=1;
            info2=cell(2,1);
            info2{1}=tem_info(tem_comb);
            info2{2}=tem_info(setdiff(1:col,tem_comb));
            tem_rank=nan(length(info2),1);
            for kk=1:length(info2)
                tem_rank(kk)=info2{kk}(1);
            end
            [~,pos]=sort(tem_rank,'ascend');
            info2=info2(pos);
            for i=1:length(info2)
                tem_info(info2{i})=i;
            end
            count_seg=2;
            while count_seg<col%%%
                chi20=0;
                %             chi2_pValue0=1;
                for i=1:length(info2)%
                    for j=1:length(info2{i})%
                        tem_table2=tem_table;
                        tem_table2(:,i)=sum(table(:,info2{i}(1:j)),2);
                        tem_table2=[tem_table2(:,1:i),sum(table(:,info2{i}(j+1:end)),2),tem_table2(:,i+1:end)];%%%%%%
                        [chi2,~] = chi2test(tem_table2);
                        [~,chi2_pValue] = chi2test(tem_table2(:,[i,i+1]));
                        if chi2>chi20 && chi2_pValue<=0.01
                            chi20=chi2;
                            %                         chi2_pValue0=chi2_pValue;
                            seg=i;
                            tem_comb=1:j;
                            tem_table3=tem_table2;
                        end
                    end
                end
                if chi20~=0
                    tem_table=tem_table3;
                    teminfo=info2;
                    info2{seg}=teminfo{seg}(tem_comb);
                    info2=[info2(1:seg);teminfo{seg}(tem_comb(end)+1:end);info2(seg+1:end)];
                    tem_rank=nan(length(info2),1);
                    for kk=1:length(info2)
                        tem_rank(kk)=info2{kk}(1);
                    end
                    [~,pos]=sort(tem_rank,'ascend');
                    info2=info2(pos);
                    for i=1:length(info2)
                        tem_info(info2{i})=i;
                    end
                else
                    break;
                end
                count_seg=count_seg+1;
            end
        else
            useful=0;
        end
    else
        seg_tem=1;
        chi20=0;
        chi2_pValue0=1;
        while col-seg_tem>=seg_tem% 至少有2列以上
            comb=nchoosek(1:col,seg_tem); % 好像多算
            if col-seg_tem==seg_tem
                comb=comb(1:size(comb,1)/2,:);
            end
            temtable=cell(size(comb,1),1);
            chi2=nan(size(comb,1),1);
            chi2_pValue=nan(size(comb,1),1);
            for i=1:size(comb,1)
                temtable{i}=zeros(size(table,1),2);
                temtable{i}(:,1)=sum(table(:,comb(i,:)),2);
                temtable{i}(:,2)=sum(table(:,setdiff(1:col,comb(i,:))),2);
                [chi2(i),chi2_pValue(i)] = chi2test(temtable{i});
            end
            [m_v,pos]=max(chi2);
            if m_v>chi20
                chi20=m_v;
                chi2_pValue0=chi2_pValue(pos);
                tem_comb=comb(pos,:);% 要用的
                tem_table=temtable{pos};% 要用的
            end
            seg_tem=seg_tem+1;
        end
        if chi2_pValue0<=0.01
            useful=1;
            count=1;
            info2=cell(2,1);
            info2{1}=tem_info(tem_comb);
            info2{2}=tem_info(setdiff(1:col,tem_comb));
            tem_rank=nan(length(info2),1);
            for kk=1:length(info2)
                tem_rank(kk)=info2{kk}(1);
            end
            [~,pos]=sort(tem_rank,'ascend');
            info2=info2(pos);
            for i=1:length(info2)
                tem_info(info2{i})=i;
            end
            count_seg=2;
            while count_seg<col%%%
                chi20=0;
                %             chi2_pValue0=1;
                for i=1:length(info2)%已划分的段数
                    if length(info2{i})==2%某一段只有1个分段点
                        tem_table2=tem_table;
                        tem_table2(:,i)=table(:,info2{i}(1));
                        tem_table2=[tem_table2(:,1:i),table(:,info2{i}(2)),tem_table2(:,i+1:end)];%%%%%%
                        [chi2,~] = chi2test(tem_table2);
                        [~,chi2_pValue] = chi2test(table(:,info2{i}));
                        if chi2>chi20 && chi2_pValue<=0.01
                            chi20=chi2;
                            %                         chi2_pValue0=chi2_pValue;
                            seg=i;
                            tem_comb=1;
                            tem_table3=tem_table2;
                        end
                    elseif length(info2{i})>2
                        seg_tem=1;
                        while length(info2{i})-seg_tem>=seg_tem% 至少有2列以上
                            comb=nchoosek(1:length(info2{i}),seg_tem);
                            if length(info2{i})-seg_tem==seg_tem
                                comb=comb(1:size(comb,1)/2,:);
                            end
                            for j=1:size(comb,1)
                                tem_table2=tem_table;
                                tem_table2(:,i)=sum(table(:,info2{i}(comb(j,:))),2);
                                tem_table2=[tem_table2(:,1:i),sum(table(:,info2{i}(setdiff(1:length(info2{i}),comb(j,:)))),2),tem_table2(:,i+1:end)];
                                [chi2,~] = chi2test(tem_table2);
                                [~,chi2_pValue] = chi2test(tem_table2(:,[i,i+1]));
                                if chi2>chi20 && chi2_pValue<=0.01
                                    chi20=chi2;
                                    %                                 chi2_pValue0=chi2_pValue;
                                    seg=i;
                                    tem_comb=comb(j,:);
                                    tem_table3=tem_table2;
                                end
                            end
                            seg_tem=seg_tem+1;
                        end
                    end
                end
                if chi20~=0
                    tem_table=tem_table3;
                    teminfo=info2;
                    info2{seg}=teminfo{seg}(tem_comb);
                    info2=[info2(1:seg);teminfo{seg}(setdiff(1:length(teminfo{seg}),tem_comb));info2(seg+1:end)];%
                    tem_rank=nan(length(info2),1);
                    for kk=1:length(info2)
                        tem_rank(kk)=info2{kk}(1);
                    end
                    [~,pos]=sort(tem_rank,'ascend');
                    info2=info2(pos);
                    for i=1:length(info2)
                        tem_info(info2{i})=i;
                    end
                else
                    break;
                end
                count_seg=count_seg+1;
            end
        else
            useful=0;
        end
    end
end
if useful==1 && count>0
    infolabel=unique(tem_info0);
    info3=cell(1,length(infolabel));
    for i=1:length(infolabel)
        info3{i}=find(tem_info0==infolabel(i));
    end
    for i=1:length(info2)
        info4{i}=cell2mat(info3(info2{i}));
    end
    info2=info4;
    tem_rank=nan(length(info2),1);
    for kk=1:length(info2)
        tem_rank(kk)=info2{kk}(1);
    end
    [~,pos]=sort(tem_rank,'ascend');
    info2=info2(pos);
    tem_info=nan(size(tem_info0));
    for i=1:length(info2)
        tem_info(info2{i})=i;
    end
    xx=cell(size(sub_xx,1),length(unique(tem_info)));
    table=nan(size(sub_xx,1),length(unique(tem_info)));
    for j=1:size(sub_xx,1)
        for i=1:length(unique(tem_info))
            tem_bases=find(tem_info==i);
            for k=1:length(tem_bases)
                xx{j,i}=[xx{j,i};sub_xx{j}(temx(sub_xx{j})==tem_bases(k))];
            end
            table(j,i)=length(xx{j,i});
        end
    end
end
if useful==0
    tem_info=zeros(1,length(tem_info0));
end
end
