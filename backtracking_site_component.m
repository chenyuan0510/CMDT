function [site_info,test_x,x,myendpoint,randnum,del_position]=backtracking_site_component(y,x,test_x,if_continuous,randnum)
% 调用方法：site_info=backtracking_site(trainy,train);
% 前向搜索
% 基于信息增益率确定引入哪个特征
% 用于连续性的组分特征
% 如果第j个x为连续型变量，则if_continuous(j)=1；否则if_continuous(j)=0；
% 结果site_info: the row is the each site position, the column is the segment of each position
% y label must 1 ,2 ,3 .....。y为分类类标，二分类的话，y取值为1,2
% x即train（不含Y）
% After this program,you need run:x(:,del_position==1)=[];site_info(del_position==1')=[];if_continuous(del_position==1)=[];test_x(:,del_position==1)=[];
%
if nargin<5
    num=length(y);
    randnum=randperm(num);
end
site_num=size(x,2);
del_position=zeros(1,site_num);
ylab=unique(y);
% xlab=unique(x(:,1));%%%%%
% site_info=repmat(1:length(xlab),site_num,1);%%%%
site_info=cell(site_num,1);
% num=length(y);
myendpoint=cell(1,site_num);
for myi=1:site_num
    if if_continuous(myi)
        [~,mybestc]=MIC_OIC_chi_1_1_class([y(randnum'),x(randnum',myi)],0.55,5);
        mybestc=sort(mybestc);
        temx=sort(x(:,myi));
        endpoint=nan(length(mybestc),1);
        %         endpoint(1)=min(temx);
        %         endpoint(end)=max(temx);
        for k=1:length(mybestc)
            endpoint(k)=temx(mybestc(k));
        end
        myendpoint{myi}=unique(endpoint);
        %         myendpoint{myi}(1)=myendpoint{myi}(1)-0.1;
        x2=x(:,myi);
        testx2=test_x(:,myi);
        for k=1:length(myendpoint{myi})
            if k==1
                x2(x(:,myi)<=myendpoint{myi}(k),1)=k;
                testx2(test_x(:,myi)<=myendpoint{myi}(k),1)=k;
                if k==length(myendpoint{myi})
                    if ~isempty(x2(x(:,myi)>myendpoint{myi}(k),1))
                        x2(x(:,myi)>myendpoint{myi}(k),1)=k+1;
                        testx2(test_x(:,myi)>myendpoint{myi}(k),1)=k+1;
                    else
                        testx2(test_x(:,myi)>myendpoint{myi}(k),1)=k;
                    end
                end
            else
                x2(x(:,myi)>myendpoint{myi}(k-1) & x(:,myi)<=myendpoint{myi}(k),1)=k;
                testx2(test_x(:,myi)>myendpoint{myi}(k-1) & test_x(:,myi)<=myendpoint{myi}(k),1)=k;
                if k==length(myendpoint{myi})
                    if ~isempty(x2(x(:,myi)>myendpoint{myi}(k),1))
                        x2(x(:,myi)>myendpoint{myi}(k),1)=k+1;
                        testx2(test_x(:,myi)>myendpoint{myi}(k),1)=k+1;
                    else
                        testx2(test_x(:,myi)>myendpoint{myi}(k),1)=k;
                    end
                end
            end
        end
        x(:,myi)=x2;
        test_x(:,myi)=testx2;
    end
    xlab=unique(x(:,myi));
    if length(xlab)<=1
        site_info{myi}=0;
    else
        site_info{myi}=1:length(xlab);
        table=zeros(length(ylab),length(xlab));
        xx=cell(length(ylab),length(xlab));
        for j=1:length(ylab)
            for i=1:length(xlab)
                temindex=find(y==j);
                xx{j,i}=temindex(x(y==j,myi)==i);
                table(j,i)=length(xx{j,i});
            end
        end
                site_info{myi}=backtracking_sub(table,site_info{myi},if_continuous(myi));
%         site_info{myi}=backtracking_sub(table,site_info{myi},0);
    end
    if sum(site_info{myi})==0
        del_position(myi)=1;
    end
end
end
function tem_info=backtracking_sub(table,tem_info,if_continuous)
col=size(table,2);
if col==2
    [~,pValue] = chi2test(table);
    if pValue>0.01% && all(sum(table,1)>=5)
        useful=0;
    else
        useful=1;
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
            info2=cell(2,1);
            info2{1}=tem_info(tem_comb);
            info2{2}=tem_info(setdiff(1:col,tem_comb));
            tem_rank=nan(length(info2),1);
            for kk=1:length(info2)
                tem_rank(kk)=info2{kk}(1);
            end
            [~,pos]=sort(tem_rank,'ascend');
            info2=info2(pos);
            tem_table=tem_table(:,pos);
            for i=1:length(info2)
                tem_info(info2{i})=i;
            end
            count_seg=2;
            while count_seg<col%%%
                chi20=0;
                for i=1:length(info2)%
                    for j=1:length(info2{i})%
                        tem_table2=tem_table;
                        tem_table2(:,i)=sum(table(:,info2{i}(1:j)),2);
                        tem_table2=[tem_table2(:,1:i),sum(table(:,info2{i}(j+1:end)),2),tem_table2(:,i+1:end)];%%%%%%
                        [chi2,~] = chi2test(tem_table2);
                        [~,chi2_pValue] = chi2test(tem_table2(:,[i,i+1]));
                        if chi2>chi20 && chi2_pValue<=0.01
                            chi20=chi2;
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
                    tem_table=tem_table(:,pos);
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
            info2=cell(2,1);
            info2{1}=tem_info(tem_comb);
            info2{2}=tem_info(setdiff(1:col,tem_comb));
            tem_rank=nan(length(info2),1);
            for kk=1:length(info2)
                tem_rank(kk)=info2{kk}(1);
            end
            [~,pos]=sort(tem_rank,'ascend');
            info2=info2(pos);
            tem_table=tem_table(:,pos);
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
                    tem_table=tem_table(:,pos);
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
if useful==0
    tem_info=zeros(1,length(tem_info));
end
end
