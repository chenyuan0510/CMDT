function [predy,myscore]=DR_pred_backtrack(trainy,site_table,site_info,rank_site,testx,if_adjust)
% direct prediction based on "introduce_site"
%调用方法：[predy,score]=DR_pred_backtrack(trainy,subset.table,subset.site_info,subset.introduce_site,test);

ylab=unique(trainy);
[row,col]=size(testx);
predy=nan(row,1);
% mylab=unique(site_info{1});
class_sitinfo=zeros(length(unique(site_info{1})),length(site_info));
% for i=1:length(mylab)
%     class_sitinfo(i,1)=mylab(i);
% end
class_sitinfo(:,1)=unique(site_info{1});
seg_info=class_sitinfo;
% seg_info=nan(length(site_info)-1,);
if if_adjust
    %     sum_row=sum(site_table,2);
    %     site_table(2,:)=site_table(2,:)*(sum_row(1)/sum_row(2));
    sum_sample=length(trainy);
    site_table=site_table.*repmat((sum_sample/size(site_table,1))./sum(site_table,2),1,size(site_table,2));
end
for i=2:length(site_info)
    count=0;
    for j=1:length(site_info{i})
        if sum(site_info{i}{j})==0% && isempty(site_info{i}{j})
            count=count+1;
        else
            temlab=unique(site_info{i}{j});
            tem_class=repmat(class_sitinfo(count+1,:),length(temlab),1);
            seg_info_tem=repmat(seg_info(count+1,:),length(temlab),1);
            seg_info_tem(:,i)=j;
%             for k=1:length(temlab)
%                 tem_class(k,i)=temlab(k);
%             end
            tem_class(:,i)=temlab;
            if count==0
                class_sitinfo=[tem_class;class_sitinfo(2:end,:)];
                seg_info=[seg_info_tem;seg_info(2:end,:)];
                count=count+length(temlab);
            else
                class_sitinfo=[class_sitinfo(1:count,:);tem_class;class_sitinfo(count+1:end,:)];
                seg_info=[seg_info(1:count,:);seg_info_tem;seg_info(count+1:end,:)];
                class_sitinfo(count+length(temlab)+1,:)=[];
                seg_info(count+length(temlab)+1,:)=[];
                count=count+length(temlab);
            end
        end
    end
end
% chi2value=nan(1,length(ylab));
myscore=nan(row,length(ylab));
for i=1:row
    pos_row=1:size(class_sitinfo,1);
    for j=1:col
        if j==1
            temk=site_info{1}(testx(i,rank_site(j)));
            tem_pos=pos_row(class_sitinfo(:,j)==temk);
        else
            if length(tem_pos)==1
                break;
            else
                if sum(class_sitinfo(tem_pos,j))~=0
                    temseg=unique(seg_info(tem_pos,j));
                    temk=site_info{j}{temseg}(testx(i,rank_site(j)));
                    tem_pos=tem_pos(class_sitinfo(tem_pos,j)==temk);
                end
            end
            
        end
    end
    %################
%     for myi=1:length(ylab)
%         temtalbe=site_table;
%         temtalbe(myi,tem_pos)=temtalbe(myi,tem_pos)+1;
%         [chi2value(myi),~] = chi2test(temtalbe);
%     end
%     [~,pos]=max(chi2value);
%     predy(i)=ylab(pos);
%     myscore(i,:)=site_table(:,tem_pos)./sum(site_table(:,tem_pos));
    %     myscore(i,:)=chi2value;%myscore(i,:)=chi2value./sum(chi2value);
    %################
    row_s=sum(sum(site_table));
    myI=(-sum((sum(site_table,2)/row_s).*mylog(sum(site_table,2)/row_s))-sum((sum(site_table)/row_s).*...
        mylog(sum(site_table)/row_s))+sum(sum((site_table/row_s).*mylog(site_table/row_s))));
    [~,pos]=max(site_table(:,tem_pos));
    predy(i)=ylab(pos);
    myscore(i,:)=(site_table(:,tem_pos)/sum(site_table(:,tem_pos)));%*(sum(site_table(:,tem_pos))/row_s)*myI;
end