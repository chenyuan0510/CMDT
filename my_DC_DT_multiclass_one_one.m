function [predy,rlt_matrix,myrandnum]=my_DC_DT_multiclass_one_one(train,test,if_continuous,randnum)
% y label must 1 ,2 ,3 .....; y为分类类标
% predy: 1th col: PK, 2th col : vote, 3th col: weighted vote
%  one vs one of PK 

ylab=unique(train(:,1));
comb=nchoosek(1:length(ylab),2);
rlt=nan(size(test,1),size(comb,1));
score=nan(size(test,1),size(comb,1)*2);
myrandnum=cell(size(comb,1),1);
for i=1:size(comb,1)
    mytrain=[train(train(:,1)==comb(i,1),2:end);train(train(:,1)==comb(i,2),2:end)];
    mytest=test(:,2:end);
    trainy=ones(size(mytrain,1),1);
    trainy(sum(train(:,1)==comb(i,1))+1:end)=2;
    if nargin<4 
        temrlt=myDC_DT(trainy,mytrain,mytest,if_continuous);
    else
        temrlt=myDC_DT(trainy,mytrain,mytest,if_continuous,randnum{i});
    end
    myrandnum{i}=temrlt.randnum;
    index1=temrlt.predy==1;
    index2=temrlt.predy==2;
    rlt(index1,i)=comb(i,1);
    rlt(index2,i)=comb(i,2);
    score(:,i*2-1:i*2)=temrlt.myscore;
end
rlt_matrix=rlt;
nindex=rlt(:,1);
rowx=length(nindex);
temcol=1:rowx;
for i=3:length(ylab)
    temindex=[nindex,zeros(rowx,1)+i];
    col_ind=index2tril(temindex',length(ylab));
    nindex = rlt(sub2ind(size(rlt),temcol',col_ind'));
end
predy=nindex;
for i=1:size(test,1)
    tem_count_c=tabulate([rlt(i,:),1:length(ylab)]);
    [~,predy(i,2)]=max(tem_count_c(:,2));
end
index_tem=reshape(comb',1,[]);
lab_score=nan(size(test,1),length(ylab));
for i=1:length(ylab)
    lab_score(:,i)=sum(score(:,index_tem==i),2);
end
[~,predy(:,3)]=max(lab_score,[],2);
% predy(:,3)=ylab(pos);
