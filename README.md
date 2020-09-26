# CMDT
A Chi-MIC based adaptive multi-branch decision tree

## 1） MATLAB is the tool of CMDT;
## 2)  The Chi-MIC should be installed:https://github.com/chenyuan0510/Chi-MIC
## 3)  Usage
 For Binary classification:
 
     rlt=myDC_DT(trainy,train_x,test_x,if_continuous,randnum)
 For Multi-class problem:
 
    [predy,rlt_matrix,myrandnum]=my_DC_DT_multiclass_one_one(train,test,if_continuous,randnum)
Input:
1. trainy is the labels of training dataset; train_x and test_x is the feature matrix for training dataset and test dataset, respectively, in which each row is an observation/sample and each column is a variable/feature/attribute.
2. if_continuous is a row vector to indicate the feature types of categorical(a nominal attribute) or numerical(a numeric attribute) data. such as if_continuous=[0,0,1,1,1,0]
means the 1th, 2th and 6th features are nominal attributes, and the remaining variables are numeric attributes.
3. randnum is a column vector to scramble the order of training samples, and it is required. We can generate it by :randnum=randperm(length(trainy))'; or randnum=randperm(size(train,1))';

Output:
