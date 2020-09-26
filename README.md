# CMDT
A Chi-MIC based adaptive multi-branch decision tree

## 1） MATLAB is the tool of CMDT;
## 2)  The Chi-MIC should be installed:https://github.com/chenyuan0510/Chi-MIC
## 3)  Usage
 For Binary classification:
 
     rlt=myDC_DT(trainy,train_x,test_x,if_continuous,randnum)
 For Multi-class problem:
 
    [predy,rlt_matrix,myrandnum]=my_DC_DT_multiclass_one_one(train,test,if_continuous,randnum)
1. trainy is the labels of training dataset; train_x and test_x is the feature matrix for training dataset and test dataset, respectively, in which each row is an observation/sample and each column is a variable/feature/attribute.
2. 
