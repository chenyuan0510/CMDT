# CMDT
A Chi-MIC based adaptive multi-branch decision tree

## 1） MATLAB is the tool of NDC;
## 2)  The Chi-MIC should be installed:https://github.com/chenyuan0510/Chi-MIC
## 3)  Usage
 For Binary classification:
     rlt=myDC_DT(trainy,train,test,if_continuous,randnum)
 For Multi-class problem:
    predy,rlt_matrix,myrandnum]=my_DC_DT_multiclass_one_one(train,test,if_continuous,randnum)
