## CMDT

A Chi-MIC based adaptive multi-branch decision tree

**1） MATLAB is the tool of CMDT;**

**2)  The Chi-MIC should be installed:https://github.com/chenyuan0510/Chi-MIC**, and add the search path by ```addpath('D:\your-path-for-Chi-MIC');```, or using the graphical user interface: Home -> Environment -> Set Path -> Add with subfolders -> Save

**3)  Usage**

 For Binary classification:
 
    rlt=myDC_DT(trainy,train_x,test_x,if_continuous,randnum)
 For Multi-class problem:
 
    [predy,rlt_matrix,myrandnum]=my_DC_DT_multiclass_one_one(train,test,if_continuous,randnum)
Input:
1. *trainy* is the labels of training dataset and presented as a column vector, the entries for it can only be 1 or 2; *train_x* and *test_x* are the feature matrixs for training dataset and test dataset, respectively, in which each row is an observation/sample and each column is a variable/feature/attribute. 
2. *train* and *test* are the feature matrixs for training dataset and test dataset, but include the labels at the first column, the entries of labels should be positive integer (1, 2, 3, ...).
3. *if_continuous* is a row vector to indicate the feature types of categorical(a nominal attribute) or numerical(a numeric attribute) data. such as if_continuous=[0,0,1,1,1,0]
means the 1th, 2th and 6th features are nominal attributes, and the remaining variables are numeric attributes.
4. *randnum* is a column vector to scramble the order of training samples, and it is required. We can generate it by:
```randnum=randperm(length(trainy))'; ```
   or 
```randnum=randperm(size(train,1))';```

Output:
1. *rlt* and *rlt_matrix* are the structure variable inclued *predy*, *myscore*, *introduce_site*, *subset*, *site_info*, *randnum*, *del_position*.
we can use ```rlt.predy ``` get the predicted results of test dataset; and use ```rlt.del_position ``` get the features that has been removed by CMDT;

Contact me: chenyuan0510@126.com
