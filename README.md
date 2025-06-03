# TNNLS-2024-P-34330: Near-Optimal Algorithms for Instance-level Constrained k-Center Clustering 

  ## Dependencies
  - Java version: 17.0.1
  - CPLEX: IBM CPLEX2211 [https://www.ibm.com/products/ilog-cplex-optimization-studio](https://www.ibm.com/products/ilog-cplex-optimization-studio)

   ## Install CPLEX
   download CPLEX by the above link
   
    - pip install cplex
    
  ## Data sets
  Download [zip](https://drive.google.com/drive/folders/1lA6-7GxS5VwpiBIqZAZ2tfHnIBkkwG9X) file here and unzip.
  - wine
  - cnae
  - kdd
  - skin
  - wide09
  - covertype

 ## Usage
 > Run the main code:

    chmod +x run.sh
    ./run.sh inputFileName

  - you have to use the above dataset name to replace the _inputFileName_.

  #### Plot the output

  -  Purity, NMI and RI calculate the degree of agreement between an algorithm's clustering result and its labels. In addition, we also report the Cost of $k$-center and runtime of our algorithms.
