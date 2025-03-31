# TNNLS-2024-P-34330: Near-Optimal Algorithms for Instance-level Constrained k-Center Clustering 

  ## Dependencies
  - Java version: 17.0.1
  - CPLEX: IBM CPLEX2211 [https://www.ibm.com/products/ilog-cplex-optimization-studio](https://www.ibm.com/products/ilog-cplex-optimization-studio)

   ## Install CPLEX
   download CPLEX by above link
    - pip install cplex
    
  ## Data sets
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

  - you need use the above datasets name replace the _inputFileName_.

  #### Plot the output

  - Use the Cost, Purity, NMI, RI and runtime to calculate the agreement degree between an algorithm's clustering result and its labels. 
