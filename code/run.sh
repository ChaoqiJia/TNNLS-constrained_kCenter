#!/bin/bash

#check
if [ $# -lt 1 ]; then
  echo "check: $0 <inputFilename>"
  exit 1
fi

inputFilename="$1"

# set CPLEX PATH
export CPLEX_HOME="/Applications/CPLEX_Studio2211/cplex"

# set CLASSPATH
export CLASSPATH=".:$CPLEX_HOME/lib/cplex.jar"

javac -cp "$CLASSPATH" addConstraints/OutPut.java

# set dataset
java -cp "$CLASSPATH" \
  -Djava.library.path="$CPLEX_HOME/bin/arm64_osx" \
  addConstraints.OutPut "$inputFilename"
