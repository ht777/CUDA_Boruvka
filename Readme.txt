->Extract Group_02_PBA.tar.gz

Executing Parallel Code
->Go to Parallel directory
->Run command $nvcc Boruvka.cu OutGPU
->Execute the exeutable file with the input file
->$./OutGPU graph_distinct_1000.txt
->(for 500) $./OutGPU graph_distinct_500.txt


Executing Serial Code
->Go to Serial directory
->Run command $nvcc -o Boruvka.cpp OutCPU
->Execute the exeutable file with the input file
->$./OutCPU graph8_distinct_1000.txt 
->(for 500) $./OutCPU graph_distinct_500.txt
