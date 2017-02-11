
#include <stdio.h>
#include <string>
#include <iostream>
//#include <string.h>
//#include <assert.h>
//#include <stdlib.h>
#include <cuda_runtime.h>
//#include <time.h>
#include <stdlib.h>



#include <fstream>
#include <fstream>
#include <iomanip>
// includes, project
////////////////////////////////////////////////////////////////////////////////
// declarations, forward

#define WIDTH 32 
int total_weight=0;
//extern "C"
// a structure to represent a weighted edge in graph
struct Edge
{
    int src, dest, weight;
};
 
// a structure to represent a connected, undirected
// and weighted graph as a collection of edges.
struct Graph
{
    // V-> Number of vertices, E-> Number of edges
    int V, E;
 
    // graph is represented as an array of edges.
    // Since the graph is undirected, the edge
    // from src to dest is also edge from dest
    // to src. Both are counted as 1 edge here.
    Edge* edge;
};
 
// A structure to represent a subset for union-find
struct subset
{
    int parent;
    int rank;
};
struct Graph* createGraph(int V, int E)
{
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->V = V;
    graph->E = E;
    graph->edge = (Edge *)malloc(E*sizeof(Edge));
    return graph;
}
struct Graph* create_graph(char *filename);
__device__ int find(struct subset subsets[], int i);
__device__ void Union(struct subset subsets[], int x, int y);


__device__ int find(struct subset subsets[], int i)
{
	

    // find root and make root as parent of i
    // (path compression)
    //printf("subsets[i].parent is = %d,i = %d",subsets[i].parent,i);

	for (int k=i;k<100;k++)
	{
	 if (subsets[i].parent == i)
	    return subsets[i].parent;
	else
		continue;
	}
 		//printf("inside kernel\n");
		    return subsets[i].parent;


}
 
// A function that does union of two sets of x and y
// (uses union by rank)
__device__ void Union(struct subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);
 
    // Attach smaller rank tree under root of high
    // rank tree (Union by Rank)
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;
 
    // If ranks are same, then make one as root and
    // increment its rank by one
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}
 
__global__ void find_subsets(struct subset* subsets, Edge* edge, int* cheapest_in, int* cheapest_out, unsigned long long* runtime)
{
	    int tid=threadIdx.x;
	unsigned long long start_time=clock64();
            // Find components (or sets) of two corners
            // of current edge
	   // printf("Kernel Address %u  and tid:%d\n",edge,tid);
		//printf ("edgesrc %d\n", edge[tid].src);
            int set1 = find(subsets, edge[tid].src);
 		//printf("inside kernel\n");
            int set2 = find(subsets, edge[tid].dest);

            // If two corners of current edge belong to
            // same set, ignore current edge
            if (set1 == set2)
               return;

            // Else check if current edge is closer to previous
            // cheapest edges of set1 and set2
           
               if (cheapest_in[set1] == -1 ||
                   edge[cheapest_in[set1]].weight > edge[tid].weight)
                 cheapest_out[set1] = tid;
 
               if (cheapest_in[set1] == -1 ||
                   edge[cheapest_in[set2]].weight > edge[tid].weight)
                 cheapest_out[set2] = tid;

		//printf ("thread id %d\n",tid);

  	unsigned long long stop_time=clock64();
	runtime[tid]=(unsigned long long)(stop_time-start_time);//runtime for each thread

}

 
__global__ void subsets_weight(struct subset* subsets, Edge* edge, int* cheapest_in, int* MSTweight, int* numTrees, unsigned long long* runtime)
{
	    int tid=threadIdx.x;

	unsigned long long start_time=clock64();
            // Check if cheapest for current set exists
            if (cheapest_in[tid] != -1)
            {

                int set1 = find(subsets, edge[cheapest_in[tid]].src);

                int set2 = find(subsets, edge[cheapest_in[tid]].dest);

               // if (set1 != set2)
                 //   continue;

                MSTweight[tid] = edge[cheapest_in[tid]].weight;
	//	printf ("mst in kernel %d\n", MSTweight[tid]);
              //  printf("Edge %d-%d-%d included in MST\n",
                //       edge[cheapest[i]].src, edge[cheapest[i]].dest,
                  //     edge[cheapest[i]].weight);
 
                // Do a union of set1 and set2 and decrease number
                // of trees

                Union(subsets, set1, set2);//--
                numTrees[tid]=0;
		//printf ("inside ker\n");
            }
	else
	{
		MSTweight[tid]=0;
		numTrees[tid]=1;
	}
  
	unsigned long long stop_time=clock64();
	runtime[tid]=(unsigned long long)(stop_time-start_time);//runtime for each thread

}

/*

void boruvkaMST(struct Graph* graph)
{


}	


*/
/**
 * Host main routine
 */
//#define V 100
//#define E 100
   #define EDGE E
   #define VERTEX V
int main(int argc, char *argv[]) 
{

   if (argc < 2) {
    printf("Error: usage: %s <program_file_1> <program_file_2> ...\n",
           argv[0]);
    exit(1);
   }

    printf("boruvka Simulator\n\n");

  //  struct Graph* graph = create_graph(argv[1]);
	

    char *file = argv[1];
    char *line_arr = (char *)malloc(100*sizeof(char));
    std::string line;
    std::ifstream myfile (file);
   int V, E;
    getline(myfile, line);
    strcpy(line_arr, line.c_str());
   // printf("string is %s\n",line_arr);
    char *graph_dims = strtok(line_arr, " ");
    V = atoi(graph_dims);
    //printf("Vertices %d\n",V);
    graph_dims = strtok(NULL, " ");
    E = atoi(graph_dims);
  //  printf("Edges %d\n",E);
    Edge *edge = (Edge *)malloc(2*E*sizeof(Edge));

    struct Graph* g = createGraph(V, 2*E);

 //   Graph g = instantiate_graph(V, 2 * E);

    int offset_count = 0;
    int temp;
    
    while(getline (myfile, line))
    {
		int src=0, dest=0;

        strcpy(line_arr, line.c_str());
      //  printf("string is %s\n",line_arr);
        char* s = strtok(line_arr, " ");
        src = atoi(s);
        s = strtok(NULL, " ");
        dest = atoi(s);
        //printf("dest is %d\n",dest);
        s = strtok(NULL, " ");
        int weight = atoi(s);
        total_weight+=weight;
        //printf("weight is %d\n",weight);
        g->edge[offset_count].src = src;
        g->edge[offset_count].dest = dest;
        g->edge[offset_count].weight = weight;
		
        offset_count++;
       // printf ("after ifstream %d\n", offset_count);

    }
    
   // printf ("after graph generation\n");
        printf ("total weight %d\n", total_weight);

    myfile.close();
    free(line_arr);
    	//printf ("eof\n");



	cudaError_t err = cudaSuccess;
    // Get data of given graph
   // int V = graph->V, E = graph->E;
    edge = g->edge;	
 
	//printf ("Edge details %d\n", edge[5].weight);
    // Allocate memory for creating V subsets.
  //  struct subset *subsets = new subset[V];
    
    struct subset *subsets = (subset *)malloc(V*sizeof(subset));
 
    // An array to store index of the cheapest edge of
    // subset.  The stored index for indexing array 'edge[]'
  //  int *cheapest = new int[V];
    int *cheapest = (int *)malloc(V*sizeof(int));
 

    // Create V subsets with single elements
    for (int v = 0; v < V; ++v)//--
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
        cheapest[v] = -1;
    }
    
//printf ("inside\n");
 
    // Initially there are V different trees.
    // Finally there will be one tree that will be MST
    int numTrees = V;
    int MSTweight = 0;
    unsigned long long net_runtime=0;//stores the total execution time
   subset *d_subsets;//=(subset *)malloc(V*sizeof(subset));
	//printf ("inside1\n");	
   err = cudaMalloc((void**)&d_subsets, V*sizeof(subset)); // TODO: Allocate context in GPU device memory
//printf ("inside2\n"); 
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device subset data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }


//printf ("inside1\n");	
    Edge *d_edge;//=(Edge *)malloc(sizeof(Edge));

//printf ("inside2\n"); 
    err = cudaMalloc((void**)&d_edge, 2*E*sizeof(Edge)); // TODO: Allocate context in GPU device memory 
//printf ("inside2\n"); 
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device edge data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    Edge *d_graph;//=(Edge *)malloc(sizeof(Edge));

//printf ("inside2\n"); 
/*
    err = cudaMalloc((void**)&d_graph, E*sizeof(Graph)); // TODO: Allocate context in GPU device memory 
//printf ("inside2\n"); 
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device graph data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
*/


    int *d_cheapest;// = (int *)malloc(100*sizeof(int));

    err = cudaMalloc((void**)&d_cheapest, V*sizeof(int)); // TODO: Allocate context in GPU device memory 
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device cheapest data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
// printf ("inside2\n");
   int* d_cheapest_out;// = (int*)malloc(100*sizeof(int));

    err = cudaMalloc((void**)&d_cheapest_out, V*sizeof(int)); // TODO: Allocate context in GPU device memory 
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device cheapest out data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }


//printf("Here\n");
/*-------------runtime related memory allocation----------*/
    unsigned long long* d_runtime;
    int r_size = EDGE*VERTEX*sizeof(unsigned long long);
    unsigned long long* runtime = (unsigned long long*)malloc(r_size);
    memset(runtime, 0, r_size);
    cudaMalloc((void**)&d_runtime, r_size);
/*-------------------------xxxxxxxxx-----------------------*/


    int* d_MSTweight;
    int mst_size = VERTEX*sizeof(int);
    int* mst_weight = (int*)malloc(mst_size);
    memset(mst_weight, 0, mst_size);
    cudaMalloc((void**)&d_MSTweight, mst_size);
/*---------------------------------xxxxxxxxxx----------------*/

    int* d_numTrees;
    int trees_size = VERTEX*sizeof(int);
    int* num_of_trees = (int*)malloc(trees_size);
    memset(num_of_trees, 0, trees_size);
    cudaMalloc((void**)&d_numTrees, trees_size);

/*------------------xxxxxxxxxxxxxx--------------------*/

// int EDGE = E;
// int VERTEX = V;

    // Keep combining components (or sets) until all
    // compnentes are not combined into single MST.
 //  printf("Number of trees %d\n",numTrees);
    while (numTrees > 1)//--
    {

        // Traverse through all edges and update
        // cheapest of every component
    //    printf("First for loop1");

    err = cudaMemcpy(d_subsets, subsets, V*sizeof(subset), cudaMemcpyHostToDevice);// TODO: Copy the input/updated context to GPU
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy subset data from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(d_edge, edge, 2*E*sizeof(Edge), cudaMemcpyHostToDevice);// TODO: Copy the input/updated context to GPU
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy edge data from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }



    err = cudaMemcpy(d_cheapest, cheapest, V*sizeof(int), cudaMemcpyHostToDevice);// TODO: Copy the input/updated context to GPU
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy cheapest data from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
/*
    err = cudaMemcpy(d_numTrees, numTrees, sizeof(int), cudaMemcpyHostToDevice);// TODO: Copy the input/updated context to GPU
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy numTrees data from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
*/



    dim3 dimGrid1(EDGE,1, 1);
    dim3 dimBlock1(1, 1, 1);
	//printf("Address %u\n",d_edge);

// Call the kernel function
    find_subsets<<<dimBlock1,dimGrid1>>>(d_subsets, d_edge, d_cheapest, d_cheapest_out, d_runtime);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Kernel-1 execution failed (error code %s)\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    cudaThreadSynchronize();

    //printf("Copy between kernel data from the CUDA device to the host memory\n");//copying the updated context from GPU to CPU
    err = cudaMemcpy(cheapest,d_cheapest_out, V*sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy between kernel data from device to host (error code %s)!\n", cudaGetErrorString(err));
     }

    err = cudaMemcpy(d_cheapest, cheapest, V*sizeof(int), cudaMemcpyHostToDevice);// TODO: Copy the input/updated context to GPU
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy cheapest data second time from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    cudaMemcpy(runtime, d_runtime, r_size, cudaMemcpyDeviceToHost);
    cudaThreadSynchronize();
    
    unsigned long long elapsed_time_EDGE = 0;
    for(int i = 0; i < EDGE; i++)
        if(elapsed_time_EDGE < runtime[i])
            elapsed_time_EDGE = runtime[i];//highest execution time among all the simultaneously running threads
    net_runtime += elapsed_time_EDGE;// calculates the total execution time, each time when the kernel is executed


    dim3 dimGrid2(VERTEX,1, 1);
    dim3 dimBlock2(1, 1, 1);
    // Call the kernel function
    subsets_weight<<<dimBlock2,dimGrid2>>>(d_subsets, d_edge, d_cheapest, d_MSTweight, d_numTrees, d_runtime);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Kernel-2 execution failed (error code %s)\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    cudaThreadSynchronize();

    err = cudaMemcpy(subsets, d_subsets, V*sizeof(subset), cudaMemcpyDeviceToHost);// TODO: Copy the input/updated context to GPU
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy subset data from device to host device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
   cudaMemcpy(runtime, d_runtime, r_size, cudaMemcpyDeviceToHost);//--
    cudaThreadSynchronize();
    
    unsigned long long elapsed_time_VER = 0;
    for(int j = 0; j < VERTEX; j++)
        if(elapsed_time_VER < runtime[j])
            elapsed_time_VER = runtime[j];//highest execution time among all the simultaneously running threads
    net_runtime += elapsed_time_VER;// calculates the total execution time, each time when the kernel is executed

    cudaMemcpy(num_of_trees, d_numTrees, trees_size, cudaMemcpyDeviceToHost);
    cudaThreadSynchronize();
    

    for(int k = 0; k < VERTEX; k++)
	numTrees += num_of_trees[k];
    numTrees = numTrees/VERTEX;
    //printf("num trees %d\n", numTrees);
//numTrees--;
   //printf ("num of trees %d\n", numTrees);
  // printf ("mst weight %d\n", d_MSTweight);
    cudaMemcpy(mst_weight, d_MSTweight, mst_size, cudaMemcpyDeviceToHost);
    cudaThreadSynchronize();
    
   int temp_weight=0;
    for(int n = 0; n < VERTEX; n++)
	//if (temp_weight<mst_weight[n])
	//	temp_weight=mst_weight[n];
	MSTweight += mst_weight[n];
    }
 cudaFree(d_subsets);
 cudaFree(d_edge);
 cudaFree(d_cheapest);
 cudaFree(d_cheapest_out);
   // printf("Weight of MST is %d\n", MSTweight);
    printf("total run time in nanoseconds on GPU %d\n",int(net_runtime*8.17)/64);
    return;

	return 0;
}







///--------working function--------///

/*
struct Graph* create_graph(char *filename)
{
    char *file = filename;
//(Edge *)malloc(sizeof(Edge))
    char *line_arr = (char *)malloc(100*sizeof(char));
    std::string line;
    std::ifstream myfile (file);
   // int V, E;
    getline(myfile, line);
    strcpy(line_arr, line.c_str());
    printf("string is %s\n",line_arr);
    char *graph_dims = strtok(line_arr, " ");
    //V = atoi(graph_dims);
    printf("Vertices %d\n",V);
    graph_dims = strtok(NULL, " ");
    //E = atoi(graph_dims);
    printf("Edges %d\n",E);

    struct Graph* g = createGraph(V, 2*E);

 //   Graph g = instantiate_graph(V, 2 * E);

    int offset_count = 0;
    int temp;
    
    while(getline (myfile, line))
    {
		int src=0, dest=0;

        strcpy(line_arr, line.c_str());
        printf("string is %s\n",line_arr);
        char* s = strtok(line_arr, " ");
        src = atoi(s);
        s = strtok(NULL, " ");
        dest = atoi(s);
        printf("dest is %d\n",dest);
        s = strtok(NULL, " ");
        int weight = atoi(s);
        total_weight+=weight;
        printf("weight is %d\n",weight);
        g->edge[offset_count].src = src;
        g->edge[offset_count].dest = dest;
        g->edge[offset_count].weight = weight;
		
        offset_count++;
        printf ("after ifstream %d\n", offset_count);

    }
    
    printf ("after graph generation\n");
        printf ("total weight %d\n", total_weight);

    myfile.close();
    free(line_arr);
    	printf ("eof\n");

    return g;
}*/
