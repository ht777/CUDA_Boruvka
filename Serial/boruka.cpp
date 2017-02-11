// Boruvka's algorithm to find Minimum Spanning
// Tree of a given connected, undirected and
// weighted graph
#include <stdio.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdbool.h>
#include <sstream>
using namespace std;
#include <list>
#include <time.h>
//#include <helper_cuda.h>
//#include <helper_functions.h> 
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
 
struct Graph* create_graph(char *filename);

// Function prototypes for union-find (These functions are defined
// after boruvkaMST() )
int find(struct subset subsets[], int i);
void Union(struct subset subsets[], int x, int y);
 
 int cnt;
// The main function for MST using Boruvka's algorithm
void boruvkaMST(struct Graph* graph)
{
	 	int start_time=clock();

    // Get data of given graph
    int V = graph->V, E = graph->E;
    Edge *edge = graph->edge;
 
    // Allocate memory for creating V subsets.
    struct subset *subsets = new subset[V];
 
    // An array to store index of the cheapest edge of
    // subset.  The stored index for indexing array 'edge[]'
    int *cheapest = new int[V];
 
    // Create V subsets with single elements
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
        cheapest[v] = -1;
    }
 
    // Initially there are V different trees.
    // Finally there will be one tree that will be MST
    int numTrees = V;
    int MSTweight = 0;
    // Keep combining components (or sets) until all
    // compnentes are not combined into single MST.
 //  printf("Number of trees %d\n",numTrees);
    while (numTrees > 1)
    {

        // Traverse through all edges and update
        // cheapest of every component
    //    printf("First for loop1");
        for (int i=0; i<E; i++)
        {
            // Find components (or sets) of two corners
            // of current edge
            int set1 = find(subsets, edge[i].src);
            int set2 = find(subsets, edge[i].dest);

            // If two corners of current edge belong to
            // same set, ignore current edge
            if (set1 == set2)
                continue;
 
            // Else check if current edge is closer to previous
            // cheapest edges of set1 and set2
            else
            {
               if (cheapest[set1] == -1 ||
                   edge[cheapest[set1]].weight > edge[i].weight)
                 cheapest[set1] = i;
 
               if (cheapest[set1] == -1 ||
                   edge[cheapest[set2]].weight > edge[i].weight)
                 cheapest[set2] = i;
            }
        }

	//	printf("second for loop1");
        // Consider the above picked cheapest edges and add them
        // to MST
        for (int i=0; i<V; i++)
        {

            // Check if cheapest for current set exists
            if (cheapest[i] != -1)
            {

                int set1 = find(subsets, edge[cheapest[i]].src);

                int set2 = find(subsets, edge[cheapest[i]].dest);

              // if (set1 != set2)
                  //continue;

                MSTweight += edge[cheapest[i]].weight;
                //printf("Edge %d-%d-%d included in MST\n",
                  //     edge[cheapest[i]].src, edge[cheapest[i]].dest,
                    // edge[cheapest[i]].weight);
 
                // Do a union of set1 and set2 and decrease number
                // of trees

                Union(subsets, set1, set2);
                numTrees--;
            }
        }
    }
    int stop_time=clock();
	int runtime=(int)(stop_time-start_time);
 
  //  printf("Weight of MST is %d\n", MSTweight);
    printf("total runtime in nanoseconds on CPU %d\n", int(runtime*8.17));
    return;
}
 
// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E)
{
    Graph* graph = new Graph;
    graph->V = V;
    graph->E = E;
    graph->edge = new Edge[E];
    return graph;
}
int total_weight=0;
struct Graph* create_graph(char *filename){
    char *file = filename;

    char *line_arr = new char[100];
    string line;
    ifstream myfile (file);
    int V, E;
    getline(myfile, line);
    strcpy(line_arr, line.c_str());
 //   printf("string is %s\n",line_arr);
    char *graph_dims = strtok(line_arr, " ");
    V = atoi(graph_dims);
   // printf("Vertices %d\n",V);
    graph_dims = strtok(NULL, " ");
    E = atoi(graph_dims);
   // printf("Edges %d\n",E);

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
      //  printf("dest is %d\n",dest);
        s = strtok(NULL, " ");
        int weight = atoi(s);
        total_weight+=weight;
      //  printf("weight is %d\n",weight);
        g->edge[offset_count].src = src;
        g->edge[offset_count].dest = dest;
        g->edge[offset_count].weight = weight;
		
        offset_count++;
      //  printf ("after ifstream %d\n", offset_count);

    }
    
   // printf ("after graph generation\n");
        printf ("total weight %d\n", total_weight);

    myfile.close();
    delete[] line_arr;
//    	printf ("eof\n");

    return g;
}
 
 int count=0;
// A utility function to find set of an element i
// (uses path compression technique)
int find(struct subset subsets[], int i)
{
	count++;
		if (count>=10)
		return 0;

    // find root and make root as parent of i
    // (path compression)
    //printf("subsets[i].parent is = %d,i = %d",subsets[i].parent,i);
    if (subsets[i].parent != i)
      subsets[i].parent =
             find(subsets, subsets[i].parent);
    return subsets[i].parent;
}
 
// A function that does union of two sets of x and y
// (uses union by rank)
void Union(struct subset subsets[], int x, int y)
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
 
// Driver program to test above functions
int main(int argc, char *argv[])
{
   if (argc < 2) {
    printf("Error: usage: %s <program_file_1> <program_file_2> ...\n",
           argv[0]);
    exit(1);
  }

  printf("boruvka Simulator\n\n");


  //  struct Graph* graph = createGraph(V, E);


	struct Graph* graph = create_graph(argv[1]);
     boruvkaMST(graph);
 
    return 0;
}
