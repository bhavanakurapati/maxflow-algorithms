## Program Execution:

A make file complies all the algorithms and to run the program we use:  
`./MaxFlow`

Then it shows:
Please input file name: input.dat (example)
Then the output will be generated in input.out file.

## Motivation and Background:

Maximum Flow is the maximum amount of flow that the network would allow to flow from source to sink. All we concern about Max Flow Problem is How to Maximize the flow in a Network from source to sink. Many Algorithms exists to find the Maximum Flow. Some of them are Dinic’s Algorithm, Ford-Fulkerson Algorithm, Edmonds-Karp Algorithm, Push Relabel Algorithm, Binary blocking flow Algorithm, MPM Algorithm etc., and the applications include Maximum cardinality bipartite Matching, Maximum Flow with vertex Capacities, Maximum number of paths from source to sink, etc. Real World Applications include Airline Scheduling, Baseball elimination, Circulation demand problem.

# Pseudocode and Implementation:

## Ford-Fulkerson Algorithm:

### Input:

`CapacityMatrix (n*n); Adjacency_list (n) // created by DFS node contains the nodes that are on path src //Source_node target //destination_node`

### Output:

`maxflow// Max flow FordFulkerson(): ResMatrix= CapacityMatrix maxflow=0; while(true): shortestPath = DFS (ResMatrix, src, target) // DFS finds shortest path if(shortestPath==NULL): break; else: cur_flow= minflow (shortest_path) flow+= cur_flow; current =t; END IF While (current! =s ): p \_vertex = adjaceny_matrix(current)-> previous ; ResMatrix[p, current] =Res_matrix[p,current]-cur_flow ResMatrix [current, p] =ResMatrix[current, p]+cur_flow current= p; END; END; return maxflow;`  
Dinic’s Algorithm:
Step 1: Read the values from file and copy them into a 2D Array of size n\*n where n is total number of vertices.
Step2: First Vertex that is read from file is considered as source and last vertex that is read is our sink.
Step3: Augment the flow until there is a path from source to sink.
While(bfs() returns false)
{
//…….
}
Step4: Find if more flow is possible or send flows until there are no flows traversing from source to sink. In this case, for every non-zero flow keep adding current path flow to maxflow.
Step4: In the Send Flow,
• Traverse all edges,
• Find min flow from adjacent vertex(u) to sink.,
• For a non-zero flows, update the residual graph flow.
Step6: At the end of send flow loop and bfs loop we will have Maximum Flow.
Edmond-Karp Algorithm:
Pseudocode: Input: capacityMatrix [n x n]: Capacity Matrix adjMatrix [n x n]: Adjacency Matrix src: Source target: sink or target Output: maxFlow: Maximum Flow Rate The Edmonds-Karp: maxFlow = 0 // Initialize the flow to 0 residualMatrix [n x n] // The residual capacity array while true: // Finding the shortest path mininum, augmentPath = BFS(capacityMatrix, adjMatrix, source, target,residualMatrix); if m = 0: break; //if shortest path is zero maxFlow = maxFlow + min // Walk through the augmenting path v = target while v != src: u = P[v]; // Reduce the residual capacity residualMatrix [u, v] = residualMatrix [u, v] – min // Increase the residual capacity of reverse edges
residualMatrix [v, u] = residualMatrix [v ,u] + min v = u return maxFlow;
General Push Re-Label:
Step 1: Read the values from file and copy them into a 2D Array of size n\*n where n is total number of vertices.
Step2: First Vertex that is read from file is considered as source and last vertex that is read is our sink.
Step3: Once the values are read into an array keep adding the edges and weights into an capacity array.
Step4: Initialize height, flow to zero and flow of every edge to 0, Also initialize all vertices adjacent to source to capacity.
for (int i = 0; i < n; i++) {
f[s][i] = cap[s][i];
f[i][s] = -f[s][i];
e[i] = cap[s][i];
}
Step5: while you have a chance to push or relabel do the loop, until excess flow do push or relabel.
While(vertex has over flow)
{
if( you cannot push)
{
//Relabel it;
}
}
Step6: In order to relabel, you have to find the adjacent vertex with min. height
While(you find min height)
{
If(edge flow==edge capacity)
{
//update min height
}
}
Step7: At the end of excess flow loop, we will have max flow.
