#include <iostream>
#include <chrono>
#include <fstream>
#include <ctime>
#include <sstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <list>
#include <climits>
#include <cstring>

using namespace std;

// 2D array 
typedef vector<vector<int> > TwoDimArray;

#define MAXN  1001
#define MAXM  20002

TwoDimArray weightGraph;
int vertexNum = 0;
int source = 0;
int sink = 0;
double elapsedTime = 0.0;

// Ford-Fulkerson algorithm  
int FordFlow[MAXN][MAXN];
bool FordVisited[MAXN];
int FordDis[MAXN][MAXN];
int FordFulkersonDFS(int s, int t, int minimum)
{
	FordVisited[s] = true;
	// if source and sink is same
	if (s == t)
		return minimum;
	for (int i = 0; i < vertexNum; i++)
	{
		int flow_capacity = FordDis[s][i] - FordFlow[s][i];
		if (!FordVisited[i] && flow_capacity > 0)
		{
			// find min capacity in dfs path
			if (int sent = FordFulkersonDFS(i, t, min(minimum, flow_capacity)))
			{
				// adjust the capacity
				FordFlow[s][i] += sent;
				FordFlow[i][s] -= sent;
				return sent;
			}
		}
	}
	return false;
}

// Returns the maximum flow from s to t in the given graph
int FordFulkersonMaxFlow()
{
	for (int i = 0; i < vertexNum; i++)
	{
		for (int j = 0; j < vertexNum; j++)
			FordDis[i][j] = weightGraph[i][j];
	}

	auto start = chrono::steady_clock::now();

	int maxFlow = 0;
	// While ther is augmenting path , from s and t, with positive flow capacity
	while (int sent = FordFulkersonDFS(source, sink, INT_MAX))
	{
		maxFlow += sent;
		// Reset visited array, for searching next path
		memset(FordVisited, 0, sizeof(FordVisited));
	}

	auto end = chrono::steady_clock::now();
	elapsedTime = chrono::duration<double>(end - start).count();
	return maxFlow;
}

// Dinic's algorithm 
int DinicH[MAXN + 1], DinicNex[MAXM], DinicNa[MAXM], DinicW[MAXM];
int DinicDep[MAXN];
int DinicDFS(int tmpS, int tmpT, int f)
{
	if (tmpS == tmpT)
		return f;
	int d;
	for (int i = DinicH[tmpS]; i; i = DinicNex[i])
	{
		int v = DinicNa[i];
		if (DinicW[i]>0 && DinicDep[v] == DinicDep[tmpS] + 1 && (d = DinicDFS(v, tmpT, min(f, DinicW[i]))))
		{
			DinicW[i] -= d;
			DinicW[i ^ 1] += d;
			return d;
		}
	}
	return 0;
}

bool DinicBFS()  
{
	queue<int> q; 
	q.push(source + 1);
	memset(DinicDep, -1, sizeof(DinicDep));
	DinicDep[source + 1] = 0;
	while (!q.empty()) 
	{
		int u = q.front(); 
		q.pop();
		for (int i = DinicH[u]; i; i = DinicNex[i])
		{
			int v = DinicNa[i];
			if (DinicW[i]>0 && DinicDep[v] == -1)
			{
				DinicDep[v] = DinicDep[u] + 1;
				q.push(v);
			}
		}
	}
	return DinicDep[sink + 1] != -1;
}

int DinicMaxFlow()
{
	int index = 1;
	for (int i = 0; i < vertexNum; i++)
	{
		for (int j = 0; j < vertexNum; j++)
		{
			if (weightGraph[i][j] <= 0)
				continue;
			int u = i + 1;
			int v = j + 1;
			int weight = weightGraph[i][j];
			DinicNa[++index] = v; DinicNex[index] = DinicH[u]; DinicH[u] = index; DinicW[index] = weight;
			DinicNa[++index] = u; DinicNex[index] = DinicH[v]; DinicH[v] = index;
		}
	}

	auto start = chrono::steady_clock::now();

	int maxFlow = 0;
	while (DinicBFS())
	{
		while (1)
		{
			int f = DinicDFS(source + 1, sink + 1, INT_MAX);
			if (!f)
				break;
			maxFlow += f;
		}
	}
	auto end = chrono::steady_clock::now();
	elapsedTime = chrono::duration<double>(end - start).count();
	return maxFlow;
}

// Edmond-Karp algorithm 
int EdmondDis[MAXN][MAXN];
int EdmondVis[MAXN], EdmondPre[MAXN];
bool EdmondBFS()
{
	memset(EdmondVis, 0, sizeof(EdmondVis));
	memset(EdmondPre, -1, sizeof(EdmondPre));
	queue<int> q;
	q.push(source + 1);
	EdmondVis[source + 1] = 1;
	while (!q.empty())
	{
		int p = q.front();
		q.pop();
		for (int i = 1; i <= vertexNum; i++)
		{
			if (!EdmondVis[i] && EdmondDis[p][i] > 0)
			{
				EdmondPre[i] = p;
				EdmondVis[i] = 1;
				if (i == sink + 1)
					return true;
				q.push(i);
			}
		}
	}
	return false;
}

int EdmondKarpMaxFlow()
{
	for (int i = 0; i < vertexNum; i++)
	{
		for (int j = 0; j < vertexNum; j++)
		{
			if (weightGraph[i][j] <= 0)
				continue;
			int u = i + 1;
			int v = j + 1;
			int weight = weightGraph[i][j];
			EdmondDis[u][v] += weight;
		}
	}

	auto start = chrono::steady_clock::now();
	int maxFlow = 0;
	while (EdmondBFS())
	{
		int d = INT_MAX;
		int u = sink + 1;
		while (EdmondPre[u] != -1)
		{
			d = min(d, EdmondDis[EdmondPre[u]][u]);
			u = EdmondPre[u];
		}
		u = sink + 1;
		while (EdmondPre[u] != -1)
		{
			EdmondDis[EdmondPre[u]][u] -= d;
			EdmondDis[u][EdmondPre[u]] += d;
			u = EdmondPre[u];
		}
		maxFlow += d;
	}
	auto end = chrono::steady_clock::now();
	elapsedTime = chrono::duration<double>(end - start).count();
	return maxFlow;
}

// Generic Push Relabel algorithm 
int PushRelabelDis[MAXN][MAXN], PushRelabelFlow[MAXN][MAXN];
int PushRelabelH[MAXN], PushRelabelEss[MAXN];

void PreFlow() 
{
	for (int i = source + 1; i <= sink + 1; ++ i) 
	{
		if (PushRelabelDis[source + 1][i] > 0) 
		{
			const int flow = PushRelabelDis[source + 1][i];
			PushRelabelFlow[source + 1][i] += flow;
			PushRelabelFlow[i][source + 1] = - PushRelabelFlow[source + 1][i];
			PushRelabelEss[source + 1] -= flow;
			PushRelabelEss[i] += flow;
		}
	}
}

void Push(int start, int end) 
{
	int flow = PushRelabelEss[start] > (PushRelabelDis[start][end] - PushRelabelFlow[start][end]) ? (PushRelabelDis[start][end] - PushRelabelFlow[start][end]) : PushRelabelEss[start];
	PushRelabelFlow[start][end] += flow;
	PushRelabelFlow[end][start] = - PushRelabelFlow[start][end];
	PushRelabelEss[start] -= flow;
	PushRelabelEss[end] += flow;
}

bool ReLabel(int index) 
{
	int minH = INT_MAX;
	for (int i = source + 1; i <= sink + 1; ++ i)
		if (PushRelabelDis[index][i] - PushRelabelFlow[index][i] > 0)
			minH = minH > PushRelabelH[i] ? PushRelabelH[i] : minH;
	if (minH == INT_MAX)
		return false;
	PushRelabelH[index] = minH + 1;

	for (int i = source + 1; i <= sink + 1; ++ i) 
	{
		if (PushRelabelEss[index] == 0)
			break;
		if (PushRelabelH[i] == minH && PushRelabelDis[index][i] > PushRelabelFlow[index][i])
			Push(index, i);
	}
	return true;
}

int PushReLabelMaxFlow()
{
	for (int i = 0; i < vertexNum; i++)
	{
		for (int j = 0; j < vertexNum; j++)
		{
			if (weightGraph[i][j] > 0)
				PushRelabelDis[i + 1][j + 1] = weightGraph[i][j];
		}
	}
	PushRelabelH[source + 1] = sink + 1;

	auto start = chrono::steady_clock::now();

	PreFlow();
	bool flag = true;
	while (true) 
	{
		if (flag == false)
			break;
		flag = false;
		for (int i = source + 1; i <= sink; ++ i)
			if (PushRelabelEss[i] > 0) 
				flag = ReLabel(i);
	}
	int maxFlow = PushRelabelEss[sink + 1];

	auto end = chrono::steady_clock::now();
	elapsedTime = chrono::duration<double>(end - start).count();
	return maxFlow;
}

struct EdgeInfo
{
	int startVertex;
	int endVertex;
	int weight;
	EdgeInfo(int _start, int _end, int _w) : startVertex(_start), endVertex(_end), weight(_w) {}
};

// Read input file to matrix
void ReadFile(string fileName)
{
	ifstream ifs(fileName);
	if (!ifs)
	{
		cerr << "Can not open file: " << fileName << "!" << endl;
		exit(0);
	}
	vector<EdgeInfo> allEdges;
	vector<int> allNodes;
	string line;
	while (getline(ifs, line))
	{
		stringstream tmpLine(line);
		int tmpNode = 0, endNode, weight;
		if (!(tmpLine >> tmpNode))
			continue;
		allNodes.push_back(tmpNode);
		char delimiter;
		while (tmpLine >> endNode >> delimiter >> weight)
			allEdges.push_back(EdgeInfo(tmpNode, endNode, weight));
	}
	vertexNum = (int)allNodes.size();
	if (vertexNum < 2)
	{
		cerr << "File: " << fileName << " content is wrong!" << endl;
		exit(0);
	}
	source = allNodes[0];
	sink = allNodes[vertexNum - 1];
	for (int i = 0; i < vertexNum; i++)
		weightGraph.push_back(vector<int>(vertexNum, 0));
	for (const auto& tmpEdge : allEdges)
		weightGraph[tmpEdge.startVertex][tmpEdge.endVertex] = tmpEdge.weight;
}

int main(int argc, char* argv[])
{

	string fileName;
	cout << "Please input file name: ";
	cin >> fileName;
	ReadFile(fileName);

	string::size_type idx = fileName.find(".dat");
	string outFileName;
	if (idx == string::npos)
		outFileName = fileName + ".out";
	else
		outFileName = fileName.substr(0, idx) + ".out";
	ofstream ofs(outFileName);

	double t1, t2, t3, t4;
	int result1 = FordFulkersonMaxFlow();
	t1 = elapsedTime;
	// cout << "Ford-Fulkerson result: " << result1 << endl;

	int result2 = DinicMaxFlow();
	// cout << "Dinic's Algorithm result: " << result2 << endl;
	t2 = elapsedTime;

	int result3 = EdmondKarpMaxFlow();
	// cout << "Edmond-Karp result: " << result3 << endl;
	t3 = elapsedTime;

	int result4 = PushReLabelMaxFlow();
	// cout << "General Push-Relabel: " << result4 << endl;
	t4 = elapsedTime;

	ofs << "Ford-Fulkerson: " << t1 << " sec" << endl
		<< "Dinic's Algorithm: " << t2 << " sec" << endl
		<< "Edmond-Karp: " << t3 << " sec" << endl
		<< "General Push-Relabel: " << t4 << " sec" << endl;
	ofs << endl << "Max Flow: " << result1 << endl;
	return 0;
}
