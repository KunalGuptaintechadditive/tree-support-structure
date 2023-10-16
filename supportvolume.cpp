#include <vtkAutoInit.h> 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);
#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCubeSource.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOctreePointLocator.h>
#include <vtkOutlineFilter.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include<vtkPolyLine.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOctreePointLocator.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkSliderWidget.h>
#include <vtkSphereSource.h>
#include <vtkSTLReader.h>
#include <vtkTextProperty.h>
#include <vtkWidgetEvent.h>
#include <vtkWidgetEventTranslator.h>
#include<vtkPointSource.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vtkTriangleFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkOctreePointLocator.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <iostream>
#include <set>
#include <vtkCursor3D.h>
#include <vtkMath.h>
#include <vtkMath.h>
#include<vtkSelectEnclosedPoints.h>
#include <vtkCubeSource.h>
#include <vtkSmartPointer.h>
#include <vtkExtractVOI.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <new.h>
#include<stack>
#include <chrono>
#include <unordered_set>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkSTLWriter.h>
#include <vtkDirectory.h>
#include<vtkAxesActor.h>
#include <vtkTransformPolyDataFilter.h>
#include<vtkTransform.h>



using namespace std;

//GLOBAL VARIABLES

long long int maxoct;
//array that stores all bottom octants voxel id..
std::vector<int> bottom;
//map that map voxelid to voxeltype and its center cordinate..
std::map<int, std::pair<char, std::vector<double>>> bmap;
double minbox[3], maxbox[3];



int totalchlild = 0;





void GetOctantVertices(double minBounds[3], double maxBounds[3], double vertices[8][3])
{
	// Extract the minimum and maximum coordinates of the bounding box
	double xmin = minBounds[0];
	double ymin = minBounds[1];
	double zmin = minBounds[2];
	double xmax = maxBounds[0];
	double ymax = maxBounds[1];
	double zmax = maxBounds[2];

	// Calculate the eight vertices of the octant
	vertices[0][0] = xmin;
	vertices[0][1] = ymin;
	vertices[0][2] = zmin;

	vertices[1][0] = xmax;
	vertices[1][1] = ymin;
	vertices[1][2] = zmin;

	vertices[2][0] = xmax;
	vertices[2][1] = ymax;
	vertices[2][2] = zmin;

	vertices[3][0] = xmin;
	vertices[3][1] = ymax;
	vertices[3][2] = zmin;

	vertices[4][0] = xmin;
	vertices[4][1] = ymin;
	vertices[4][2] = zmax;

	vertices[5][0] = xmax;
	vertices[5][1] = ymin;
	vertices[5][2] = zmax;

	vertices[6][0] = xmax;
	vertices[6][1] = ymax;
	vertices[6][2] = zmax;

	vertices[7][0] = xmin;
	vertices[7][1] = ymax;
	vertices[7][2] = zmax;
}

void filltest(double x, double y, double z, double test[3]) {
	test[0] = x;
	test[1] = y;
	test[2] = z;
}

//left right front back  down
//current voxel id and bmap

//created for left voxel return type
struct VoxelData {
	int voxelId;
	double centerX;
	double centerY;
	double centerZ;
};


VoxelData leftvoxel(int voxelid, double center[3], double voxelsize) {
	VoxelData ans;
	ans.voxelId = voxelid - 1;
	ans.centerX = center[0] - voxelsize;
	ans.centerY = center[1];
	ans.centerZ = center[2];

	//if left voxel does not exist
	if (ans.centerX < minbox[0]) {
		ans.voxelId = -100;
	}

	return ans;

}


//right
//right
VoxelData rightvoxel(int voxelid, double center[3], double voxelsize) {
	VoxelData ans;
	ans.voxelId = voxelid + 1;
	ans.centerX = center[0] + voxelsize;
	ans.centerY = center[1];
	ans.centerZ = center[2];

	//if right voxel does not exist
	if (ans.centerX > maxbox[0]) {
		ans.voxelId = -100;
	}

	return ans;

}

//bottom

VoxelData bottomvoxel(int voxelid, double center[3], double voxelsize) {
	VoxelData ans;
	int minus = std::cbrt(maxoct);
	ans.voxelId = voxelid - minus;
	ans.centerX = center[0];
	ans.centerY = center[1] - voxelsize;
	ans.centerZ = center[2];

	//if bottom voxel does not exist
	if (ans.centerY < minbox[1]) {
		ans.voxelId = -100;
	}

	return ans;
}

//top

VoxelData topvoxel(int voxelid, double center[3], double voxelsize) {
	VoxelData ans;
	int plus = std::cbrt(maxoct);
	ans.voxelId = voxelid + plus;
	ans.centerX = center[0];
	ans.centerY = center[1] + voxelsize;
	ans.centerZ = center[2];

	//if bottom voxel does not exist
	if (ans.centerY > maxbox[1]) {
		ans.voxelId = -100;
	}

	return ans;
}



//front

VoxelData frontvoxel(int voxelid, double center[3], double voxelsize) {
	VoxelData ans;
	int minus = std::cbrt(maxoct);
	minus = minus * minus;
	ans.voxelId = voxelid - minus;
	ans.centerX = center[0];
	ans.centerY = center[1];
	ans.centerZ = center[2] - voxelsize;

	//if bottom voxel does not exist
	if (ans.centerZ < minbox[2]) {
		ans.voxelId = -100;
	}

	return ans;
}

//back
VoxelData backvoxel(int voxelid, double center[3], double voxelsize) {
	VoxelData ans;
	int plus = std::cbrt(maxoct);
	plus = plus * plus;
	ans.voxelId = voxelid + plus;
	ans.centerX = center[0];
	ans.centerY = center[1];
	ans.centerZ = center[2] + voxelsize;

	//if bottom voxel does not exist
	if (ans.centerZ > maxbox[2]) {
		ans.voxelId = -100;
	}

	return ans;
}


//tells the position of voxel inside outside or partial
//voxelid = node
int position(double cord[3], int voxelid, double voxelsize, vtkPolyData* shape) {
	//ans =1  inside
	//ans=2  partial
	//ans=3   outside

	int ans = 0;
	double center[3];
	filltest(cord[0], cord[1], cord[2], center);
	//drawing for 60th voxel
	//insdie,outside,partial
	double minb[3], maxb[3];
	minb[0] = cord[0] - ((0.5) * voxelsize);
	minb[1] = cord[1] - ((0.5) * voxelsize);
	minb[2] = cord[2] - ((0.5) * voxelsize);
	//cout << "minbbbb " << minb[0] << " " << minb[1] << " " << minb[2] << endl;
	maxb[0] = cord[0] + ((0.5) * voxelsize);
	maxb[1] = cord[1] + ((0.5) * voxelsize);
	maxb[2] = cord[2] + ((0.5) * voxelsize);

	double vertices[8][3];
	GetOctantVertices(minb, maxb, vertices);

	//inside and outside vertex count....
	int inside = 0;
	int outside = 0;


	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	double test0[3];
	double test1[3];
	double test2[3];
	double test3[3];
	double test4[3];
	double test5[3];
	double test6[3];
	double test7[3];


	filltest(vertices[0][0], vertices[0][1], vertices[0][2], test0);
	points->InsertNextPoint(test0);
	filltest(vertices[1][0], vertices[1][1], vertices[1][2], test1);
	points->InsertNextPoint(test1);
	filltest(vertices[2][0], vertices[2][1], vertices[2][2], test2);
	points->InsertNextPoint(test2);
	filltest(vertices[3][0], vertices[3][1], vertices[3][2], test3);
	points->InsertNextPoint(test3);
	filltest(vertices[4][0], vertices[4][1], vertices[4][2], test4);
	points->InsertNextPoint(test4);
	filltest(vertices[5][0], vertices[5][1], vertices[5][2], test5);
	points->InsertNextPoint(test5);
	filltest(vertices[6][0], vertices[6][1], vertices[6][2], test6);
	points->InsertNextPoint(test6);
	filltest(vertices[7][0], vertices[7][1], vertices[7][2], test7);
	points->InsertNextPoint(test7);

	//cout << "test0" << test0[0] << " " << test0[1] << " " << test0[2] << endl;


	/*vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
	double test[3] = {2.5,5.5,5.5};
	points1->InsertNextPoint(test);
	vtkSmartPointer<vtkPolyData> pointsPolydata1 = vtkSmartPointer<vtkPolyData>::New();;
	pointsPolydata1->SetPoints(points1);
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints1 = vtkSmartPointer<vtkSelectEnclosedPoints>::New();;
	selectEnclosedPoints1->SetInputData(pointsPolydata1);
	selectEnclosedPoints1->SetSurfaceData(shape);
	selectEnclosedPoints1->Update();
	if (selectEnclosedPoints1->IsInside(0) == 1) cout << "iiiiiiinside";
	else cout << "outsideeeee";*/

	vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
	pointsPolydata->SetPoints(points);
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();;
	selectEnclosedPoints->SetInputData(pointsPolydata);
	selectEnclosedPoints->SetSurfaceData(shape);
	selectEnclosedPoints->Update();

	for (unsigned int xx = 0; xx < 8; xx++) {
		if (selectEnclosedPoints->IsInside(xx) == 1)
		{
			inside++;
		}
		else
		{
			outside++;
		}
	}

	std::vector<double> c{ cord[0],cord[1],cord[2] };

	if (inside == 8) {
		ans = 1;
		//cout << "insideoctant" <<voxelid << endl;
		bmap[voxelid] = std::make_pair('I', c);
		//cout << "did this one inside voxelid===== " << voxelid << endl;
	}
	else {
		if (outside == 8) {
			ans = 3;
			//cout << "outsideocatnat"<<voxelid << endl;
			bmap[voxelid] = std::make_pair('O', c);
			//cout << "did this one outside voxelid===== " << voxelid << endl;
		}
		else {
			ans = 2;
			//cout << "partialoctant                        "<<voxelid<< endl;
			bmap[voxelid] = std::make_pair('P', c);
			//cout << "did this one partial voxelid===== " << voxelid << endl;

			// 1 posistion call of bottomvoxel extra here so commenting it..
			/*int xxx = bottomvoxel(voxelid, voxelsize, bmap, maxoct);
			if (xxx != -1) {
				char c = bmap[xxx].first;
				if (c == 'O') bottom.push_back(voxelid);
			}*/
		}
	}

	return ans;
}

vector<VoxelData> adjacent(int node, double center[3], double voxelsize) {
	vector<VoxelData> adj;

	VoxelData left = leftvoxel(node, center, voxelsize);
	if (left.voxelId != -100) {
		adj.push_back(left);
	}

	//right
	VoxelData right = rightvoxel(node, center, voxelsize);
	if (right.voxelId != -100) {
		adj.push_back(right);
	}

	VoxelData top = topvoxel(node, center, voxelsize);
	if (top.voxelId != -100) {
		adj.push_back(top);
	}

	VoxelData bottom = bottomvoxel(node, center, voxelsize);
	if (bottom.voxelId != -100) {
		adj.push_back(bottom);
	}

	VoxelData front = frontvoxel(node, center, voxelsize);
	if (front.voxelId != -100) {
		adj.push_back(front);
	}

	VoxelData back = backvoxel(node, center, voxelsize);
	if (back.voxelId != -100) {
		adj.push_back(back);
	}

	return adj;
}



double lowest_y = INT_MAX;
//dfs(start, adj, vis, insidels, partialls);
//dfs function returns list of nodes that are internal or partial... 
void dfs(int node, vector<int>& vis, vector<int>& insidels, vector<int>& partialls, double voxelsize, vtkPolyData* shape) {

	stack<int> s;
	s.push(node);

	while (!s.empty()) {
		int var = s.top();
		s.pop();

		if (!vis[var])
		{

			vis[var] = 1;
			std::vector<double> center = bmap[var].second;

			double cen[3] = { center[0],center[1],center[2] };
			int ans = position(cen, var, voxelsize, shape);
			//	cout << "posistion is =======================================" << ans << endl;
			if (ans == 1) { insidels.push_back(var); 
				if (lowest_y > center[1]) lowest_y = center[1];
			}
			if (ans == 2) {		
				partialls.push_back(var);
				if (lowest_y > center[1]) lowest_y = center[1];
			}
			if (ans == 3) continue;
			else {
				vector<VoxelData> adj = adjacent(var, cen, voxelsize);

				for (const auto& it : adj) {
					if (!vis[it.voxelId]) {

						s.push(it.voxelId);

					}
				}
			}
		}
	}
}


unordered_set<int> bottompartialvoxel;
//function to get all the bottom partial voxels from partial voxels
vector<int> getbottompartial(vector<int> partialls) {
	vector<int> bottom;
	int minus = std::cbrt(maxoct);

	for (auto it : partialls) {
		//no need to check for existence as bottom voxel of every partial voxel is added into adjacent list so with posistion function its bmap is added if it exist and existence of bottom voxel check in adjacent list
		if (bmap[it - minus].first == 'I') {
			bottom.push_back(it);
			bottompartialvoxel.insert(it);
		}
	}

	return bottom;

}
//bfs
bool find_shortest_path(int vid, int bottom_y, double voxelsize, map<int, int>& pred, int& dest) {

	std::deque<int> queue;
	unordered_set<int> visited;

	queue.push_back(vid);
	visited.insert(vid);

	while (!queue.empty()) {

		int voxel_id = queue.front();

		queue.pop_front();

		//cout << "queue.size = " << queue.size() << endl;
		vector<double> c = bmap[voxel_id].second;
		double cen[3] = { c[0],c[1],c[2] };
		// Check adjacent voxels and add them to the queue if not visited and not obstructed by mesh..
		vector<VoxelData> adj = adjacent(voxel_id, cen, voxelsize);

		for (int i = 0; i < adj.size(); i++) {
			//cout << "adjacentes are" << adj[i].voxelId << endl;
			vector<double> ccord = bmap[adj[i].voxelId].second;
			double cc[3] = { ccord[0],ccord[1],ccord[2] };

			char ans = bmap[adj[i].voxelId].first;
			if (visited.find(adj[i].voxelId) != visited.end() || ans == 'I' || ans == 'P') continue;
			else
			{
				visited.insert(adj[i].voxelId);
				//pred[adj[i].voxelId] = voxel_id;
				if (pred.find(adj[i].voxelId) == pred.end()) pred.insert({ adj[i].voxelId ,voxel_id });
				queue.push_back(adj[i].voxelId);

				//stopping bfs when we reach bottom floor....

				//vector<double> center = bmap[adj[i].voxelId].second;
				int ycord_now = adj[i].centerY;

				if (ycord_now == bottom_y) {
					dest = adj[i].voxelId;
					//cout << "destination is             " << dest<<endl;
				//	cout << "ycord_now =   " << ycord_now << endl;
					visited.clear();
					//cout << "pred ka size is" << pred.size();
					return true;
				}
			}
		}
	}

	return false;  // No path found
}



//v to pata hi nahi hai
//vector<int> BFS(int node, int bottom_y, double voxelsize)
//{
//	list<int> queue;
//	vector<bool> visited, pred, dist;
//
//
//	for (int i = 0; i < v; i++) {
//		//visited[i] = false;
//		dist[i] = INT_MAX;
//		pred[i] = -1;
//	}
//
//	visited[src] = true;
//	dist[src] = 0;
//	queue.push_back(src);
//
//	while (!queue.empty()) {
//		int u = queue.front();
//		queue.pop_front();
//		vector<double> cen = bmap[u].second;
//		double ucen[3] = {cen[0],cen[2],cen[3]};
//		vector<VoxelData> adj = adjacent(u, ucen, voxelsize);
//		
//		for (int i = 0; i < adj.size(); i++) {
//			if (visited[adj[i].voxelId] == false) {
//				visited[adj[i].voxelId] = true;
//				dist[adj[i].voxelId] = dist[u] + 1;
//				pred[adj[i].voxelId] = u;
//				queue.push_back(adj[i].voxelId);
//
//				//stopping bfs when we reach bottom floor....
//
//				vector<double> center = bmap[adj[i].voxelId].second;
//				int ycord_now = center[1];
//				if (ycord_now == bottom_y)	return true;
//			}
//		}
//	}
//
//	return false;
//}

//stores all the centers of 
vector<double*> grppoints;

unordered_set<int> bottom_visited;

// 1 layer
// this would be faster than 2 layer code 
//return joining center and vector of points that we want to join.....
pair<double*, vector<int>> centerpoint(int node, double center[3], double voxelsize) {
	//center would be the grp bottom not 
	double* joincenter = new double[3];
	vector<int> joinee;
	int no_of_points = 0;
	double x = 0, y = 0, z = 0;

	/*bool minus = 0;*/
	// minus = 0   minus the voxelsize from center 
	// minus = 1   minus the voxelsize/2 from center


	//itself ko add
	vector<double> c = bmap[node].second;
	x += c[0];
	y += c[1];
	z += c[2];
	no_of_points++;
	joinee.push_back(node);


	//right
	VoxelData right = rightvoxel(node, center, voxelsize);
	if (right.voxelId != -100 && bottompartialvoxel.find(right.voxelId) != bottompartialvoxel.end() && bottom_visited.find(right.voxelId) == bottom_visited.end()) {
		x += right.centerX;
		y += right.centerY;
		z += right.centerZ;
		no_of_points++;
		joinee.push_back(right.voxelId);
	}

	//right ka top
	double cr[3] = { right.centerX,right.centerY,right.centerZ };
	VoxelData righttop = topvoxel(right.voxelId, cr, voxelsize);
	if (righttop.voxelId != -100 && bottompartialvoxel.find(righttop.voxelId) != bottompartialvoxel.end() && bottom_visited.find(righttop.voxelId) == bottom_visited.end()) {
		x += righttop.centerX;
		y += righttop.centerY;
		z += righttop.centerZ;
		no_of_points++;
		joinee.push_back(righttop.voxelId);
	}

	//right ka top ka back
	double rtbcenter[3] = { righttop.centerX,righttop.centerY,righttop.centerZ };
	VoxelData rtb = topvoxel(righttop.voxelId, rtbcenter, voxelsize);
	if (rtb.voxelId != -100 && bottompartialvoxel.find(rtb.voxelId) != bottompartialvoxel.end() && bottom_visited.find(rtb.voxelId) == bottom_visited.end()) {
		x += rtb.centerX;
		y += rtb.centerY;
		z += rtb.centerZ;
		no_of_points++;
		joinee.push_back(rtb.voxelId);
	}


	VoxelData top = topvoxel(node, center, voxelsize);
	if (top.voxelId != -100 && bottompartialvoxel.find(top.voxelId) != bottompartialvoxel.end() && bottom_visited.find(top.voxelId) == bottom_visited.end()) {
		x += top.centerX;
		y += top.centerY;
		z += top.centerZ;
		no_of_points++;
		joinee.push_back(top.voxelId);
	}


	VoxelData back = backvoxel(node, center, voxelsize);
	if (back.voxelId != -100 && bottompartialvoxel.find(back.voxelId) != bottompartialvoxel.end() && bottom_visited.find(back.voxelId) == bottom_visited.end()) {
		x += back.centerX;
		y += back.centerY;
		z += back.centerZ;
		no_of_points++;
		joinee.push_back(back.voxelId);
	}

	//backkatop 
	double btcenter[3] = { back.centerX,back.centerY,back.centerZ };

	VoxelData bt = backvoxel(back.voxelId, btcenter, voxelsize);
	if (bt.voxelId != -100 && bottompartialvoxel.find(bt.voxelId) != bottompartialvoxel.end() && bottom_visited.find(bt.voxelId) == bottom_visited.end()) {
		x += bt.centerX;
		y += bt.centerY;
		z += bt.centerZ;
		no_of_points++;
		joinee.push_back(bt.voxelId);
	}

	//backkaright
	double brcenter[3] = { back.centerX,back.centerY,back.centerZ };

	VoxelData br = rightvoxel(back.voxelId, brcenter, voxelsize);
	if (br.voxelId != -100 && bottompartialvoxel.find(br.voxelId) != bottompartialvoxel.end() && bottom_visited.find(br.voxelId) == bottom_visited.end()) {
		x += br.centerX;
		y += br.centerY;
		z += br.centerZ;
		no_of_points++;
		joinee.push_back(br.voxelId);
	}

	x = x / no_of_points;
	y = y / no_of_points;
	z = z / no_of_points;

	x = x;
	y = y - voxelsize;
	z = z;

	joincenter[0] = x;
	joincenter[1] = y;
	joincenter[2] = z;

	double pushback[3] = { x,y,z };
	grppoints.push_back(joincenter);

	pair<double*, vector<int>> m = { joincenter,joinee };
	return m;

}

struct DoubleArray {
	double data[3];

	bool operator==(const DoubleArray& other) const {
		return data[0] == other.data[0] && data[1] == other.data[1] && data[2] == other.data[2];
	}
};

// Define a custom hash function for DoubleArray
struct DoubleArrayHash {
	std::size_t operator()(const DoubleArray& arr) const {
		std::size_t seed = 0;
		for (int i = 0; i < 3; ++i) {
			// Combine the hash of each double in the array
			seed ^= std::hash<double>{}(arr.data[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

// Define a custom equality function for DoubleArray
struct DoubleArrayEqual {
	bool operator()(const DoubleArray& arr1, const DoubleArray& arr2) const {
		return arr1 == arr2;
	}
};

//std::unordered_set<double*, DoubleArrayPointerHash, DoubleArrayPointerEqual> grpvisited2;

vector<double*> grpnextpoints;
std::unordered_set<DoubleArray, DoubleArrayHash, DoubleArrayEqual> grpvisited;

//grp size 4 join center
//problem in find here
pair<double*, vector<double*>> grpcenterpoint(double node[3], double voxelsize, int depth) {
	//center to where we need to join
	double* joincenter = new double[3];
	//all the points which we need to join
	vector<double*> joinee;

	double x, y, z;
	x = node[0];
	y = node[1];
	z = node[2];
	int no_of_points = 1;
	//finding the range of node according to 4*4*4 box
	//min point and max point of range in x,y,z

	//min point 
	double xmin, ymin, zmin;
	xmin = node[0] - voxelsize / 2;
	ymin = node[1] - voxelsize / 2;
	zmin = node[2] - voxelsize / 2;

	//max point
	double xmax, ymax, zmax;
	//group box range
	int plus = 4 * (pow(2, depth - 2)) * voxelsize;
	xmax = xmin + plus;
	ymax = ymin + plus;
	zmax = zmin + plus;

	double* nodde = new double[3];
	nodde = node;

	joinee.push_back(nodde);

	//grpvisited.insert(nodde);


	//cout << " inseted 111111 gprvisisted size is" << grpvisited.size() << endl;
	for (int i = 0; i < grppoints.size(); i++) {
		cout << "i here " << i << endl;
		double* cord = new double[3];
		cord = grppoints[i];

		//	double cordd[3] = { cord[0],cord[1],cord[2] };
		DoubleArray arr = { {cord[0], cord[1], cord[2]} };
		// check here
		if (grpvisited.find(arr) == grpvisited.end()) {
			cout << "runned " << i << endl;
			if (xmin <= cord[0] && xmax >= cord[0]) {
				if (ymin <= cord[1] && ymax >= cord[1]) {
					if (zmin <= cord[2] && zmax >= cord[2]) {
						x += cord[0];
						y += cord[1];
						z += cord[2];
						joinee.push_back(cord);
						//grpvisited.insert(cord);
						//cout << " inseted 222222222 gprvisisted size is" << grpvisited.size() << endl;
						no_of_points++;

						/*if (no_of_points == 4) {
							break;
						}*/
					}
				}
			}
		}
		else cout << "not runned " << i << endl;
	}


	x = x / no_of_points;
	y = y / no_of_points;
	z = z / no_of_points;

	int minus = pow(2, depth - 1) * voxelsize;

	y = y - (minus);

	joincenter[0] = x;
	joincenter[1] = y;
	joincenter[2] = z;

	double pushback[3] = { joincenter[0],joincenter[1],joincenter[2] };
	grpnextpoints.push_back(joincenter);

	pair<double*, vector<double*>> m = { joincenter,joinee };
	return m;
}









int main()
{

	vtkNew<vtkNamedColors> colors;
	auto start_time = std::chrono::high_resolution_clock::now();

	std::string stl_file = "C:\\Users\\kunal.g\\Desktop\\STL_Files\\supportvolume.stl";

	// "D:\\DOWNLOADS\\cat.stl"  "D:\\DOWNLOADS\\cylinder.stl" "C:\\Users\\dell\\Documents\\Sample.stl"
	//  C:\Users\dell\Documents\Horn.stl
	// "D:\\DOWNLOADS\\Modular_Bar_and_Snap_Covid_Mask_Ear_Saver_4569503\\files\\14pbar_03.stl"
	// Value                    1         2       4.5       4         8      16
	// HeatSink1               1hr      5min              52 sec    nn     5sec
	// Engine_Mount_Bracket               nn              10min            30 sec
	// cat                      doint         10min           2min
	// Aerospace_AlSi10Mg_Waveguide block 01      4hrs
	// 
	// starter casing -3 arm.stl   value 1 time



	//tree time             1,2.5          2,2             4,1.5
	//heatsink1              (42 min)     (5min)         1min
	// 
	//						 4,4          8,2       16,1.8
	//Engine_Mount_Bracket    8min					  36sec
	// 
	//						 2,4				4,3			8,2
	//cat					9min				1 min		 30sec
	// 
	// 
	// Aerospace_AlSi10Mg_Waveguide block 01
	// 
	//
	//starter casing -3 arm.stl
	// Read the STL file

	vtkNew<vtkSTLReader> reader;
	reader->SetFileName(stl_file.c_str());
	reader->Update();

	vtkPolyData* shape = reader->GetOutput();



	double c[3];
	shape->GetCenter(c);
	cout << "center   " << c[0] << " " << c[1] << " " << c[2] << endl;





	// Step 2: Create a bounding box using vtkBoundingBox
	vtkBoundingBox boundingBox;
	boundingBox.SetBounds(shape->GetBounds());
	//maximum length in bounding box
	double length[3];
	boundingBox.GetLengths(length);
	cout << length[0] << " " << length[1] << " " << length[2] << endl;
	// Find the maximum side length among x, y, and z dimensions
	double maxlength = std::max(length[0], std::max(length[1], length[2]));
	//maxlength = maxlength * 1.5;
	cout << "maxlength" << maxlength << endl;
	boundingBox.Inflate((maxlength - length[0]) / 2, (maxlength - length[1]) / 2, (maxlength - length[2]) / 2);

	double ll[3];
	boundingBox.GetLengths(ll);
	cout << ll[0] << " " << ll[1] << " " << ll[2] << endl;
	double center[3];
	boundingBox.GetCenter(center);




	//transform before messing other terms
	//
	double tbounds[6];
	boundingBox.GetBounds(tbounds);
	//double translation[3] = { -tbounds[0] , -tbounds[2]  , -tbounds[4] };
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	//transform->Translate(translation);
	transform->Translate(center); // Translate to the center
	transform->RotateWXYZ(-90.0, 1.0, 0.0, 0.0); // Rotate 90 degrees around the Z-axis
	transform->Translate(-center[0], -center[1], -center[2]); // Translate back to the original position



	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetInputData(shape);
	transformFilter->SetTransform(transform);
	transformFilter->Update();

	shape = transformFilter->GetOutput();


	vtkSmartPointer<vtkPolyData> polyData = transformFilter->GetOutput();


	reader->Update();

	// Convert polygons to triangles
	vtkNew<vtkTriangleFilter> triangleFilter;
	triangleFilter->SetInputData(transformFilter->GetOutput());
	triangleFilter->Update();

	// Create a mapper for the triangles
	vtkNew<vtkPolyDataMapper> triangleMapper;
	triangleMapper->SetInputData(triangleFilter->GetOutput());

	// Create an actor for the triangles
	vtkNew<vtkActor> triangleActor;
	triangleActor->SetMapper(triangleMapper);
	triangleActor->GetProperty()->SetInterpolationToFlat();
	triangleActor->GetProperty()->SetColor(colors->GetColor4d("Red").GetData());
	triangleActor->GetProperty()->SetOpacity(0.2);
	// previous opacity is 0.7






	vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();

	//double bound[6];
	//boundingBox.GetBounds(bound);
	//cubeSource->SetBounds(bound);


	//enter the size of voxel
	cout << "enter the size of voxel" << endl;
	cout << "FOR DEFAULT VALUE(i:e maximum side length of bounding box/10) ENTER 0" << endl;
	double voxelsize;
	cin >> voxelsize;
	//divide by 15 bcoz done maxlenght*1.5 above
	if (voxelsize == 0) voxelsize = maxlength / 15;
	int level;
	cout << "Enter the level of tree support structure (FOR DEFAULT VALUE(i:e 2) ENTER 0 )" << endl;
	cin >> level;
	if (level == 0) level = 2;
	//bounding box size increase function
	// non factor and float value included in voxelsize
	float increase = maxlength / voxelsize;
	increase = ceil(increase);
	maxlength = increase * voxelsize;

	double newtranslatedcenter[3] = { maxlength / 2,maxlength / 2,maxlength / 2 };
	cubeSource->SetCenter(center);

	//Previous one where we make the bounding box as cubic with all sides equal to the largest side of orignal bounding box
	// But if the stl file is also cubic box then all voxels will be inside then it will show problem so making the default size acc to machine 
	//so changing the maxlength to default value bcoz maxlength is used to calculate other terms also

	cubeSource->SetXLength(maxlength);
	cubeSource->SetYLength(maxlength);
	cubeSource->SetZLength(maxlength);

	/*double defaultvalue = 350;
	cubeSource->SetXLength(defaultvalue);
	cubeSource->SetYLength(defaultvalue);
	cubeSource->SetZLength(defaultvalue);*/
	cubeSource->Update();

	double center1[3];
	cubeSource->GetCenter(center1);
	cout << "center of cubesource " << center1[0] << " " << center1[1] << " " << center1[2] << endl;

	double cubebound[6] = { 0 };
	cubeSource->GetBounds(cubebound);

	//// Access the coordinates of the minimum and maximum points
	double xmin = cubebound[0];
	double ymin = cubebound[2];
	double zmin = cubebound[4];

	double xmax = cubebound[1];
	double ymax = cubebound[3];
	double zmax = cubebound[5];



	//if out of bound of machine plate of 350*350*350 whose min point is origin
	/*if (xmin >= 0 && xmin <= 350 && ymin >= 0 && ymin <= 350 && zmin >= 0 && zmin <= 350) {
		cout << "ERROR WARNING STL OUT OF BOUND OF MACHINE PLATE" << endl;
		return 0;
	}

	if (xmax >= 0 && xmax <= 350 && ymax >= 0 && ymax <= 350 && zmax >= 0 && zmax <= 350) {
		cout << "ERROR WARNING STL OUT OF BOUND OF MACHINE PLATE" << endl;
		return 0;
	}*/

	//// Use the coordinates (xmin, ymin, zmin) and (xmax, ymax, zmax) as needed
	std::cout << "Minimum point: (" << xmin << ", " << ymin << ", " << zmin << ")" << std::endl;
	std::cout << "Maximum point: (" << xmax << ", " << ymax << ", " << zmax << ")" << std::endl;


	////machinebox min pt is origin and side lengthis 350
	//vtkSmartPointer<vtkCubeSource> machinebox = vtkSmartPointer<vtkCubeSource>::New();
	//double origin[3] = { 175, 175, 175 };
	//machinebox->SetCenter(origin);
	//machinebox->SetXLength(350);
	//machinebox->SetYLength(350);
	//machinebox->SetZLength(350);
	//machinebox->Update();

	//vtkSmartPointer<vtkPolyDataMapper> mapper =
	//	vtkSmartPointer<vtkPolyDataMapper>::New();
	//mapper->SetInputConnection(machinebox->GetOutputPort());

	//// 3. Create an actor
	//vtkSmartPointer<vtkActor> machineactor =
	//	vtkSmartPointer<vtkActor>::New();
	//machineactor->SetMapper(mapper);
	//machineactor->GetProperty()->SetOpacity(0.5);





	//// Create an axes actor
	//vtkSmartPointer<vtkAxesActor> axes =
	//	vtkSmartPointer<vtkAxesActor>::New();



	//global variable
	minbox[0] = xmin;
	minbox[1] = ymin;
	minbox[2] = zmin;
	maxbox[0] = xmax;
	maxbox[1] = ymax;
	maxbox[2] = zmax;

	int voxelid = 1;
	double cord[3] = { xmin + (voxelsize / 2) ,ymin + (voxelsize / 2),zmin + (voxelsize / 2) };
	double duplicate[3] = { xmin + (voxelsize / 2) ,ymin + (voxelsize / 2),zmin + (voxelsize / 2) };

	int dimension = maxlength / voxelsize;
	int insideoct = 0;
	int outsideoct = 0;
	int partialoct = 0;
	cout << "center" << cord[0] << " " << cord[1] << " " << cord[2] << endl;

	int no_of_points = 0;
	//checking for octants which are in bottom
	int bottomcheck = 0;
	//when to stop building the line when bottom reaached
	int stop = 0;
	//bottom floor y value (lowest y axis)..
	//int bottom_y = duplicate[1] + (voxelsize * dimension) - voxelsize;
	double bottom_y = duplicate[1];
	cout << "bottom y is" << bottom_y << endl;

	int bool_bottomcheck = 0;

	double center_y = ymin + (maxlength / 2);
	cout << "center y is " << center_y << endl;

	//max number of octants formed
	maxoct = (long long)dimension * (long long)dimension * (long long)dimension;
	cout << "maxoct  " << maxoct << endl;



	//setting default value of bmap
	//std::pair<char, std::vector<double>> defaultValue = { 'O', {0, 0, 0} };

	//// Inserting the default value for each key
	//for (int i = 0; i <= maxoct+1; ++i) {
	//	bmap[i] = defaultValue;
	//}

	double cordd[3] = { xmin + (voxelsize / 2) ,ymin + (voxelsize / 2),zmin + (voxelsize / 2) };


	for (int z = 0; z < dimension; z++) {
		cordd[2] = duplicate[2] + z * voxelsize;
		for (int y = 0; y < dimension; y++) {
			cordd[1] = duplicate[1] + y * voxelsize;
			for (int x = 0; x < dimension; x++) {
				cordd[0] = duplicate[0] + x * voxelsize;
				//cout << "voxelid " << voxelid << " cordinates  (" << cordd[0] << " " << cordd[1] << " " << cordd[2] << ") " << endl;
				// bmap default value set
				bmap[voxelid] = { 'O',{cordd[0],cordd[1],cordd[2]} };
				voxelid++;
			}
		}
	}



	double testcase[3];
	vtkIdType i = 0;
	//a single point which is inside the polydata... it get filled in testcase..
	polyData->GetPoint(i, testcase);
	cout << "testcase review dubug         mmmmmmmmmmmmmmmmmmmmmmmmmm " << endl;

	cout << testcase[0] << " " << testcase[1] << " " << testcase[2] << endl;
	//power problem
	//int dist = (sqrt((testcase[0] - xmin)* (testcase[0] - xmin) + (testcase[1] - ymin)* (testcase[1] - ymin) + (testcase[2] - zmin)* (testcase[2] - zmin)))/voxelsize;
	//check 2.......
	//finding the center process start
	double ix = (testcase[0] - duplicate[0]) / voxelsize;
	double iy = (testcase[1] - duplicate[1]) / voxelsize;
	double iz = (testcase[2] - duplicate[2]) / voxelsize;

	/*calculate the remainder
		rem = cord(x)-duplic(x) % v.size
		rem > vsize /2              ---->>>>    ceil(ix)
		rem < vsize/2				---->>>>    floor(ix)
		in case if rem == v.size/2 then point is on boundary take a lower
		or higher point in the boundary....
	*/

	cout << "testcase[0]           " << testcase[0] << endl;
	cout << "duplicate[0]          " << duplicate[0] << endl;
	cout << " voxelsize            " << voxelsize << endl;
	double remx = remainder((testcase[0] - duplicate[0]), voxelsize);
	//if (remx < 0) remx = voxelsize + remx;
	cout << "remainder of xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << remx << endl;
	double remy = remainder((testcase[1] - duplicate[1]), voxelsize);
	cout << "remy is " << remy << endl;
	//if (remy < 0) remy = voxelsize + remy;
	cout << "remainder of xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << remy << endl;
	double remz = remainder((testcase[2] - duplicate[2]), voxelsize);
	//if (remz < 0) remz = voxelsize + remz;
	cout << "remainder of xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << remz << endl;
	/*if (remx > (voxelsize / 2)) ix = ceil(ix);
	else ix = floor(ix);
	if (remy > (voxelsize / 2)) iy = ceil(iy);
	else iy = floor(iy);
	if (remz > (voxelsize / 2)) iz = ceil(iz);
	else iz = floor(iz);*/
	if (remx > 0) ix = floor(ix);
	else ix = ceil(ix);
	if (remy > 0) iy = floor(iy);
	else iy = ceil(iy);
	if (remz > 0) iz = floor(iz);
	else iz = ceil(iz);
	double cen[3] = { duplicate[0] + ix * voxelsize,duplicate[1] + iy * voxelsize, duplicate[2] + iz * voxelsize };
	cout << "z cordinate is                  " << duplicate[2] << endl;
	cout << "ix   = " << ix << endl;
	cout << "iy     = " << iy << endl;
	cout << "iz        =" << iz << endl;
	cout << "the center is                   vvvvvvvvvvvvv" << endl;
	cout << cen[0] << "     " << cen[1] << "      " << cen[2] << endl;
	cout << "maxoctttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt" << maxoct << endl;
	//finding the start value minus the center..
	cout << "cen[0]-cord[0]                     =================" << cen[0] << endl;
	double start = (abs(cen[0] - cord[0]) * (1 / voxelsize)) + (abs(cen[1] - cord[1]) * (dimension / voxelsize)) + (abs(cen[2] - cord[2]) * (dimension / voxelsize) * (dimension)) + 1;
	cout << "startttttttttttttttttttttttttttttttttttttttttt" << start << endl;



	// dfs
	//n= h* l * b / v.size;
	cout << length[0] << " " << length[1] << " " << length[2] << endl;
	//visited[n]
	long long int len = ceil((length[0] * length[1] * length[2]) / (voxelsize * voxelsize * voxelsize));
	cout << "lennnnnnnnnnnnnnnnnnnnnnnn" << len << endl;
	len = maxoct + 1;
	//length of vis is problem
	vector<int> vis(len, 0);
	//start=center voxelid
	//adjacent nodes will be inserted through left right bottom top front back algo
	//will find the adjecent in the dfs loop only...
	vector<int> adj;
	//list of all inside voxels
	vector<int> insidels;
	//list of all partial voxels
	vector<int> partialls;

	dfs(start, vis, insidels, partialls, voxelsize, shape);


	// int maxoct = insideoct + outsideoct + partialoct;

	int minus = std::cbrt(maxoct);
	//cout << "minus" << minus << endl;

	//array that stores all bottom octants voxel id..
	bottom = getbottompartial(partialls);
	// making support for all bottom voxels....

	// A renderer and render window
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);


	// An interactor
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	//now only fix this
	//support generation for all bottom voxels + deviation of shortest path
	cout << "bottom size" << bottom.size() << endl;

	//number of octants in one edge
	int oct_in_edge = dimension;
	cout << "maxoct here are" << maxoct << endl;
	cout << "oct_in_edge here are " << oct_in_edge << endl;

	// 0th layer
	for (int i = 0; i < bottom.size(); i++) {
		//4 top vertex to center
		vtkSmartPointer<vtkPoints> centerline = vtkSmartPointer<vtkPoints>::New();

		int vid = bottom[i];

		std::vector<double> cord = bmap[vid].second;
		double linep[3] = { cord[0],cord[1],cord[2] };

		//filltest(cord[0], cord[1], cord[2], linep);
		//traverse to add points in linepoints
		// j is voxelid
		//going down from single bottom partial voxel center

		int m = voxelsize / 2;
		double vertex1[3] = { cord[0] - m,cord[1] + m,cord[2] - m };
		double vertex2[3] = { cord[0] + m,cord[1] + m,cord[2] - m };
		double vertex3[3] = { cord[0] + m,cord[1] + m,cord[2] + m };
		double vertex4[3] = { cord[0] - m,cord[1] + m,cord[2] + m };

		centerline->InsertNextPoint(linep);    //point 0 
		centerline->InsertNextPoint(vertex1);  // point 1
		centerline->InsertNextPoint(vertex2);  // point 2
		centerline->InsertNextPoint(vertex3);  // point 3
		centerline->InsertNextPoint(vertex4);  // point 4


		//center line top 4 vertex to center
		// Create a cell array to store the lines in and add the lines to it...				
		vtkNew<vtkCellArray> centercells;
		vtkIdType lineIndices01[2] = { 0, 1 }; // Connect Point 4 to Point 0
		centercells->InsertNextCell(2, lineIndices01);

		vtkIdType lineIndices02[2] = { 0, 2 };
		centercells->InsertNextCell(2, lineIndices02);

		vtkIdType lineIndices03[2] = { 0 ,3 };
		centercells->InsertNextCell(2, lineIndices03);

		vtkIdType lineIndices04[2] = { 0 , 4 };
		centercells->InsertNextCell(2, lineIndices04);

		// Create a polydata to store everything in
		vtkNew<vtkPolyData> centerlinepolyData;

		// Add the points to the dataset
		centerlinepolyData->SetPoints(centerline);

		// Add the lines to the dataset
		centerlinepolyData->SetLines(centercells);

		// Setup actor and mapper
		vtkNew<vtkPolyDataMapper> centerlinemapper;
		centerlinemapper->SetInputData(centerlinepolyData);

		vtkNew<vtkActor> centerlineactor;
		centerlineactor->SetMapper(centerlinemapper);
		centerlineactor->GetProperty()->SetColor(colors->GetColor3d("black").GetData());
		int centerline_width = 10;  //# Change this value to the desired line width
		centerlineactor->GetProperty()->SetLineWidth(centerline_width);
		renderer->AddActor(centerlineactor);

	}

	sort(bottom.begin(), bottom.end());

	//@  1st layer
	for (int i = 0; i < bottom.size(); i++) {
		//straight bottome line
		if (bottom_visited.find(bottom[i]) == bottom_visited.end()) {

			vtkSmartPointer<vtkPoints> linepoints = vtkSmartPointer<vtkPoints>::New();

			int vid = bottom[i];
			//cout << "vid is" << vid << endl;
			int no_of_points = 0;

			std::vector<double> cord = bmap[vid].second;
			double linep[3] = { cord[0],cord[1],cord[2] };


			auto result = centerpoint(vid, linep, voxelsize);
			double* joincenter;
			joincenter = result.first;
			vector<int> joinee = result.second;



			// find grp joining point
			/*

			front     back
			3 4       7 8
			1 2       5 6
			*/

			//if(pos_in_grp == 1)

			// if grp is exception(last left out voxels)
			// 
			//path for finding the linepoint  ###@

			//map<int, int> pred;
			//int dest = 0;

			//find_shortest_path(vid, bottom_y, voxelsize, pred, dest);

			//vector<int> path;
			//int crawl = dest;

			//linepoints->InsertNextPoint(linep);
			//no_of_points++;
			//path.push_back(dest);
			//no_of_points++;

			//while (pred.find(crawl) != pred.end()) {
			//	path.push_back(pred[crawl]);
			//	crawl = pred[crawl];
			//	no_of_points++;
			//	cout << "no of points" << no_of_points << endl;
			//}

			//pred.clear();
			//reverse(path.begin(), path.end());
			// Create a cell array to store the lines in and add the lines to it
			vtkNew<vtkCellArray> cells;

			linepoints->InsertNextPoint(joincenter); // main point 0

			for (int i = 0; i < joinee.size(); i++) {

				bottom_visited.insert(joinee[i]);
				vector<double> d = bmap[joinee[i]].second;

				double center[3] = { d[0],d[1],d[2] };
				linepoints->InsertNextPoint(center);
				int z = i + 1;
				vtkIdType lineIndices[2] = { 0, z }; // Connect Point 4 to Point 0
				cells->InsertNextCell(2, lineIndices);
			}

			/*	vtkNew<vtkPolyLine> polyLine;
				polyLine->GetPointIds()->SetNumberOfIds(no_of_points);*/

				// Create a polydata to store everything in
			vtkNew<vtkPolyData> linepolyData;

			// Add the points to the dataset
			linepolyData->SetPoints(linepoints);

			// Add the lines to the dataset
			linepolyData->SetLines(cells);

			// Setup actor and mapper
			vtkNew<vtkPolyDataMapper> linemapper;
			linemapper->SetInputData(linepolyData);

			vtkNew<vtkActor> lineactor;
			lineactor->SetMapper(linemapper);
			lineactor->GetProperty()->SetColor(colors->GetColor3d("black").GetData());
			int line_width = 5;  //# Change this value to the desired line width
			lineactor->GetProperty()->SetLineWidth(line_width);
			renderer->AddActor(lineactor);
		}
	}


	cout << "grppoints are   ====================" << endl;
	for (int i = 0; i < grppoints.size(); i++) {
		cout << grppoints[i][0] << "  " << grppoints[i][1] << "  " << grppoints[i][2] << endl;
	}
	cout << "finished" << endl;
	//@@
	int depth = 1;
	//grppoints from after level 1
	//from 2 layer
	while (++depth <= level) {

		for (int i = 0; i < grppoints.size(); i++) {

			cout << "grppoints.size() is    " << grppoints.size() << endl;
			cout << "grpvisited.size() is     " << grpvisited.size() << endl;
			double fvertex[3] = { grppoints[i][0] ,grppoints[i][1] ,grppoints[i][2] };
			DoubleArray arrfind = { {grppoints[i][0], grppoints[i][1], grppoints[i][2]} };
			cout << "the grppoint is  " << grppoints[i][0] << " " << grppoints[i][1] << " " << grppoints[i][2] << endl;
			if (grpvisited.find(arrfind) == grpvisited.end()) {
				cout << "test1 =====================================" << i << endl;
				vtkSmartPointer<vtkPoints> linepoints = vtkSmartPointer<vtkPoints>::New();

				double node[3] = { grppoints[i][0],grppoints[i][1],grppoints[i][2] };
				auto x = grpcenterpoint(node, voxelsize, depth);

				double* joincenter = new double[3];
				joincenter = x.first;
				vector<double*> joinee = x.second;

				//double joincenterdata[3] = { joincenter[0],joincenter[1],joincenter[2] };

				cout << "joincenter is " << joincenter[0] << "   " << joincenter[1] << " " << joincenter[2] << endl;
				vtkNew<vtkCellArray> cells;
				linepoints->InsertNextPoint(joincenter); // main point 0

				for (int i = 0; i < joinee.size(); i++) {

					double vertex[3] = { joinee[i][0],joinee[i][1],joinee[i][2] };
					DoubleArray arrinsert = { {joinee[i][0], joinee[i][1], joinee[i][2]} };
					grpvisited.insert(arrinsert);
					linepoints->InsertNextPoint(vertex); // point i

					cout << "joinee isssssss " << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl;
					int z = i + 1;
					vtkIdType lineIndices[2] = { 0, z }; // Connect Point i to Point 0
					cells->InsertNextCell(2, lineIndices);

				}

				// Create a polydata to store everything in
				vtkNew<vtkPolyData> linepolyData;

				// Add the points to the dataset
				linepolyData->SetPoints(linepoints);

				// Add the lines to the dataset
				linepolyData->SetLines(cells);

				// Setup actor and mapper
				vtkNew<vtkPolyDataMapper> linemapper;
				linemapper->SetInputData(linepolyData);

				vtkNew<vtkActor> lineactor;
				lineactor->SetMapper(linemapper);
				lineactor->GetProperty()->SetColor(colors->GetColor3d("black").GetData());
				int line_width = 1;  //# Change this value to the desired line width
				lineactor->GetProperty()->SetLineWidth(line_width);
				renderer->AddActor(lineactor);
			}
		}


		//cout << "printing grpvisited ===================" << endl;
		/*for (double* ptr : grpvisited) {
			std::cout << "Array elements: ";
			for (int i = 0; i < 3; ++i) {
				std::cout << ptr[i] << " ";
			}
			std::cout << std::endl;
		}*/

		//default
		grppoints.clear();
		//for (const double* arr : grpnextpoints) {
		//	// Create a new double array and copy the values
		//	double* newarr = new double[3];
		//	for (int i = 0; i < 3; ++i) {
		//		newarr[i] = arr[i];
		//	}
		//	grppoints.push_back(newarr);
		//}
		grppoints = grpnextpoints;
		grpnextpoints.clear();
		grpvisited.clear();
	}

	//bottom_y changed so support fall till the machinecube floor
	//bottom_y = 0;
	//final layer of going straight down
	bottom_y = lowest_y-(0.5*voxelsize);

	for (int i = 0; i < grppoints.size(); i++) {
		vtkSmartPointer<vtkPoints> linepoints = vtkSmartPointer<vtkPoints>::New();

		double node[3] = { grppoints[i][0],grppoints[i][1],grppoints[i][2] }; // main point

		linepoints->InsertNextPoint(node);
		vtkNew<vtkCellArray> cells;

		if (node[1] > bottom_y) {
			double connect[3] = { node[0],bottom_y,node[2] };
			linepoints->InsertNextPoint(connect);
			vtkIdType lineIndices[2] = { 0, 1 }; // Connect Point 4 to Point 0
			cells->InsertNextCell(2, lineIndices);
			// Create a polydata to store everything in
			vtkNew<vtkPolyData> linepolyData;

			// Add the points to the dataset
			linepolyData->SetPoints(linepoints);

			// Add the lines to the dataset
			linepolyData->SetLines(cells);

			// Setup actor and mapper
			vtkNew<vtkPolyDataMapper> linemapper;
			linemapper->SetInputData(linepolyData);

			vtkNew<vtkActor> lineactor;
			lineactor->SetMapper(linemapper);
			lineactor->GetProperty()->SetColor(colors->GetColor3d("black").GetData());
			int line_width = 5;  //# Change this value to the desired line width
			lineactor->GetProperty()->SetLineWidth(line_width);
			renderer->AddActor(lineactor);

			/*	stlWriter->SetInputData(linepolyData);*/
		}
		else continue;
	}



	auto end_time = std::chrono::high_resolution_clock::now();
	auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	cout << "execution time in ms " << execution_time << endl;
	double time = execution_time / 60000;
	std::cout << "Execution time: " << time << "minutes" << std::endl;

	// Create mapper and actor for the bounding box
	vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputData(cubeSource->GetOutput());
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->GetProperty()->SetOpacity(0.3);
	actor1->SetMapper(mapper1);

	cout << "a comp" << endl;
	// Add the actors to the scene 
	// ###$
	renderer->AddActor(triangleActor);
	renderer->AddActor(actor1);
	//renderer->AddActor(machineactor);
	//renderer->AddActor(axes);
	renderer->SetBackground(colors->GetColor3d("white").GetData());
	//previous midnight blue
	cout << "b" << endl;

	// Render an image (lights and cameras are created automatically)
	renderWindow->SetWindowName("OctreeVisualize");
	renderWindow->SetSize(600, 600);
	renderWindow->Render();

	renderWindowInteractor->Initialize();
	renderWindow->Render();
	renderWindowInteractor->Start();
	return EXIT_SUCCESS;
}
