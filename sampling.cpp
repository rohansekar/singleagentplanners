/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <string>
#include <stdexcept>
#include <regex>	// For regex and split logic
#include <iostream> // cout, endl
#include <fstream>	// For reading/writing files
#include <assert.h>

#include "planner.h"

/* Input Arguments */
#define MAP_IN prhs[0]
#define ARMSTART_IN prhs[1]
#define ARMGOAL_IN prhs[2]
#define PLANNER_ID_IN prhs[3]

/* Planner Ids */
#define RRT 0
#define RRTCONNECT 1
#define RRTSTAR 2
#define PRM 3
#define Epsilon 2

/* Output Arguments */
#define PLAN_OUT plhs[0]
#define PLANLENGTH_OUT plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y * XSIZE + X)

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

// the length of each link in the arm
#define LINKLENGTH_CELLS 10

// Some potentially helpful imports
using std::array;
using std::cout;
using std::endl;
using std::make_tuple;
using std::runtime_error;
using std::string;
using std::tie;
using std::tuple;
using std::vector;
using namespace std;
// RRTclass *rrt = new RRTclass;
//*******************************************************************************************************************//
//                                                                                                                   //
//                                                GIVEN FUNCTIONS                                                    //
//                                                                                                                   //
//*******************************************************************************************************************//

/// @brief
/// @param filepath
/// @return map, x_size, y_size
tuple<double *, int, int> loadMap(string filepath)
{
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f)
	{
	}
	else
	{
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2)
	{
		throw runtime_error("Invalid loadMap parsing map metadata");
	}

	////// Go through file and add to m_occupancy
	double *map = new double[height * width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			char c;
			do
			{
				if (fscanf(f, "%c", &c) != 1)
				{
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0'))
			{
				map[y + x * width] = 1; // Note transposed from visual
			}
			else
			{
				map[y + x * width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string &str, const string &delim)
{
	// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
	const std::regex ws_re(delim);
	return {std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator()};
}

double *doubleArrayFromString(string str)
{
	vector<string> vals = split(str, ",");
	double *ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i)
	{
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double *v1, double *v2, int size)
{
	for (int i = 0; i < size; ++i)
	{
		// cout << "v1 and v2" << v1[i] << v2[i] << endl;
		if (abs(v1[i] - v2[i]) > 1e-3)
		{
			cout << endl;
			return false;
		}
	}
	return true;
}

typedef struct
{
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;

void ContXY2Cell(double x, double y, short unsigned int *pX, short unsigned int *pY, int x_size, int y_size)
{
	double cellsize = 1.0;
	// take the nearest cell
	*pX = (int)(x / (double)(cellsize));
	if (x < 0)
		*pX = 0;
	if (*pX >= x_size)
		*pX = x_size - 1;

	*pY = (int)(y / (double)(cellsize));
	if (y < 0)
		*pY = 0;
	if (*pY >= y_size)
		*pY = y_size - 1;
}
double dist_joints(double *angle1, double *angle2, int numofDOFs)
{
	double dist = 0;
	for (int i = 0; i < numofDOFs; i++)
	{
		dist = dist + ((angle1[i] - angle2[i]) * (angle1[i] - angle2[i]));
	}
	return sqrt(dist);
}

void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
	params->UsingYIndex = 0;

	if (fabs((double)(p2y - p1y) / (double)(p2x - p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
	{
		params->Y1 = p1x;
		params->X1 = p1y;
		params->Y2 = p2x;
		params->X2 = p2y;
	}
	else
	{
		params->X1 = p1x;
		params->Y1 = p1y;
		params->X2 = p2x;
		params->Y2 = p2y;
	}

	if ((p2x - p1x) * (p2y - p1y) < 0)
	{
		params->Flipped = 1;
		params->Y1 = -params->Y1;
		params->Y2 = -params->Y2;
	}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX = params->X2 - params->X1;
	params->DeltaY = params->Y2 - params->Y1;

	params->IncrE = 2 * params->DeltaY * params->Increment;
	params->IncrNE = 2 * (params->DeltaY - params->DeltaX) * params->Increment;
	params->DTerm = (2 * params->DeltaY - params->DeltaX) * params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
	if (params->UsingYIndex)
	{
		*y = params->XIndex;
		*x = params->YIndex;
		if (params->Flipped)
			*x = -*x;
	}
	else
	{
		*x = params->XIndex;
		*y = params->YIndex;
		if (params->Flipped)
			*y = -*y;
	}
}

int get_next_point(bresenham_param_t *params)
{
	if (params->XIndex == params->X2)
	{
		return 0;
	}
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else
	{
		params->DTerm += params->IncrNE;
		params->YIndex += params->Increment;
	}
	return 1;
}

int IsValidLineSegment(double x0, double y0, double x1, double y1, double *map,
					   int x_size, int y_size)
{
	bresenham_param_t params;
	int nX, nY;
	short unsigned int nX0, nY0, nX1, nY1;

	// printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);

	// make sure the line segment is inside the environment
	if (x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	// printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	// iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do
	{
		get_current_point(&params, &nX, &nY);
		if (map[GETMAPINDEX(nX, nY, x_size, y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double *angles, int numofDOFs, double *map,
							int x_size, int y_size)
{
	double x0, y0, x1, y1;
	int i;

	// iterate through all the links starting with the base
	x1 = ((double)x_size) / 2.0;
	y1 = 0;
	for (i = 0; i < numofDOFs; i++)
	{
		// compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS * cos(2 * PI - angles[i]);
		y1 = y0 - LINKLENGTH_CELLS * sin(2 * PI - angles[i]);

		// check the validity of the corresponding line segment
		if (!IsValidLineSegment(x0, y0, x1, y1, map, x_size, y_size))
			return 0;
	}
	return 1;
}
double *extend(double *joint, double *next_joint, double *map, int numofDOFs, int x_size, int y_size, double epsilon)
{
	double distance = 0;
	int i, j;
	for (j = 0; j < numofDOFs; j++)
	{
		if (distance < fabs(joint[j] - next_joint[j]))
			distance = fabs(joint[j] - next_joint[j]);
	}
	int numofsamples = (int)(distance / (PI / 20));
	double *temp_joints = (double *)malloc(numofDOFs * sizeof(double));
	double *final_joints = NULL;
	for (i = 1; i < numofsamples; i++)
	{
		for (j = 0; j < numofDOFs; j++)
		{
			temp_joints[j] = joint[j] + ((double)(i) / (numofsamples - 1)) * (next_joint[j] - joint[j]);
		}
		if (IsValidArmConfiguration(temp_joints, numofDOFs, map, x_size, y_size) &&
			dist_joints(temp_joints, joint, numofDOFs) < epsilon)
		{
			if (i == 1)
			{
				final_joints = (double *)malloc(numofDOFs * sizeof(double));
			}
			memcpy(final_joints, temp_joints, numofDOFs * sizeof(double));
		}
		else
		{
			break;
		}
	}
	free(temp_joints);
	return final_joints;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                          DEFAULT PLANNER FUNCTION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void planner(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	// no plan by default
	*plan = NULL;
	*planlength = 0;

	// for now just do straight interpolation between start and goal checking for the validity of samples

	double distance = 0;
	int i, j;
	for (j = 0; j < numofDOFs; j++)
	{
		if (distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
			distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
	}
	int numofsamples = (int)(distance / (PI / 20));
	if (numofsamples < 2)
	{
		printf("the arm is already at the goal\n");
		return;
	}
	*plan = (double **)malloc(numofsamples * sizeof(double *));
	int firstinvalidconf = 1;
	for (i = 0; i < numofsamples; i++)
	{
		(*plan)[i] = (double *)malloc(numofDOFs * sizeof(double));
		for (j = 0; j < numofDOFs; j++)
		{
			(*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i) / (numofsamples - 1)) * (armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
		}
		if (!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
		{
			firstinvalidconf = 1;
			printf("ERROR: Invalid arm configuration!!!\n");
		}
	}
	*planlength = numofsamples;

	return;
}
int RRTclass::nearestnode(double *random_joints)
{
	double min_dist = 100000;
	int min_index = -1;
	for (int i = 0; i < node.size(); i++)
	{
		double dist = 0;
		for (int j = 0; j < numDOFs; j++)
		{
			dist += pow(random_joints[j] - node[i]->joint[j], 2);
		}
		if (dist < min_dist)
		{
			min_dist = dist;
			min_index = i;
		}
	}
	return min_index;
}

RRTclass::RRTclass(int numDOFs_, double *start_angles_, double *goal_angles_) : numDOFs(numDOFs_), start_angles(start_angles_), goal_angles(goal_angles_)
{
	node.clear();
	Node *start_node = new Node;
	start_node->joint = start_angles;
	start_node->cost = 0;
	start_node->parent = NULL;
	node.push_back(start_node);
}
void randomnode(double *random_angles, double *map, int numofDOFs, int x_size, int y_size)
{
	while (true)
	{
		for (int i = 0; i < numofDOFs; i++)
		{
			int temp = rand() % 360;
			random_angles[i] = (double)temp / 360 * 2 * PI;
		}
		if (IsValidArmConfiguration(random_angles, numofDOFs, map, x_size, y_size))
		{
			break;
		}
	}
}
void randomnode1(double *random_angles, double *map, int numofDOFs, int x_size, int y_size)
{
	srand(time(0));
	while (true)
	{
		for (int i = 0; i < numofDOFs; i++)
		{
			int temp = rand() % 360;
			random_angles[i] = (double)temp / 360 * 2 * PI;
		}
		if (IsValidArmConfiguration(random_angles, numofDOFs, map, x_size, y_size))
		{
			break;
		}
	}
}

// Check if the two nodes angles are equal
bool is_angles_equal(double *angle1, double *angle2, int numofDOFs)
{
	if (angle1 == NULL || angle2 == NULL)
		return false;
	for (int i = 0; i < numofDOFs; i++)
	{
		if (angle1[i] != angle2[i])
		{
			return false;
			break;
		}
	}
	return true;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              RRT IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRT(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	*plan = NULL;
	*planlength = 0;

	RRTclass tree(numofDOFs, armstart_anglesV_rad, armgoal_anglesV_rad);
	int samples = 0;
	bool reached = false;
	double *random_joints = new double[numofDOFs];
	while (!reached && samples < 50000)
	{

		samples++;
		randomnode(random_joints, map, numofDOFs, x_size, y_size);
		int nearest_node_index = tree.nearestnode(random_joints);
		double *extended_joints = extend(tree.node[nearest_node_index]->joint, random_joints, map, numofDOFs, x_size, y_size, Epsilon);
		if (extended_joints != NULL && (IsValidArmConfiguration(extended_joints, numofDOFs, map, x_size, y_size)))
		{
			Node *new_node = new Node;
			new_node->joint = extended_joints;
			new_node->parent = tree.node[nearest_node_index];

			tree.node.push_back(new_node);
			if (dist_joints(extended_joints, armgoal_anglesV_rad, numofDOFs) < Epsilon)
			{
				reached = true;
				Node *goal_node = new Node;
				goal_node->joint = armgoal_anglesV_rad;
				goal_node->parent = tree.node[tree.node.size() - 1];

				tree.node.push_back(goal_node);
				cout << "Path Found" << endl;
				break;
			}
		}
	}
	if (reached)
	{
		vector<Node *> path;
		Node *current_node = tree.node[tree.node.size() - 1];
		while (current_node->parent != NULL)
		{
			path.push_back(current_node);
			current_node = current_node->parent;
		}
		path.push_back(current_node);
		*plan = (double **)malloc(path.size() * sizeof(double *));
		for (int i = 0; i < path.size(); ++i)
		{
			(*plan)[i] = (double *)malloc(numofDOFs * sizeof(double));

			(*plan)[i] = path[path.size() - 1 - i]->joint;
		}
		*planlength = path.size();
	}
	else
	{
		cout << "Path Not Found" << endl;
	}
	return;
	// plannerRRT(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                         RRT CONNECT IMPLEMENTATION                                                //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRTConnect(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	*plan = NULL;
	*planlength = 0;
	RRTclass start_tree(numofDOFs, armstart_anglesV_rad, armgoal_anglesV_rad);
	RRTclass end_tree(numofDOFs, armgoal_anglesV_rad, armstart_anglesV_rad);
	int samples = 0;
	bool reached = false;
	double *random_joints = new double[numofDOFs];
	double *extended_joints = new double[numofDOFs];
	double *extended_joints_end = new double[numofDOFs];
	bool front = true;
	bool back = true;
	while (!reached && samples < 50000)
	{
		samples++;
		if (front)
		{
			randomnode(random_joints, map, numofDOFs, x_size, y_size);
			int nearest_node_index = start_tree.nearestnode(random_joints);
			double *extended_joints = extend(start_tree.node[nearest_node_index]->joint, random_joints, map, numofDOFs, x_size, y_size, Epsilon);
			if (extended_joints == NULL)
				continue;
			Node *new_node = new Node;
			new_node->joint = extended_joints;
			new_node->parent = start_tree.node[nearest_node_index];
			start_tree.node.push_back(new_node);
			int nearest_node_index_end = end_tree.nearestnode(extended_joints);
			double *extended_joints_end = extend(end_tree.node[nearest_node_index_end]->joint, extended_joints, map, numofDOFs, x_size, y_size, Epsilon);
			if (extended_joints_end == NULL)
				continue;
			Node *new_node_end = new Node;
			new_node_end->joint = extended_joints_end;
			new_node_end->parent = end_tree.node[nearest_node_index_end];
			end_tree.node.push_back(new_node_end);
			if (dist_joints(extended_joints, extended_joints_end, numofDOFs) < Epsilon)
			{
				reached = true;

				cout << "Path Found" << endl;

				break;
			}
		}
		else
		{
			randomnode(random_joints, map, numofDOFs, x_size, y_size);
			int nearest_node_index = end_tree.nearestnode(random_joints);
			double *extended_joints_end = extend(end_tree.node[nearest_node_index]->joint, random_joints, map, numofDOFs, x_size, y_size, Epsilon);
			if (extended_joints_end == NULL)
				continue;
			Node *new_node = new Node;
			new_node->joint = extended_joints_end;
			new_node->parent = end_tree.node[nearest_node_index];
			end_tree.node.push_back(new_node);

			int nearest_node_index_start = start_tree.nearestnode(extended_joints_end);
			double *extended_angles_rad_start = extend(start_tree.node[nearest_node_index_start]->joint, extended_joints_end, map, numofDOFs, x_size, y_size, Epsilon);
			if (extended_angles_rad_start == NULL)
				continue;
			Node *new_node_start = new Node;
			new_node_start->joint = extended_angles_rad_start;
			new_node_start->parent = start_tree.node[nearest_node_index_start];
			start_tree.node.push_back(new_node_start);

			if (dist_joints(extended_joints_end, extended_angles_rad_start, numofDOFs) < Epsilon)
			{
				reached = true;

				cout << "Path Found" << endl;
				break;
			}
		}
		front = !front;
	}
	if (reached)
	{
		vector<Node *> path;
		Node *current_node = start_tree.node[0];

		Node *current_node_end = end_tree.node[end_tree.node.size() - 2];
		// }
		static int j = 0;
		static int k = end_tree.node.size() - 2;
		while (j <= start_tree.node.size() - 1)
		{
			path.push_back(current_node);
			j += 1;
			current_node = start_tree.node[j];
		}
		while (k >= 0)
		{
			path.push_back(current_node_end);
			k -= 1;
			current_node_end = end_tree.node[k];
		}

		*plan = (double **)malloc(path.size() * sizeof(double *));
		int x = 0;
		for (int i = 0; i < path.size(); i++)
		{
			x += 1;
			(*plan)[i] = (double *)malloc(numofDOFs * sizeof(double));
			(*plan)[i] = path[i]->joint;
		}

		*planlength = path.size();
	}
	else
	{
		cout << "Path Not Found" << endl;
	}
	// plannerRRTConnect(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                           RRT STAR IMPLEMENTATION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRTStar(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	*plan = NULL;
	*planlength = 0;
	int samples = 0;

	RRTStarTree tree(numofDOFs, armstart_anglesV_rad);

	bool found = false;

	double *random_joints = new double[numofDOFs];

	while (!found && samples < 50000)
	{
		samples++;
		randomnode(random_joints, map, numofDOFs, x_size, y_size);

		int nn_id = tree.nearestnode(random_joints);
		double *extended_joints = extend(tree.GetVertex(nn_id), random_joints, map, numofDOFs, x_size, y_size, Epsilon);
		if (extended_joints != NULL)
		{
			int extend_id = tree.Vertex(extended_joints);
			tree.Edge(nn_id, extend_id);
			extended_joints = extend(extended_joints, armgoal_anglesV_rad, map, numofDOFs, x_size, y_size, Epsilon);
			if (is_angles_equal(extended_joints, armgoal_anglesV_rad, numofDOFs))
			{
				int goal_id = tree.Vertex(armgoal_anglesV_rad);
				tree.Edge(extend_id, goal_id);
				found = true;
				break;
			}
		}
	}
	if (found)
	{
		cout << "Path Found" << endl;
		vector<int> path;
		int next_id = tree.LastElement();
		while (next_id != 0)
		{
			path.insert(path.begin(), next_id);
			next_id = tree.GetParent(next_id);
		}
		path.insert(path.begin(), 0);

		*planlength = path.size();

		*plan = (double **)malloc(path.size() * sizeof(double *));

		for (int i = 0; i < path.size(); i++)
		{
			(*plan)[i] = (double *)malloc(numofDOFs * sizeof(double));
			(*plan)[i] = tree.GetVertex(path[i]);
		}
	}

	else
	{
		cout << "path not found" << endl;
	}

	return;
	// planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              PRM IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerPRM(
	double *map,
	int x_size,
	int y_size,
	double *armstart_anglesV_rad,
	double *armgoal_anglesV_rad,
	int numofDOFs,
	double ***plan,
	int *planlength)
{
	/* TODO: Replace with your implementation */
	*plan = NULL;
	*planlength = 0;
	PRMclass graph(numofDOFs, armstart_anglesV_rad, armgoal_anglesV_rad);
	graph.preprocess(map, armgoal_anglesV_rad, 7000, Epsilon, 5, x_size, y_size);
	unordered_set<int> closed, open;
	unordered_map<int, int> parent;

	bool found = false;
	open.insert(0);
	graph.g[0] = 0;
	int goal_id = 1;

	int vertices_expanded = 0;
	while (!open.empty() && !found)
	{
		unordered_set<int>::iterator close_index;
		double min_f = 10000;

		vertices_expanded++;

		for (auto itr = open.begin(); itr != open.end(); ++itr)
		{
			double f = graph.g[*itr] + dist_joints(graph.GetVertex(*itr), graph.GetVertex(goal_id), numofDOFs);
			if (f < min_f)
			{
				min_f = f;
				close_index = itr;
			}
		}
		int current_id = *close_index;

		closed.insert(current_id);
		open.erase(close_index);

		vector<int> successors = graph.successors(current_id);
		for (auto id : successors)
		{
			if (id == goal_id)
			{
				found = true;
				parent[goal_id] = current_id;
				break;
			}

			if (closed.count(id) == 0)
			{
				if (open.count(id) == 0)
				{
					open.insert(id);
					graph.g[id] = graph.g[current_id] + dist_joints(graph.GetVertex(id), graph.GetVertex(current_id), numofDOFs);
					parent[id] = current_id;
				}
				if (open.count(id) != 0)
				{
					if (graph.g[id] > graph.g[current_id] + dist_joints(graph.GetVertex(id), graph.GetVertex(current_id), numofDOFs))
					{
						graph.g[id] = graph.g[current_id] + dist_joints(graph.GetVertex(id), graph.GetVertex(current_id), numofDOFs);
						parent[id] = current_id;
					}
				}
			}
		}
	}
	if (found)
	{
		cout << "Path Found" << endl;
		vector<int> path;
		int follow_id = goal_id;
		while (follow_id != 0)
		{
			path.insert(path.begin(), follow_id);
			follow_id = parent[follow_id];
		}
		path.insert(path.begin(), 0);

		*planlength = path.size();

		*plan = (double **)malloc(path.size() * sizeof(double *));

		for (int i = 0; i < path.size(); i++)
		{
			(*plan)[i] = (double *)malloc(numofDOFs * sizeof(double));
			(*plan)[i] = graph.GetVertex(path[i]);
		}
	}

	else
	{
		cout << "No Path Found" << endl;
	}
	return;
	// planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                MAIN FUNCTION                                                      //
//                                                                                                                   //
//*******************************************************************************************************************//

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos,
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char **argv)
{
	double *map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double *startPos = doubleArrayFromString(argv[3]);
	double *goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if (!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size) ||
		!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size))
	{
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double **plan = NULL;
	int planlength = 0;
	double *random_joints = new double[3];
	randomnode1(random_joints, map, 3, x_size, y_size);

	// Call the corresponding planner function
	if (whichPlanner == PRM)
	{
		plannerPRM(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else if (whichPlanner == RRT)
	{
		plannerRRT(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else if (whichPlanner == RRTCONNECT)
	{
		plannerRRTConnect(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else if (whichPlanner == RRTSTAR)
	{
		plannerRRTStar(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}
	else
	{
		planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as my
	//// grading script will not work and you will recieve a 0.
	///////////////////////////////////////

	// Your solution's path should start with startPos and end with goalPos
	if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) ||
		!equalDoubleArrays(plan[planlength - 1], goalPos, numOfDOFs))
	{
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open())
	{
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i)
	{
		for (int k = 0; k < numOfDOFs; ++k)
		{
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
PRMclass::PRMclass(int num, double *armstart_anglesV_rad, double *armgoal_anglesV_rad)
{
	numofDOFs = num;
	vertices.clear();
	edges.clear();
	g.clear();
	int start_id = Vertex(armstart_anglesV_rad);
	start_angles_rad = armstart_anglesV_rad;
	goal_angles_rad = armgoal_anglesV_rad;
}
void PRMclass::preprocess(double *map, double *armgoal_anglesV_rad, int N, double R, int k, int x_size, int y_size)
{

	int i = 0;
	int goal_id = Vertex(armgoal_anglesV_rad);

	int edges = 0;

	while (i < N)
	{

		i++;
		int connected = 0;

		double *random_joints = (double *)malloc(numofDOFs * sizeof(double));

		randomnode(random_joints, map, numofDOFs, x_size, y_size);

		double *temp_joints = (double *)malloc(numofDOFs * sizeof(double));
		temp_joints = random_joints;

		int random_id = Vertex(temp_joints);

		for (int j = 0; j < vertices.size() - 1; j++)
		{
			if (dist_joints(vertices[j], temp_joints, numofDOFs) < R)
			{
				double *extended_joints = extend(vertices[j], temp_joints, map, numofDOFs, x_size, y_size, R);
				if (is_angles_equal(extended_joints, temp_joints, numofDOFs))
				{
					Edge(j, random_id);
					connected++;
					edges++;
					if (connected >= k)
						break;
				}
			}
		}
	}
}
int PRMclass::Vertex(double *sample_angelsV_rad)
{
	int vid = vertices.size();
	vector<int> first;
	vertices.push_back(sample_angelsV_rad);
	edges.push_back(first);

	return vid;
}

void PRMclass::Edge(int id1, int id2)
{
	edges[id2].push_back(id1);
	edges[id1].push_back(id2);
}

int PRMclass::numofVertices()
{
	return vertices.size();
}

double *PRMclass::GetVertex(int id)
{
	return vertices[id];
}

vector<int> PRMclass::successors(int id)
{
	return edges[id];
}
RRTStarTree::RRTStarTree(int num, double *armstart_anglesV_rad)
{
	numofDOFs = num;
	vertices.clear();
	edges.clear();
	vertices.push_back(armstart_anglesV_rad);
}
double *RRTStarTree::GetVertex(int id)
{
	return vertices[id];
}

int RRTStarTree::numofVertices()
{
	return vertices.size();
}

int RRTStarTree::GetParent(int id)
{
	return edges[id];
}
int RRTStarTree::Vertex(double *sample_angelsV_rad)
{
	int vid = vertices.size();
	vertices.push_back(sample_angelsV_rad);

	return vid;
}

int RRTStarTree::LastElement()
{
	return vertices.size() - 1;
}

void RRTStarTree::Edge(int sid, int eid)
{
	edges[eid] = sid;
}

int RRTStarTree::nearestnode(double *joint_angles)
{
	int nearest_id = 0;
	double nearest_dist = 100000;
	for (int i = 0; i < vertices.size(); i++)
	{
		int tmp_dist = 0;
		for (int j = 0; j < numofDOFs; j++)
		{
			tmp_dist += pow((vertices[i][j] - joint_angles[j]),2);
		}
		tmp_dist = sqrt(tmp_dist);
		if (tmp_dist < nearest_dist)
		{
			nearest_dist = tmp_dist;
			nearest_id = i;
		}
	}

	return nearest_id;
}