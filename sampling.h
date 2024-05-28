#include <vector>
using namespace std;
struct Node
{
    double *joint;
    Node *parent = nullptr;
    double cost = 0;
};

class RRTclass
{

public:
    RRTclass(int numDOFs, double *start_angles, double *goal_angles);
    int nearestnode(double *random_sample_angles);
    vector<Node *> node;
    int numDOFs;
    double *start_angles;
    double *goal_angles;
};

class PRMclass
{
private:
    vector<double *> vertices;
    vector<vector<int>> edges;
    int numofDOFs;
    double *start_angles_rad;
    double *goal_angles_rad;

public:
    PRMclass(int num, double *armstart_anglesV_rad, double *armgoal_anglesV_rad);
    unordered_map<int, double> g;
    void preprocess(double *map, double *armgoal_anglesV_rad, int N, double R, int k, int x_size, int y_size);
    int Vertex(double *sample_angelsV_rad);

    void Edge(int id1, int id2);

    int numofVertices();

    double *GetVertex(int id);

    vector<int> successors(int id);
};
class RRTStarTree
{
private:
    vector<double *> vertices;
    unordered_map<int, int> edges;
    int numofDOFs;

public:
    RRTStarTree(int num, double *armstart_anglesV_rad);

    int Vertex(double *sample_angelsV_rad);

    int LastElement();

    void Edge(int sid, int eid);

    int nearestnode(double *joint_angles);

    double *GetVertex(int id);

    int numofVertices();

    int GetParent(int id);
};
