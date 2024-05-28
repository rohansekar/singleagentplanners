#ifndef PLANNER_H
#define PLANNER_H
#include <vector>
struct cell
{
    float g{100000000}, h, f;
    int cell_x, cell_y, cell_t;
    cell *parent = nullptr;
    bool collision;
};
struct bestval
{
    int x, y, cost;
};
// Declare the plan function
void planner(
    int *map,
    int collision_thresh,
    int x_size,
    int y_size,
    int robotposeX,
    int robotposeY,
    int target_steps,
    int *target_traj,
    int targetposeX,
    int targetposeY,
    int curr_time,
    int *action_ptr);
class Planner
{
public:
    std::vector<std::vector<cell *>> grid;
    std::vector<int> optimalx, optimaly, optimalcost;
    cell *child;
    int robot_x, robot_y, target_x, target_y, map_x, map_y, collision_threshold, curr_time, time_elapsed, target_steps, buffer, trajpoint, cost, finaltargetx, finaltargety;
    int finalGoalX, finalGoalY;
    int *target_traj;
    int maxcost = 100000000;
    std::vector<std::pair<int, int>> start;
    int *map;
    int dx[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int dy[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    std::vector<std::pair<int, int>> open1;
    float heuristic(int x, int y);
    void make_grid(int *map);
    void find_target(bool heuristicflag);
    std::vector<std::pair<int, int>> path(int target_x, int target_y);
    void reset_values();
};

#endif
// PLANNER_H