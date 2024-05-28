/*=================================================================
 *
 * planner.cpp
 *
 *=================================================================*/
#include "planner.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <queue>
#include <chrono>
#include <thread>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) ((Y - 1) * XSIZE + (X - 1))

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#define NUMOFDIRS 8
Planner *p = new Planner;
int cost = 0;
int finaltargetx = 0;
int finaltargety = 0;
bool plannerflag = false;
bool plannerflag2 = false;
bool heuristicflag = false;
bool breakflag = false;
bool flag3 = false;
int i = 0;
using namespace std;
class lowestf
{
public:
    bool operator()(const cell *a, const cell *b)
    {
        return a->f > b->f;
    }
};
class lowestcost
{
public:
    bool operator()(const bestval *a, const bestval *b)
    {
        return a->cost > b->cost;
    }
};
float Planner::heuristic(int x, int y)
{
    return sqrt(pow(x - target_x, 2) + pow(y - target_y, 2));
}
void Planner::make_grid(int *map)
{

    for (int i = 0; i < map_x; i++)
    {
        vector<cell *> seg;

        for (int j = 0; j < map_y; j++)
        {
            cell *square = new cell;
            square->cell_x = i;
            square->cell_y = j;
            if ((i > 0 && j > 0 && i < map_x - 1 && j < map_y - 1))
            {
                if (int(map[GETMAPINDEX(i, j, map_x, map_y)]) >= collision_threshold)
                {
                    square->collision = true;
                }
            }

            seg.push_back(square);
        }
        grid.push_back(seg);
    }
}

void Planner::find_target(bool heuristicflag)
{
    std::priority_queue<cell *, vector<cell *>, lowestf> open;
    reset_values();
    int x = robot_x;
    int y = robot_y;
    int i_t = 0;
    int new_t = 0;
    cell *curr = grid[robot_x][robot_y];
    curr->cell_x = robot_x;
    curr->cell_y = robot_y;
    int newcost = 0;
    curr->g = 0;
    curr->h = heuristic(curr->cell_x, curr->cell_y);
    curr->f = curr->g + curr->h;
    open.push(curr);
    open1.clear();
    while (open.size() > 0)
    {
        curr = open.top();
        open.pop();
        x = curr->cell_x;
        y = curr->cell_y;
        i_t = curr_time;

        for (int i = 0; i < NUMOFDIRS; i++)
        {
            int new_x = x + dx[i];
            int new_y = y + dy[i];
            new_t = i_t + 1;

            if ((new_x > 1 && new_y > 1 && new_x < map_x - 1 && new_y < map_y - 1))
            {
                child = grid[new_x][new_y];
                child->cell_x = new_x;
                child->cell_y = new_y;
                child->cell_t = new_t;
                if (!heuristicflag)
                {
                    newcost = int(map[GETMAPINDEX(child->cell_x, child->cell_y, map_x, map_y)]);
                }
                else
                {
                    newcost = 1;
                }
                if ((!child->collision) && (child->g > curr->g + newcost))
                {
                    child->g = curr->g + newcost;
                    child->h = heuristic(child->cell_x, child->cell_y);
                    child->f = child->g + child->h;
                    child->parent = curr;
                    open.push(child);
                }
            }
            if (curr->cell_x == target_x && curr->cell_y == target_y)
            {
                path(target_x, target_y);
                break;
            }
        }
    }
}
vector<std::pair<int, int>> Planner::path(int target_x, int target_y)
{
    cell *start;
    start = grid[target_x][target_y];
    start->cell_x = target_x;
    start->cell_y = target_y;
    while (start->parent != nullptr)
    {
        open1.push_back(make_pair(start->cell_x, start->cell_y));
        start = start->parent;
    }
    return open1;
}
void Planner::reset_values()
{
    for (int i = 1; i < p->map_x - 1; i++)
    {
        for (int j = 1; j < p->map_y - 1; j++)
        {
            p->grid[i][j]->g = 1000000;
            p->grid[i][j]->f = 1000000;
            p->grid[i][j]->parent = nullptr;
        }
    }
}
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
    int *action_ptr)
{
    // 8-connected grid
    int dX[NUMOFDIRS] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int dY[NUMOFDIRS] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int buffer;
    int time_elapsed;
    int goalposeX = target_traj[target_steps - 1];
    int goalposeY = target_traj[target_steps - 1 + target_steps];
    clock_t tStart = clock();
    buffer = 20;
    p->buffer = buffer;
    p->cost = cost;
    p->finaltargetx = finaltargetx;
    p->finaltargety = finaltargety;
    p->map_x = x_size;
    p->map_y = y_size;
    p->target_steps = target_steps;
    p->collision_threshold = collision_thresh;
    p->map = map;
    int t = 15;
    static int count = 0;
    time_elapsed = buffer + (int)((clock() - tStart) / CLOCKS_PER_SEC);
    p->time_elapsed = time_elapsed;
    if (!flag3)
    {
        for (static bool first = true; first; first = false)
        {
            p->make_grid(map);
        }
        flag3 = true;
    }
    if (!plannerflag)
    {
        if (MAX(p->map_x, p->map_y) > 1500)
        {
            t = 100;
            buffer = 20;
            heuristicflag = true;
        }
        for (int trajPoint = 0; trajPoint < target_steps; (trajPoint++))
        {
            if (trajPoint % t != 0 || trajPoint == target_steps - 1)
            {
                continue;
            }
            static int leastCost = 1000000;
            p->robot_x = robotposeX;
            p->robot_y = robotposeY;
            goalposeX = (int)target_traj[trajPoint];
            goalposeY = (int)target_traj[trajPoint + target_steps];
            int CurrentCost = 0;
            p->target_x = goalposeX;
            p->target_y = goalposeY;
            p->curr_time = curr_time;
            p->target_traj = target_traj;
            p->find_target(heuristicflag);
            float thresh = 0.5;
            bool caught = false;
            targetposeX = target_traj[trajPoint];
            targetposeY = target_traj[trajPoint + target_steps];
            if (p->open1.size() != 0)
            {
                if (p->open1.size() + buffer <= trajPoint)
                {
                    CurrentCost = p->open1.size() + (int)(trajPoint - p->open1.size()) * int(map[GETMAPINDEX(goalposeX, goalposeY, p->map_x, p->map_y)]);
                    if (CurrentCost < leastCost)
                    {
                        leastCost = CurrentCost;
                        p->finalGoalX = goalposeX;
                        p->finalGoalY = goalposeY;
                        count = count + 1;
                    }
                }
            }
        }
        if (count == 0)
        {
            heuristicflag = true;
            p->finalGoalX = target_traj[target_steps - 1];
            p->finalGoalY = target_traj[target_steps - 1 + target_steps];
        }
        plannerflag = true;
    }

    p->robot_x = robotposeX;
    p->robot_y = robotposeY;
    goalposeX = p->finalGoalX;
    goalposeY = p->finalGoalY;
    p->target_x = p->finalGoalX;
    p->target_y = p->finalGoalY;
    if (goalposeX == robotposeX && goalposeY == robotposeY)
    {
        action_ptr[0] = robotposeX;
        action_ptr[1] = robotposeY;
        return;
    }
    if (!breakflag)
    {
        p->find_target(heuristicflag);
        reverse(p->open1.begin(), p->open1.end());
        breakflag = true;
    }

    plannerflag2 = false;
    action_ptr[0] = p->open1[i].first;
    action_ptr[1] = p->open1[i].second;
    i += 1;
    return;
}
