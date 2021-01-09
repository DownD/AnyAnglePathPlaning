/*
This code was developed by Tansel Uras (turas@usc.edu) at USC.
The code is hosted at 'http://idm-lab.org/anyangle'.
If you use this code in your research, please  cite our SoCS paper:

T. Uras and S. Koenig,  2015. An Empirical Comparison of Any-Angle Path-Planning Algorithms. In: Proceedings of the 8th Annual Symposium on Combinatorial
Search. Code available at: http://idm-lab.org/anyangle

Bibtex:
@inproceedings{uras:15,
  author = "T. Uras and S. Koenig",
  title = "An Empirical Comparison of Any-Angle Path-Planning Algorithms",
  booktitle = {Proceedings of the 8th Annual Symposium on Combinatorial Search},
  year = "2015",
  note = "Code available at: http://idm-lab.org/anyangle",
}
*/

#include "ANYA.h"

#ifdef ANY_ANGLE_RUNNING_IN_HOG
ANYA::ANYA(MapEnvironment *env) :
    AnyAngleAlgorithm::AnyAngleAlgorithm(env)
{
    ANYAInitialize();
}

//#define ANYA_DEBUG
#endif

ANYA::ANYA(std::vector<bool> &bits, int _width, int _height) :
    AnyAngleAlgorithm::AnyAngleAlgorithm(bits, _width, _height)
{
    ANYAInitialize();
}
void ANYA::ANYAInitialize()
{
    ComputeClearances();

    states_.reserve(10000);
    open_list_.reserve(10000);

#ifdef ANY_ANGLE_STATISTICS
    SetupStatistics();
#endif
}
ANYA::~ANYA()
{
    for (unsigned int x = 0; x < width_ - 1; x++) {
        delete[] corners_[x];
    }
    delete[] corners_;
}

void ANYA::ComputeClearances()
{
    corners_ = new Corner*[width_ - 1];
    for (unsigned int x = 0; x < width_ - 1; x++) {
        corners_[x] = new Corner[height_ - 1];
    }

    // Left clearances
    for (int y = 0; y < (int) height_ - 1; y++) {
        int clearance = 0;
        int topClearance = 0;
        int botClearance = 0;

        bool topBlocked = true;
        bool botBlocked = true;

        for (int x = 0; x < (int) width_ - 1; x++) {
            bool nextTopBlocked = !IsTraversable(NorthEastCell(x, y));
            bool nextBotBlocked = !IsTraversable(SouthEastCell(x, y));

            if (topBlocked && botBlocked)
                clearance = 0;

            corners_[x][y].interval_clearance[LEFT] = clearance;
            corners_[x][y].interval_blocked[TOP][LEFT] = topBlocked;
            corners_[x][y].interval_blocked[BOT][LEFT] = botBlocked;

            if (topBlocked)
                topClearance = 0;

            if (botBlocked)
                botClearance = 0;

            corners_[x][y].clearance[TOP][LEFT] = topClearance;
            corners_[x][y].clearance[BOT][LEFT] = botClearance;

            if (topBlocked != nextTopBlocked || botBlocked != nextBotBlocked) {
                clearance = 0;
                topBlocked = nextTopBlocked;
                botBlocked = nextBotBlocked;
            }

            clearance++;
            topClearance++;
            botClearance++;
        }
    }

    // Right clearances
    for (int y = 0; y < (int) height_ - 1; y++) {
        int clearance = 0;
        int topClearance = 0;
        int botClearance = 0;

        bool topBlocked = true;
        bool botBlocked = true;

        for (int x = (int) width_ - 2; x >= 0; x--) {
            bool nextTopBlocked = !IsTraversable(NorthWestCell(x, y));
            bool nextBotBlocked = !IsTraversable(SouthWestCell(x, y));

            if (topBlocked && botBlocked)
                clearance = 0;

            corners_[x][y].interval_clearance[RIGHT] = clearance;
            corners_[x][y].interval_blocked[TOP][RIGHT] = topBlocked;
            corners_[x][y].interval_blocked[BOT][RIGHT] = botBlocked;

            if (topBlocked)
                topClearance = 0;

            if (botBlocked)
                botClearance = 0;

            corners_[x][y].clearance[TOP][RIGHT] = topClearance;
            corners_[x][y].clearance[BOT][RIGHT] = botClearance;

            if (topBlocked != nextTopBlocked || botBlocked != nextBotBlocked) {
                clearance = 0;
                topBlocked = nextTopBlocked;
                botBlocked = nextBotBlocked;
                corners_[x][y].is_transition_point = true;
            }
            else
                corners_[x][y].is_transition_point = false;

            corners_[x][y].is_convex_corner = IsConvexCorner(x, y);

            clearance++;
            topClearance++;
            botClearance++;
        }
    }
}

const ANYAStateKey ANYA::GenerateKey(const xyLoc root, const Interval interval) const
{
    uint64_t hash = interval.GetHash();
    hash = (hash << 24) | (root.x << 12) | root.y;

    return ANYAStateKey(hash, interval.f_left, interval.f_right);
}

const stateId ANYA::GenerateState(const xyLoc root, const Interval interval, const xyLoc goal)
{
    // Generate the key for root, interval.
    ANYAStateKey key = GenerateKey(root, interval);
    stateId ind = anya_index_table_.Get(key);

    // If the state already exists, return its index.
    if (ind != ANYA_KEY_NOT_FOUND)
        return ind;

    // Otherwise, generate the state and add it to the hash table.
    ind = states_.size();
    states_.push_back(State(root, interval));
    anya_index_table_.Add(key, ind);

    // Initialize its values.
    states_.back().h_val = HeuristicDistance(root, interval, goal);
    states_.back().g_val = INFINITE_COST;
    states_.back().list = NO_LIST;

#ifdef ANY_ANGLE_STATISTICS
    num_generated_++;
#endif
    return ind;
}

void ANYA::GetSubIntervals(const int row, const double left, const double right, std::vector<Interval> & sub) const
{
    // Stretch the end-points to (x1, x2) so that they are integers (will be addressed later on).
    int x_left = floor(left + EPSILON);
    int x_right = ceil(right - EPSILON);

    // Start from the left end-point and find intervals towards right.
    int x = x_left;

    while (x < x_right) {
        int delta_x = corners_[x][row].interval_clearance[RIGHT];  // How far the next interval extends.

        Interval I(row, x, x + delta_x);
        I.blocked[TOP] = corners_[x][row].interval_blocked[TOP][RIGHT];
        I.blocked[BOT] = corners_[x][row].interval_blocked[BOT][RIGHT];
        sub.push_back(I);

        x += delta_x;
        if (delta_x == 0)
            break;
    }

    if (sub.empty())
        return;

    // Adjust the left end-point of the first interval and the right end-point of the second interval, in case the initial end-points are float values.
    sub.front().f_left = left;
    sub.back().f_right = right;
    sub.back().i_right = ceil(right - EPSILON);
}
void ANYA::GenerateSubIntervals(const int row, const double left, const double right, const xyLoc root, const xyLoc goal, std::vector<stateId> & ids)
{
    std::vector<Interval> intervals;

    GetSubIntervals(row, left, right, intervals);

    for (unsigned i = 0; i < intervals.size(); i++) {
        ids.push_back(GenerateState(root, intervals[i], goal));
    }
}

ANYA::Interval ANYA::GetLeftInterval(const int x, const int y) const
{
    Interval I(y, x - corners_[x][y].interval_clearance[LEFT], x);
    I.blocked[TOP] = corners_[x][y].interval_blocked[TOP][LEFT];
    I.blocked[BOT] = corners_[x][y].interval_blocked[BOT][LEFT];
    return I;
}
ANYA::Interval ANYA::GetRightInterval(const int x, const int y) const
{
    Interval I(y, x, x + corners_[x][y].interval_clearance[RIGHT]);
    I.blocked[TOP] = corners_[x][y].interval_blocked[TOP][RIGHT];
    I.blocked[BOT] = corners_[x][y].interval_blocked[BOT][RIGHT];
    return I;
}

void ANYA::GenerateIntervalsTowardsLeft(const double fx, const int y, const int unblockedDir, const xyLoc root, const xyLoc goal,
        std::vector<stateId> & successors, const double bound)
{
    // Determine how far the set of intervals extend and then call GenerateSubIntervals(..)
    int x = GetLeftColumn(fx);
    double left_bound = x - corners_[x][y].clearance[unblockedDir][LEFT];

    if (left_bound < bound)
        left_bound = bound;

    if (left_bound < fx) // TODO: maybe <= (same for intervals toward right)
    {
        GenerateSubIntervals(y, left_bound, fx, root, goal, successors);
    }
}
void ANYA::GenerateIntervalsTowardsRight(const double fx, const int y, const int unblockedDir, const xyLoc root, const xyLoc goal,
        std::vector<stateId> & successors, const double bound)
{
    // Determine how far the set of intervals extend and then call GenerateSubIntervals(..)
    int x = GetRightColumn(fx);
    double right_bound = x + corners_[x][y].clearance[unblockedDir][RIGHT];

    if (right_bound > bound)
        right_bound = bound;

    if (right_bound > fx) {
        GenerateSubIntervals(y, fx, right_bound, root, goal, successors);
    }
}

void ANYA::GenerateSuccessors_SameRow(const xyLoc root, const Interval interval, const xyLoc goal, std::vector<stateId> & successors)
{
    // We can assume that the end points of this interval are integers. // TODO: Why?

    // If the root is towards the right side of the interval, extend towards left.
    if (root.x >= interval.i_right) {
        xyLoc left_end_point(interval.i_left, interval.row);

        if (corners_[left_end_point.x][left_end_point.y].interval_clearance[LEFT] > 0) {    // If there is room to extend left
            // Observable successor	(on the same row).
            successors.push_back(GenerateState(root, GetLeftInterval(left_end_point.x, left_end_point.y), goal));

            // Unobservable successors (Row above).
            if (interval.blocked[TOP]) // If the top of the interval is not blocked, the string would not be taut
                GenerateIntervalsTowardsLeft(left_end_point.x, left_end_point.y - 1, BOT, left_end_point, goal, successors);

            // Unobservable successors (Row below).
            if (interval.blocked[BOT]) // If the bot of the interval is not blocked, the string would not be taut
                GenerateIntervalsTowardsLeft(left_end_point.x, left_end_point.y + 1, TOP, left_end_point, goal, successors);
        }
    }

    // If the root is towards the left side of the interval, extend towards right.
    if (root.x <= interval.i_left) {
        xyLoc right_end_point(interval.i_right, interval.row);

        if (corners_[right_end_point.x][right_end_point.y].interval_clearance[RIGHT] > 0) { // If there is room to extend right
            // Observable successor	(on the same row)
            successors.push_back(GenerateState(root, GetRightInterval(right_end_point.x, right_end_point.y), goal));

            // Unobservable successors (Row above)
            if (interval.blocked[TOP]) // If the top of the interval is not blocked, the string would not be taut
                GenerateIntervalsTowardsRight(right_end_point.x, right_end_point.y - 1, BOT, right_end_point, goal, successors);

            // Unobservable successors (Row below)
            if (interval.blocked[BOT]) // If the bot of the interval is not blocked, the string would not be taut
                GenerateIntervalsTowardsRight(right_end_point.x, right_end_point.y + 1, TOP, right_end_point, goal, successors);

        }
    }
}
void ANYA::GenerateSuccessors_DifferentRow(const xyLoc root, const Interval interval, const xyLoc goal, std::vector<stateId> & successors)
{
    // Notice that this function is called iff interval.row != root.y.

    // Determine the direction of the next interval.
    int new_row;    // In addition to the interval's row, we are looking to generate intervals on 'new_row".
    int dir;        // The 'new_row's direction relative to the interval's row (TOP or BOT).
    int rev_dir;    // Reverse of 'dir'.

    // If the root is to the bottom of the interval's row, the next row should be at towards the top of the interval's row.
    if (interval.row < root.y) {
        dir = TOP;
        rev_dir = BOT;
        new_row = interval.row - 1;
    }
    // If the root is to the top of the interval's row, the next row should be at towards the bottom of the interval's row.
    else {
        dir = BOT;
        rev_dir = TOP;
        new_row = interval.row + 1;
    }

    // Figure out the interval in the 'new_row' that is observable from root through 'interval'.

    // First, draw lines from the root towards the end-points of the 'interval' and figure out where they cross the next_row.
    double intersecting_left = GetIntersectingX(root.x, root.y, interval.f_left, interval.row, new_row);
    double intersecting_right = GetIntersectingX(root.x, root.y, interval.f_right, interval.row, new_row);

    // Adjust the bounds based on blocked cells (ignore the blocked cells that share an edge with the 'interval' for now).
    int left_bound = interval.i_left - corners_[interval.i_left][interval.row].clearance[dir][LEFT];
    int right_bound = interval.i_right + corners_[interval.i_right][interval.row].clearance[dir][RIGHT];

    // Any observable intervals generated on the next row should be contained within these end-points.
    double new_left = (intersecting_left < left_bound) ? left_bound : intersecting_left;
    double new_right = (intersecting_right > right_bound) ? right_bound : intersecting_right;

    // Generate the observable successors (if the 'interval' is blocked, there aren't any observable successors).
    if (!interval.blocked[dir])
        GenerateSubIntervals(new_row, new_left, new_right, root, goal, successors);

    // Generate the unobservable successors

    // TOWARDS LEFT:
    // We only generate unobservable successors towards left if the left end-point of the interval is an integer, and is at a convex corner of an obstacle.
    // This is so, because, unobservable successors are generated by a taut turn, which can only happen at a convex corner of an obstacle.

    // Assume that dir = TOP, that is, the root is towards the bottom of the 'interval'. Let, TL (Top-Left), TR (Top-Right), BL (Bottom-Left), BR (Bottom-Right)
    // be the cells surrounding the left end-point of the 'interval'.

    // TODO: Finish comments
    // - BR cannot be blocked because otherwise 'interval' cannot be observable by the root.
    // - If TR is blocked (and, therefore, TL is unblocked, because, otherwise, the left end-point will not be at a convex corner of an obstacle),

    // Check if we can generate an interval towards the left side of the 'interval'. For this, the left end-point of the interval should be at a convex corner
    // of an obstacle. If so, the blocked cell that produces the convex corner can either be towards 'dir' or 'rev_dir'
    // (there can also be two diagonally adjacent blocked cells that produce a single corner).

    if (fabs(interval.f_left - (double) interval.i_left) < EPSILON && corners_[interval.i_left][interval.row].is_convex_corner) {

        // If the blocked cell is towards 'rev_dir' (and, therefore, towards left, since 'interval' is visible from the root),
        // then, turning that corner produces a taut path.
        // If this is the case, then generate a new interval to the left of the 'interval', on the same row.
        if (corners_[interval.i_left][interval.row].interval_blocked[rev_dir][LEFT])
            successors.push_back(GenerateState(xyLoc(interval.i_left, interval.row), GetLeftInterval(interval.i_left, interval.row), goal));

        // Case for two diagonally adjacent obstacles.
        if (corners_[interval.i_left][interval.row].interval_blocked[dir][RIGHT] && corners_[interval.i_left][interval.row].interval_blocked[rev_dir][LEFT])
            GenerateIntervalsTowardsLeft(interval.i_left, new_row, rev_dir, xyLoc(interval.i_left, interval.row), goal, successors);

        else {
            if (corners_[interval.i_left][interval.row].interval_blocked[dir][LEFT]) // If the corner is from an obstacle between I's row and the next interval's row
                if (interval.i_left > root.x) // The root of the state should be towards the left of I for the path to be taut
                    GenerateIntervalsTowardsRight(interval.i_left, new_row, rev_dir, xyLoc(interval.i_left, interval.row), goal, successors, new_left);

            if (corners_[interval.i_left][interval.row].interval_blocked[rev_dir][LEFT]) // If the corner is from an obstacle that is away from the next interval's row
                if (new_left == intersecting_left)
                    GenerateIntervalsTowardsLeft(new_left, new_row, rev_dir, xyLoc(interval.i_left, interval.row), goal, successors);
        }

    }

    // Check if we can generate an interval towards the right of the expanded interval

    if (fabs(interval.f_right - (double) interval.i_right) < EPSILON && // If the right end of the interval is an integer
            corners_[interval.i_right][interval.row].is_convex_corner) // If the corresponding corner is a convex corner of an obstacle
    {
        if (corners_[interval.i_right][interval.row].interval_blocked[rev_dir][RIGHT]) // If turning that corner produces a taut path)
            successors.push_back(GenerateState(xyLoc(interval.i_right, interval.row), GetRightInterval(interval.i_right, interval.row), goal));

        // Case for two diagonally adjacent obstacles
        if (corners_[interval.i_right][interval.row].interval_blocked[dir][LEFT] && corners_[interval.i_right][interval.row].interval_blocked[rev_dir][RIGHT])
            GenerateIntervalsTowardsRight(interval.i_right, new_row, rev_dir, xyLoc(interval.i_right, interval.row), goal, successors);

        else {
            if (corners_[interval.i_right][interval.row].interval_blocked[dir][RIGHT]) // If the corner is from an obstacle between I's row and the next interval's row
                if (interval.i_right < root.x) // The root of the state should be towards the right of I for the path to be taut
                    GenerateIntervalsTowardsLeft(interval.i_right, new_row, rev_dir, xyLoc(interval.i_right, interval.row), goal, successors, new_right);

            if (corners_[interval.i_right][interval.row].interval_blocked[rev_dir][RIGHT]) // If the corner is from an obstacle that is away from the next interval's row
                if (new_right == intersecting_right)
                    GenerateIntervalsTowardsRight(new_right, new_row, rev_dir, xyLoc(interval.i_right, interval.row), goal, successors);
        }
    }
}
void ANYA::GenerateSuccessors(const xyLoc l, const Interval I, xyLoc const goal, std::vector<stateId> & successors)
{
    if (I.f_right - I.f_left <= EPSILON)
        return;

    bool sameRow = (l.y == I.row); // The point is on the same row as the line

    if (!sameRow) {
        GenerateSuccessors_DifferentRow(l, I, goal, successors);
    }

    else // The root is at the same row as the interval
    {
        GenerateSuccessors_SameRow(l, I, goal, successors);
    }
}

void ANYA::InitializeSearch(const xyLoc from, const xyLoc to)
{
    ResetSearch();

    std::vector<stateId> initial_states;

    // Generate the intervals to the top and bottom of the start location.
    for (int dir = 0; dir < 2; dir++) // 0 = TOP and 1 = BOT
    {
        int leftBound = from.x - corners_[from.x][from.y].clearance[dir][LEFT];
        int rightBound = from.x + corners_[from.x][from.y].clearance[dir][RIGHT];

        if (leftBound != rightBound) {
            int delta_y = dir * 2 - 1; // Gives -1 for 0 (TOP), and +1 for 1 (BOT)
            GenerateSubIntervals(from.y + delta_y, leftBound, rightBound, from, to, initial_states);
        }
    }

    // Add the interval to the left.
    if (corners_[from.x][from.y].interval_clearance[LEFT] > 0) {
        Interval I = GetLeftInterval(from.x, from.y);
        initial_states.push_back(GenerateState(from, I, to));
    }

    // Add the interval to the right.
    if (corners_[from.x][from.y].interval_clearance[RIGHT] > 0) {
        Interval I = GetRightInterval(from.x, from.y);
        initial_states.push_back(GenerateState(from, I, to));
    }

    // Put all the generated states into the open list.
    for (unsigned int i = 0; i < initial_states.size(); i++) {
        stateId s = initial_states[i];
        states_[s].g_val = 0;
        AddToOpen(s);
    }
}

cost ANYA::FindXYLocPath(xyLoc from, xyLoc to, std::vector<xyLoc> &path)
{
    path.clear();

    // This special case is necessary because the start state is automatically expanded in the InitializeSearch function.
    if (from == to) {
        path.push_back(from);
        return 0;
    }

    InitializeSearch(from, to);

    bool path_found = false;
    cost solution_cost = INFINITE_COST;
    stateId goal;

    while (!open_list_.empty()) {
        // Get the state with the minimum f-value.
        stateId curr = GetMin().id;
        PopMin();

        // Check if it is a goal state (that is, if its interval contains the goal location).
        if (states_[curr].interval.row == to.y && states_[curr].interval.f_left - EPSILON <= to.x && states_[curr].interval.f_right + EPSILON >= to.x) {
            path_found = true;
            solution_cost = states_[curr].g_val;
            goal = curr;
            break;
        }

#ifdef ANY_ANGLE_STATISTICS
        num_expansions_++;
#endif

        std::vector<stateId> successors;
        GenerateSuccessors(states_[curr].l, states_[curr].interval, to, successors);

        for (unsigned int i = 0; i < successors.size(); i++) {
            stateId succ = successors[i];

            if (states_[succ].list != CLOSED_LIST) {
                cost newGVal = states_[curr].g_val + EuclideanDistance(states_[curr].l, states_[succ].l);

                if (newGVal + EPSILON < states_[succ].g_val) {
                    states_[succ].g_val = newGVal;
                    states_[succ].parent = curr;
                    AddToOpen(succ);
                }
            }
        }
    }

    // Extract path
    if (path_found) {
#ifdef ANY_ANGLE_RUNNING_IN_HOG
        state_path_.clear(); // Used for visualizing the states on the path in HOG.
#endif
        path.push_back(to);
        stateId curr = goal;
        stateId next = curr;

        while (!(states_[curr].l == from)) {
            if (!(states_[curr].l == states_[next].l))
                path.push_back(states_[curr].l);

            curr = next;
            next = states_[curr].parent;

#ifdef ANY_ANGLE_RUNNING_IN_HOG
            state_path_.push_back(curr);
#endif
        }
#ifdef ANY_ANGLE_RUNNING_IN_HOG
        state_path_.push_back(curr);
#endif
        path.push_back(from);
        std::reverse(path.begin(), path.end());
    }

#ifdef ANYA_DEBUG
    for (int i = 0; i < state_path_.size(); i++) {
        //PrintState(state_path_[i]);
    }

    for (int i = 0; i < states_.size(); i++) {
        State* s = &states_[i];
        if (s->l.y == s->interval.row && s->l.y == to.y)
            PrintState(i);
    }
#endif

#ifdef ANY_ANGLE_ASSERTIONS
    // Make sure that the heap property is satisfied for the open list
    for (int i = 0; i < open_list_.size(); i++) {
        int child1 = (i << 1) + 1;
        int child2 = (i << 1) + 2;

        if (child2 < open_list_.size()) {
            assert(open_list_[i].f_val < open_list_[child1].f_val + 100*EPSILON);
            assert(open_list_[i].f_val < open_list_[child2].f_val + 100*EPSILON);
        }
    }
#endif

    return solution_cost;
}

cost ANYA::HeuristicDistance(const xyLoc root, const Interval interval, const xyLoc goal) const
{
    // If the root, goal, and the interval are all on the same row, there is a special case.
    if (interval.row == root.y && interval.row == goal.y) {
        // If both the root and the goal are on the left side of the interval, return the cost of the path root -> left end-point of the interval -> goal.
        if (root.x < interval.f_left && goal.x < interval.f_left)
            return EuclideanDistance(root.x, root.y, interval.f_left, interval.row) + EuclideanDistance(goal.x, goal.y, interval.f_left, interval.row);

        // If both the root and the goal are on the right side of the interval, return the cost of the path root -> right end-point of the interval -> goal.
        else if (root.x > interval.f_right && goal.x > interval.f_right)
            return EuclideanDistance(root.x, root.y, interval.f_right, interval.row) + EuclideanDistance(goal.x, goal.y, interval.f_right, interval.row);

        // Otherwise (the root and the goal are both inside the interval, or on the opposite sides of the interval), return the cost of the path root -> goal.
        else
            return EuclideanDistance(root, goal);
    }

    int goal_y = goal.y;

    // If the root and the goal are both on the same side of the interval, take the mirror image of the goal wrt the interval
    int delta_root = (int) root.y - (int) interval.row; // Normalize the y-coordinates so that the interval's y coordinate corresponds to 0.
    int delta_goal = (int) goal.y - (int) interval.row;

    if (delta_root * delta_goal > 0) {
        goal_y = (int) interval.row - delta_goal;
    }

    // Now, the root and the goal (or its mirror image) are on different sides of the interval's row.
    // Find the x-coordinate of the intersection between the line (root, goal) and the interval's row.
    double x = GetIntersectingX(root.x, root.y, goal.x, goal_y, interval.row);

    // If x does not lie in the interval, set it to one of its endpoints.
    if (x < interval.f_left)
        x = interval.f_left;
    if (x > interval.f_right)
        x = interval.f_right;

    return EuclideanDistance(root.x, root.y, x, interval.row) + EuclideanDistance(goal.x, goal.y, x, interval.row);
}

void ANYA::ResetSearch()
{
    states_.clear();
    anya_index_table_.Reset();
    open_list_.clear();
}

// Heap operations.
void ANYA::AddToOpen(const stateId id)
{
    // If it is already in the open list, update its location.
    if (states_[id].list == OPEN_LIST) {
        int i = states_[id].heap_index;
        open_list_[i].g_val = states_[id].g_val;
        open_list_[i].f_val = states_[id].g_val + states_[id].h_val;
        PercolateUp(i);
    }

    // Otherwise, add it to the open list.
    else {
        states_[id].list = OPEN_LIST;
        states_[id].heap_index = open_list_.size();
        open_list_.push_back(HeapElement(id, states_[id].g_val, states_[id].g_val + states_[id].h_val));
        PercolateUp(open_list_.size() - 1);
    }
}
void ANYA::PopMin()
{
    states_[open_list_[0].id].list = CLOSED_LIST;
    open_list_[0] = open_list_.back();
    states_[open_list_[0].id].heap_index = 0;
    open_list_.pop_back();
    PercolateDown(0);
}
void ANYA::PercolateUp(int index)
{
    HeapElement elem = open_list_[index];
    int parent = (index - 1) >> 1;

    while (index > 0 && elem < open_list_[parent]) {
        open_list_[index] = open_list_[parent];
        states_[open_list_[index].id].heap_index = index;

        index = parent;
        parent = (index - 1) >> 1;

#ifdef ANY_ANGLE_STATISTICS
        num_percolations_++;
#endif
    }

    open_list_[index] = elem;
    states_[elem.id].heap_index = index;
}
void ANYA::PercolateDown(int index)
{
    HeapElement elem = open_list_[index];
    int maxSize = open_list_.size();

    while (true) {
        int child1 = (index << 1) + 1;
        if (child1 >= maxSize)
            break;

        int child2 = child1 + 1;

        // If the first child has smaller key.
        if (child2 == maxSize || open_list_[child1] < open_list_[child2]) {
            if (open_list_[child1] < elem) {
                open_list_[index] = open_list_[child1];
                states_[open_list_[index].id].heap_index = index;
                index = child1;

#ifdef ANY_ANGLE_STATISTICS
                num_percolations_++;
#endif
            }
            else
                break;
        }

        else if (open_list_[child2] < elem) {
            open_list_[index] = open_list_[child2];
            states_[open_list_[index].id].heap_index = index;
            index = child2;

#ifdef ANY_ANGLE_STATISTICS
            num_percolations_++;
#endif
        }
        else
            break;
    }

    open_list_[index] = elem;
    states_[elem.id].heap_index = index;
}

#ifdef ANY_ANGLE_RUNNING_IN_HOG
void ANYA::ProcessDebugLoc(xyLoc l)
{
    printf("\nDebug point: (%d, %d)\n", l.x, l.y);

    printf("Left interval clearance = %u (top blocked: %u, bot blocked: %u)\nLeft clearance (Top = %u, Bot = %u)\n",
            corners_[l.x][l.y].interval_clearance[LEFT], corners_[l.x][l.y].interval_blocked[TOP][LEFT], corners_[l.x][l.y].interval_blocked[BOT][LEFT],
            corners_[l.x][l.y].clearance[TOP][LEFT], corners_[l.x][l.y].clearance[BOT][LEFT]);

    printf("Right interval clearance = %u (top blocked: %u, bot blocked: %u)\nRight clearance (Top = %u, Bot = %u)\n",
            corners_[l.x][l.y].interval_clearance[RIGHT], corners_[l.x][l.y].interval_blocked[TOP][RIGHT], corners_[l.x][l.y].interval_blocked[BOT][RIGHT],
            corners_[l.x][l.y].clearance[TOP][RIGHT], corners_[l.x][l.y].clearance[BOT][RIGHT]);

    if (corners_[l.x][l.y].is_transition_point)
        printf("Transition point!\n");

    if (corners_[l.x][l.y].is_convex_corner)
        printf("Convex corner!\n");
}
void ANYA::VisualizeAlgorithm(const MapEnvironment *env)
{
    bool displayGeneratedStates = 0;

#ifdef ANYA_DEBUG
    displayGeneratedStates = 1;
#endif

//    env->SetColor(1, 0, 0);
//    DrawPoint(env, debug_loc_);

    if (displayGeneratedStates) {

        for (unsigned int j = 0; j < state_path_.size(); j++) {
            int i = state_path_[j];

            xyLoc l = states_[i].l;
            Interval I = states_[i].interval;
            xyLocCont left = xyLocCont(I.f_left, I.row);
            xyLocCont right = xyLocCont(I.f_right, I.row);
            xyLocCont lCont = xyLocCont(l.x, l.y);

            if (left.x < 1) // Otherwise, there is an issue with the display that forces a crash
                left.x = 1;

            if (right.x < 1)
                right.x = 1;

            if (l.x == 0)
                continue;

            if (I.blocked[BOT] && I.blocked[TOP])
                env->SetColor(1, 0, 0);

            if (I.blocked[BOT] && !I.blocked[TOP])
                env->SetColor(0, 1, 0);

            if (!I.blocked[BOT] && I.blocked[TOP])
                env->SetColor(0, 0, 1);

            if (!I.blocked[BOT] && !I.blocked[TOP])
                env->SetColor(0, 0, 0);

            //if (I.blocked[BOT] && I.blocked[TOP])
            {
                DrawLine(env, left, right);
                DrawLine(env, left, lCont);
                DrawLine(env, right, lCont);

                env->SetColor(1, 0, 0);
                DrawPoint(env, l);
            }

        }
    }
}
#endif

