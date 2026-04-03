#include "field_d_star_planner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace {

constexpr double kSqrt2 = 1.4142135623730951;

bool keyLess(const field_d_star::Key& a, const field_d_star::Key& b)
{
    if (a.k1 < b.k1 - field_d_star::kEps) {
        return true;
    }
    if (a.k1 > b.k1 + field_d_star::kEps) {
        return false;
    }
    return a.k2 < b.k2 - field_d_star::kEps;
}

bool keyEqual(const field_d_star::Key& a, const field_d_star::Key& b)
{
    return std::abs(a.k1 - b.k1) <= field_d_star::kEps &&
           std::abs(a.k2 - b.k2) <= field_d_star::kEps;
}

bool samePoint(const field_d_star::Point2D& a, const field_d_star::Point2D& b)
{
    return std::abs(a.x - b.x) <= field_d_star::kTolerance &&
           std::abs(a.y - b.y) <= field_d_star::kTolerance;
}

}  // namespace

namespace field_d_star {

const double kInf = std::numeric_limits<double>::infinity();

bool FieldDStarPlanner::PQElement::operator>(const PQElement& other) const
{
    if (key.k1 == other.key.k1) {
        if (key.k2 == other.key.k2) {
            return version > other.version;
        }
        return key.k2 > other.key.k2;
    }
    return key.k1 > other.key.k1;
}

FieldDStarPlanner::FieldDStarPlanner(
    const CostGrid& grid,
    Point2D start,
    GridCoord goal,
    double penalty_weight)
{
    // Configuration is split across setters so the same validation rules can be
    // reused by callers that build the planner step by step
    setGrid(grid);
    setPenaltyWeight(penalty_weight);
    setStart(start);
    setGoal(goal);
}

void FieldDStarPlanner::validateGrid(const CostGrid& grid) const
{
    if (grid.empty()) {
        throw std::invalid_argument("FieldDStarPlanner grid must not be empty");
    }
    if (grid.front().empty()) {
        throw std::invalid_argument("FieldDStarPlanner grid columns must not be empty");
    }

    const std::size_t expected_height = grid.front().size();
    for (const auto& column : grid) {
        if (column.size() != expected_height) {
            throw std::invalid_argument("FieldDStarPlanner grid must be rectangular and indexed as grid[x][y]");
        }
    }
}

void FieldDStarPlanner::validateGoal(GridCoord goal) const
{
    if (goal.x < 0 || goal.x > width_ || goal.y < 0 || goal.y > height_) {
        throw std::out_of_range("FieldDStarPlanner goal must be a valid corner node");
    }
}

void FieldDStarPlanner::validateStart(Point2D start) const
{
    if (start.x < 0.0 || start.x > static_cast<double>(width_) ||
        start.y < 0.0 || start.y > static_cast<double>(height_)) {
        throw std::out_of_range("FieldDStarPlanner start must lie inside the grid footprint");
    }
}

void FieldDStarPlanner::ensureConfigured() const
{
    if (!configured_) {
        throw std::logic_error("FieldDStarPlanner requires grid, start and goal before planning");
    }
}

void FieldDStarPlanner::setGrid(const CostGrid& grid)
{
    validateGrid(grid);
    grid_ = grid;
    width_ = static_cast<int>(grid_.size());
    height_ = static_cast<int>(grid_.front().size());

    // Field D* searches over cell corners, so a WxH cell grid induces
    // (W + 1) * (H + 1) corner nodes
    nodes_.assign(static_cast<std::size_t>((width_ + 1) * (height_ + 1)), {});

    for (int y = 0; y <= height_; ++y) {
        for (int x = 0; x <= width_; ++x) {
            Node& node = nodes_[getIndex(x, y)];
            node.x = x;
            node.y = y;
        }
    }

    goal_index_ = -1;
    configured_ = false;
    initialized_ = false;
    is_optimal_ = false;
    path_.clear();
    open_entries_.clear();
    std::priority_queue<PQElement, std::vector<PQElement>, std::greater<PQElement>> empty_queue;
    std::swap(open_list_, empty_queue);
}

void FieldDStarPlanner::setStart(Point2D start)
{
    if (width_ == 0 || height_ == 0) {
        throw std::logic_error("Set the grid before setting the start");
    }
    validateStart(start);
    start_pos_ = start;
    configured_ = (goal_index_ >= 0);
    is_optimal_ = false;
    path_.clear();
}

void FieldDStarPlanner::moveStart(Point2D new_start)
{
    ensureConfigured();
    validateStart(new_start);

    // Incremental replanning keeps the existing search tree and only updates the
    // D* Lite key modifier by the distance traveled since the previous start
    k_m_ += std::hypot(new_start.x - start_pos_.x, new_start.y - start_pos_.y);
    start_pos_ = new_start;
    path_.clear();
    if (initialized_) {
        extractPath();
    }
}

void FieldDStarPlanner::setGoal(GridCoord goal)
{
    if (width_ == 0 || height_ == 0) {
        throw std::logic_error("Set the grid before setting the goal");
    }
    validateGoal(goal);
    goal_coord_ = goal;
    goal_index_ = getIndex(goal.x, goal.y);
    configured_ = true;
    initialized_ = false;
    is_optimal_ = false;
    path_.clear();
}

void FieldDStarPlanner::setPenaltyWeight(double penalty_weight)
{
    if (penalty_weight <= 0.0) {
        throw std::invalid_argument("Penalty weight must be positive");
    }
    penalty_weight_ = penalty_weight;
    is_optimal_ = false;
    path_.clear();
}

void FieldDStarPlanner::updateCellCost(int x, int y, int value)
{
    ensureConfigured();
    if (x < 0 || x >= width_ || y < 0 || y >= height_) {
        throw std::out_of_range("Cell update is outside the grid");
    }
    if (value != kUndefinedValue && (value < 0 || value > kObstacleValue)) {
        throw std::invalid_argument("Cell cost must be -1 or in the range [0, 254]");
    }
    if (grid_[x][y] == value) {
        return;
    }

    grid_[x][y] = value;
    is_optimal_ = false;
    path_.clear();

    if (!initialized_) {
        return;
    }

    // Only the four corners touching this cell can have their local interpolation geometry changed directly by a cell-cost update
    for (int nx : {x, x + 1}) {
        for (int ny : {y, y + 1}) {
            updateNode(getIndex(nx, ny));
        }
    }
}

void FieldDStarPlanner::initialize()
{
    ensureConfigured();

    std::priority_queue<PQElement, std::vector<PQElement>, std::greater<PQElement>> empty_queue;
    std::swap(open_list_, empty_queue);
    open_entries_.clear();
    path_.clear();

    k_m_ = 0.0;
    is_optimal_ = false;
    initialized_ = true;

    for (Node& node : nodes_) {
        node.g = kInf;
        node.rhs = kInf;
        node.key = {};
        node.open_version = 0;
    }

    // The search is rooted at the goal
    nodes_[goal_index_].rhs = 0.0;
    pushOpen(goal_index_);
}

double FieldDStarPlanner::getCellCost(int x, int y) const
{
    if (x < 0 || x >= width_ || y < 0 || y >= height_) {
        return kInf;
    }

    const int value = grid_[x][y];
    // obstacle or unknown cells are treated as non-traversable, otherwise integer costs are scaled into cost per meter
    if (value >= kObstacleValue || value == kUndefinedValue) {
        return kInf;
    }
    return 1.0 + static_cast<double>(value) / penalty_weight_;
}

Key FieldDStarPlanner::calculateKeyValue(int index) const
{
    const double best = std::min(nodes_[index].g, nodes_[index].rhs);
    const double h = std::hypot(
        static_cast<double>(nodes_[index].x) - start_pos_.x,
        static_cast<double>(nodes_[index].y) - start_pos_.y);
    return {best + h + k_m_, best};
}

bool FieldDStarPlanner::isConsistent(int index) const
{
    return (std::abs(nodes_[index].g - nodes_[index].rhs) <= kEps) ||
           (!std::isfinite(nodes_[index].g) && !std::isfinite(nodes_[index].rhs));
}

void FieldDStarPlanner::pushOpen(int index)
{
    // The heap uses lazy deletion: every refresh increments the node version and
    // stale entries are discarded when they reach the top
    nodes_[index].key = calculateKeyValue(index);
    nodes_[index].open_version += 1;
    open_entries_[index] = {nodes_[index].key, nodes_[index].open_version};
    open_list_.push({nodes_[index].key, nodes_[index].open_version, index});
}

void FieldDStarPlanner::removeFromOpen(int index)
{
    open_entries_.erase(index);
}

Key FieldDStarPlanner::topKey()
{
    while (!open_list_.empty()) {
        const PQElement& top = open_list_.top();
        const auto active = open_entries_.find(top.node_index);
        if (active != open_entries_.end() &&
            active->second.version == top.version &&
            keyEqual(active->second.key, top.key)) {
            return top.key;
        }
        open_list_.pop();
    }
    return {};
}

bool FieldDStarPlanner::popValid(int& out_index, Key& out_key)
{
    while (!open_list_.empty()) {
        PQElement top = open_list_.top();
        open_list_.pop();

        const auto active = open_entries_.find(top.node_index);
        if (active != open_entries_.end() &&
            active->second.version == top.version &&
            keyEqual(active->second.key, top.key)) {
            out_index = top.node_index;
            out_key = top.key;
            open_entries_.erase(active);
            return true;
        }
    }
    return false;
}

void FieldDStarPlanner::computeShortestPath()
{
    ensureConfigured();
    if (!initialized_) {
        initialize();
    }

    const int start_ref = closestNodeToStart();
    int iterations = 0;

    while ((keyLess(topKey(), calculateKeyValue(start_ref)) || !isConsistent(start_ref)) &&
           iterations < kMaxIterations) {
        int u = -1;
        Key old_key;
        if (!popValid(u, old_key)) {
            break;
        }

        const Key current_key = calculateKeyValue(u);
        if (keyLess(old_key, current_key)) {
            // The stale key was too optimistic -> Reinsert with the refreshed key
            pushOpen(u);
        } else if (nodes_[u].g > nodes_[u].rhs) {
            // Overconsistent node: accept the improved rhs estimate as the new g
            nodes_[u].g = nodes_[u].rhs;
            for (const int predecessor : predecessors(u)) {
                updateNode(predecessor);
            }
        } else {
            // Underconsistent node: invalidate its g value and let the neighborhood rebuild the best successor relation
            nodes_[u].g = kInf;
            updateNode(u);
            for (const int predecessor : predecessors(u)) {
                updateNode(predecessor);
            }
        }

        ++iterations;
    }

    if (iterations >= kMaxIterations) {
        throw std::runtime_error("Field D* exceeded the iteration limit");
    }

    is_optimal_ = true;
}

void FieldDStarPlanner::updateNode(int index)
{
    if (index != goal_index_) {
        nodes_[index].rhs = computeMinRHS(index);
    }

    if (isConsistent(index)) {
        removeFromOpen(index);
    } else {
        pushOpen(index);
    }
}

double FieldDStarPlanner::computeMinRHS(int index)
{
    if (index == goal_index_) {
        return 0.0;
    }
    return bestSuccessorChoice(index).cost;
}

FieldDStarPlanner::InterpChoice FieldDStarPlanner::computeCostChoice(
    double c,
    double b,
    double g1,
    double g2) const
{
    // This is the ComputeCost routine for one consecutive neighbor pair
    // It evaluates whether the optimal path leaves the current node by going
    // straight to the perpendicular node, diagonally to the corner, or by
    // cutting across the cell boundary at an interpolated point
    if (std::isnan(c) || std::isnan(b) || std::isnan(g1) || std::isnan(g2)) {
        return {kInf, 0.0, InterpChoice::Mode::kNan, false};
    }
    if (std::min(c, b) == kInf) {
        return {kInf, 0.0, InterpChoice::Mode::kBlocked, false};
    }
    if (g1 <= g2 + kEps) {
        return {std::min(c, b) + g1, 0.0, InterpChoice::Mode::kStraight, true};
    }

    const double f = g1 - g2;
    if (f <= b + kEps) {
        if (c <= f + kEps) {
            return {c * kSqrt2 + g2, 1.0, InterpChoice::Mode::kDiag, true};
        }

        const double val = c * c - f * f;
        double y = 1.0;
        if (val > kEps) {
            y = f / std::sqrt(val);
        }
        y = std::min(y, 1.0);
        return {
            c * std::sqrt(1.0 + y * y) + f * (1.0 - y) + g2,
            y,
            InterpChoice::Mode::kCut,
            true};
    }

    if (c <= b + kEps) {
        return {c * kSqrt2 + g2, 1.0, InterpChoice::Mode::kDiag, true};
    }

    const double val = c * c - b * b;
    double x_sub = 1.0;
    if (val > kEps) {
        x_sub = b / std::sqrt(val);
    }
    const double x = 1.0 - std::min(x_sub, 1.0);
    return {
        c * std::sqrt(1.0 + (1.0 - x) * (1.0 - x)) + b * x + g2,
        x,
        InterpChoice::Mode::kBoundary,
        true};
}

FieldDStarPlanner::SuccessorChoice FieldDStarPlanner::bestSuccessorChoice(int index) const
{
    if (index == goal_index_) {
        return {0.0, -1, 0.0, InterpChoice::Mode::kGoal, index, index, {}, {}, true};
    }

    // Possible neighbors
    static constexpr int ds1x[8] = {1, 0, -1, 0, -1, 0, 1, 0};
    static constexpr int ds1y[8] = {0, 1, 0, 1, 0, -1, 0, -1};
    static constexpr int ds2x[8] = {1, 1, -1, -1, -1, -1, 1, 1};
    static constexpr int ds2y[8] = {1, 1, 1, 1, -1, -1, -1, -1};
    static constexpr int dcx[8] = {0, 0, -1, -1, -1, -1, 0, 0};
    static constexpr int dcy[8] = {0, 0, 0, 0, -1, -1, -1, -1};
    static constexpr int dbx[8] = {0, -1, -1, 0, -1, 0, 0, -1};
    static constexpr int dby[8] = {-1, 0, -1, 0, 0, -1, 0, -1};

    const int sx = nodes_[index].x;
    const int sy = nodes_[index].y;
    SuccessorChoice best;

    // The eight pairs correspond to consecutive neighbors around the node
    // Each pair defines one candidate exit edge and the two cells needed by te analytical cost formula
    for (int pair_index = 0; pair_index < 8; ++pair_index) {
        const int s1x = sx + ds1x[pair_index];
        const int s1y = sy + ds1y[pair_index];
        const int s2x = sx + ds2x[pair_index];
        const int s2y = sy + ds2y[pair_index];

        if (s1x < 0 || s1x > width_ || s1y < 0 || s1y > height_ ||
            s2x < 0 || s2x > width_ || s2y < 0 || s2y > height_) {
            continue;
        }

        const int s1_idx = getIndex(s1x, s1y);
        const int s2_idx = getIndex(s2x, s2y);
        const double g1 = nodes_[s1_idx].g;
        const double g2 = nodes_[s2_idx].g;
        const double c_cost = getCellCost(sx + dcx[pair_index], sy + dcy[pair_index]);
        const double b_cost = getCellCost(sx + dbx[pair_index], sy + dby[pair_index]);

        const InterpChoice choice = computeCostChoice(c_cost, b_cost, g1, g2);
        if (!choice.valid) {
            continue;
        }

        if (!best.valid || choice.cost + kEps < best.cost) {
            best.cost = choice.cost;
            best.pair_index = pair_index;
            best.edge_param = choice.edge_param;
            best.mode = choice.mode;
            best.s1_idx = s1_idx;
            best.s2_idx = s2_idx;
            best.c_cell = {sx + dcx[pair_index], sy + dcy[pair_index]};
            best.b_cell = {sx + dbx[pair_index], sy + dby[pair_index]};
            best.valid = true;
        }
    }

    if (!best.valid) {
        best.mode = InterpChoice::Mode::kUnreachable;
    }
    return best;
}

std::vector<int> FieldDStarPlanner::getNeighbors(int index) const
{
    std::vector<int> neighbors;
    const int cx = nodes_[index].x;
    const int cy = nodes_[index].y;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int nx = cx + dx;
            const int ny = cy + dy;
            if (nx >= 0 && nx <= width_ && ny >= 0 && ny <= height_) {
                neighbors.push_back(getIndex(nx, ny));
            }
        }
    }
    return neighbors;
}

std::vector<int> FieldDStarPlanner::predecessors(int index) const
{
    return getNeighbors(index);
}

Point2D FieldDStarPlanner::pointOnEdge(int n1_idx, int n2_idx, double t) const
{
    const Node& n1 = nodes_[n1_idx];
    const Node& n2 = nodes_[n2_idx];
    return {
        static_cast<double>(n1.x) + static_cast<double>(n2.x - n1.x) * t,
        static_cast<double>(n1.y) + static_cast<double>(n2.y - n1.y) * t};
}

double FieldDStarPlanner::interpolateEdgeCost(int n1_idx, int n2_idx, double t) const
{
    const double g1 = nodes_[n1_idx].g;
    const double g2 = nodes_[n2_idx].g;

    if (!std::isfinite(g1) && !std::isfinite(g2)) {
        return kInf;
    }
    if (!std::isfinite(g1)) {
        return g2;
    }
    if (!std::isfinite(g2)) {
        return g1;
    }
    return g1 + (g2 - g1) * t;
}

std::vector<GridCoord> FieldDStarPlanner::cellsContainingPoint(const Point2D& point) const
{
    std::vector<GridCoord> cells;
    if (width_ == 0 || height_ == 0) {
        return cells;
    }

    // A fractional point can lie on a cell edge or corner, so multiple incident
    // cells may need to be considered during extraction
    const int min_cx = std::max(0, static_cast<int>(std::floor(point.x - kTolerance)));
    const int max_cx = std::min(width_ - 1, static_cast<int>(std::floor(point.x + kTolerance)));
    const int min_cy = std::max(0, static_cast<int>(std::floor(point.y - kTolerance)));
    const int max_cy = std::min(height_ - 1, static_cast<int>(std::floor(point.y + kTolerance)));

    for (int cx = min_cx; cx <= max_cx; ++cx) {
        for (int cy = min_cy; cy <= max_cy; ++cy) {
            cells.push_back({cx, cy});
        }
    }
    return cells;
}

bool FieldDStarPlanner::pointInCell(const Point2D& point, const GridCoord& cell) const
{
    return cell.x - kTolerance <= point.x && point.x <= cell.x + 1 + kTolerance &&
           cell.y - kTolerance <= point.y && point.y <= cell.y + 1 + kTolerance;
}

std::array<int, 4> FieldDStarPlanner::cellCorners(const GridCoord& cell) const
{
    return {
        getIndex(cell.x, cell.y),
        getIndex(cell.x + 1, cell.y),
        getIndex(cell.x + 1, cell.y + 1),
        getIndex(cell.x, cell.y + 1)};
}

Point2D FieldDStarPlanner::choicePoint(const SuccessorChoice& choice) const
{
    return pointOnEdge(choice.s1_idx, choice.s2_idx, choice.edge_param);
}

double FieldDStarPlanner::choicePointCost(const SuccessorChoice& choice) const
{
    return interpolateEdgeCost(choice.s1_idx, choice.s2_idx, choice.edge_param);
}

double FieldDStarPlanner::directGoalCost(const Point2D& point) const
{
    const Point2D goal_point{static_cast<double>(goal_coord_.x), static_cast<double>(goal_coord_.y)};

    // If the current point and the goal share a traversable cell, the direct
    // segment inside that cell can be evaluated exactly
    for (const GridCoord& cell : cellsContainingPoint(point)) {
        if (!pointInCell(goal_point, cell)) {
            continue;
        }

        const double cell_cost = getCellCost(cell.x, cell.y);
        if (std::isfinite(cell_cost)) {
            return cell_cost * std::hypot(point.x - goal_point.x, point.y - goal_point.y);
        }
    }
    return kInf;
}

std::pair<int, int> FieldDStarPlanner::normalizeEdge(int n1_idx, int n2_idx) const
{
    if (n1_idx <= n2_idx) {
        return {n1_idx, n2_idx};
    }
    return {n2_idx, n1_idx};
}

std::vector<FieldDStarPlanner::BoundaryCandidate> FieldDStarPlanner::localSuccessorCandidates(
    const Point2D& point) const
{
    std::vector<BoundaryCandidate> candidates;

    // Instead of sampling arbitrary edge points, reuse the interpolation choices
    // already computed for the surrounding corner nodes
    for (const GridCoord& cell : cellsContainingPoint(point)) {
        const double cell_cost = getCellCost(cell.x, cell.y);
        if (!std::isfinite(cell_cost)) {
            continue;
        }

        for (const int corner_idx : cellCorners(cell)) {
            const SuccessorChoice choice = bestSuccessorChoice(corner_idx);
            if (!choice.valid || !std::isfinite(choice.cost) || choice.pair_index < 0) {
                continue;
            }

            const Point2D candidate_point = choicePoint(choice);
            if (!pointInCell(candidate_point, cell)) {
                continue;
            }
            if (!((choice.c_cell.x == cell.x && choice.c_cell.y == cell.y) ||
                  (choice.b_cell.x == cell.x && choice.b_cell.y == cell.y))) {
                continue;
            }

            // The segment from the current point to the chosen boundary point
            // lies inside cell, so one cell traversal cost is sufficient
            const double segment_cost =
                cell_cost * std::hypot(candidate_point.x - point.x, candidate_point.y - point.y);
            const double point_cost_value = choicePointCost(choice);
            if (!std::isfinite(point_cost_value)) {
                continue;
            }

            bool duplicate = false;
            for (const BoundaryCandidate& existing : candidates) {
                if (existing.cell.x == cell.x &&
                    existing.cell.y == cell.y &&
                    samePoint(existing.point, candidate_point)) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate) {
                continue;
            }

            BoundaryCandidate candidate;
            candidate.point = candidate_point;
            candidate.point_cost = point_cost_value;
            candidate.total_cost = segment_cost + point_cost_value;
            candidate.cell = cell;
            candidate.edge_n1_idx = choice.s1_idx;
            candidate.edge_n2_idx = choice.s2_idx;
            candidate.valid = true;
            candidates.push_back(candidate);
        }
    }

    const double goal_cost = directGoalCost(point);
    if (std::isfinite(goal_cost)) {
        // The goal itself is a valid candidate whenever it is locally visible
        BoundaryCandidate candidate;
        candidate.point = {static_cast<double>(goal_coord_.x), static_cast<double>(goal_coord_.y)};
        candidate.point_cost = 0.0;
        candidate.total_cost = goal_cost;
        candidate.cell = {-1, -1};
        candidate.edge_n1_idx = goal_index_;
        candidate.edge_n2_idx = goal_index_;
        candidate.valid = true;
        candidates.push_back(candidate);
    }

    return candidates;
}

double FieldDStarPlanner::pointCost(
    const Point2D& point,
    const std::vector<std::pair<int, int>>& forbidden_edges) const
{
    double best = kInf;

    // This helper is used by one-step lookahead and can optionally forbid the
    // edge that led to the current point to avoid immediate backtracking
    for (const BoundaryCandidate& candidate : localSuccessorCandidates(point)) {
        const auto edge = normalizeEdge(candidate.edge_n1_idx, candidate.edge_n2_idx);
        if (std::find(forbidden_edges.begin(), forbidden_edges.end(), edge) != forbidden_edges.end()) {
            continue;
        }
        best = std::min(best, candidate.total_cost);
    }

    return best;
}

std::vector<FieldDStarPlanner::BoundaryCandidate> FieldDStarPlanner::nextBoundaryCandidates(
    const Point2D& point) const
{
    std::vector<BoundaryCandidate> candidates;

    // One-step lookahead during extraction to reduce interpolation error in hard local configurations
    for (const BoundaryCandidate& candidate : localSuccessorCandidates(point)) {
        if (std::hypot(candidate.point.x - point.x, candidate.point.y - point.y) <= kTolerance) {
            continue;
        }

        const std::vector<std::pair<int, int>> forbidden = {
            normalizeEdge(candidate.edge_n1_idx, candidate.edge_n2_idx)};
        const double lookahead = pointCost(candidate.point, forbidden);

        BoundaryCandidate scored = candidate;
        if (std::isfinite(lookahead) && candidate.cell.x >= 0 && candidate.cell.y >= 0) {
            const double cell_cost = getCellCost(candidate.cell.x, candidate.cell.y);
            scored.total_cost =
                std::hypot(candidate.point.x - point.x, candidate.point.y - point.y) * cell_cost + lookahead;
        }
        candidates.push_back(scored);
    }

    return candidates;
}

int FieldDStarPlanner::closestNodeToStart() const
{
    const int nx = std::min(std::max(static_cast<int>(std::round(start_pos_.x)), 0), width_);
    const int ny = std::min(std::max(static_cast<int>(std::round(start_pos_.y)), 0), height_);
    return getIndex(nx, ny);
}

const std::vector<Point2D>& FieldDStarPlanner::extractPath()
{
    ensureConfigured();
    path_.clear();
    path_.push_back(start_pos_);

    Point2D current = start_pos_;
    const Point2D goal_point{static_cast<double>(goal_coord_.x), static_cast<double>(goal_coord_.y)};
    std::vector<std::pair<std::int64_t, std::int64_t>> visited;

    // Walk through successive boundary waypoints until the goal is reached or no numerically stable progress can be made
    for (int step = 0; step < 300; ++step) {
        if (std::hypot(current.x - goal_point.x, current.y - goal_point.y) < kTolerance) {
            break;
        }

        const std::vector<BoundaryCandidate> candidates = nextBoundaryCandidates(current);
        if (candidates.empty()) {
            break;
        }

        // Prefer the smallest total cost, breaking ties by the interpolated point cost
        const auto best_it = std::min_element(
            candidates.begin(),
            candidates.end(),
            [](const BoundaryCandidate& lhs, const BoundaryCandidate& rhs) {
                if (lhs.total_cost < rhs.total_cost - kEps) {
                    return true;
                }
                if (lhs.total_cost > rhs.total_cost + kEps) {
                    return false;
                }
                return lhs.point_cost < rhs.point_cost;
            });

        const Point2D next_point = best_it->point;
        // Loop protection for degenerate local choices or repeated waypoints
        const std::pair<std::int64_t, std::int64_t> rounded = {
            static_cast<std::int64_t>(std::llround(next_point.x * 100.0)),
            static_cast<std::int64_t>(std::llround(next_point.y * 100.0))};
        if (std::find(visited.begin(), visited.end(), rounded) != visited.end()) {
            break;
        }

        visited.push_back(rounded);
        current = next_point;
        path_.push_back(current);
    }

    if (std::hypot(current.x - goal_point.x, current.y - goal_point.y) < kTolerance) {
        if (!samePoint(current, goal_point)) {
            path_.push_back(goal_point);
        }
    } else {
        const double goal_distance = std::hypot(current.x - goal_point.x, current.y - goal_point.y);
        if (goal_distance <= kSqrt2 + kTolerance && !samePoint(current, goal_point)) {
            path_.push_back(goal_point);
        }
    }

    return path_;
}

std::vector<Point2D> FieldDStarPlanner::plan()
{
    initialize();
    computeShortestPath();
    return extractPath();
}

bool FieldDStarPlanner::advanceStartToNextWaypoint()
{
    if (path_.size() < 2) {
        return false;
    }
    moveStart(path_[1]);
    return true;
}

}  // namespace field_d_star
