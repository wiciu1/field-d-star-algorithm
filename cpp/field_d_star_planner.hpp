#ifndef FIELD_D_STAR_PLANNER_HPP_
#define FIELD_D_STAR_PLANNER_HPP_

#include <array>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace field_d_star {

/// Numerical tolerances and default grid-cost constants used by the planner
inline constexpr double kTolerance = 0.01;
inline constexpr double kEps = 1e-9;
inline constexpr int kDefaultPenalty = 20;
inline constexpr int kObstacleValue = 254;
inline constexpr int kUndefinedValue = -1;

/// Continuous point expressed in grid coordinates
struct Point2D {
    /// Horizontal coordinate in the grid reference frame
    double x = 0.0;
    /// Vertical coordinate in the grid reference frame
    double y = 0.0;
};

/// Integer grid coordinate
/// For goals this refers to a corner node
struct GridCoord {
    int x = -1;
    int y = -1;
};

/// D* Lite priority key
struct Key {
    /// Primary key component: min(g, rhs) + heuristic + key modifier
    double k1 = std::numeric_limits<double>::infinity();
    /// Secondary key component: min(g, rhs)
    double k2 = std::numeric_limits<double>::infinity();
};

/**
 * Plain C++ implementation of the 2D Field D* algorithm
 *
 * The planner operates on a uniform cost grid indexed as grid[x][y]
 * Traversal costs are attached to cells, while search states live on the
 * corners of those cells. The extracted path is a polyline in continuous
 * grid coordinates, not a sequence of grid-node hops
 */
class FieldDStarPlanner {
public:
    using CostGrid = std::vector<std::vector<int>>;

    FieldDStarPlanner() = default;

    /// Constructs and configures the planner with a grid, fractional start and corner goal
    FieldDStarPlanner(const CostGrid& grid, Point2D start, GridCoord goal, double penalty_weight = kDefaultPenalty);

    /**
     * Replaces the entire cost grid and resets internal search state
     *
     * The grid is expected to be rectangular and indexed as grid[x][y]
     * Each entry stores an integer traversal descriptor:
     * - 0..253: for traversable cells,
     * - 254: for obstacles,
     * - -1: for unknown / undefined.
     */
    void setGrid(const CostGrid& grid);

    /**
     * Sets the current agent position in continuous grid coordinates
     * The start may lie anywhere inside the grid footprint, not only on a node
     */
    void setStart(Point2D start);

    /**
     * Moves the current start and updates the key modifier k_m
     */
    void moveStart(Point2D new_start);

    /// Sets the goal corner node.
    /// The goal must be one of the grid corners
    void setGoal(GridCoord goal);

    /**
     * Sets the divisor used to map integer cell values to traversal cost
     */
    void setPenaltyWeight(double penalty_weight);

    /**
     * Updates one cell cost and marks the incident corner nodes for replanning
     *
     * If the planner has already been initialized, the four corners of the
     * modified cell are refreshed so that a later computeShortestPath()
     * incrementally repairs the solution
     */
    void updateCellCost(int x, int y, int value);

    /**
     * Initializes a fresh search from the current start to the current goal
     *
     * All node g/rhs values are reset to infinity, the goal rhs is set to
     * zero and the goal node is inserted into OPEN
     */
    void initialize();

    /**
     * Repairs inconsistent nodes until the start reference becomes locally consistent
     */
    void computeShortestPath();

    /**
     * Extracts a continuous path using boundary waypoints and one-step lookahead
     *
     * The resulting polyline is stored in path() and returned by reference
     * It starts at the current start point and ends at the goal if a path can
     * be completed or approximated robustly
     */
    const std::vector<Point2D>& extractPath();

    /// Convenience wrapper that runs initialize, shortest-path propagation and extraction.
    std::vector<Point2D> plan();

    /**
     * Advances the start to the next extracted waypoint, if one exists
     */
    bool advanceStartToNextWaypoint();

    const std::vector<Point2D>& path() const { return path_; }
    Point2D start() const { return start_pos_; }
    GridCoord goal() const { return goal_coord_; }
    int width() const { return width_; }
    int height() const { return height_; }
    bool isOptimal() const { return is_optimal_; }

private:
    /// Corner node
    struct Node {
        int x = 0;
        int y = 0;
        /// Current best-known path cost from this node to the goal
        double g = std::numeric_limits<double>::infinity();
        /// One-step lookahead estimate used by D* Lite consistency checks
        double rhs = std::numeric_limits<double>::infinity();
        /// Last computed priority key for insertion into OPEN
        Key key;
        /// Monotonic version used for lazy deletion of stale heap entries
        int open_version = 0;
    };

    /// Result of evaluating one consecutive-neighbor edge pair around a node
    struct InterpChoice {
        enum class Mode {
            kBlocked,
            kStraight,
            kDiag,
            kCut,
            kBoundary,
            kGoal,
            kUnreachable,
            kNan
        };

        /// Best total cost obtained for this edge pair
        double cost = std::numeric_limits<double>::infinity();
        /// Interpolation parameter on the exit edge in the range [0, 1]
        double edge_param = 0.0;
        /// Geometric mode selected by the analytical ComputeCost routine
        Mode mode = Mode::kUnreachable;
        /// Whether the candidate is traversable and numerically valid
        bool valid = false;
    };

    /// Best successor geometry chosen for a corner node
    struct SuccessorChoice {
        /// Total cost-to-go through the selected exit edge
        double cost = std::numeric_limits<double>::infinity();
        /// Index of the consecutive neighbor pair that produced the optimum
        int pair_index = -1;
        /// Interpolation parameter on the chosen edge
        double edge_param = 0.0;
        /// Mode reported by computeCostChoice()
        InterpChoice::Mode mode = InterpChoice::Mode::kUnreachable;
        /// First endpoint of the chosen exit edge: the perpendicular neighbor
        int s1_idx = -1;
        /// Second endpoint of the chosen exit edge: the diagonal neighbor
        int s2_idx = -1;
        /// Primary cell cut through by the path
        GridCoord c_cell;
        /// Secondary boundary cell sharing the chosen edge
        GridCoord b_cell;
        /// Whether this node has a valid finite successor geometry
        bool valid = false;
    };

    /// Candidate waypoint lying on a cell boundary during path extraction
    struct BoundaryCandidate {
        /// Fractional waypoint on a cell edge or the goal point itself
        Point2D point;
        /// Interpolated cost-to-go value of the waypoint
        double point_cost = std::numeric_limits<double>::infinity();
        /// Total local score used for selecting the next waypoint
        double total_cost = std::numeric_limits<double>::infinity();
        /// Cell through which the current segment travels
        GridCoord cell;
        /// First endpoint of the supporting edge for loop-avoidance bookkeeping
        int edge_n1_idx = -1;
        /// Second endpoint of the supporting edge for loop-avoidance bookkeeping
        int edge_n2_idx = -1;
        /// Whether the candidate is valid
        bool valid = false;
    };

    /// Active OPEN-list metadata used for lazy deletion of stale heap entries
    struct OpenRecord {
        /// Last valid key associated with the node
        Key key;
        /// Last valid OPEN insertion version
        int version = 0;
    };

    /// Heap element stored in OPEN. The node itself is addressed by index
    struct PQElement {
        /// Lexicographic priority of the node
        Key key;
        /// Version number associated with this heap insertion
        int version = 0;
        /// Index into nodes_
        int node_index = -1;

        bool operator>(const PQElement& other) const;
    };

    /// Converts a node coordinate into the flat nodes_ storage index
    inline int getIndex(int x, int y) const { return y * (width_ + 1) + x; }

    /// Validation helpers
    void validateGrid(const CostGrid& grid) const;
    void validateGoal(GridCoord goal) const;
    void validateStart(Point2D start) const;
    void ensureConfigured() const;

    /// Returns the traversal cost per unit length of one cell
    double getCellCost(int x, int y) const;

    /// Computes the priority key for one node
    Key calculateKeyValue(int index) const;

    /// Returns whether a node currently satisfies the local consistency condition
    bool isConsistent(int index) const;

    /// Recomputes rhs(index) and synchronizes the node with OPEN
    void updateNode(int index);

    /// Computes the minimum rhs value from all local interpolation successors
    double computeMinRHS(int index);

    /**
     * Implements the analytical Field D* ComputeCost routine.
     *
     * Parameters:
     * - c: traversal cost of the primary cell cut by the path,
     * - b: traversal cost of the boundary-adjacent cell sharing the exit edge
     * - g1: cost-to-go of the perpendicular neighbor
     * - g2: cost-to-go of the diagonal neighbor
     */
    InterpChoice computeCostChoice(double c, double b, double g1, double g2) const;

    /// Returns the best successor geometry for a corner node
    SuccessorChoice bestSuccessorChoice(int index) const;

    /// OPEN-list maintenance with logical deletion of stale entries
    void pushOpen(int index);
    void removeFromOpen(int index);
    Key topKey();
    bool popValid(int& out_index, Key& out_key);

    /// Returns the 8-connected corner neighbors of a node
    std::vector<int> getNeighbors(int index) const;

    /// In this uniform grid formulation, predecessors and successors coincide
    std::vector<int> predecessors(int index) const;

    /// Interpolates a fractional point on the edge connecting two corner nodes
    Point2D pointOnEdge(int n1_idx, int n2_idx, double t) const;

    /// Interpolates the edge cost-to-go from the g values of its endpoints
    double interpolateEdgeCost(int n1_idx, int n2_idx, double t) const;

    /// Returns all cells incident to a fractional point, accounting for tolerance
    std::vector<GridCoord> cellsContainingPoint(const Point2D& point) const;

    /// Checks whether a point lies inside the closed footprint of a cell
    bool pointInCell(const Point2D& point, const GridCoord& cell) const;

    /// Returns the four corner-node indices of one cell
    std::array<int, 4> cellCorners(const GridCoord& cell) const;

    /// Converts a node successor choice into a concrete fractional waypoint
    Point2D choicePoint(const SuccessorChoice& choice) const;

    /// Returns the interpolated cost-to-go at the chosen waypoint
    double choicePointCost(const SuccessorChoice& choice) const;

    /// Returns the direct in-cell cost to the goal if the goal is locally visible
    double directGoalCost(const Point2D& point) const;

    /// Builds local boundary-waypoint candidates around a fractional point.
    std::vector<BoundaryCandidate> localSuccessorCandidates(const Point2D& point) const;

    /// Returns the best local cost-to-go for a fractional point, optionally forbidding one edge.
    double pointCost(const Point2D& point, const std::vector<std::pair<int, int>>& forbidden_edges = {}) const;

    /// Scores local waypoint candidates using a one-step lookahead
    std::vector<BoundaryCandidate> nextBoundaryCandidates(const Point2D& point) const;

    /// Returns the corner node used as the D* Lite start reference
    int closestNodeToStart() const;

    /// Canonicalizes an edge endpoint pair so it can be compared or stored in sets
    std::pair<int, int> normalizeEdge(int n1_idx, int n2_idx) const;

    CostGrid grid_;
    int width_ = 0;
    int height_ = 0;
    double penalty_weight_ = static_cast<double>(kDefaultPenalty);

    std::vector<Node> nodes_;
    /// OPEN priority queue containing current and stale entries
    std::priority_queue<PQElement, std::vector<PQElement>, std::greater<PQElement>> open_list_;
    /// Metadata describing which OPEN entries are currently valid
    std::unordered_map<int, OpenRecord> open_entries_;

    Point2D start_pos_;
    GridCoord goal_coord_;
    std::vector<Point2D> path_;

    int goal_index_ = -1;
    /// key modifier tracking how far the start has moved
    double k_m_ = 0.0;
    bool configured_ = false;
    bool initialized_ = false;
    bool is_optimal_ = false;

    /// Safety bound protecting against accidental infinite propagation loops
    static constexpr int kMaxIterations = 1000000;
};

}  // namespace field_d_star

#endif
