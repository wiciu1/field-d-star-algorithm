#include "field_d_star_planner.hpp"

#include <SDL2/SDL.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using field_d_star::FieldDStarPlanner;
using field_d_star::GridCoord;
using field_d_star::Point2D;

constexpr int kCols = 20;
constexpr int kRows = 15;
constexpr int kCellSize = 50;
constexpr int kWindowWidth = kCols * kCellSize;
constexpr int kWindowHeight = kRows * kCellSize;
constexpr int kTargetFrameMs = 16;

enum class BrushMode {
    kCost,
    kObstacle,
    kErase
};

using CostGrid = FieldDStarPlanner::CostGrid;

CostGrid makeGrid(int cols, int rows)
{
    return CostGrid(static_cast<std::size_t>(cols), std::vector<int>(static_cast<std::size_t>(rows), 0));
}

void fillRect(CostGrid& grid, int min_x, int min_y, int max_x, int max_y, int value)
{
    for (int x = std::max(0, min_x); x <= std::min(max_x, static_cast<int>(grid.size()) - 1); ++x) {
        for (int y = std::max(0, min_y); y <= std::min(max_y, static_cast<int>(grid[x].size()) - 1); ++y) {
            grid[x][y] = value;
        }
    }
}

CostGrid buildInitialGrid()
{
    CostGrid grid = makeGrid(kCols, kRows);

    fillRect(grid, 12, 2, 16, 4, 90);
    fillRect(grid, 4, 5, 7, 7, 50);
    fillRect(grid, 14, 9, 16, 11, field_d_star::kObstacleValue);

    return grid;
}

SDL_Color cellColor(int cost_value)
{
    if (cost_value >= field_d_star::kObstacleValue) {
        return SDL_Color{101, 67, 33, 255};
    }

    const double ratio = std::clamp(static_cast<double>(cost_value) / 253.0, 0.0, 1.0);
    const Uint8 red = 255;
    const Uint8 green = static_cast<Uint8>(255.0 - 135.0 * ratio);
    const Uint8 blue = static_cast<Uint8>(255.0 - 255.0 * ratio);
    return SDL_Color{red, green, blue, 255};
}

void fillCircle(SDL_Renderer* renderer, int cx, int cy, int radius, SDL_Color color)
{
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
    for (int dy = -radius; dy <= radius; ++dy) {
        const int dx = static_cast<int>(std::sqrt(radius * radius - dy * dy));
        SDL_RenderDrawLine(renderer, cx - dx, cy + dy, cx + dx, cy + dy);
    }
}

class FieldDStarApp {
public:
    FieldDStarApp()
        : grid_(buildInitialGrid()),
          start_pos_{1.5, kRows - 1.5},
          goal_node_{kCols - 2, 2},
          planner_(grid_, start_pos_, goal_node_)
    {
        if (SDL_Init(SDL_INIT_VIDEO) != 0) {
            throw std::runtime_error(std::string("SDL_Init failed: ") + SDL_GetError());
        }

        window_ = SDL_CreateWindow(
            "Field D*",
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED,
            kWindowWidth,
            kWindowHeight,
            SDL_WINDOW_SHOWN);
        if (window_ == nullptr) {
            throw std::runtime_error(std::string("SDL_CreateWindow failed: ") + SDL_GetError());
        }

        renderer_ = SDL_CreateRenderer(window_, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
        if (renderer_ == nullptr) {
            renderer_ = SDL_CreateRenderer(window_, -1, SDL_RENDERER_SOFTWARE);
        }
        if (renderer_ == nullptr) {
            throw std::runtime_error(std::string("SDL_CreateRenderer failed: ") + SDL_GetError());
        }

        SDL_SetRenderDrawBlendMode(renderer_, SDL_BLENDMODE_BLEND);

        planner_.initialize();
        updateWindowTitle();
        printControls();
    }

    ~FieldDStarApp()
    {
        if (renderer_ != nullptr) {
            SDL_DestroyRenderer(renderer_);
        }
        if (window_ != nullptr) {
            SDL_DestroyWindow(window_);
        }
        SDL_Quit();
    }

    void run()
    {
        const int headless_frames = readHeadlessFrames();
        int frame_count = 0;

        while (running_) {
            handleEvents();
            render();

            if (headless_frames > 0) {
                ++frame_count;
                if (frame_count >= headless_frames) {
                    running_ = false;
                }
            }

            SDL_Delay(kTargetFrameMs);
        }
    }

private:

    static int readHeadlessFrames()
    {
        const char* env_value = std::getenv("FIELD_DSTAR_HEADLESS_FRAMES");
        if (env_value == nullptr) {
            return 0;
        }
        return std::max(0, std::atoi(env_value));
    }

    void printControls() const
    {
        std::cout
            << "Field D* GUI controls\n"
            << "  Space: compute path / move agent by one waypoint\n"
            << "  1: paint cost mode\n"
            << "  2: obstacle mode\n"
            << "  3: erase mode\n"
            << "  Up / Down or mouse wheel: change brush value in cost mode\n"
            << "  Left drag: paint using current brush\n"
            << "  Right drag: erase to zero\n"
            << "  Esc or window close: quit\n";
    }


    void updateWindowTitle() const
    {
        const char* mode_text = "Cost";
        if (brush_mode_ == BrushMode::kObstacle) {
            mode_text = "Obstacle";
        } else if (brush_mode_ == BrushMode::kErase) {
            mode_text = "Erase";
        }

        const std::string title =
            "Field D* | Mode: " + std::string(mode_text) +
            " | Brush: " + std::to_string(brush_value_) +
            " | Space=plan/move, 1/2/3=mode, wheel=brush";
        SDL_SetWindowTitle(window_, title.c_str());
    }

    int currentBrushValue() const
    {
        if (brush_mode_ == BrushMode::kObstacle) {
            return field_d_star::kObstacleValue;
        }
        if (brush_mode_ == BrushMode::kErase) {
            return 0;
        }
        return brush_value_;
    }

    void paintFromMouse(int mouse_x, int mouse_y, bool erase_override)
    {
        const int cx = mouse_x / kCellSize;
        const int cy = mouse_y / kCellSize;
        applyCellValue(cx, cy, erase_override ? 0 : currentBrushValue());
    }

    void applyCellValue(int cx, int cy, int value)
    {
        if (cx < 0 || cx >= kCols || cy < 0 || cy >= kRows) {
            return;
        }

        value = std::clamp(value, 0, field_d_star::kObstacleValue);

        const Point2D agent = planner_.start();
        if (value >= field_d_star::kObstacleValue &&
            cx <= agent.x && agent.x <= cx + 1 &&
            cy <= agent.y && agent.y <= cy + 1) {
            return;
        }

        if (grid_[cx][cy] == value) {
            return;
        }

        grid_[cx][cy] = value;
        planner_.updateCellCost(cx, cy, value);
    }

    void handleSpacePress()
    {
        if (!planner_.isOptimal()) {
            planner_.computeShortestPath();
            planner_.extractPath();
            return;
        }

        const Point2D start = planner_.start();
        if (std::hypot(start.x - goal_node_.x, start.y - goal_node_.y) < 0.1) {
            std::cout << "Goal achieved\n";
            return;
        }

        if (!planner_.advanceStartToNextWaypoint()) {
            std::cout << "No path to the goal\n";
        }
    }

    void handleEvents()
    {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_QUIT:
                running_ = false;
                break;
            case SDL_KEYDOWN:
                handleKeyDown(event.key.keysym.sym);
                break;
            case SDL_MOUSEBUTTONDOWN:
                if (event.button.button == SDL_BUTTON_LEFT) {
                    left_mouse_down_ = true;
                    paintFromMouse(event.button.x, event.button.y, false);
                } else if (event.button.button == SDL_BUTTON_RIGHT) {
                    right_mouse_down_ = true;
                    paintFromMouse(event.button.x, event.button.y, true);
                }
                break;
            case SDL_MOUSEBUTTONUP:
                if (event.button.button == SDL_BUTTON_LEFT) {
                    left_mouse_down_ = false;
                } else if (event.button.button == SDL_BUTTON_RIGHT) {
                    right_mouse_down_ = false;
                }
                break;
            case SDL_MOUSEMOTION:
                if (left_mouse_down_) {
                    paintFromMouse(event.motion.x, event.motion.y, false);
                }
                if (right_mouse_down_) {
                    paintFromMouse(event.motion.x, event.motion.y, true);
                }
                break;
            case SDL_MOUSEWHEEL:
                adjustBrushValue(event.wheel.y);
                break;
            default:
                break;
            }
        }
    }

    void handleKeyDown(SDL_Keycode key)
    {
        switch (key) {
        case SDLK_ESCAPE:
            running_ = false;
            break;
        case SDLK_SPACE:
            handleSpacePress();
            break;
        case SDLK_1:
            brush_mode_ = BrushMode::kCost;
            updateWindowTitle();
            break;
        case SDLK_2:
            brush_mode_ = BrushMode::kObstacle;
            updateWindowTitle();
            break;
        case SDLK_3:
            brush_mode_ = BrushMode::kErase;
            updateWindowTitle();
            break;
        case SDLK_UP:
            adjustBrushValue(1);
            break;
        case SDLK_DOWN:
            adjustBrushValue(-1);
            break;
        default:
            break;
        }
    }

    void adjustBrushValue(int delta)
    {
        brush_value_ = std::clamp(brush_value_ + delta * 5, 0, 253);
        updateWindowTitle();
    }

    void render()
    {
        SDL_SetRenderDrawColor(renderer_, 255, 255, 255, 255);
        SDL_RenderClear(renderer_);

        drawGrid();
        drawPath();
        drawGoal();
        drawAgent();

        SDL_RenderPresent(renderer_);
    }

    void drawGrid() const
    {
        for (int x = 0; x < kCols; ++x) {
            for (int y = 0; y < kRows; ++y) {
                const SDL_Color fill = cellColor(grid_[x][y]);
                SDL_Rect rect{x * kCellSize, y * kCellSize, kCellSize, kCellSize};

                SDL_SetRenderDrawColor(renderer_, fill.r, fill.g, fill.b, fill.a);
                SDL_RenderFillRect(renderer_, &rect);

                SDL_SetRenderDrawColor(renderer_, 221, 221, 221, 255);
                SDL_RenderDrawRect(renderer_, &rect);
            }
        }
    }

    void drawPath() const
    {
        const auto& path = planner_.path();
        if (path.size() < 2) {
            return;
        }

        SDL_SetRenderDrawColor(renderer_, 0, 102, 255, 255);
        for (std::size_t i = 1; i < path.size(); ++i) {
            const Point2D& a = path[i - 1];
            const Point2D& b = path[i];
            SDL_RenderDrawLine(
                renderer_,
                static_cast<int>(std::lround(a.x * kCellSize)),
                static_cast<int>(std::lround(a.y * kCellSize)),
                static_cast<int>(std::lround(b.x * kCellSize)),
                static_cast<int>(std::lround(b.y * kCellSize)));
        }
    }

    void drawGoal() const
    {
        fillCircle(
            renderer_,
            goal_node_.x * kCellSize,
            goal_node_.y * kCellSize,
            8,
            SDL_Color{128, 0, 128, 255});
    }

    void drawAgent() const
    {
        const Point2D start = planner_.start();
        fillCircle(
            renderer_,
            static_cast<int>(std::lround(start.x * kCellSize)),
            static_cast<int>(std::lround(start.y * kCellSize)),
            6,
            SDL_Color{0, 170, 0, 255});
    }

    CostGrid grid_;
    Point2D start_pos_;
    GridCoord goal_node_;
    FieldDStarPlanner planner_;

    SDL_Window* window_ = nullptr;
    SDL_Renderer* renderer_ = nullptr;

    BrushMode brush_mode_ = BrushMode::kCost;
    int brush_value_ = 80;
    bool running_ = true;
    bool left_mouse_down_ = false;
    bool right_mouse_down_ = false;
};

}  // namespace

int main()
{
    try {
        FieldDStarApp app;
        app.run();
    } catch (const std::exception& error) {
        std::cerr << "Application error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
