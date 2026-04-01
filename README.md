# Field D* Path Planning Algorithm

## Visualization
![img.png](img.png)

## Overview
Classical grid-based planners usually restrict motion to headings spaced by 45 degrees. That makes the resulting paths unnaturally jagged, even when the planner is optimal on the underlying grid graph.

Field D* improves this by using linear interpolation between corner-node values. Instead of forcing the path to move only along grid edges and diagonals, it can leave a cell through interpolated boundary points and produce smoother, lower-cost trajectories.

This repository contains:
- a reference prototype in Python,
- a plain C++ port of the algorithm,
- an SDL2 GUI simulation in C++ that mirrors the Python app.

## Repository Structure
- [field_d_star.py](/home/wiktor/C/Studia/Stellar/field_d_star_algorithm/field_d_star.py) - reference Python implementation of the Field D* core.
- [simulation.py](/home/wiktor/C/Studia/Stellar/field_d_star_algorithm/simulation.py) - Tkinter-based Python simulator.
- [cpp/field_d_star_planner.hpp](/home/wiktor/C/Studia/Stellar/field_d_star_algorithm/cpp/field_d_star_planner.hpp) - plain C++ API and documented data structures.
- [cpp/field_d_star_planner.cpp](/home/wiktor/C/Studia/Stellar/field_d_star_algorithm/cpp/field_d_star_planner.cpp) - C++ implementation of the planner.
- [cpp/main.cpp](/home/wiktor/C/Studia/Stellar/field_d_star_algorithm/cpp/main.cpp) - SDL2 GUI simulation for the C++ version.
- [CMakeLists.txt](/home/wiktor/C/Studia/Stellar/field_d_star_algorithm/CMakeLists.txt) - CMake build for the C++ planner and GUI demo.

## Python Version
The original prototype requires only Python 3 and the standard library.

Run:
```bash
python3 simulation.py
```

## C++ Version
The C++ port is a plain implementation of Field D*. It operates on:
- a cost grid indexed as `grid[x][y]`,
- a fractional start position `Point2D`,
- a goal corner `GridCoord`.

The GUI simulator in C++ uses SDL2 and exposes the same style of interaction as the Python version.

### Build With CMake
Requirements:
- C++17 compiler
- CMake 3.16+
- SDL2
- `pkg-config`

Build:
```bash
cmake -S . -B build
cmake --build build
```

Run:
```bash
./build/field_d_star_gui
```

### Manual Build
If needed, the GUI can also be built directly:
```bash
g++ -std=c++17 -Wall -Wextra -pedantic -Icpp $(pkg-config --cflags sdl2) \
    cpp/main.cpp cpp/field_d_star_planner.cpp \
    $(pkg-config --libs sdl2) \
    -o field_d_star_gui
```

## C++ GUI Controls
- `Space` - first compute the path, then move the agent by one waypoint on each next press
- `1` - paint traversal cost
- `2` - paint obstacle (`254`)
- `3` - erase cell to `0`
- `Up` / `Down` - change brush value in cost-paint mode
- mouse wheel - change brush value in cost-paint mode
- left mouse drag - paint using the current brush
- right mouse drag - erase to `0`
- `Esc` - quit

## Literature
- https://www.ri.cmu.edu/pub_files/pub4/ferguson_david_2006_3/ferguson_david_2006_3.pdf
- https://www-robotics.jpl.nasa.gov/media/documents/fdstar3d.pdf
