# Project Overview: Teste Motor (Non-Linear Structural Analysis)

This project is a C++ application designed for the non-linear structural analysis of beams using the Finite Element Method (FEM). It provides a graphical interface for modeling, solving, and visualizing structural responses.

## Main Technologies
- **Language:** C++ (Standard C++14)
- **Graphics & Windowing:** [Raylib](https://www.raylib.com/)
- **User Interface:** [Dear ImGui](https://github.com/ocornut/imgui)
- **Plotting:** [ImPlot](https://github.com/epezent/implot)
- **Linear Algebra:** [Eigen](https://eigen.tuxfamily.org/)
- **Integration:** `rlImGui` for Raylib-ImGui integration.

## Architecture
The application is currently structured with its primary logic contained within `main.cpp`, supported by several libraries integrated directly into the source tree.
- **`main.cpp`**: Contains the core Finite Element Analysis (FEA) logic, including classes for nodes (`No`), finite elements (`ElementoFinito`, `Viga2DLinear`), and the main simulation loop.
- **Root Directory**: Houses the source files for ImGui, ImPlot, and the main application.
- **`include/`**: Contains headers for Eigen, ImGui, Raylib, and project-specific headers like `eigenpch.h`.
- **`resources/`**: Stores assets such as fonts (`segoeuisl.ttf`) and images.
- **`build/`**: Target directory for compiled executables.

## Building and Running

### Prerequisites
- **Compiler:** GCC (MinGW-w64 recommended on Windows).
- **Environment:** Raylib must be installed or accessible in the path specified in the `Makefile`.

### Desktop Build
To build the project on a desktop platform, use the provided `Makefile`.

```bash
# Using mingw32-make on Windows
mingw32-make RAYLIB_PATH=C:/raylib/raylib PROJECT_NAME=main OBJS=*.cpp

# Using standard make on Linux/macOS
make PROJECT_NAME=main OBJS=*.cpp
```

### VS Code Integration
The project includes a `.vscode/tasks.json` with pre-configured build tasks:
- **build debug**: Compiles with debug symbols using `mingw32-make`.
- **build release**: Compiles with optimizations.

### Running
After building, execute the generated binary:
```bash
./main.exe
```
*Note: A pre-built version may exist at `build/analise-nao-linear.exe`.*

### Android Build
A `Makefile.Android` is provided for building Android APKs, requiring the Android NDK and SDK to be configured.

## Development Conventions
- **Source Organization:** While the `obj/` directory suggests a modular history (`app/`, `editor/`, etc.), the current active development is centered in the root directory and `main.cpp`.
- **Naming Conventions:** The codebase uses a mix of Portuguese (e.g., `normalizaAngulo`, `EsforcosLocais`) and English (e.g., `getMatrizRigidezGlobal`). Maintain consistency with the existing mathematical and structural terminology.
- **Dependencies:** Avoid adding external dependencies that are not already present in the `include/` or root directory without updating the `Makefile`.
- **Headers:** Prefer including `eigenpch.h` for Eigen-related operations to potentially benefit from precompiled headers if configured.

## TODO / Future Improvements
- [ ] Refactor `main.cpp` into modular components (e.g., `model/`, `solver/`, `ui/`) to match the intended architecture seen in `obj/`.
- [ ] Standardize language usage in code (Portuguese vs. English).
- [ ] Update `Makefile` to automatically discover all `.cpp` files in a modular structure.
