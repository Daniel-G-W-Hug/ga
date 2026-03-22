# VTK 3D Visualization Migration - Conversation Context

## User Background & Requirements

### Current Setup

- **Existing Framework**: 2D visualization framework written in C++ based on Qt6
- **Current Capabilities**:
  - Display coordinate systems
  - Render 2D primitives
  - User interaction via mouse
- **Platform**: Windows (primary), macOS (secondary)
- **VR Hardware**: Pimax headset on Windows (SteamVR/OpenXR); macOS headset TBD

### Project Goals

- Migrate to 3D visualization using C++ framework
- Support both on-screen (desktop) and VR visualization
- Preferably use the same codebase for both modes
- Focus on scientific visualization with interaction capabilities
- Support **Python-using coworkers** without requiring full Python wrapper of GA lib

### Initial Requirements

- Draw 3D coordinate systems
- Visualize intersecting planes with visible intersection lines
- Interactive capabilities (mouse for desktop, VR controllers for headset)
- Integrate with existing Qt6 application
- C++ primary, Python secondary (via PyVista — see below)

## Toolkit Decision

### Alternatives Considered

| Toolkit | C++ | Python | VR | Verdict |
| --- | --- | --- | --- | --- |
| **VTK** | First-class | First-class (PyVista) | OpenXR (Windows) | **Selected** |
| Unreal Engine | Native | Plugin only | Excellent | Overkill, heavy build system |
| Open3D | Yes | Yes | None | No VR support |
| Godot | Via GDExtension | Community bindings | OpenXR | C++ boundary layer awkward |
| WebXR (Three.js) | Via Emscripten/WASM | Via JS | Most portable | Paradigm shift from Qt/C++ |
| OpenSceneGraph | Yes | PyOSG (unmaintained) | Limited | Aging |

### Python Wrapper: Build First or After?

**Decision: After.** Do not block 3D toolkit work on writing a Python wrapper for the GA
library. PyVista (a high-level Python wrapper around VTK) gives Python coworkers access to
the visualization layer immediately. A GA Python binding can be added later when there is
concrete demand.

## Recommended Solution: VTK (Visualization Toolkit)

### Why VTK?

- Purpose-built for scientific/medical/engineering visualization
- Excellent C++ API with comprehensive interaction support
- Native mouse and VR controller input handling
- **Same visualization pipeline works for both desktop and VR**
- Direct Qt6 integration via `QVTKOpenGLNativeWidget`
- Native VR support through OpenXR (no external frameworks like ParaView needed)
- **PyVista** provides a first-class Python interface on top of VTK without custom wrapper
  work

### Key VTK Components Used

#### Desktop Mode

- `vtkGenericOpenGLRenderWindow` - Standard rendering to screen
- `QVTKOpenGLNativeWidget` - Qt widget for VTK integration
- `vtkInteractorStyleTrackballCamera` - Mouse interaction

#### VR Mode

- `vtkOpenXRRenderWindow` - VR headset rendering via OpenXR
- `vtkOpenXRRenderer` - VR-specific renderer
- `vtkOpenXRRenderWindowInteractor` - VR controller interaction

#### Scene Components

- `vtkRenderer` - Main rendering engine
- `vtkAxesActor` - 3D coordinate system visualization
- `vtkPlaneSource` - Geometric plane creation
- `vtkCutter` - Calculate plane intersections
- `vtkTubeFilter` - Visualize intersection lines

## Migration Strategy

### Phase 1: Hybrid Approach (Week 1-2)

- Keep existing Qt 2D rendering
- Add VTK widget alongside for 3D experimentation
- Learn VTK basics with coordinate systems and planes

### Phase 2: Gradual Migration (Week 3-4)

- Optionally migrate 2D primitives to VTK
- VTK handles 2D well with `vtkActor2D`, `vtkPolyLine`
- Benefit: Same interaction code for 2D and 3D

### Phase 3: Full 3D + VR

- All visualization in VTK
- Easy VR support (just swap render window type)

## Complete Working Implementation

### Project Structure

```text
project/
├── CMakeLists.txt
├── main.cpp
└── Scene3D.h
```

### CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.12)
project(VTK3DVisualization)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_AUTOMOC ON)

# Find VTK with Qt and OpenXR support
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersSources
    GUISupportQt
    InteractionStyle
    RenderingCore
    RenderingOpenGL2
    RenderingOpenXR  # For VR support (Windows/OpenXR)
)

find_package(Qt6 COMPONENTS Core Widgets REQUIRED)

add_executable(${PROJECT_NAME} 
    main.cpp
    Scene3D.h
)

target_link_libraries(${PROJECT_NAME}
    ${VTK_LIBRARIES}
    Qt6::Core
    Qt6::Widgets
)

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)
```

### Scene3D.h (Reusable Scene - Works for Both Desktop & VR)

```cpp
#ifndef SCENE3D_H
#define SCENE3D_H

#include <vtkRenderer.h>
#include <vtkAxesActor.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkTubeFilter.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>

class Scene3D {
public:
    static void setupScene(vtkRenderer* renderer) {
        // Background color
        renderer->SetBackground(0.1, 0.1, 0.15);
        
        // 1. Add 3D coordinate system with labels
        vtkNew<vtkAxesActor> axes;
        axes->SetTotalLength(2.5, 2.5, 2.5);
        axes->SetShaftTypeToCylinder();
        axes->SetCylinderRadius(0.02);
        axes->SetAxisLabels(1);
        axes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
        axes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
        axes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
        renderer->AddActor(axes);
        
        // 2. First plane (tilted XY plane) - Red
        vtkNew<vtkPlaneSource> planeXY;
        planeXY->SetOrigin(-2, -2, 0);
        planeXY->SetPoint1(2, -2, 0.6);
        planeXY->SetPoint2(-2, 2, 0.4);
        planeXY->SetResolution(20, 20);
        
        vtkNew<vtkPolyDataMapper> mapperXY;
        mapperXY->SetInputConnection(planeXY->GetOutputPort());
        
        vtkNew<vtkActor> actorXY;
        actorXY->SetMapper(mapperXY);
        actorXY->GetProperty()->SetColor(0.9, 0.3, 0.3);
        actorXY->GetProperty()->SetOpacity(0.7);
        actorXY->GetProperty()->SetInterpolationToPhong();
        actorXY->GetProperty()->BackfaceCullingOff();
        renderer->AddActor(actorXY);
        
        // 3. Second plane (tilted YZ plane) - Green
        vtkNew<vtkPlaneSource> planeYZ;
        planeYZ->SetOrigin(0, -2, -2);
        planeYZ->SetPoint1(0.6, 2, -2);
        planeYZ->SetPoint2(0.4, -2, 2);
        planeYZ->SetResolution(20, 20);
        
        vtkNew<vtkPolyDataMapper> mapperYZ;
        mapperYZ->SetInputConnection(planeYZ->GetOutputPort());
        
        vtkNew<vtkActor> actorYZ;
        actorYZ->SetMapper(mapperYZ);
        actorYZ->GetProperty()->SetColor(0.3, 0.9, 0.3);
        actorYZ->GetProperty()->SetOpacity(0.7);
        actorYZ->GetProperty()->SetInterpolationToPhong();
        actorYZ->GetProperty()->BackfaceCullingOff();
        renderer->AddActor(actorYZ);
        
        // 4. Calculate intersection line
        // Create a plane from the second plane source
        vtkNew<vtkPlane> cuttingPlane;
        cuttingPlane->SetOrigin(0, 0, 0);
        cuttingPlane->SetNormal(0.6, 2, -2);  // Direction from plane setup
        
        vtkNew<vtkCutter> cutter;
        cutter->SetInputConnection(planeXY->GetOutputPort());
        cutter->SetCutFunction(cuttingPlane);
        
        // Make the intersection line visible as a tube
        vtkNew<vtkTubeFilter> tubeFilter;
        tubeFilter->SetInputConnection(cutter->GetOutputPort());
        tubeFilter->SetRadius(0.03);
        tubeFilter->SetNumberOfSides(16);
        
        vtkNew<vtkPolyDataMapper> lineMapper;
        lineMapper->SetInputConnection(tubeFilter->GetOutputPort());
        
        vtkNew<vtkActor> lineActor;
        lineActor->SetMapper(lineMapper);
        lineActor->GetProperty()->SetColor(1.0, 1.0, 0.0);  // Yellow
        lineActor->GetProperty()->SetSpecular(0.8);
        lineActor->GetProperty()->SetSpecularPower(50);
        renderer->AddActor(lineActor);
        
        // 5. Set nice initial camera position
        renderer->ResetCamera();
        vtkCamera* camera = renderer->GetActiveCamera();
        camera->Azimuth(45);
        camera->Elevation(30);
        camera->Dolly(1.2);
        renderer->ResetCameraClippingRange();
    }
};

#endif // SCENE3D_H
```

### main.cpp (Desktop & VR Mode Switching)

```cpp
#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>
#include <QHBoxLayout>
#include <QLabel>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

// For VR support
#include <vtkOpenXRRenderWindow.h>
#include <vtkOpenXRRenderWindowInteractor.h>
#include <vtkOpenXRRenderer.h>

#include "Scene3D.h"

class MainWindow : public QMainWindow {
    Q_OBJECT
    
public:
    MainWindow() {
        setupUI();
        setupDesktopMode();
    }
    
private slots:
    void switchToDesktop() {
        if (currentMode == Desktop) return;
        
        qDebug() << "Switching to Desktop mode...";
        currentMode = Desktop;
        setupDesktopMode();
        statusLabel->setText("Mode: Desktop - Use mouse to rotate/pan/zoom");
    }
    
    void switchToVR() {
        if (currentMode == VR) return;
        
        qDebug() << "Switching to VR mode...";
        currentMode = VR;
        setupVRMode();
        statusLabel->setText("Mode: VR - Put on your headset!");
    }
    
private:
    enum Mode { Desktop, VR };
    Mode currentMode = Desktop;
    
    QVTKOpenGLNativeWidget* vtkWidget;
    QLabel* statusLabel;
    vtkNew<vtkRenderer> renderer;
    
    void setupUI() {
        QWidget* centralWidget = new QWidget(this);
        QVBoxLayout* mainLayout = new QVBoxLayout(centralWidget);
        
        // Control buttons
        QHBoxLayout* buttonLayout = new QHBoxLayout();
        QPushButton* desktopBtn = new QPushButton("Desktop Mode");
        QPushButton* vrBtn = new QPushButton("VR Mode");
        
        connect(desktopBtn, &QPushButton::clicked, this, &MainWindow::switchToDesktop);
        connect(vrBtn, &QPushButton::clicked, this, &MainWindow::switchToVR);
        
        buttonLayout->addWidget(desktopBtn);
        buttonLayout->addWidget(vrBtn);
        buttonLayout->addStretch();
        
        // Status label
        statusLabel = new QLabel("Mode: Desktop");
        
        // VTK Widget
        vtkWidget = new QVTKOpenGLNativeWidget(this);
        
        mainLayout->addLayout(buttonLayout);
        mainLayout->addWidget(statusLabel);
        mainLayout->addWidget(vtkWidget);
        
        setCentralWidget(centralWidget);
        setWindowTitle("VTK 3D Visualization - Desktop & VR");
        resize(1024, 768);
    }
    
    void setupDesktopMode() {
        // Create standard OpenGL render window
        vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
        vtkWidget->setRenderWindow(renderWindow);
        
        // Clear and re-add renderer
        renderWindow->AddRenderer(renderer);
        
        // Setup the 3D scene (same scene for both modes!)
        renderer->RemoveAllViewProps();
        Scene3D::setupScene(renderer);
        
        // Enable mouse interaction
        vtkNew<vtkInteractorStyleTrackballCamera> style;
        renderWindow->GetInteractor()->SetInteractorStyle(style);
        renderWindow->GetInteractor()->Initialize();
        
        renderWindow->Render();
    }
    
    void setupVRMode() {
        // Create OpenXR render window for VR
        vtkNew<vtkOpenXRRenderWindow> vrRenderWindow;
        vrRenderWindow->Initialize();
        
        // Create VR renderer
        vtkNew<vtkOpenXRRenderer> vrRenderer;
        vrRenderWindow->AddRenderer(vrRenderer);
        
        // Copy the scene to VR renderer (same setup!)
        Scene3D::setupScene(vrRenderer);
        
        // Create VR interactor
        vtkNew<vtkOpenXRRenderWindowInteractor> vrInteractor;
        vrInteractor->SetRenderWindow(vrRenderWindow);
        vrInteractor->Initialize();
        
        // Hide Qt window (VR will show in headset)
        vtkWidget->hide();
        
        // Start VR rendering loop
        vrInteractor->Start();
        
        // When VR exits, show desktop again
        vtkWidget->show();
        switchToDesktop();
    }
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    
    // Needed for VTK 9.0+
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    
    MainWindow mainWindow;
    mainWindow.show();
    
    return app.exec();
}

#include "main.cpp.moc"
```

## Desktop vs VR Switching - Technical Details

### Native VTK Support (No ParaView Required)

VTK natively supports both rendering modes through different render window types:

#### Desktop Mode part 2

- **Render Window**: `vtkGenericOpenGLRenderWindow`
- **Display**: Qt widget embedded in application
- **Interaction**: Mouse/keyboard via `vtkInteractorStyleTrackballCamera`
- **Output**: Standard monitor

#### VR Mode part 2

- **Render Window**: `vtkOpenXRRenderWindow`
- **Display**: VR headset via OpenXR API
- **Interaction**: VR controllers automatically handled
- **Output**: Stereoscopic rendering to headset

### Key Principle: Same Scene, Different Window

```cpp
// The SAME scene setup function works for both!
Scene3D::setupScene(renderer);

// Only the render window type changes:
// Desktop: vtkGenericOpenGLRenderWindow
// VR:      vtkOpenXRRenderWindow
```

### Interaction Capabilities

#### Desktop (Automatic)

- **Left Mouse Drag**: Rotate camera around scene
- **Middle Mouse Drag**: Pan (translate) camera
- **Right Mouse Drag / Scroll**: Zoom in/out
- **R Key**: Reset camera to initial view

#### VR (Automatic)

- **Controller Triggers**: Grab and manipulate scene
- **Joysticks/Touchpads**: Navigate/teleport
- **Grip Buttons**: Scale scene

## Building with VR Support

### VTK Build Requirements

```bash
cmake -DVTK_MODULE_ENABLE_VTK_RenderingOpenXR=YES \
      -DVTK_GROUP_ENABLE_Qt=YES \
      /path/to/vtk/source
```

### Runtime Requirements for VR

- **OpenXR Runtime** installed:
  - SteamVR (most common, works with most headsets including Pimax)
  - Oculus Runtime (for Meta Quest Link)
  - Windows Mixed Reality (for WMR headsets)
- **VR Headset** connected to PC
- Compatible headsets: Pimax, Meta Quest, Valve Index, HP Reverb, HTC Vive, etc.

### macOS VR — Important Limitations

**macOS VR support is fundamentally fragmented from Windows:**

| Scenario | Status |
| --- | --- |
| Windows + Pimax via SteamVR/OpenXR | Works with VTK OpenXR |
| macOS + Apple Vision Pro | Requires visionOS/RealityKit — separate stack, no VTK path |
| macOS + Meta Quest (standalone browser) | WebXR only — different architecture |
| macOS + SteamVR | Not available (Valve dropped macOS support in 2020) |
| macOS + OpenXR runtime | No mainstream consumer runtime available |

**Recommendation**: Implement screen-based 3D and Windows VR first. Treat macOS VR as a
separate future concern once the target headset is confirmed. No single toolkit eliminates
this platform fragmentation — they all have platform-specific backends underneath.

### Project Build

```bash
mkdir build && cd build
cmake ..
cmake --build .
./VTK3DVisualization
```

## Integration with Existing Qt Framework

### Recommended Approach

1. **Keep existing Qt UI structure**
2. **Add VTK widget to layout**:

   ```cpp
   QVTKOpenGLNativeWidget* vtkWidget = new QVTKOpenGLNativeWidget();
   yourExistingLayout->addWidget(vtkWidget);
   ```

3. **Maintain existing 2D views** (or migrate gradually to VTK 2D)
4. **Use Qt signals/slots** to connect UI controls to VTK scene updates

### VTK Can Handle 2D Too

- `vtkActor2D` for 2D overlays
- `vtkPolyLine` for 2D lines
- `vtkCoordinate` for 2D coordinate systems
- Benefit: Unified codebase, same interaction model

## Common Usage Patterns

### Pattern 1: Runtime Mode Selection (Example Above)

User clicks button to switch between desktop and VR

### Pattern 2: Startup Configuration

```cpp
int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "--vr") {
        runVRMode();
    } else {
        runDesktopMode();
    }
}
```

### Pattern 3: Configuration File

```json
{
  "visualization_mode": "desktop",
  "vr_headset": "steamvr"
}
```

## What Software Does the User Actually Run?

**Not ParaView.** ParaView is an end-user application built on top of VTK for loading and
analysing large scientific datasets (CFD, FEA, etc.). It is the wrong tool here.

The user builds and runs their **own Qt6 application** — an extended `ga_view` — that
embeds VTK as a rendering library, exactly as `ga_view` embeds Qt's 2D canvas today.

```text
ParaView  = an application that uses VTK internally (you load data into it)
Your app  = your own Qt6 application that uses VTK as a rendering library
```

For Python coworkers: **PyVista** opens its own interactive window (`pyvista.Plotter`) —
they run a Python script, no separate application needed.

## Interactive Elements: VTK Widgets

VTK has a full widget system mapping directly to the draggable interactive points in
`ga_view`. These are `vtkInteractiveWidget` subclasses:

| Widget | Purpose |
| --- | --- |
| `vtkHandleWidget` | Draggable 3D point — direct equivalent of ga_view interactive points |
| `vtkSphereWidget2` | Draggable sphere with position callback |
| `vtkLineWidget2` | Interactive line between two handle points |
| `vtkPlaneWidget` | Interactive plane (drag corners/normal) |
| `vtkSliderWidget` | On-screen slider for scalar parameters |
| `vtkOrientationMarkerWidget` | Always-visible axes in viewport corner |
| `vtkBoxWidget2` | Interactive bounding box |

Each widget fires an `InteractionEvent` connected to a callback — update GA computation
there and call `Render()`.

### Handle Widget Pattern (equivalent to ga_view draggable point)

```cpp
vtkNew<vtkPointHandleRepresentation3D> handleRep;
handleRep->SetWorldPosition(point.data());

vtkNew<vtkHandleWidget> handle;
handle->SetInteractor(renderWindow->GetInteractor());
handle->SetRepresentation(handleRep);
handle->EnabledOn();

handle->AddObserver(vtkCommand::InteractionEvent, [&](
    vtkObject*, unsigned long, void*) {
    double pos[3];
    handleRep->GetWorldPosition(pos);
    // pos → update GA vector → update vtkActor geometry
    renderWindow->Render();
});
```

## Typical Development Workflow

### One iteration cycle

```text
1. Define GA objects (vectors, bivectors, motors, ...) as C++ state
2. Map each GA object to a VTK actor (arrow, disc, sphere, line, ...)
3. Add vtkHandleWidget for each user-movable point
4. In widget callback: recompute GA → update actor geometry → Render()
5. Embed in QVTKOpenGLNativeWidget inside your Qt6 window
```

### Animations: Qt timer approach

```cpp
QTimer* timer = new QTimer(this);
connect(timer, &QTimer::timeout, [&]() {
    t += dt;                       // advance GA simulation state
    updateActorsFromGAState();     // push new geometry to VTK actors
    vtkWidget->renderWindow()->Render();
});
timer->start(16); // ~60 fps
```

VTK also has `vtkAnimationScene`/`vtkAnimationCue` for keyframe-style animations, but the
Qt timer approach is simpler and sufficient for interactive GA demonstrations.

## Summary: What the Application Looks Like

One executable — the GA viewer application — that:

- Has a Qt6 window with existing UI controls
- Contains a `QVTKOpenGLNativeWidget` showing the 3D scene
- Has draggable `vtkHandleWidget` instances for user-controlled GA elements
- Optionally has a VR button that swaps the render window to `vtkOpenXRRenderWindow`

Python coworkers use a PyVista script with `add_sphere_widget` / `add_slider_widget` /
picking callbacks for their interactive needs — no shared application binary required.

## Next Steps for User

The user will want to:

1. **Share their existing Qt C++ 2D framework code** for integration advice
2. **Discuss specific 2D primitives** they need to migrate
3. **Plan the migration strategy** based on their codebase
4. **Potentially add custom interaction** (object picking, measurement tools)
5. **Optimize the scene setup** for their specific scientific visualization needs

## Python Architecture

C++ users and Python users share the same VTK scene infrastructure:

```text
GA Library (C++ header-only)
         │
         ├── C++ users ──► VTK C++ API ──► OpenXR (VR on Windows/Pimax)
         │
         └── Python users ──► PyVista (VTK Python wrapper)
                                    └── GA Python bindings (add later, when needed)
```

**PyVista** (`pip install pyvista`) is a high-level Python wrapper around VTK. Python
coworkers can use it directly without requiring a custom GA Python binding. The GA binding
can be added as a separate effort driven by concrete demand.

## Key Advantages of This Approach

✅ **Single Codebase**: Same scene definition for desktop and VR
✅ **Qt6 Integration**: Works seamlessly with existing Qt6 applications
✅ **Native VTK**: No external dependencies like ParaView needed
✅ **OpenXR Standard**: Compatible with all major VR headsets on Windows
✅ **Scientific Focus**: VTK is designed for scientific visualization
✅ **Interaction Built-in**: Mouse and VR controller support automatic
✅ **Gradual Migration**: Can keep existing 2D code while adding 3D
✅ **Python via PyVista**: No custom wrapper needed for Python coworkers
⚠️ **macOS VR**: Requires separate stack — handle after Windows path is working

## Resources

- **VTK Documentation**: <https://vtk.org/documentation/>
- **VTK Examples**: <https://examples.vtk.org/site/Cxx/>
- **Qt + VTK**: Look for examples under "GUI/Qt" section
- **OpenXR**: <https://www.khronos.org/openxr/>
- **VTK Discourse Forum**: <https://discourse.vtk.org/>

---

## ga_view Architecture Reference

Understanding the existing 2D code is essential for the 3D migration.

### Module Structure

```text
ga_view/src/
├── coordsys.hpp/cpp          - Axis + Coordsys: coordinate system with pan/zoom/tick logic
├── coordsys_model.hpp/cpp    - Coordsys_model: data model (all primitives + active items)
├── w_coordsys.hpp/cpp        - QGraphicsView: owns mouse events, maps them to active points
├── w_mainwindow.hpp/cpp      - QMainWindow: model switcher (slider/spinbox), hosts w_Coordsys
├── w_statusbar.hpp/cpp       - Status bar (mouse position, mode display)
├── w_common.hpp/cpp          - Shared interaction enums (pz_action, move_mode, etc.)
├── active_common.hpp         - Shared primitive types (pt2d, vt2d, marks, active structs)
├── item_*.hpp/cpp            - QGraphicsItem subclasses for passive primitives
└── active_*.hpp/cpp          - QGraphicsItem subclasses for interactive primitives
```

### Primitive Inventory

**Passive primitives** (drawn, not interactive):

| 2D Type | Description |
| --- | --- |
| `pt2d` | Euclidean 2D point |
| `pt2dp` | Projective 2D point |
| `ln2d` | Polyline (vector of `pt2d`) |
| `cln2dp` | Projective polyline (vector of `pt2dp`) |
| `vt2d` | Vector (beg→end with arrowhead) |
| `bivt2dp` | Projective bivector (= a line in PGA2D) |

**Active primitives** (index-based, interactive — driven by draggable `apt` points):

| 2D Type | Description |
| --- | --- |
| `apt` | Active point (draggable handle) |
| `avt2d` | Active vector (two apt indices) |
| `abivt2d` | Active bivector (two apt indices) |
| `aproj2d` | Active projection (three apt indices) |
| `arefl2d` | Active reflection (two apt indices) |
| `abivt2dp` | Active projective bivector (two apt indices) |
| `arefl2dp` | Active projective reflection (four apt indices) |
| `akinematic2dp` | Active kinematic chain (N apt indices) |
| `aode_spring2d` | ODE spring-mass on active fixation point |
| `aode_plate_pga2dp` | ODE rigid plate pendulum on active pivot point |

**Key architectural pattern:** Active items hold indices into the `apt` vector, not geometry.
`w_Coordsys` finds the nearest `apt` on mouse-press and drags it. All dependent items recompute
from current `apt` positions on every paint event. This pattern maps cleanly to VTK handle
widget callbacks.

### 2D → 3D Primitive Mapping

| 2D (ga_view) | 3D (VTK) | Notes |
| --- | --- | --- |
| `Coordsys` (axes, grid, ticks) | `vtkAxesActor` + `vtkOrientationMarkerWidget` | |
| `pt2d` / `pt2dp` | `vtkSphereSource` + `vtkActor` | Small sphere |
| `ln2d` (polyline) | `vtkPolyData` + `vtkActor` | |
| `vt2d` (arrow) | `vtkArrowSource` + `vtkTransform` + `vtkActor` | Scale+orient beg→end |
| `bivt2dp` (projective line) | `vtkPlaneSource` + `vtkActor` | Semi-transparent plane |
| `apt` (draggable point) | `vtkHandleWidget` + `vtkPointHandleRepresentation3D` | |
| Active vector | Recomputed `vtkArrowSource` from two handle positions | Same index pattern |
| Active bivector | Recomputed `vtkPlaneSource` from handle positions | |
| Active reflection | Recomputed actor from handle positions | |
| Kinematic chain | Chain of `vtkArrowSource` actors from N handles | |
| ODE simulation | Qt timer + state update → actor position/orientation | |
| Model switcher (slider) | Keep Qt slider → hide/show actor sets per model | |

### Coordsys_model → 3D Scene Mapping

`Coordsys_model` stays as-is for data storage. A new `Scene3D` class owns VTK actors and
maps model data to them, analogous to how `w_Coordsys` maps model data to `QGraphicsItem`s.

## Architectural Assessment: Is ga_view Suitable for 3D Extension?

### What the current structure does well

The core design is sound and maps naturally to VTK:

- `Coordsys_model` as a data store separate from rendering is correct
- Index-based active items (`avt2d` stores `beg_idx`/`end_idx` into `apt`) maps directly to
  VTK handle widget callbacks
- The model-switching pattern (vector of models, slider) is reusable as-is
- Signal/slot notification chain maps to VTK observer callbacks

### Three structural problems that will hurt the 3D extension

#### Problem 1: Qt rendering types embedded in the data model (biggest barrier)

In `active_common.hpp`, the mark structs carry `QPen` and `QColor`:

```cpp
struct pt2d_mark {
    QPen pen{...};   // Qt type inside what should be pure data
    Symbol symbol;
    int nsize;
};
```

These flow into `Coordsys_model`. A 3D VTK scene has no `QPen` — it uses RGB floats,
opacity, and line width. As long as Qt types live in the model, the model cannot be
shared between a 2D Qt renderer and a 3D VTK renderer.

#### Problem 2: Latent dangling reference in `active_pt2d`

`active_pt2d` holds `pt2d& m_pos` — a reference into `Coordsys_model::apt` vector. Works
today because `apt` is never resized after items are created. Any future code that adds an
`apt` after item construction silently invalidates all existing handles. For a 3D system
with dynamic scene construction this is dangerous.

#### Problem 3: Physics fused with Qt rendering in `active_ode`

`active_ode` mixes three concerns in one class: `QTimer*` (Qt event loop), RK4 integration
(pure physics, no Qt needed), and `QGraphicsItem` (Qt rendering). The physics cannot be
reused in a VTK context without extracting it.

### Recommended Refactoring: Three Focused Changes

#### Change 1 — Separate data styles from Qt styles (do before starting 3D)

Split `active_common.hpp` into two files:

```text
ga_view_types.hpp   ← pure data: pt2d aliases, vt2d, avt2d, aode structs, Symbol enum
                       NO Qt types — reusable by both 2D and 3D
ga_view_mark.hpp    ← Qt-specific rendering styles: pt2d_mark, ln2d_mark, vt2d_mark
                       includes QPen, QColor — stays in ga_view only
```

This change is purely additive and rename-based — no logic changes.

#### Change 2 — Fix active point ownership (do together with Change 1)

Replace the reference-into-vector with a safer pattern:

```cpp
// current (fragile):
pt2d& m_pos;   // reference into apt vector — dangling if vector reallocates

// proposed (safe):
std::vector<pt2d>* m_apt_vec;
size_t m_idx;
pt2d& pos() { return (*m_apt_vec)[m_idx]; }
```

Behavior is identical. The pointer-to-vector is stable across vector reallocations.

#### Change 3 — Extract ODE physics from Qt rendering (defer to when 3D ODE is needed)

Split `active_ode` into:

```cpp
// Pure physics — no Qt, no VTK, steppable by any timer source
class ODE_spring_system {
    void step(double dt);
    vec2dp mass_position() const;
    // ... state, RK4, params
};

// Qt wrapper — thin QGraphicsItem + QTimer driving the physics object
class active_ode : public QObject, public QGraphicsItem {
    ODE_spring_system m_physics;
    QTimer* m_timer;
};
```

The 3D equivalent reuses the same `ODE_spring_system` with a VTK/Qt timer.

### Do NOT attempt

A universal `Renderer` abstract base class spanning both Qt and VTK would be a significant
project in itself and provides marginal value. The 2D and 3D sides have separate `active_*`
implementations. The only shared layer that matters is the pure data model (Change 1).

### Proposed file structure after refactoring

```text
ga_view/src/
├── ga_view_types.hpp        ← NEW: pure GA+data types (no Qt), shared with ga_view3d
├── ga_view_mark.hpp         ← NEW: Qt-specific mark/style structs (QPen, QColor)
├── coordsys_model.hpp/cpp   ← adjusted: mark types come from ga_view_mark.hpp
├── active_pt2d.hpp/cpp      ← adjusted: index-based ownership (Change 2)
├── active_ode.hpp/cpp       ← unchanged for now (Change 3 deferred)
└── ... rest unchanged

ga_view3d/src/
├── ga_view3d_mark.hpp       ← NEW: 3D style structs (RGB float, opacity, width)
├── scene3d.hpp/cpp          ← NEW: VTK actor manager, analogous to w_Coordsys
└── active_pt3d.hpp/cpp      ← NEW: vtkHandleWidget wrapper (uses ga_view_types.hpp)
```

Both `ga_view` and `ga_view3d` include `ga_view_types.hpp` — the shared, Qt-free data layer.

## Preparatory Work Before Starting

### 1. VTK Dependency Decision (Critical — Do First)

VTK is not yet in `vcpkg.json`. Options on Windows:

| Option | Effort | Notes |
| --- | --- | --- |
| vcpkg `vtk` package | Medium | Large build (~30 min); needs Qt+OpenXR feature flags |
| Pre-built binary from vtk.org | Low | Fastest start; set `VTK_DIR` in CMake |
| System install (like Qt6/Lua pattern) | Low | Add detection to `cmake/vcpkg_dependencies.cmake` |

**Recommended**: Pre-built VTK binary for initial development, using same system-install
detection pattern as Lua in `cmake/vcpkg_dependencies.cmake`.

Required VTK modules:

```cmake
find_package(VTK REQUIRED COMPONENTS
    CommonCore CommonDataModel
    FiltersCore FiltersSources
    GUISupportQt
    InteractionStyle InteractionWidgets   # vtkHandleWidget etc.
    RenderingCore RenderingOpenGL2
    RenderingAnnotation                   # vtkAxesActor
    RenderingOpenXR                       # VR (Windows only, optional)
)
```

### 2. Build System Stub (Validate Before Implementation)

Create `ga_view3d/` as a new module alongside `ga_view/`:

```text
ga_view3d/
├── CMakeLists.txt    - finds VTK, links Qt6+VTK, compiles hello-world
└── src/
    └── main.cpp      - minimal QVTKOpenGLNativeWidget + vtkAxesActor
```

Add `add_subdirectory(ga_view3d)` to the root `CMakeLists.txt`, guarded with
`if(VTK_FOUND)` so the existing build is unaffected if VTK is not installed yet.

### 3. VTK + vcpkg.json (if vcpkg path is chosen instead of pre-built binary)

```json
{ "name": "vtk", "features": ["qt", "opengl"] }
```

Note: OpenXR support in vcpkg VTK may require a custom overlay or manual build flag.
Decide this only after the pre-built binary path has been validated.

---

## Instructions for Next Claude Instance

When the user continues this conversation:

1. They have an **existing Qt6 C++ 2D visualization framework** (`ga_view`)
2. They want to **migrate to 3D using VTK** — building their own application, not using
   ParaView
3. They need **desktop + VR support**: Windows/Pimax via OpenXR is the priority; macOS VR
   is deferred
4. They have **Python-using coworkers** — PyVista is the answer, no custom GA Python
   wrapper needed yet
5. Interactive elements (draggable points from `ga_view`) map to **`vtkHandleWidget`** in
   3D
6. Animations use a **Qt timer** calling `Render()` each frame
7. The toolkit decision is settled: **VTK**. Do not re-open this unless explicitly asked.
8. The `ga_view` architecture has been fully analysed — see "ga_view Architecture Reference"
   section above. The 2D→3D primitive mapping is documented there.
9. Preparatory work is documented — see "Preparatory Work Before Starting" section above.
   First task: decide VTK install method and validate it builds with Qt6+MSVC.

The user may ask about:

- Implementing the `ga_view3d` build system stub
- Choosing between pre-built VTK binary vs vcpkg for the dependency
- Implementing `Scene3D` class (VTK actor manager analogous to `w_Coordsys`)
- Migrating specific primitives (start with `vtkAxesActor` + `vtkHandleWidget`)
- Mapping GA types (Vec3d, BiVec3d, Motor) to VTK visual representations
- Animation of GA operations (rotors, motors, interpolation) via Qt timer
- VR mode activation (swap render window type)
- PyVista usage patterns for Python coworkers
