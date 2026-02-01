# VTK 3D Visualization Migration - Conversation Context

## User Background & Requirements

### Current Setup

- **Existing Framework**: 2D visualization framework written in C++ based on Qt
- **Current Capabilities**:
  - Display coordinate systems
  - Render 2D primitives
  - User interaction via mouse

### Project Goals

- Migrate to 3D visualization using C++ framework
- Support both on-screen (desktop) and VR visualization
- Preferably use the same codebase for both modes
- Focus on scientific visualization with interaction capabilities

### Initial Requirements

- Draw 3D coordinate systems
- Visualize intersecting planes with visible intersection lines
- Interactive capabilities (mouse for desktop, VR controllers for headset)
- Integrate with existing Qt application

## Recommended Solution: VTK (Visualization Toolkit)

### Why VTK?

- Purpose-built for scientific/medical/engineering visualization
- Excellent C++ API with comprehensive interaction support
- Native mouse and VR controller input handling
- **Same visualization pipeline works for both desktop and VR**
- Direct Qt integration via `QVTKOpenGLNativeWidget`
- Native VR support through OpenXR (no external frameworks like ParaView needed)

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

set(CMAKE_CXX_STANDARD 17)
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
    RenderingOpenXR  # For VR support
)

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)

add_executable(${PROJECT_NAME} 
    main.cpp
    Scene3D.h
)

target_link_libraries(${PROJECT_NAME} 
    ${VTK_LIBRARIES}
    Qt5::Core
    Qt5::Widgets
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
  - SteamVR (most common, works with most headsets)
  - Oculus Runtime (for Meta Quest Link)
  - Windows Mixed Reality (for WMR headsets)
- **VR Headset** connected to PC
- Compatible headsets: Meta Quest, Valve Index, HP Reverb, HTC Vive, etc.

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

## Next Steps for User

The user will want to:

1. **Share their existing Qt C++ 2D framework code** for integration advice
2. **Discuss specific 2D primitives** they need to migrate
3. **Plan the migration strategy** based on their codebase
4. **Potentially add custom interaction** (object picking, measurement tools)
5. **Optimize the scene setup** for their specific scientific visualization needs

## Key Advantages of This Approach

✅ **Single Codebase**: Same scene definition for desktop and VR  
✅ **Qt Integration**: Works seamlessly with existing Qt applications  
✅ **Native VTK**: No external dependencies like ParaView needed  
✅ **OpenXR Standard**: Compatible with all major VR headsets  
✅ **Scientific Focus**: VTK is designed for scientific visualization  
✅ **Interaction Built-in**: Mouse and VR controller support automatic  
✅ **Gradual Migration**: Can keep existing 2D code while adding 3D  

## Resources

- **VTK Documentation**: <https://vtk.org/documentation/>
- **VTK Examples**: <https://examples.vtk.org/site/Cxx/>
- **Qt + VTK**: Look for examples under "GUI/Qt" section
- **OpenXR**: <https://www.khronos.org/openxr/>
- **VTK Discourse Forum**: <https://discourse.vtk.org/>

---

## Instructions for Next Claude Instance

When the user continues this conversation:

1. They have an **existing Qt C++ 2D visualization framework**
2. They want to **migrate to 3D using VTK**
3. They need **both desktop and VR support**
4. Focus is on **scientific visualization** with **interaction**
5. The complete working example above is the starting point
6. Next step: Review their existing code and plan integration strategy
7. Be prepared to discuss specific Qt integration patterns, custom interactions, and scene
   optimization

The user may ask about:

- Integrating VTK with their specific Qt widgets/layouts
- Migrating specific 2D primitives (coordinate systems, lines, etc.)
- Custom interaction handlers (picking, measurements)
- Performance optimization
- VR controller customization
- Multi-viewport setups
- Data visualization pipelines specific to their domain
