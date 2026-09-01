#include <ttkSierpinskiKnoppWassersteinDistance.h>

#include <ttkPersistenceDiagramUtils.h>
#include <ttkMacros.h>

#include <vtkCellType.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

vtkStandardNewMacro(ttkSierpinskiKnoppWassersteinDistance);

namespace {

static void collectUnstructuredGrids(vtkDataObject *obj,
                                     std::vector<vtkUnstructuredGrid *> &out) {
  if(obj == nullptr) {
    return;
  }

  if(auto ug = vtkUnstructuredGrid::SafeDownCast(obj)) {
    out.emplace_back(ug);
    return;
  }

  if(auto mb = vtkMultiBlockDataSet::SafeDownCast(obj)) {
    const unsigned int nBlocks = mb->GetNumberOfBlocks();
    for(unsigned int i = 0; i < nBlocks; ++i) {
      collectUnstructuredGrids(mb->GetBlock(i), out);
    }
  }
}

} // namespace

ttkSierpinskiKnoppWassersteinDistance::ttkSierpinskiKnoppWassersteinDistance() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkSierpinskiKnoppWassersteinDistance::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    // Accept both vtkUnstructuredGrid and vtkMultiBlockDataSet. We declare the
    // generic vtkDataObject type here and restrict the ParaView UI in the XML.
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");

    // The second input is optional. If only the first input is connected and it
    // is a vtkMultiBlockDataSet containing at least two diagrams, the first two
    // vtkUnstructuredGrid blocks found in it are used.
    if(port == 1) {
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }

    return 1;
  }
  return 0;
}

int ttkSierpinskiKnoppWassersteinDistance::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkSierpinskiKnoppWassersteinDistance::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer tm{};

  std::vector<vtkUnstructuredGrid *> inputDiagrams{};

  vtkDataObject *input0 = vtkDataObject::GetData(inputVector[0], 0);
  collectUnstructuredGrids(input0, inputDiagrams);

  if(inputVector[1] != nullptr
     && inputVector[1]->GetNumberOfInformationObjects() > 0) {
    vtkDataObject *input1 = vtkDataObject::GetData(inputVector[1], 0);
    collectUnstructuredGrids(input1, inputDiagrams);
  }

  if(inputDiagrams.size() < 2) {
    this->printErr("SK_PROJECTION_ROBUST_BUILD: Expected either two "
                   "vtkUnstructuredGrid persistence diagrams on two input "
                   "ports, or one vtkMultiBlockDataSet containing at least "
                   "two vtkUnstructuredGrid diagrams. Found "
                   + std::to_string(inputDiagrams.size())
                   + " vtkUnstructuredGrid object(s).");
    return 0;
  }

  if(inputDiagrams.size() > 2) {
    this->printWrn("SK_PROJECTION_ROBUST_BUILD: More than two diagrams were "
                   "provided. Only the first two vtkUnstructuredGrid diagrams "
                   "will be projected.");
  }

  auto inputDiagram0 = inputDiagrams[0];
  auto inputDiagram1 = inputDiagrams[1];

  ttk::DiagramType diagram0{};
  ttk::DiagramType diagram1{};

  if(VTUToDiagram(diagram0, inputDiagram0, *this) < 0) {
    this->printErr("Could not read persistence diagram 0.");
    return 0;
  }

  if(VTUToDiagram(diagram1, inputDiagram1, *this) < 0) {
    this->printErr("Could not read persistence diagram 1.");
    return 0;
  }

  std::vector<ttk::SierpinskiKnoppWassersteinDistance::ProjectionPoint>
    projectedPoints{};

  if(this->execute(diagram0, diagram1, projectedPoints, this->L) != 1) {
    this->printErr("SK projection failed.");
    return 0;
  }

  auto output = vtkUnstructuredGrid::GetData(outputVector, 0);
  if(output == nullptr) {
    this->printErr("Invalid output.");
    return 0;
  }

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> grid{};

  const auto prefix = this->OutputArrayPrefix;
  const double nan = std::numeric_limits<double>::quiet_NaN();

  vtkNew<vtkDoubleArray> skCoordinate{};
  skCoordinate->SetName((prefix + "Coordinate").c_str());
  skCoordinate->SetNumberOfComponents(1);

  vtkNew<vtkIntArray> diagramId{};
  diagramId->SetName((prefix + "DiagramId").c_str());
  diagramId->SetNumberOfComponents(1);

  vtkNew<vtkIntArray> isDiagonal{};
  isDiagonal->SetName((prefix + "IsDiagonal").c_str());
  isDiagonal->SetNumberOfComponents(1);

  vtkNew<vtkIntArray> category{};
  category->SetName((prefix + "Category").c_str());
  category->SetNumberOfComponents(1);

  vtkNew<vtkIntArray> pairId{};
  pairId->SetName((prefix + "PairId").c_str());
  pairId->SetNumberOfComponents(1);

  vtkNew<vtkDoubleArray> birth{};
  birth->SetName((prefix + "Birth").c_str());
  birth->SetNumberOfComponents(1);

  vtkNew<vtkDoubleArray> death{};
  death->SetName((prefix + "Death").c_str());
  death->SetNumberOfComponents(1);

  vtkNew<vtkDoubleArray> persistence{};
  persistence->SetName((prefix + "Persistence").c_str());
  persistence->SetNumberOfComponents(1);

  vtkNew<vtkDoubleArray> trianglePoint{};
  trianglePoint->SetName((prefix + "TrianglePoint").c_str());
  trianglePoint->SetNumberOfComponents(2);

  const auto insertPointData = [&](const double t,
                                   const int did,
                                   const int diag,
                                   const int cat,
                                   const int pid,
                                   const double b,
                                   const double d,
                                   const double pers,
                                   const double tx,
                                   const double ty) {
    skCoordinate->InsertNextTuple1(t);
    diagramId->InsertNextTuple1(did);
    isDiagonal->InsertNextTuple1(diag);
    category->InsertNextTuple1(cat);
    pairId->InsertNextTuple1(pid);
    birth->InsertNextTuple1(b);
    death->InsertNextTuple1(d);
    persistence->InsertNextTuple1(pers);
    trianglePoint->InsertNextTuple2(tx, ty);
  };

  grid->SetPoints(points);

  if(this->AddLineSegment) {
    const vtkIdType p0 = points->InsertNextPoint(0.0, 0.0, 0.0);
    insertPointData(0.0, -1, -1, -1, -1, nan, nan, nan, nan, nan);

    const vtkIdType p1 = points->InsertNextPoint(1.0, 0.0, 0.0);
    insertPointData(1.0, -1, -1, -1, -1, nan, nan, nan, nan, nan);

    vtkIdType lineIds[2] = {p0, p1};
    grid->Allocate(static_cast<vtkIdType>(projectedPoints.size()) + 1);
    grid->InsertNextCell(VTK_LINE, 2, lineIds);
  } else {
    grid->Allocate(static_cast<vtkIdType>(projectedPoints.size()));
  }

  for(const auto &p : projectedPoints) {
    const vtkIdType id = points->InsertNextPoint(p.t, 0.0, 0.0);
    vtkIdType vertexId[1] = {id};
    grid->InsertNextCell(VTK_VERTEX, 1, vertexId);

    insertPointData(p.t,
                    p.diagramId,
                    p.isDiagonal,
                    p.category,
                    p.pairId,
                    p.birth,
                    p.death,
                    p.persistence,
                    p.x,
                    p.y);
  }

  grid->GetPointData()->AddArray(skCoordinate);
  grid->GetPointData()->AddArray(diagramId);
  grid->GetPointData()->AddArray(isDiagonal);
  grid->GetPointData()->AddArray(category);
  grid->GetPointData()->AddArray(pairId);
  grid->GetPointData()->AddArray(birth);
  grid->GetPointData()->AddArray(death);
  grid->GetPointData()->AddArray(persistence);
  grid->GetPointData()->AddArray(trianglePoint);
  grid->GetPointData()->SetActiveScalars((prefix + "Category").c_str());

  output->ShallowCopy(grid);

  this->printMsg("SK_PROJECTION_ROBUST_BUILD: Projected "
                   + std::to_string(projectedPoints.size())
                   + " SK atoms from two diagrams",
                 1.0, tm.getElapsedTime(), this->threadNumber_);
  return 1;
}
