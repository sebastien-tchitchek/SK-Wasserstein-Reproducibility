/// \ingroup vtk
/// \class ttkSierpinskiKnoppWassersteinDistance
/// \brief VTK wrapper for projecting two normalized persistence diagrams onto the 
/// unit interval through the SK curve.
///
/// Input 0 and Input 1 must be persistence diagrams encoded as
/// vtkUnstructuredGrid objects. The output is a vtkUnstructuredGrid containing
/// one vertex for every off-diagonal point and one vertex for its diagonal
/// projection. Point-data arrays identify the source diagram and whether the
/// output point is off-diagonal or diagonal.

#pragma once

// VTK Module
#include <ttkSierpinskiKnoppWassersteinDistanceModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <SierpinskiKnoppWassersteinDistance.h>

#include <string>

class TTKSIERPINSKIKNOPPWASSERSTEINDISTANCE_EXPORT ttkSierpinskiKnoppWassersteinDistance
  : public ttkAlgorithm,
    protected ttk::SierpinskiKnoppWassersteinDistance {

private:
  int L{30};
  bool AddLineSegment{true};
  std::string OutputArrayPrefix{"SK"};

public:
  static ttkSierpinskiKnoppWassersteinDistance *New();
  vtkTypeMacro(ttkSierpinskiKnoppWassersteinDistance, ttkAlgorithm);

  vtkSetMacro(L, int);
  vtkGetMacro(L, int);

  vtkSetMacro(AddLineSegment, bool);
  vtkGetMacro(AddLineSegment, bool);
  vtkBooleanMacro(AddLineSegment, bool);

  vtkSetMacro(OutputArrayPrefix, const std::string &);
  vtkGetMacro(OutputArrayPrefix, std::string);

protected:
  ttkSierpinskiKnoppWassersteinDistance();
  ~ttkSierpinskiKnoppWassersteinDistance() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
