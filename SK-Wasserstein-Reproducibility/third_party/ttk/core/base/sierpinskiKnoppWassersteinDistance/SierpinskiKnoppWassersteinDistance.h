/// \ingroup base
/// \class ttk::SierpinskiKnoppWassersteinDistance
/// \brief Projection of two normalized persistence diagrams onto the unit interval
/// through the SK curve.
///
/// This module projects the off-diagonal points of two normalized persistence
/// diagrams, together with their diagonal projections, onto the unit interval
/// through a finite-level Sierpinski--Knopp first-hit selector.

#pragma once

#include <Debug.h>
#include <PersistenceDiagramUtils.h>

#include <string>
#include <vector>

namespace ttk {

  class SierpinskiKnoppWassersteinDistance : virtual public Debug {

  public:
    struct ProjectionPoint {
      double t{}; // SK coordinate in [0,1]
      double x{}; // point represented in the persistence triangle
      double y{};
      double birth{}; // original pair birth
      double death{}; // original pair death
      double persistence{};
      int diagramId{}; // 0 or 1
      int isDiagonal{}; // 0: off-diagonal atom, 1: diagonal projection
      int pairId{}; // index of the original persistence pair in its diagram
      int category{}; // 2*diagramId + isDiagonal
    };

    SierpinskiKnoppWassersteinDistance();

    int execute(const DiagramType &diagram0,
                const DiagramType &diagram1,
                std::vector<ProjectionPoint> &output,
                const int L) const;
  };

} // namespace ttk
