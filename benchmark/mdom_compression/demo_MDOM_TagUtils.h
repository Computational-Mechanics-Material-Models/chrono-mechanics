// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Utility helpers for multidomain demos in branches where ChArchiveSetUniqueTags
// is not available.
// =============================================================================

#ifndef DEMO_MDOM_TAGUTILS_H
#define DEMO_MDOM_TAGUTILS_H

#include "chrono/fea/ChMesh.h"
#include "chrono/physics/ChSystem.h"

namespace chrono {
namespace multidomain {

inline void SetUniqueTags(ChSystem& sys, int start_tag = 1) {
    int tag = start_tag;

    const_cast<ChAssembly&>(sys.GetAssembly()).SetTag(tag++);

    for (auto& body : sys.GetBodies())
        body->SetTag(tag++);

    for (auto& link : sys.GetLinks())
        link->SetTag(tag++);

    for (auto& mesh : sys.GetMeshes()) {
        mesh->SetTag(tag++);
        for (auto& node : mesh->GetNodes())
            node->SetTag(tag++);
        for (auto& element : mesh->GetElements())
            element->SetTag(tag++);
    }

    for (auto& item : sys.GetOtherPhysicsItems())
        item->SetTag(tag++);
}

}  // namespace multidomain
}  // namespace chrono

#endif
