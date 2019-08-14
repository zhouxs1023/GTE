// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <Applications/Window3.h>
#include <Mathematics/NURBSSphere.h>
using namespace gte;

class NURBSSphereWindow3 : public Window3
{
public:
    NURBSSphereWindow3(Parameters& parameters);

    virtual void OnIdle() override;
    virtual bool OnCharPress(unsigned char key, int x, int y) override;

private:
    void CreateScene();
    void CreateEighthSphere();
    void CreateHalfSphere();
    void CreateFullSphere();

    std::shared_ptr<RasterizerState> mNoCullSolidState;
    std::shared_ptr<RasterizerState> mNoCullWireState;
    NURBSEighthSphereDegree4<float> mEighthSphere;
    std::shared_ptr<Visual> mEighthSphereVisual;
    NURBSHalfSphereDegree3<float> mHalfSphere;
    std::shared_ptr<Visual> mHalfSphereVisual;
    NURBSFullSphereDegree3<float> mFullSphere;
    std::shared_ptr<Visual> mFullSphereVisual;

    std::shared_ptr<Visual> mCurrentVisual;
};
