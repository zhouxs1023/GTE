// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <Applications/Window3.h>
#include <Graphics/BillboardNode.h>
using namespace gte;

#define DEMONSTRATE_VIEWPORT_BOUNDING_RECTANGLE
#define DEMONSTRATE_POST_PROJECTION_REFLECTION

class BillboardNodesWindow3 : public Window3
{
public:
    BillboardNodesWindow3(Parameters& parameters);

    virtual void OnIdle() override;
    virtual bool OnCharPress(unsigned char key, int x, int y) override;
    virtual bool OnMouseMotion(MouseButton button, int x, int y, unsigned int modifiers) override;

private:
    bool SetEnvironment();
    void CreateScene();

    // All triangle meshes have this common vertex format.
    struct Vertex
    {
        Vector3<float> position;
        Vector2<float> tcoord;
    };

    Culler mCuller;
    std::shared_ptr<Node> mScene;
    std::shared_ptr<Texture2> mGroundTexture;
    std::shared_ptr<Texture2> mSkyTexture;

    // Billboard 0 has a rectangle attached.  Billboard 1 has a torus
    // attached.
    std::shared_ptr<Visual> mGround, mRectangle, mTorus;
    std::shared_ptr<BillboardNode> mBillboard0, mBillboard1;

#if defined(DEMONSTRATE_VIEWPORT_BOUNDING_RECTANGLE)
    // Compute the bounding rectangle in normalized display coordinates
    // [-1,1]^2 for the torus.
    void ComputeTorusBoundingRectangle();

    std::shared_ptr<BlendState> mBlendState;
    std::shared_ptr<OverlayEffect> mOverlay;
    std::shared_ptr<RasterizerState> mNoCullState;
#endif

#if defined(DEMONSTRATE_POST_PROJECTION_REFLECTION)
    std::shared_ptr<RasterizerState> mCullCWState;
#endif
};
