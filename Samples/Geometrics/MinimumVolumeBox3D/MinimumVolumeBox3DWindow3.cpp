// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2020.09.05

#include "MinimumVolumeBox3DWindow3.h"
#include <Applications/Timer.h>
#include <Graphics/MeshFactory.h>
#include <Graphics/VertexColorEffect.h>
#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/MinimumVolumeBox3.h>
#include <iostream>
#include <random>

MinimumVolumeBox3DWindow3::MinimumVolumeBox3DWindow3(Parameters& parameters)
    :
    Window3(parameters),
    mVertices(NUM_POINTS)
{
    mWireState = std::make_shared<RasterizerState>();
    mWireState->cullMode = RasterizerState::CULL_NONE;
    mWireState->fillMode = RasterizerState::FILL_WIREFRAME;
    mEngine->SetRasterizerState(mWireState);

    CreateScene();
    InitializeCamera(60.0f, GetAspectRatio(), 0.1f, 100.0f, 0.001f, 0.1f,
        { 0.0f, 0.0f, -2.0f }, { 0.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 0.0f });
    mPVWMatrices.Update();
}

void MinimumVolumeBox3DWindow3::OnIdle()
{
    mTimer.Measure();

    if (mCameraRig.Move())
    {
        mPVWMatrices.Update();
    }

    mEngine->ClearBuffers();
    mEngine->Draw(mPoints);
    mEngine->Draw(mPolytope);
    mEngine->Draw(mBoxMesh);
    mEngine->Draw(8, mYSize - 8, { 0.0f, 0.0f, 0.0f, 1.0f }, mTimer.GetFPS());
    mEngine->DisplayColorBuffer(0);

    mTimer.UpdateFrameCount();
}

void MinimumVolumeBox3DWindow3::CreateScene()
{
    mScene = std::make_shared<Node>();

    std::mt19937 mte;
    std::uniform_real_distribution<float> rnd(-1.0f, 1.0f);
    Vector3<float> center{ 0.0f, 0.0f, 0.0f };
    Vector3<float> extent{ 1.0f, 0.25f, 0.125f };
    Vector3<float> axis[3] = {
        { 1.0f, 1.0f, 0.0f },
        { -1.0f, 1.0f, 0.0f },
        { 0.0f, 0.0f, 1.0f }
    };
    Normalize(axis[0]);
    Normalize(axis[1]);
    Normalize(axis[2]);
    for (auto& v : mVertices)
    {
        float theta = rnd(mte) * (float)GTE_C_TWO_PI;
        float phi = rnd(mte) * (float)GTE_C_PI;
        float radius = 0.5f * (rnd(mte) + 1.0f);
        float x = extent[0] * std::cos(theta) * std::sin(phi);
        float y = extent[1] * std::sin(theta) * std::sin(phi);
        float z = extent[2] * std::cos(phi);
        v = center + radius * (x * axis[0] + y * axis[1] + z * axis[2]);
    }

    struct Vertex
    {
        Vector3<float> position;
        Vector4<float> color;
    };
    VertexFormat vformat;
    vformat.Bind(VA_POSITION, DF_R32G32B32_FLOAT, 0);
    vformat.Bind(VA_COLOR, DF_R32G32B32A32_FLOAT, 0);
    auto vbuffer = std::make_shared<VertexBuffer>(vformat, NUM_POINTS);
    Vertex* vertex = vbuffer->Get<Vertex>();
    for (int i = 0; i < NUM_POINTS; ++i)
    {
        vertex[i].position[0] = (float)mVertices[i][0];
        vertex[i].position[1] = (float)mVertices[i][1];
        vertex[i].position[2] = (float)mVertices[i][2];
        vertex[i].color[0] = 0.5f * (rnd(mte) + 1.0f);
        vertex[i].color[1] = 0.5f * (rnd(mte) + 1.0f);
        vertex[i].color[2] = 0.5f * (rnd(mte) + 1.0f);
        vertex[i].color[3] = 1.0f;
    }

    auto ibuffer = std::make_shared<IndexBuffer>(IP_POLYPOINT, NUM_POINTS);
    auto effect = std::make_shared<VertexColorEffect>(mProgramFactory);

    mPoints = std::make_shared<Visual>(vbuffer, ibuffer, effect);
    mPVWMatrices.Subscribe(mPoints->worldTransform, effect->GetPVWMatrixConstant());
    mScene->AttachChild(mPoints);

    // Compute the minimum-volume box. The convex hull of the vertices is
    // computed internally.
    uint32_t const numThreads = 4;
    size_t const lgMaxSample = 5;
    MinimumVolumeBox3<float, false> mvb3(numThreads);
    OrientedBox3<float> minBox;
    float volume;
    mvb3(NUM_POINTS, &mVertices[0], lgMaxSample, minBox, volume);

    // Recompute the convex hull for visualization.
    ConvexHull3<float> ch3;
    ch3(NUM_POINTS, &mVertices[0], 0.0f);
    std::vector<TriangleKey<true>> const& triangles = ch3.GetHullUnordered();
    int const* indices = static_cast<int const*>(&triangles[0].V[0]);
    ibuffer = std::make_shared<IndexBuffer>(IP_TRIMESH,
        static_cast<uint32_t>(triangles.size()), sizeof(int));
    std::memcpy(ibuffer->GetData(), indices, ibuffer->GetNumBytes());
    mPolytope = std::make_shared<Visual>(vbuffer, ibuffer, effect);
    mScene->AttachChild(mPolytope);

    MeshFactory mf;
    mf.SetVertexFormat(vformat);
    mBoxMesh = mf.CreateBox(1.0f, 1.0f, 1.0f);
    vbuffer = mBoxMesh->GetVertexBuffer();
    vertex = vbuffer->Get<Vertex>();
    std::array<Vector3<float>, 8> corner;
    minBox.GetVertices(corner);
    for (int i = 0; i < 8; ++i)
    {
        vertex[i].position[0] = corner[i][0];
        vertex[i].position[1] = corner[i][1];
        vertex[i].position[2] = corner[i][2];
        vertex[i].color[0] = 0.5f * (rnd(mte) + 1.0f);
        vertex[i].color[1] = 0.5f * (rnd(mte) + 1.0f);
        vertex[i].color[2] = 0.5f * (rnd(mte) + 1.0f);
        vertex[i].color[3] = 1.0f;
    }
    mBoxMesh->SetEffect(effect);
    mScene->AttachChild(mBoxMesh);

    mTrackBall.Attach(mScene);
    mTrackBall.Update();
}

