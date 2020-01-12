// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

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

    // Choose the number of threads to use.  The default constructor for
    // MinimumVolumeBox3 uses a default of 1, in which case all computations
    // are on the main thread.  The timings below are for a 64-bit release
    // build (no debugger attached) on Intel Core i7-3930K CPUs running at
    // 3.20 GHz.
    unsigned int numThreads = 1;

#if 0
    // Compute the convex hull internally using arbitrary precision
    // arithmetic.  This is slower than computing the hull explicitly using
    // the maximum fixed precision; see the other conditional block of code.
    Timer timer;
    typedef BSRational<UIntegerAP32> MVBRational;
    MinimumVolumeBox3<float, MVBRational> mvb3(numThreads);
    OrientedBox3<float> minBox = mvb3(NUM_POINTS, &mVertices[0]);
    std::cout << "mvb3 seconds = " << timer.GetSeconds() << std::endl;
    // numThreads = 1, seconds = 7.09
    // numThreads = 2, seconds = 6.22
#else
    // If mVertices were to use 'double', you would need the template type
    // UIntegerFP32<167> to compute the convex hull.
    Timer timer;
    typedef BSNumber<UIntegerFP32<27>> CHRational;
    ConvexHull3<float, CHRational> ch3(numThreads);
    ch3(NUM_POINTS, &mVertices[0], 0.0f);
    std::vector<TriangleKey<true>> const& triangles = ch3.GetHullUnordered();
    int const numIndices = static_cast<int>(3 * triangles.size());
    int const* indices = static_cast<int const*>(&triangles[0].V[0]);
    typedef BSRational<UIntegerAP32> MVBRational;
    MinimumVolumeBox3<float, MVBRational> mvb3(numThreads);
    OrientedBox3<float> minBox = mvb3(NUM_POINTS, &mVertices[0], numIndices, indices);
    std::cout << "mvb3 seconds = " << timer.GetSeconds() << std::endl;
    // numThreads = 1, seconds = 2.69
    // numThreads = 2, seconds = 2.01
#endif

    std::vector<int> const& hull = mvb3.GetHull();
    ibuffer = std::make_shared<IndexBuffer>(IP_TRIMESH, static_cast<int>(hull.size() / 3), sizeof(int));
    std::memcpy(ibuffer->GetData(), &hull[0], ibuffer->GetNumBytes());
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

