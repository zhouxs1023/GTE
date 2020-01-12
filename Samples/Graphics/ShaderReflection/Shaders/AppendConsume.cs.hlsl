// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

struct Particle
{
    float3 position;
    float3 color;
    float size;
};

ConsumeStructuredBuffer<Particle> particlesIn;
AppendStructuredBuffer<Particle> particlesOut;

[numthreads(1, 1, 1)]
void CSMain(int t : SV_DispatchThreadID)
{
    Particle p = particlesIn.Consume();
    particlesOut.Append(p);
}
