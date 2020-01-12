// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

Texture2D<float4> stateTexture;
SamplerState stateSampler;

struct PS_INPUT
{
    float2 vertexTCoord : TEXCOORD0;
};

struct PS_OUTPUT
{
    float4 pixelColor : SV_TARGET0;
};

PS_OUTPUT PSMain(PS_INPUT input)
{
    // Map velocity channels to colors and modulate by density.
    PS_OUTPUT output;
    float4 current = stateTexture.Sample(stateSampler, input.vertexTCoord);
    float3 color = 0.5f + 0.5f * current.xyz / (1.0f + abs(current.xyz));
    output.pixelColor = float4(current.w * color, 1.0f);
    return output;
}
