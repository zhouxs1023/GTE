// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

struct PS_INPUT
{
    float3 vertexColor : COLOR0;
    noperspective float2 planeConstant : TEXCOORD0;
};

struct PS_OUTPUT
{
    float4 pixelColor : SV_TARGET0;
    float2 planeConstant : SV_TARGET1;
};

PS_OUTPUT PSMain(PS_INPUT input)
{
    PS_OUTPUT output;
    output.pixelColor = float4(input.vertexColor, 1.0f);
    output.planeConstant = input.planeConstant;
    return output;
}
