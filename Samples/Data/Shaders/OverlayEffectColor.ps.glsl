// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

uniform sampler2D imageSampler;

layout(location = 0) in vec2 vertexTCoord;
layout(location = 0) out vec4 pixelColor;

void main()
{
    pixelColor = texture(imageSampler, vertexTCoord);
}
