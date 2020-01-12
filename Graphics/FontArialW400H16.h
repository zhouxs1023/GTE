// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <Graphics/Font.h>

namespace gte
{
    class FontArialW400H16 : public Font
    {
    public:
        virtual ~FontArialW400H16() = default;
        FontArialW400H16(std::shared_ptr<ProgramFactory> const& factory, int maxMessageLength);

    private:
        static int msWidth;
        static int msHeight;
        static unsigned char msTexels[];
        static float msCharacterData[];
    };
}
