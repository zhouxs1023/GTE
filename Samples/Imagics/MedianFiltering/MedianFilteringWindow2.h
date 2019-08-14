// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <Applications/Window2.h>
using namespace gte;

class MedianFilteringWindow2 : public Window2
{
public:
    MedianFilteringWindow2(Parameters& parameters);

    virtual void OnIdle() override;
    virtual bool OnCharPress(unsigned char key, int x, int y) override;

private:
    bool SetEnvironment();
    bool CreatePrograms(unsigned int txWidth, unsigned int txHeight);

    std::shared_ptr<Texture2> mOriginal;
    std::shared_ptr<Texture2> mImage[2];
    std::shared_ptr<OverlayEffect> mOverlay[2];

    // 0 = median 3x3 by insertion sort
    // 1 = median 3x3 by min-max
    // 2 = median 5x5 by insertion sort
    // 3 = median 5x5 by min-max
    unsigned int mSelection;
    std::shared_ptr<ComputeProgram> mMedianProgram[4];
    std::shared_ptr<ComputeProgram> mCProgram;
    unsigned int mNumXGroups, mNumYGroups;
    static std::string msName[4];
};
