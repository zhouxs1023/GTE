// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2020.01.10

#include "GpuGaussianBlur2Window2.h"
#include <Applications/LogReporter.h>
#include <Applications/Command.h>

int main(int numArguments, char* arguments[])
{
#if defined(_DEBUG)
    LogReporter reporter(
        "LogReport.txt",
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL);
#endif

    Command command(numArguments, arguments);
    bool useDirichlet = (command.GetBoolean("d") > 0 ? true : false);

    // The window size is that of the Head_U16_X256_Y256.binary image.
    GpuGaussianBlur2Window2::Parameters parameters(L"GpuGaussianBlur2Window2",
        0, 0, 256, 256, useDirichlet);
    auto window = TheWindowSystem.Create<GpuGaussianBlur2Window2>(parameters);
    TheWindowSystem.MessagePump(window, TheWindowSystem.DEFAULT_ACTION);
    TheWindowSystem.Destroy(window);
    return 0;
}
