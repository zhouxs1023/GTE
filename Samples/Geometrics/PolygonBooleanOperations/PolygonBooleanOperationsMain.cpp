// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2020.01.10

#include "PolygonBooleanOperationsWindow2.h"
#include <Applications/LogReporter.h>

int main()
{
#if defined(_DEBUG)
    LogReporter reporter(
        "LogReport.txt",
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL);
#endif

    Window::Parameters parameters(L"PolygonBooleanOperationsWindow2", 0, 0, 512, 512);
    auto window = TheWindowSystem.Create<PolygonBooleanOperationsWindow2>(parameters);
    TheWindowSystem.MessagePump(window, TheWindowSystem.NO_IDLE_LOOP);
    TheWindowSystem.Destroy(window);
    return 0;
}
