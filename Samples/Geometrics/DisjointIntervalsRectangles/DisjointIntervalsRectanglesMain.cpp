// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2020.01.10

#include "DisjointIntervalsRectanglesConsole.h"
#include <Applications/LogReporter.h>
using namespace gte;

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

    Console::Parameters parameters(L"DisjointIntervalsRectanglesConsole");
    auto console = TheConsoleSystem.Create<DisjointIntervalsRectanglesConsole>(parameters);
    TheConsoleSystem.Execute(console);
    TheConsoleSystem.Destroy(console);
    return 0;
}
