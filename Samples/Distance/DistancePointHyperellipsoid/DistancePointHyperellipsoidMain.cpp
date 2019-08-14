// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#include "DistancePointHyperellipsoidConsole.h"
#include <Applications/LogReporter.h>
using namespace gte;

int main(int, char const*[])
{
#if defined(_DEBUG)
    LogReporter reporter(
        "LogReport.txt",
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL);
#endif

    Console::Parameters parameters(L"DistancePointHyperellipsoidConsole");
    auto console = TheConsoleSystem.Create<DistancePointHyperellipsoidConsole>(parameters);
    TheConsoleSystem.Execute(console);
    TheConsoleSystem.Destroy(console);
    return 0;
}
