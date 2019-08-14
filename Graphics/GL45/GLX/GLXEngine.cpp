// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#include <Graphics/GL45/GTGraphicsGL45PCH.h>
#include <Applications/GLX/WindowSystem.h>
#include <Graphics/GL45/GLX/GLXEngine.h>

// TODO: This appears to be here to allow access to TheWindowSystem.
// Can this #include be avoided (and the explicit application
// dependency avoided) by passing a functor?
//#include <Applications/GLX/GteWindowSystem.h>

#include <X11/Xlib.h>
#include <GL/glx.h>
using namespace gte;

GLXEngine::~GLXEngine()
{
    Terminate();
}

GLXEngine::GLXEngine(Display* display, unsigned long window, GLXContext context,
    int xSize, int ySize, bool useDepth24Stencil8, bool saveDriverInfo, int requiredMajor, int requiredMinor)
    :
    GL45Engine(),
    mDisplay(display),
    mWindow(window),
    mImmediate(context),
    mIsComputeWindow(false)
{
    mXSize = xSize;
    mYSize = ySize;
    Initialize(requiredMajor, requiredMinor, useDepth24Stencil8, saveDriverInfo);
}

GLXEngine::GLXEngine(bool useDepth24Stencil8, bool saveDriverInfo, int requiredMajor, int requiredMinor)
    :
    GL45Engine(),
    mDisplay(nullptr),
    mWindow(0),
    mImmediate(nullptr),
    mIsComputeWindow(false)
{
    if (TheWindowSystem.Create(mDisplay, mImmediate, mWindow, useDepth24Stencil8))
    {
        mIsComputeWindow = true;
        Initialize(requiredMajor, requiredMinor, useDepth24Stencil8, saveDriverInfo);
    }
}

bool GLXEngine::IsActive() const
{
    return mImmediate == glXGetCurrentContext();
}

void GLXEngine::MakeActive()
{
    if (mImmediate != glXGetCurrentContext())
    {
        glXMakeCurrent(mDisplay, mWindow, mImmediate);
    }
}

void GLXEngine::DisplayColorBuffer(unsigned int syncInterval)
{
    // TODO: Disable vertical sync if possible.
    (void)syncInterval;

    glXSwapBuffers(mDisplay, mWindow);
}

bool GLXEngine::Initialize(int requiredMajor, int requiredMinor, bool useDepth24Stencil8, bool saveDriverInfo)
{
    if (!glXMakeCurrent(mDisplay, mWindow, mImmediate))
    {
        LogError("glXMakeCurrent failed.");
        glXDestroyContext(mDisplay, mImmediate);
        mDisplay = nullptr;
        mWindow = 0;
        mImmediate = nullptr;
        return false;
    }

    // Get the function pointers for OpenGL; initialize the viewport,
    // default global state, and default font.
    return GL45Engine::Initialize(requiredMajor, requiredMinor, useDepth24Stencil8, saveDriverInfo);
}

void GLXEngine::Terminate()
{
    GL45Engine::Terminate();

    if (mDisplay && mImmediate)
    {
        glXDestroyContext(mDisplay, mImmediate);
    }

    if (mIsComputeWindow)
    {
        XDestroyWindow(mDisplay, mWindow);
    }
}
