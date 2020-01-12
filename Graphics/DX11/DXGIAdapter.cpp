// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#include <Graphics/DX11/GTGraphicsDX11PCH.h>
#include <Graphics/DX11/DXGIAdapter.h>
using namespace gte;

DXGIAdapter::~DXGIAdapter()
{
    if (mAdapter)
    {
        mAdapter->Release();
    }
}

DXGIAdapter::DXGIAdapter(DXGIAdapter const& object)
    :
    mAdapter(nullptr)
{
    *this = object;
}

DXGIAdapter::DXGIAdapter(IDXGIAdapter1* adapter)
    :
    mAdapter(adapter)
{
    ZeroMemory(&mDescription, sizeof(DXGI_ADAPTER_DESC1));
    if (mAdapter)
    {
        DX11Log(mAdapter->GetDesc1(&mDescription));
        for (UINT i = 0; /**/; ++i)
        {
            IDXGIOutput* output = nullptr;
            HRESULT hr = mAdapter->EnumOutputs(i, &output);
            if (hr != DXGI_ERROR_NOT_FOUND)
            {
                mOutputs.push_back(DXGIOutput(output));
            }
            else
            {
                break;
            }
        }
    }
}

DXGIAdapter& DXGIAdapter::operator=(DXGIAdapter const& object)
{
    if (object.mAdapter)
    {
        object.mAdapter->AddRef();
    }

    DX11::SafeRelease(mAdapter);

    mAdapter = object.mAdapter;
    mDescription = object.mDescription;
    mOutputs = object.mOutputs;
    return *this;
}

void DXGIAdapter::Enumerate(std::vector<DXGIAdapter>& adapters)
{
    adapters.clear();

    IDXGIFactory1* factory = nullptr;
    DX11Log(CreateDXGIFactory1(__uuidof(IDXGIFactory1), (void**)& factory));

    for (UINT i = 0; /**/; ++i)
    {
        IDXGIAdapter1* adapter = nullptr;
        HRESULT hr = factory->EnumAdapters1(i, &adapter);
        if (hr != DXGI_ERROR_NOT_FOUND)
        {
            adapters.push_back(DXGIAdapter(adapter));
        }
        else
        {
            break;
        }
    }

    DX11::SafeRelease(factory);
}
