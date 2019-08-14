// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <vector>

// An unordered set of objects stored in contiguous memory.  The type T must
// have the following member functions:
//   T::T();
//   T::~T();
//   T& operator=(const T&);
//   bool operator==(const T&) const;
// This class exists only to for the ported Wild Magic 5 sample application;
// std::unordered_set did not exist at the time of creation of that sample.
// UnorderedSet has specialized properties, mainly that the storage is
// contiguous.  At a later date the sample will be revised to use GTEngine
// manifold data structures and (if necessary) std::unordered_set.

template <typename T>
class UnorderedSet
{
public:
    enum
    {
        DEFAULT_GROW = 8
    };

    UnorderedSet(int maxNumElements = 0, int grow = 0)
        :
        mMaxNumElements(maxNumElements > 0 ? maxNumElements : DEFAULT_GROW),
        mGrow(grow > 0 ? grow : DEFAULT_GROW),
        mNumElements(0),
        mElements(mMaxNumElements)
    {
    }

    ~UnorderedSet() = default;

    void Reset(int maxNumElements = 0, int grow = 0)
    {
        mMaxNumElements = (maxNumElements > 0 ? maxNumElements : DEFAULT_GROW);
        mGrow = (grow > 0 ? grow : DEFAULT_GROW);
        mNumElements = 0;
        mElements.resize(mMaxNumElements);
    }

    void Clear()
    {
        mNumElements = 0;
    }

    inline int GetMaxNumElements() const
    {
        return mMaxNumElements;
    }

    inline int GetGrow() const
    {
        return mGrow;
    }

    inline int GetNumElements() const
    {
        return mNumElements;
    }

    T const& operator[](int i) const
    {
        return mElements[i];
    }

    T& operator[] (int i)
    {
        return mElements[i];
    }

    int Find(T const& element) const
    {
        for (int i = 0; i < mNumElements; ++i)
        {
            if (element == mElements[i])
            {
                return i;
            }
        }
        return -1;
    }

    bool Insert(T const& element)
    {
        for (int i = 0; i < mNumElements; ++i)
        {
            if (element == mElements[i])
            {
                return false;
            }
        }
        GrowArray();
        mElements[mNumElements++] = element;
        return true;
    }

    int Append(T const& element)
    {
        GrowArray();
        int location = mNumElements++;
        mElements[location] = element;
        return location;
    }

    bool Remove(T const& element, int* oldIndex = nullptr, int* newIndex = nullptr)
    {
        for (int i = 0; i < mNumElements; i++)
        {
            if (element == mElements[i])
            {
                RemoveElement(i, oldIndex, newIndex);
                return true;
            }
        }

        if (oldIndex)
        {
            *oldIndex = -1;
        }
        if (newIndex)
        {
            *newIndex = -1;
        }
        return false;
    }

    bool RemoveAt(int i, int* oldIndex = nullptr, int* newIndex = nullptr)
    {
        if (0 <= i && i < mNumElements)
        {
            RemoveElement(i, oldIndex, newIndex);
            return true;
        }

        if (oldIndex)
        {
            *oldIndex = -1;
        }
        if (newIndex)
        {
            *newIndex = -1;
        }
        return false;
    }

protected:
    void GrowArray()
    {
        if (mNumElements == mMaxNumElements)
        {
            int newMaxNumElements = mMaxNumElements + mGrow;
            mElements.resize(newMaxNumElements);
            mMaxNumElements = newMaxNumElements;
        }
    }

    // This function is called only when mNumElements is positive, so it is
    // valid to decrement mNumElements.  The members of mElements must remain
    // contiguous, so on removal, the last member of mElements is copied to
    // the location vacated at index i.  On return, oldIndex stores the index
    // for the last member of mElements before the removal.  If i does not
    // point to the last member of mElements before the removal, newIndex
    // stores the index i to which the last member has been copied.  If i does
    // point to the last member, no copy occurs and newIndex is set to -1.
    void RemoveElement(int i, int* oldIndex, int* newIndex)
    {
        --mNumElements;
        if (oldIndex)
        {
            *oldIndex = mNumElements;
        }

        if (i != mNumElements)
        {
            mElements[i] = mElements[mNumElements];
            if (newIndex)
            {
                *newIndex = i;
            }
        }
        else
        {
            if (newIndex)
            {
                *newIndex = -1;
            }
        }
    }

    int mMaxNumElements, mGrow, mNumElements;
    std::vector<T> mElements;
};
