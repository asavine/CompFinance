#pragma once

#define EPS 1.0e-08
//  Utility for filling data
template<class T>
//  Returns filled vector 
//      has all original points
//      plus additional ones if requested
//      plus additional ones so maxDx is not exceeded
//  Original vector and addPOints must be sorted
//  Returned vector is sorted
inline vector<T>
fillData(
    //  The original data
    const vector<T>&                original,
    //  The maximum spacing allowed
    const T&                        maxDx,
    //  Specific points to add or nullptr if none
    const vector<T>*                addPoints = nullptr,
    //  Minimum distance for equality
    const T&                        minDx = T(0.0))
{
    //  Results
    vector<T> filled;

    //  Add points?
    vector<T> added;
    if (addPoints)
    {
        added.reserve(original.size() + addPoints->size());
        set_union(
            original.begin(),
            original.end(),
            addPoints->begin(),
            addPoints->end(),
            back_inserter(added),
            [minDx](const T x, const T y) { return x < y - minDx; });
    }
    const vector<T>& sequence = addPoints ? added : original;

    //  Position on the start, add it
    auto it = sequence.begin();
    filled.push_back(*it);
    ++it;

    while (it != sequence.end())
    {
        Time current = filled.back();
        Time next = *it;
        //  Must supplement?
        if (next - current > maxDx)
        {
            //  Number of points to add
            int addPoints = int((next - current) / maxDx - EPS) + 1;
            //  Spacing between supplementary points
            Time spacing = (next - current) / addPoints;
            //  Add the steps
            Time t = current + spacing;
            while (t < next)
            {
                filled.push_back(t);
                t += spacing;
            }
        }
        //  Push the next step on the product timeline and advance
        filled.push_back(*it);
        ++it;
    }

    return filled;
}
