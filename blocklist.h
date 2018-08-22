
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

#pragma once
#pragma warning(disable : 4018)

//  Blocklist data structure for AAD memory management
//  See chapter 10, unchanged with expression templates of chapter 15

#include <array>
#include <list>
#include <iterator>
using namespace std;

template <class T, size_t block_size>
class blocklist
{
    //  Container = list of blocks
    list<array<T, block_size>>  data;

    using list_iter = decltype(data.begin());
    using block_iter = decltype(data.back().begin());

    //  Current block
    list_iter           cur_block;

	//  Last block
	list_iter			last_block;
	
	//  Next free space in current block
    block_iter          next_space;

    //  Last free space (+1) in current block
    block_iter          last_space;

    //  Mark
    list_iter           marked_block;
    block_iter          marked_space;

    //  Create new array
    void newblock()
    {
        data.emplace_back();
        cur_block = last_block = prev(data.end());
        next_space = cur_block->begin();
        last_space = cur_block->end();
    }

    //  Move on to next array
    void nextblock()
    {
        //  This is the last array: create new
        if (cur_block == last_block)
        {
            newblock();
        }
        else
        {
            ++cur_block;
            next_space = cur_block->begin();
            last_space = cur_block->end();
        }
    }

public:

    //  Create first block on construction
    blocklist()
    {
        newblock();
    }

    //  Factory reset
    void clear()
    {
        data.clear();
        newblock();
    }

    //  Rewind but keep all blocks
    void rewind()
    {
        cur_block = data.begin();
        next_space = cur_block->begin();
        last_space = cur_block->end();
    }

	//	Memset
	void memset(unsigned char value = 0)
	{
		for (auto& arr : data)
		{
			std::memset(&arr[0], value, block_size * sizeof(T));
		}
	}

    //  Construct object of type T in place
    //      in the next free space and return a pointer on it
    //  Implements perfect forwarding of constructor arguments
    template<typename ...Args>
    T* emplace_back(Args&& ...args)
    {
        //  No more space in current array
        if (next_space == last_space)
        {
            nextblock();
        }
        //  Placement new, construct in memory pointed by next
        T* emplaced = new (&*next_space)    //  memory pointed by next as T*
            T(forward<Args>(args)...);      //  perfect forwarding of ctor arguments

        //  Advance next
        ++next_space;

        //  Return
        return emplaced;
    }

    //  Overload for default constructed
    T* emplace_back()
    {
        //  No more space in current array
        if (next_space == last_space)
        {
            nextblock();
        }

        //  Current space
        auto old_next = next_space;

        //  Advance next
        ++next_space;

        //  Return
        return &*old_next;
    }

	//  Stores n default constructed elements 
    //      and returns a pointer on the first

	//	Version 1: n known at compile time
	template <size_t n>
	T* emplace_back_multi()
	{
		//  No more space in current array
		if (distance(next_space, last_space) < n)
		{
			nextblock();
		}

		//  Current space
		auto old_next = next_space;

		//  Advance next
		next_space += n;

		//  Return
		return &*old_next;
	}

	//	Version 2: n unknown at compile time
	T* emplace_back_multi(const size_t n)
	{
		//  No more space in current array
		if (distance(next_space, last_space) < n)
		{
			nextblock();
		}

		//  Current space
		auto old_next = next_space;

		//  Advance next
		next_space += n;

		//  Return
		return &*old_next;
	}

	//	Marks

    //  Set mark
    void setmark()
    {
        marked_block = cur_block;
        marked_space = next_space;
    }

    //  Rewind to mark
    void rewind_to_mark()
    {
        cur_block = marked_block;
        next_space = marked_space;
		last_space = cur_block->end();
    }

    //  Iterator

    class iterator 
    {
        //  List and block
        list_iter        cur_block;		//  current block
        block_iter       cur_space;		//  current space
        block_iter       first_space;	//  first space in block
        block_iter       last_space;	//  last (+1) space in block

    public:

        //  iterator traits
        using difference_type = ptrdiff_t;
        using reference = T&;
        using pointer = T*;
        using value_type = T;
        using iterator_category = bidirectional_iterator_tag;

        //	Default constructor 
        iterator() {}

        //	Constructor 
        iterator(list_iter cb, block_iter cs, block_iter fs, block_iter ls) :
            cur_block(cb), cur_space(cs), first_space(fs), last_space(ls) {}

        //	Pre-increment (we do not provide post)
        iterator& operator++()
        {
            ++cur_space;
            if (cur_space == last_space)
            {
                ++cur_block;
                first_space = cur_block->begin();
                last_space = cur_block->end();
				cur_space = first_space;
            }

            return *this;
        }

        //	Pre-decrement 
        iterator& operator--()
        {
            if (cur_space == first_space)
            {
                --cur_block;
                first_space = cur_block->begin();
                last_space = cur_block->end();
				cur_space = last_space;
            }

            --cur_space;

            return *this;
        }

        //	Access to contained elements
        T& operator*()
        {
            return *cur_space;
        }
        const T& operator*() const
        {
            return *cur_space;
        }
        T* operator->()
        {
            return &*cur_space;
        }
        const T* operator->() const
        {
            return &*cur_space;
        }

        //	Check equality
        bool operator ==(const iterator& rhs) const
        {
            return (cur_block == rhs.cur_block && cur_space == rhs.cur_space);
        }
        bool operator !=(const iterator& rhs) const
        {
            return (cur_block != rhs.cur_block|| cur_space != rhs.cur_space);
        }
    };

    //  Access to iterators

    iterator begin()
    {
        return iterator(data.begin(), data.begin()->begin(),
            data.begin()->begin(), data.begin()->end());
    }

    iterator end()
    {
        auto last_block = prev(data.end());
        return iterator(cur_block, next_space,
            cur_block->begin(), cur_block->end());
    }

    //  Iterator on mark
    iterator mark()
    {
        return iterator(marked_block, marked_space,
            marked_block->begin(), marked_block->end());
    }

    //  Find element, by pointer, searching sequentially from the end
    iterator find(const T* const element)
    {
        //	Search from the end
        iterator it = end();
        iterator b = begin();

        while (it != b)
        {
            --it;
            if (&*it == element) return it;
        }

        if (&*it == element) return it;

        return end();
    }
};