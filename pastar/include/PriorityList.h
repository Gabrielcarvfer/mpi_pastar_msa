/*!
 * \class PriorityList
 * \author Daniel Sundfeld
 * \copyright MIT License
 *
 * \brief Class that implements a Best-First Search OpenList
 */
#ifndef _PRIORITY_LIST
#define _PRIORITY_LIST
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <limits>

#include "Coord.h"
#include "Node.h"


#include <fstream>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/set.hpp>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#ifdef NO_LIB_BOOST
#error This version requires Libboost. Please checkout an older version of the program.
#endif

using boost::multi_index_container;
using namespace boost::multi_index;

struct pos{};
struct priority{};

template <int N>
struct openlist_multiindex {
    typedef boost::multi_index_container<
        Node<N>,
        indexed_by<
            ordered_unique<
                  tag<pos>,  BOOST_MULTI_INDEX_MEMBER(Node<N>,Coord<N>,pos)>,
            ordered_non_unique<
                  tag<priority>, BOOST_MULTI_INDEX_MEMBER(Node<N>,int,m_f)> >
    > type;
};

template <int N>
class PriorityList {
    public:
        bool dequeue(Node<N> &n);
        bool enqueue(const Node<N> &n);
        bool conditional_enqueue(const Node<N> &n);
        int get_highest_priority() const;
		Node<N> get_highest_priority_node(void) const;

        typename openlist_multiindex<N>::type::iterator find(const Coord<N> &c) const { return m_openlist.find(c); };
        typename openlist_multiindex<N>::type::iterator end() const { return m_openlist.end(); };
        typename openlist_multiindex<N>::type::size_type erase(const Coord<N> &c) { return m_openlist.erase(c); };
        size_t size(void) { return m_openlist.size(); };

        bool empty() { return m_openlist.empty(); };

        void merge(const PriorityList& other);

    private:
        typename openlist_multiindex<N>::type m_openlist;

        friend class boost::serialization::access;
        // When the class Archive corresponds to an output archive, the
        // & operator is defined similar to <<.  Likewise, when the class Archive
        // is a type of input archive the & operator is defined similar to >>.
        template<class Archive>
        void serialize(Archive & ar, unsigned int version)
        {
            ar & m_openlist;
        }
};

template <int N>
bool PriorityList<N>::dequeue(Node<N> &n)
{
    auto it = get<priority>(m_openlist).begin();
    if (it == get<priority>(m_openlist).end())
        return false;
    n = *it;
    get<priority>(m_openlist).erase(it);
    return true;
}

template <int N>
bool PriorityList<N>::enqueue(const Node<N> &c)
{
    auto it = m_openlist.find(c.pos);
    if (it == m_openlist.end())
        return m_openlist.insert(c).second;
    return m_openlist.modify(it, change_node<N>(c.get_f(), c.get_g(), c.get_parenti()));
}

template <int N>
bool PriorityList<N>::conditional_enqueue(const Node<N> &c)
{
    auto it = m_openlist.find(c.pos);
    if (it == m_openlist.end())
        return m_openlist.insert(c).second;
    if (c.get_f() >= it->get_f())
        return true;
    return m_openlist.modify(it, change_node<N>(c.get_f(), c.get_g(), c.get_parenti()));
}

template <int N>
int PriorityList<N>::get_highest_priority(void) const
{
    auto it = get<priority>(m_openlist).begin();
    if (it == get<priority>(m_openlist).end())
        return std::numeric_limits<int>::max();
    return it->get_f();
}

template <int N>
Node<N> PriorityList<N>::get_highest_priority_node(void) const
{
	auto it = get<priority>(m_openlist).begin();
	return *it;
}

template <int N>
void PriorityList<N>::merge(const PriorityList<N>& other)
{
    auto it = get<priority>(other.m_openlist).begin();
    while (it != get<priority>(other.m_openlist).end())
    {
        conditional_enqueue(*it);
        ++it;
    }
    return;
}
#endif
