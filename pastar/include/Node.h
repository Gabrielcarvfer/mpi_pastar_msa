/*!
 * \class Node
 * \author Daniel Sundfeld
 * \copyright MIT License
 *
 * \brief Class that hold all Nodes atributes, like the cost from the
 * origin, heuristic estimative, parent
 */
#ifndef _NODE_H
#define _NODE_H
#include <vector>
#include "Coord.h"

#include <fstream>
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

template < int N > class Node;
template < int N > std::ostream& operator<< (std::ostream &lhs, const Node<N> &rhs);
template < int N > struct change_node;

template < int N >
class Node
{
    friend struct change_node< N >;
    public:
        Coord<N> pos; //!< Multidimensional coordinate of the node
        int m_f; //!< priority
        //bool remote;
        Node();
        Node(const int g, const Coord<N> &pos, const int &parenti);
        friend std::ostream &operator<< <>(std::ostream &lhs, const Node &rhs);
        bool operator!=(const Node &rhs) const;
        void set_max();
#ifndef WIN32
        int getNeigh(std::vector<Node> a[], int vec_size = 1);
#else
        int getNeigh(std::vector<Node> a[], int vec_size = 1, int * pairwise_costs = NULL);
#endif
        int get_g() const { return m_g; };
        int get_f() const { return m_f; };
        int get_h() const { return m_f - m_g; }; //!< heuristc estimated cost to the goal
        int get_parenti() const { return parenti; };
        inline Coord<N> get_parent() const { return pos.parent(parenti); };

    private:
        int m_g; //!< exact cost of the path from the start
        int parenti; //!< Integer representing the parent
        bool borderCheck(const Coord<N> &c) const;
        inline int pairCost(const int &neigh_num, const int &mm_cost, const int &s1, const int &s2) const;

        friend class boost::serialization::access;
        // When the class Archive corresponds to an output archive, the
        // & operator is defined similar to <<.  Likewise, when the class Archive
        // is a type of input archive the & operator is defined similar to >>.
        template<class Archive>
        void serialize(Archive & ar, unsigned int version)
        {
            ar & m_g;
            ar & m_f;
            ar & pos;
            ar & parenti;
            //ar & remote;
        }
};

/*!
 * While using libboost implementation, it is possible to use a function
 * on an "update" operation. Coord always have the same value and must
 * not be updated.
 */
template < int N >
struct change_node
{
    change_node(int new_f, int new_g, int new_parenti):new_f(new_f), new_g(new_g), new_parenti(new_parenti){}

    void operator()(Node<N> &n)
    {
        n.m_f = new_f;
        n.m_g = new_g;
        n.parenti = new_parenti;
    }
    private:
        int new_f;
        int new_g;
        int new_parenti;
};
#endif //_NODE_H
