/*!
 * \author Daniel Sundfeld
 * \copyright MIT License
 */
#ifndef _BACKTRACE_H
#define _BACKTRACE_H
#include <map>

#include "Coord.h"
#include "Node.h"

template <int N> 
void backtrace(std::map< Coord<N>, Node<N> > *ClosedList, int list_size = 1);

template <int N> 
Node<N> backtrace_origin(std::map< Coord<N>, Node<N> > *ClosedList, int list_size = 1, Node<N> target_node = nullptr, int min = 0, int max = 0);
#endif
