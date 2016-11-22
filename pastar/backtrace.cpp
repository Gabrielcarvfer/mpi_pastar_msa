/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 */
#include "include/backtrace.h"

#ifndef WIN32
#include <sys/ioctl.h>
#include <unistd.h>
#endif

#include <iostream>
#include <iomanip>
#include <limits>

#include "include/Sequences.h"
#include "include/TimeCounter.h"

// Decide the best lenght size to print
int get_print_size()
{
    int size = 80;
#ifdef __linux
    struct winsize w;

    // If it is a file, we dont care about lenght
    if (!isatty(1))
        return std::numeric_limits<int>::max();

    // If it is a terminal, get the lenght
    if ((ioctl(0, TIOCGWINSZ, &w) == 0) && (w.ws_col > 1))
        size = w.ws_col - 1;
#endif
    return size;
}

/*!
 * Using the last node on \a ClosedList do a backtrace, verifing
 * gaps, matches and mismatches, saving characteres o the
 * \a backstrace_alignment.
 * The \a ClosedList is an array with size \a list_size
 *
 * At the end of the process, \a alignments contains the
 * answer for every sequence.
 */
template <int N>
void backtrace_create_alignment(std::list<char> *alignments, std::map<Coord<N>, Node<N> > *ClosedList, int list_size)
{
    Sequences *seq = Sequences::getInstance();

    int id = seq->get_final_coord<N>().get_id(list_size);
    Node<N> current = ClosedList[id][seq->get_final_coord<N>()];
    std::cout << "Final Score: " << current << std::endl;
    do
    {
        //std::cout << current.pos.get_id(list_size) << ":" << current << std::endl;
        for (int i = 0; i < N; i++)
        {
            char c;
            if (current.pos[i] != current.get_parent()[i])
                c = seq->get_seq(i)[current.pos[i] - 1];
            else
                c = '-';
            alignments[i].push_front(c);
        }
        id = current.get_parent().get_id(list_size);
        current = ClosedList[id][current.get_parent()];
    } while (current.pos != Sequences::get_initial_coord<N>());
}

/*!
* Using the last node passed, do a partial backtrace until you find an external node.
* When you find an external node, return the last local node, that is checked.
* When checking node, if its the node at origin, end of backtrace, if not, send to owner of that node
* to continue the partial backtrace.
*/
template <int N>
Node<N> partial_backtrace_alignment(std::list<char> *alignments, std::map<Coord<N>, Node<N> > *ClosedList, Node<N> currentE, int list_size, int min, int max)
{
    Sequences *seq = Sequences::getInstance();
    Node<N> current = currentE;
    int id = current.pos.get_id(list_size);

    //Backtrace until find an node from external thread
    do
    {
      
        //std::cout << current.pos.get_id(list_size) << ":" << current << std::endl;
        for (int i = 0; i < N; i++)
        {
            char c;
            if (current.pos[i] != current.get_parent()[i])
                c = seq->get_seq(i)[current.pos[i] - 1];
            else
                c = '-';
            alignments[i].push_front(c);
        }
        
        id = current.get_parent().get_id(list_size);
        // if next node is remote, stop and return the node
        if ( (id < min) | (id >= max) )
        {
            break;
        }

        current = ClosedList[id-min][current.get_parent()];
    } while (current.pos != Sequences::get_initial_coord<N>());
    return current;
}

/*!
* Using the node passed to backtrace who send that node to current rank
*/
template <int N>
void backtrace_origin(std::map<Coord<N>, Node<N> > *ClosedList, int list_size, Node<N> target_node, int min, int max)
{
	Sequences *seq = Sequences::getInstance();

	int id = seq->get_final_coord<N>().get_id(list_size);
	Node<N> current = target_node;
	do
	{
		id = current.get_parent().get_id(list_size);
		//We've to brake before reaching the parent node, that is on a remote computer and cant be loaded if tried to
		if ((id < min) | (id >= max))
			break;
		current = ClosedList[id][current.get_parent()];
	} while (current.pos != Sequences::get_initial_coord<N>());
	return current;
}

/*!
 * Print the similarity of the sequences in \a alignemnts
 */
template <int N>
void backtrace_print_similarity(std::list<char> *alignments)
{
    // all alignments[i] have same size
    int total = 0;
    int equal = 0;

    std::list<char>::iterator its[N];
    for (int i = 0; i < N; ++i)
        its[i] = alignments[i].begin();

    while (its[0] != alignments[0].end())
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = i + 1; j < N; ++j)
            {
                if (*its[i] == *its[j])
                    ++equal;
                ++total;
            }
        }

        for (int i = 0; i < N; ++i)
            ++its[i];
    }
    float percent = (equal * 100) / (float) total;
    std::cout << "Similarity: "
              << std::fixed << std::setprecision(2)
              << percent << "%" << std::endl;
}

/*!
 * Print the answer in \a alignment. Use a good lenght to print
 * considering the current terminal (linux-only).
 */
template <int N>
void backtrace_print_alignment(std::list<char> *alignments)
{
    int size = get_print_size();

    while (!alignments[0].empty())
    {
        std::cout << std::endl;
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < size; i++)
            {
                if (alignments[j].empty())
                    break;
                std::cout << alignments[j].front();
                alignments[j].pop_front();
            }
            std::cout << std::endl;
        }
    }
}

/*!
 * MSA-Node backtrace functions prints the answer. Using the
 * \a ClosedList it backtrace every node until the origin is reached
 */
template <int N>
void backtrace(std::map< Coord<N>, Node<N> > *ClosedList, int list_size)
{
    TimeCounter t("Phase 3 - backtrace: ");
    std::list<char> alignments[N];

    backtrace_create_alignment<N>(alignments, ClosedList, list_size);
    backtrace_print_similarity<N>(alignments);
    backtrace_print_alignment<N>(alignments);
}

template <int N>
void print_entire_backtrace(std::list<char> alignments[])
{
    backtrace_print_similarity<N>(alignments);
    backtrace_print_alignment<N>(alignments);
}



#define DECLARE_BACKTRACE_TEMPLATE( X ) \
template void backtrace< X >(std::map< Coord< X >, Node< X > >*ClosedList, int list_size); \

#define DECLARE_BACKTRACET_TEMPLATE( X ) \
template class Node<X> partial_backtrace_alignment< X >(std::list<char> *alignments, std::map<Coord<X>, Node<X> > *ClosedList, Node<X> currentE, int list_size,  int min, int max); \

#define DECLARE_BACKTRACETT_TEMPLATE( X ) \
template void print_entire_backtrace< X >(std::list<char> alignments[X]); \

MAX_NUM_SEQ_HELPER(DECLARE_BACKTRACE_TEMPLATE);
MAX_NUM_SEQ_HELPER(DECLARE_BACKTRACET_TEMPLATE);
MAX_NUM_SEQ_HELPER(DECLARE_BACKTRACETT_TEMPLATE);
