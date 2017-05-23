//BlockPruner.hpp

template <int N>
class BlockPruner
{
public:

bool isPruned<N>(Coord<N> coord);
bool processBlock<N>();
bool isInBlock(Coord<N> coord);
private:
void * block;

};

#include "../BlockPruner.cpp"