//BlockPruner

template <int N>
bool BlockPruner::isPruned(Coord<N> coord)
{

	Coord<N> c(*coord);
	for (int i = 0; i < N; i++)
	{
		c[i] -= 1;
		if (!currBlock[c])
			return false
	}
	return true;
}

template <int N>
bool BlockPruner::processBlock(int blockNum)
{
	currBlock = block[blockNum];
	if (!isPruned<N>(block))
	{
		//process block
	}
}

template <int N>
bool BlockPruner::isInBlock(Coord<N> coord)
{
	for (int i = 0; i < N; i++)
	{
		
	}
}