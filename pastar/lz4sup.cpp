/*!
 * \authors {Gabriel Ferreira}
 * \copyright MIT License
 *
 * \brief LZ4 generic encode and decode for this project
 */

#include "include/lz4sup.h"
#include <lz4.h> 
#include <stdint.h>


int pastar_lz4_en(const char * input, char ** buffer, int lenStr)
{
	//Experimental lz4 compression
		uint32_t len = (uint32_t) (lenStr + 1);
		uint32_t lz4len = LZ4_compressBound((int)len);
		uint32_t lz4len2 = 0;

		//Allocate memory to store lz4 compressed data
		char * tempBuff = new char[lz4len + (sizeof(uint32_t) * 3)]();

		//Compress data
		lz4len2 = LZ4_compress(input, &tempBuff[sizeof(uint32_t) * 3], (int) len);

		//Prefix LZ4 compressed data with:

			// total size
			((uint32_t*)tempBuff)[0] = lz4len2+(sizeof(uint32_t)*3);

			// lz4 maximum decompressed length
			((uint32_t*)tempBuff)[1] = lz4len;

			// original data size
			((uint32_t*)tempBuff)[2] = len;


		//Point buffer to lz4 compressed data
		*buffer = tempBuff;

		//Return size of compressed data
		return (int) (lz4len2 + ( sizeof(uint32_t) * 3 ) );
}

int pastar_lz4_dec(char ** buffer)
{
	//Experimental LZ4 decompression
	//Fetch some parameters as:

	// maximum decompressed length to decompress data
	uint32_t lz4len = ((uint32_t*)*buffer)[1];

	// original data size 
	uint32_t len    = ((uint32_t*)*buffer)[2];

	// allocate a buffer to decompress data
	char * tempBuff = new char[lz4len]();

	//Decompress data
	LZ4_decompress_fast(&((*buffer)[sizeof(uint32_t) * 3]), tempBuff, len);

	//Clear old data pointed by buffer pointer
	delete[] *buffer;

	//Point buffer to lz4 decompressed data
	*buffer = tempBuff;

	//Return size of data
	return (int)len;
}
