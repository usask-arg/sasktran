#if !defined(NXBASE_NXBITMANIP_H)
#define NXBASE_NXBITMANIP_H 1

#if !defined(MAKE_DWORD)
#define MAKE_DWORD( lword, hword )   (((nxDWORD)(hword) << 16) | ((nxDWORD)(lword)))	//!< Make a 32 bit word from two 16 bits words
#endif

#define nxHIWORD(x)	(nxWORD)(x >> 16)													//!< Get the top 16 bits as an unsigned 16 bit integer
#define nxLOWORD(x)	(nxWORD)(x & 0xFFFF)												//!< Get the bottom 16 bits as an unsigned 16 bit integer
#define SetBit(regs,pattern)     (regs) |= (pattern)									//!< Set the pattern bits in "x"
#define ClearBit(regs,pattern)   (regs) &= ~(pattern)									//!< Clear the pattern bits in "x"
#define BitIsSet(regs,pattern)   (((regs) & (pattern)) == (pattern))					//!< Return nxTRUE if the pattern bits in x are all set
#define BitIsClear(regs,pattern) (((regs) & (pattern)) == 0)							//!< Return nxTRUE if the pattern bits in x are all clear

#endif


