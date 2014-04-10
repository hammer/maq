#ifndef LH3_DWORD_H
#define LH3_DWORD_H

#include "const.h"

template<class TYPE>
struct dword_t
{
	TYPE w1, w2;
	static const int n_bits = sizeof(TYPE) * 8;
	inline dword_t() { w1 = TYPE(0); w2 = TYPE(0); }
	inline dword_t(bit16_t _w1) { w1 = TYPE(_w1); w2 = TYPE(0); }
	inline dword_t(bit32_t _w1) { w1 = TYPE(_w1); w2 = TYPE(0); }
	inline dword_t(bit64_t _w1) { w1 = TYPE(_w1); w2 = TYPE(0); }
	inline dword_t(int _w1) : w1(TYPE(_w1)), w2(TYPE(0)) {}
	inline dword_t(const TYPE &_w1, const TYPE &_w2) : w1(_w1), w2(_w2) {}
	inline dword_t(const dword_t<TYPE> &x) : w1(x.w1), w2(x.w2) {}
	inline const dword_t<TYPE> &operator >>= (int k) {
		if (k == 0) return *this;
		if (k < n_bits) {
			w1 = (w1>>k) | (TYPE(w2)<<(-k)); w2 >>= k;
		} else if (k < n_bits*2) {
			w1 = TYPE(w2)>>(k-n_bits); w2 = TYPE(0);
		} else w1 = w2 = TYPE(0);
		return *this;
	}
	inline const dword_t<TYPE> &operator <<= (int k) {
		if (k == 0) return *this;
		if (k < n_bits) {
			w2 = (w2<<k) | (w1>>(n_bits-k)); w1 <<= k;
		} else if (k < n_bits*2) {
			w2 = w1 << (k-n_bits); w1 = TYPE(0);
		} else w1 = w2 = TYPE(0);
		return *this;
	}
	inline const dword_t<TYPE> &operator &= (const dword_t<TYPE> &w) {
		w1 &= w.w1; w2 &= w.w2;
		return *this;
	}
	inline const dword_t<TYPE> &operator |= (const dword_t<TYPE> &w) {
		w1 |= w.w1; w2 |= w.w2;
		return *this;
	}
	inline dword_t<TYPE> &operator |= (int w) {
		w1 |= w;
		return *this;
	}
	inline dword_t<TYPE> &operator ^= (const dword_t<TYPE> &w) {
		w1 ^= w.w1; w2 ^= w.w2;
		return *this;
	}
	inline operator bit64_t() const { return bit64_t(w1); }
};

template<class TYPE>
inline dword_t<TYPE> operator >> (const dword_t<TYPE> &w, int k) {
	return (k == 0)? w
		: (k < dword_t<TYPE>::n_bits)? dword_t<TYPE>((w.w1>>k) | (TYPE(w.w2)<<(dword_t<TYPE>::n_bits-k)), w.w2>>k)
		: (k < dword_t<TYPE>::n_bits*2)? dword_t<TYPE>(TYPE(w.w2)>>(k-dword_t<TYPE>::n_bits), TYPE(0))
		: dword_t<TYPE>(TYPE(0), TYPE(0));
}
template<class TYPE>
inline dword_t<TYPE> operator << (const dword_t<TYPE> &w, int k) {
	return (k == 0)? w
		: (k < dword_t<TYPE>::n_bits)? dword_t<TYPE>(w.w1<<k, (w.w2<<k) | (w.w1>>(dword_t<TYPE>::n_bits-k)))
		: (k < dword_t<TYPE>::n_bits*2)? dword_t<TYPE>(TYPE(0), TYPE(w.w1<<(k-dword_t<TYPE>::n_bits)))
		: dword_t<TYPE>(TYPE(0), TYPE(0));
	
}
template<class TYPE>
inline dword_t<TYPE> operator & (const dword_t<TYPE> &w, const dword_t<TYPE> &ww) {
	return dword_t<TYPE>(w.w1 & ww.w1, w.w2 & ww.w2);
}
template<class TYPE>
inline int operator < (const dword_t<TYPE> &w, const dword_t<TYPE> &ww) {
	return (w.w2 < ww.w2 || (w.w2 == ww.w2 && w.w1 < ww.w1));
}
template<class TYPE>
inline dword_t<TYPE> operator | (const dword_t<TYPE> &w, const dword_t<TYPE> &ww) {
	return dword_t<TYPE>(w.w1 | ww.w1, w.w2 | ww.w2);
}
template<class TYPE>
inline dword_t<TYPE> operator ^ (const dword_t<TYPE> &w, const dword_t<TYPE> &ww) {
	return dword_t<TYPE>(w.w1 ^ ww.w1, w.w2 ^ ww.w2);
}
template<class TYPE>
inline int operator != (const dword_t<TYPE> &w, const dword_t<TYPE> &ww) {
	return (w.w1 != ww.w1 || w.w2 != ww.w2);
}
template<class TYPE>
inline int operator == (const dword_t<TYPE> &w, const dword_t<TYPE> &ww) {
	return (w.w1 == ww.w1 && w.w2 == ww.w2);
}
template<class TYPE>
inline dword_t<TYPE> operator ~ (const dword_t<TYPE> &w) {
	return dword_t<TYPE>(~w.w1, ~w.w2);
}
template<class TYPE, class T2>
inline T2 operator & (const dword_t<TYPE> &w, T2 k) {
	return T2(w.w1 & k);
}
template<class TYPE, class T2>
inline dword_t<TYPE> operator | (const dword_t<TYPE> &w, T2 k) {
	return dword_t<TYPE>(w.w1 | k, w.w2);
}

#endif
