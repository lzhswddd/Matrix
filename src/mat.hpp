#pragma once
#ifndef __MAT_H__
#define __MAT_H__

//#define MAT_LIB

//#ifdef MAT_LIB
//#define 
//#else
//#ifdef MAT_EXPORTS
//#define  __declspec(dllexport)
//#else 
//#define  __declspec(dllimport)
//#endif
//#endif

#include <vector>
#include <string>


#ifdef _DEBUG
#define MAT_DEBUG
#else
#define MAT_RELEASE
#endif

#ifndef LENGTH
#define LENGTH(x) (sizeof(x)/sizeof(mat_t))
#endif
#ifndef ILength
#define ILength(x) (sizeof(x)/sizeof(int))
#endif

#define MALLOC_ALIGN 16

#define FREE_PTR(ptr) if((ptr)!=nullptr) {delete (ptr); (ptr)=nullptr;}

#define ARRAY_LEN(start, end) int(((end) - (start)))

#define RANGE(low, high, v) ((low)<(v)&&(v)>(high))

	// exchange-add operation for atomic operations on reference counters
#if defined __INTEL_COMPILER && !(defined WIN32 || defined _WIN32)
	// atomic increment on the linux version of the Intel(tm) compiler
#  define MAT_XADD(addr, delta) (int)_InterlockedExchangeAdd(const_cast<void*>(reinterpret_cast<volatile void*>(addr)), delta)
#elif defined __GNUC__
#  if defined __clang__ && __clang_major__ >= 3 && !defined __ANDROID__ && !defined __EMSCRIPTEN__ && !defined(__CUDACC__)
#    ifdef __ATOMIC_ACQ_REL
#      define MAT_XADD(addr, delta) __c11_atomic_fetch_add((_Atomic(int)*)(addr), delta, __ATOMIC_ACQ_REL)
#    else
#      define MAT_XADD(addr, delta) __atomic_fetch_add((_Atomic(int)*)(addr), delta, 4)
#    endif
#  else
#    if defined __ATOMIC_ACQ_REL && !defined __clang__
// version for gcc >= 4.7
#      define MAT_XADD(addr, delta) (int)__atomic_fetch_add((unsigned*)(addr), (unsigned)(delta), __ATOMIC_ACQ_REL)
#    else
#      define MAT_XADD(addr, delta) (int)__sync_fetch_and_add((unsigned*)(addr), (unsigned)(delta))
#    endif
#  endif
#elif defined _MSC_VER && !defined RC_INVOKED
#  include <intrin.h>
#  define MAT_XADD(addr, delta) (int)_InterlockedExchangeAdd((long volatile*)addr, delta)
#else
static inline void MAT_XADD(int* addr, int delta) { int tmp = *addr; *addr += delta; return tmp; }
#endif

#define MAT_FLOAT

namespace lzh {
#ifdef MAT_FLOAT
	typedef float mat_t;
#elif defined MAT_DOUBLE
	typedef double mat_t;
#else
	typedef long double mat_t;
#endif

#define _T(v) (lzh::mat_t(v))

	static const mat_t PI = _T(3.1415926535897932384626433832795);
	//LeakyReLU的超参
	static const mat_t LReLU_alpha = _T(0.2);
	//ELU的超参
	static const mat_t ELU_alpha = _T(1.6732632423543772848170429916717);
	//SELU的超参
	static const mat_t SELU_scale = _T(1.0507009873554804934193349852946);

	static const mat_t shrink_factor = _T(1.247330950103979);
	
	typedef unsigned char uchar;
	typedef unsigned int uint;
	enum PrintType {
		FIXED = 0,
		SCIENTIFIC
	};
	enum BorderTypes {
		BORDER_CONSTANT = 0, //!< `iiiiii|abcdefgh|iiiiiii`  with some specified `i`
		BORDER_REPLICATE = 1, //!< `aaaaaa|abcdefgh|hhhhhhh`
		BORDER_REFLECT = 2, //!< `fedcba|abcdefgh|hgfedcb`
		BORDER_WRAP = 3, //!< `cdefgh|abcdefgh|abcdefg`
		BORDER_REFLECT_101 = 4, //!< `gfedcb|abcdefgh|gfedcba`
		BORDER_TRANSPARENT = 5, //!< `uvwxyz|abcdefgh|ijklmno`
		BORDER_ISOLATED = 16 //!< do not look outside of ROI
	};
	enum MatErrorInfo {
		ERR_INFO_EMPTY = 0,
		ERR_INFO_SQUARE,
		ERR_INFO_ADJ,
		ERR_INFO_INV,
		ERR_INFO_POW,
		ERR_INFO_IND,
		ERR_INFO_CON,
		ERR_INFO_EIGEN,
		ERR_INFO_LEN,
		ERR_INFO_MEMOUT,
		ERR_INFO_UNLESS,
		ERR_INFO_SIZE,
		ERR_INFO_MULT,
		ERR_INFO_NORM,
		ERR_INFO_VALUE,
		ERR_INFO_PINV,
		ERR_INFO_DET,
		ERR_INFO_DIM,
		ERR_INFO_PTR,
		ERR_INFO_NOT2D,
		ERR_INFO_FILE,
	};
	static const char *errinfo[] = {
		"error 0: 矩阵为空!\0",
		"error 1: 矩阵不是方阵!\0",
		"error 2: 矩阵不是方阵，不能设置伴随矩阵!\0",
		"error 3: 矩阵不是方阵，不能设置逆矩阵!\0",
		"error 4: 矩阵不是方阵，不能进行次幂运算!\0",
		"error 5: 矩阵不是方阵，不能设置为单位矩阵!\0",
		"error 6: 矩阵不收敛!\0",
		"error 7: 矩阵没有实数特征值!\0",
		"error 8: 矩阵维度为0!\0",
		"error 9: 矩阵索引出界!\0",
		"error 10: 矩阵索引无效!\0",
		"error 11: 两个矩阵维度不一致!\0",
		"error 12: 两个矩阵维度不满足乘法条件!\0",
		"error 13: 矩阵维度不为1，不是向量!\0",
		"error 14: 参数违法!\0",
		"error 15: 计算逆矩阵失败!\0",
		"error 16: 行列式为0!\0",
		"error 17: 不支持三维操作!\0",
		"error 18: 指针为空!\0",
		"error 19: 输入矩阵维度必须为2D!\0",
		"error 20: 没有找到文件!\0"
	};
	/**
	SPECIAL_SOLUTION	方程有特解
	GENERAL_SOLUTION	方程有通解
	NO_SOLUTION			方程无解
	*/
	enum EQUATION
	{
		SPECIAL_SOLUTION = 0,	//方程有特解
		GENERAL_SOLUTION,		//方程有通解
		NO_SOLUTION				//方程无解
	};
	/**
	MIN_TO_MAX	从小到大
	MAX_TO_MIN	从大到小
	*/
	enum ORDER
	{
		MIN_TO_MAX = 0,
		MAX_TO_MIN
	};
	/**
	ROW		行
	COL		列	
	CHANNEL	通道
	*/
	enum X_Y_Z {
		ROW = 0,
		COL,
		CHANNEL
	};
	enum Dire {
		LEFT = 0,
		RIGHT
	};
	/**
	EqualIntervalSampling 等间隔采样
	LocalMean 局部均值
	*/
	enum ReductionMothed
	{
		EqualIntervalSampling = 0,
		LocalMean
	};
	/**
	旋转方向顺时针
	ROTATE_90_ANGLE 90度
	ROTATE_180_ANGLE 180度
	ROTATE_270_ANGLE 270度
	*/
	enum RotateAngle
	{
		ROTATE_90_ANGLE = 0,
		ROTATE_180_ANGLE,
		ROTATE_270_ANGLE
	};

	template<typename _Tp>
	static inline _Tp *alignPtr(_Tp *ptr, int n = (int) sizeof(_Tp)) {
		return (_Tp *)(((size_t)ptr + n - 1) & -n);
	}

	// Aligns a buffer size to the specified number of bytes
	// The function returns the minimum number that is greater or equal to sz and is divisible by n
	// sz Buffer size to align
	// n Alignment size that must be a power of two
	static inline size_t alignSize(size_t sz, int n) {
		return (sz + n - 1) & -n;
	}

	static inline void *fastMalloc(size_t size) {
		unsigned char *udata = (unsigned char *)malloc(size + sizeof(void *) + MALLOC_ALIGN);
		if (!udata)
			return 0;
		unsigned char **adata = alignPtr((unsigned char **)udata + 1, MALLOC_ALIGN);
		adata[-1] = udata;
		return adata;
	}

	static inline void fastFree(void *ptr) {
		if (ptr) {
			unsigned char *udata = ((unsigned char **)ptr)[-1];
			free(udata);
		}
	}
	static void THROW_INFO(lzh::MatErrorInfo info);
	/**
	内存管理类
	*/
	template<class Tp_>
	class MatPtr
	{
	public:
		MatPtr<Tp_>()
			: p(nullptr), refcount(nullptr){}
		MatPtr<Tp_>(int size)
		{
			create(size);
		}
		MatPtr<Tp_>(const MatPtr<Tp_> &m)
			: p(m.p), refcount(m.refcount) 
		{
			if (refcount)
				MAT_XADD(refcount, 1);
		}
		~MatPtr<Tp_>() {
			release();
		}
		void create(uint len){
			uint totalsize = len * sizeof(Tp_);
			p = (Tp_*)fastMalloc(totalsize + sizeof(*refcount));
			refcount = (int*)(((unsigned char *)p) + totalsize);
			*refcount = 1;
		}
		void addref() {
			if (refcount)
				MAT_XADD(refcount, 1);
		}
		void release() {
			if (refcount && MAT_XADD(refcount, -1) == 1)
				fastFree(p);
			refcount = nullptr; 
			p = nullptr;
		}
		MatPtr<Tp_> &operator=(const MatPtr<Tp_> &m) {
			if (this == &m)
				return *this;

			if (m.refcount)
				MAT_XADD(m.refcount, 1);

			release();
			p = m.p;
			refcount = m.refcount;
			return *this;
		}
		MatPtr<Tp_> &operator = (Tp_ *m) {
			if (this->p == m)
				return *this;
			if (refcount && MAT_XADD(refcount, -1) == 1)
				fastFree(p);
			refcount = nullptr;
			p = m;
			return *this;
		}
		bool operator == (void *m) const {
			return (this->p == m);
		}
		bool operator != (void *m) const {
			return (this->p != m);
		}
		operator Tp_ *() {
			return p;
		}
		operator const Tp_ *() const {
			return p;
		}
		Tp_& operator [] (int idx) const {
			return p[idx];
		}
		Tp_ *p;
	private:
		int *refcount;
	};
	class Rect
	{
	public:
		Rect() : x(0), y(0), width(0), height(0) {}
		Rect(int x, int y, int width, int height) : x(x), y(y), width(width), height(height) {}
		int x;
		int y;
		int width;
		int height;
		int area()const {
			return width * height;
		}
		friend std::ostream & operator << (std::ostream &out, const Rect &t)
		{
			out << "Rect(" << t.x << "," << t.y << "," << t.width << "," << t.height << ")";
			return out;
		}
	};
	class RectF
	{
	public:
		RectF() : x(0), y(0), width(0), height(0) {}
		RectF(float x, float y, float width, float height) : x(x), y(y), width(width), height(height) {}
		float x;
		float y;
		float width;
		float height;
		float area()const {
			return width * height;
		}
		friend std::ostream & operator << (std::ostream &out, const RectF &t)
		{
			out << "Rect(" << t.x << "," << t.y << "," << t.width << "," << t.height << ")";
			return out;
		}
	};
	class Size
	{
	public:
		Size() :h(0), w(0) {}
		Size(int width, int height) :h(height), w(width) {}
		~Size() {}
		int h;
		int w;
		int area()const {
			return h * w;
		}
		bool operator == (Size size)const
		{
			return (h == size.h && w == size.w);
		}
		bool operator != (Size size)const
		{
			return !(*this == size);
		}
		friend std::ostream & operator << (std::ostream &out, const Size &t)
		{
			out << "Size(" << t.h << "," << t.w << ")";
			return out;
		}
	};
	class Size3
	{
	public:
		explicit Size3() : h(0), w(0), c(0) {}
		Size3(Size size) : h(size.h), w(size.w), c(1) {}
		Size3(int w, int h, int c = 1) : h(h), w(w), c(c) {}
		int h;
		int w;
		int c;
		int area()const {
			return h * w * c;
		}
		bool operator == (Size3 size)const
		{
			return (h == size.h && w == size.w && c == size.c);
		}
		bool operator != (Size3 size)const
		{
			return !(*this == size);
		}
		friend std::ostream & operator << (std::ostream &out, const Size3 &t)
		{
			out << "Size(" << t.h << "," << t.w << "," << t.c << ")";
			return out;
		}
	};
	template<class Tp_>
	class Point2
	{
	public:
		Point2() :x(), y() {}
		Point2(Tp_ x, Tp_ y) :x(x), y(y) {}
		template<typename T2_>
		Point2(const Point2<T2_> &p) :x(Tp_(p.x)), y(Tp_(p.y)) {}
		~Point2() {}
		bool operator == (const Point2<Tp_> &P)const
		{
			return (x == P.x) && (y == P.y);
		}
		bool operator != (const Point2<Tp_> &P)const
		{
			return !((*this) == P);
		}
		Tp_ x;
		Tp_ y;
		friend std::ostream & operator << (std::ostream &out, const Point2<Tp_> &t)
		{
			out << "Point(" << t.x << "," << t.y << ")";
			return out;
		}
	};
	template<typename Tp_>
	inline const Point2<Tp_> operator + (const Tp_& v, const Point2<Tp_> &P)
	{
		return Point2<Tp_>(P.x + v, P.y + v);
	}
	template<typename Tp_>
	inline const Point2<Tp_> operator + (const Point2<Tp_> &P1, const Point2<Tp_>& P2)
	{
		return Point2<Tp_>(P1.x + P2.x, P1.y + P2.y);
	}
	template<typename Tp_>
	inline const Point2<Tp_> operator - (const Tp_& v, const Point2<Tp_> &P)
	{
		return Point2<Tp_>(v - P.x, v - P.y);
	}
	template<typename Tp_>
	inline const Point2<Tp_> operator - (const Point2<Tp_> &P1, const Point2<Tp_>& P2)
	{
		return Point2<Tp_>(P1.x - P2.x, P1.y - P2.y);
	}
	template<typename Tp_>
	inline const Point2<Tp_> operator * (const Tp_& v, const Point2<Tp_> &P)
	{
		return Point2<Tp_>(P.x * v, P.y * v);
	}
	template<typename Tp_>
	inline const Tp_ operator * (const Point2<Tp_> &P1, const Point2<Tp_>& P2)
	{
		return P1.x * P2.x + P1.y * P2.y;
	}
	typedef Point2<char> Point2c;
	typedef Point2<uchar> Point2b;
	typedef Point2<int> Point2i;
	typedef Point2<uint> Point2Ui;
	typedef Point2<float> Point2f;
	typedef Point2<double> Point2d;
	typedef Point2<mat_t> Point2m;
	typedef Point2i Point;

	template<class Tp_>
	class Point3
	{
	public:
		Point3() :x(), y(), z() {}
		Point3(Tp_ x, Tp_ y, Tp_ z) :x(x), y(y), z(z) {}
		template<typename T2_>
		Point3(const Point3<T2_> &p) : x(Tp_(p.x)), y(Tp_(p.y)), z(Tp_(p.z)) {}
		~Point3() {}
		bool operator == (Point2<Tp_> &P)const
		{
			return (x == P.x) && (y == P.y) && (z == P.z);
		}
		bool operator != (Point2<Tp_> &P)const
		{
			return !((*this) == P);
		}
		Tp_ x;
		Tp_ y;
		Tp_ z;
		friend std::ostream & operator << (std::ostream &out, const Point3<Tp_> &t)
		{
			out << "Point(" << t.x << "," << t.y << "," << t.z << ")";
			return out;
		}
	};
	template<typename Tp_>
	inline const Point3<Tp_> operator + (const Tp_& v, const Point3<Tp_> &P)
	{
		return Point3<Tp_>(P.x + v, P.y + v, P.z + v);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator + (const Point3<Tp_> &P, const Tp_& v)
	{
		return Point3<Tp_>(P.x + v, P.y + v, P.z + v);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator + (const Point3<Tp_> &P1, const Point3<Tp_>& P2)
	{
		return Point3<Tp_>(P1.x + P2.x, P1.y + P2.y, P1.z + P2.z);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator - (const Tp_& v, const Point3<Tp_> &P)
	{
		return Point3<Tp_>(v - P.x, v - P.y, v - P.z);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator - (const Point3<Tp_> &P, const Tp_& v)
	{
		return Point3<Tp_>(P.x - v, P.y - v, P.z - v);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator - (const Point3<Tp_> &P1, const Point3<Tp_>& P2)
	{
		return Point2<Tp_>(P1.x - P2.x, P1.y - P2.y, P1.z - P2.z);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator * (const Tp_& v, const Point3<Tp_> &P)
	{
		return Point3<Tp_>(P.x * v, P.y * v, P.z * v);
	}
	template<typename Tp_>
	inline const Point3<Tp_> operator * (const Point3<Tp_> &P, const Tp_& v)
	{
		return Point3<Tp_>(P.x * v, P.y * v, P.z * v);
	}
	template<typename Tp_>
	inline const Tp_ operator * (const Point3<Tp_> &P1, const Point3<Tp_>& P2)
	{
		return P1.x * P2.x + P1.y * P2.y + P1.z * P2.z;
	}
	typedef Point3<char> Point3c;
	typedef Point3<uchar> Point3b;
	typedef Point3<int> Point3i;
	typedef Point3<uint> Point3Ui;
	typedef Point3<float> Point3f;
	typedef Point3<double> Point3d;
	typedef Point3<mat_t> Point3m;

	class Color
	{
	public:
		Color(uchar v) : r(v), g(v), b(v) {}
		Color(uchar r, uchar g, uchar b) : r(r), g(g), b(b) {}
		uchar r;
		uchar g;
		uchar b;
	};
	template<class Type>
	class Vec
	{
	public:
		explicit Vec()
			: len(0){}
		Vec(int size)
			: data(size), len(size)
		{
			memset(data, 0, sizeof(Type) * size);
		}
		Vec(Type *data, int size = 3) : len(size) {
			this->data = data;
		}
		Vec(Type d1, Type d2, Type d3) : data(3), len(3){
			data[0] = d1;
			data[1] = d2;
			data[2] = d3;
		}
		~Vec()
		{
			data.release();
			len = 0;
		}
		Type* begin()
		{
			return data;
		}
		const Type* begin()const
		{
			return data;
		}
		Type* end()
		{
			if (data == nullptr)return nullptr;
			return data + len;
		}
		const Type* end()const
		{
			if (data == nullptr)return nullptr;
			return data + len;
		}
		Type& operator [](const int index)const {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		void operator = (const Vec<Type> &vec)
		{
			len = vec.len;
			data = vec.data;
		}
		void operator = (Color &color)
		{
			if (len == 1) {
				data[0] = color.r;
			}
			else if (len == 3) {
				data[0] = color.r;
				data[1] = color.g;
				data[2] = color.b;
			}
		}
		int len;
		MatPtr<Type> data;
	};
	template<class Type>
	class Vec2
	{
	public:
		explicit Vec2() {
			data.create(len);
		}
		Vec2(const Type &v0, const Type &v1)
		{
			data.create(len);
			data[0] = v0;
			data[1] = v1;
		}
		template<typename T2_>
		Vec2(const Vec2<T2_> &v)
		{
			data.create(len);
			data[0] = Type(v[0]);
			data[1] = Type(v[1]);
		}
		Type& operator [](const int index) {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		const Type& operator [](const int index)const {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		Type& at(int index) {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		const Type& at(const int index)const {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		Type* begin()
		{
			return data;
		}
		const Type* begin()const
		{
			return data;
		}
		Type* end()
		{
			if (data == nullptr)return nullptr;
			return data + len;
		}
		const Type* end()const
		{
			if (data == nullptr)return nullptr;
			return data + len;
		}
		const Vec2<Type> operator + (const Type& v)const
		{
			return Vec2<Type>(data[0] + v, data[1] + v);
		}
		const Vec2<Type> operator - (const Type& v)const
		{
			return Vec2<Type>(data[0] - v, data[1] - v);
		}
		const Vec2<Type> operator * (const Type& v)const
		{
			return Vec2<Type>(data[0] * v, data[1] * v);
		}
		const Vec2<Type> operator / (const Type& v)const
		{
			return Vec2<Type>(data[0] / v, data[1] / v);
		}
	private:
		const int len = 2;
		MatPtr<Type> data;
	};
	typedef Vec2<char> Vec2c;
	typedef Vec2<uchar> Vec2b;
	typedef Vec2<int> Vec2i;
	typedef Vec2<uint> Vec2Ui;
	typedef Vec2<float> Vec2f;
	typedef Vec2<double> Vec2d;
	typedef Vec2<mat_t> Vec2m;
	template<class Type>
	class Vec3
	{
	public:
		explicit Vec3() {}
		Vec3(const Type &v0, const Type &v1, const Type &v2)
		{
			data[0] = v0;
			data[1] = v1;
			data[1] = v2;
		}
		template<typename T2_>
		Vec3(const Vec3<T2_> &v)
		{
			data[0] = Type(v[0]);
			data[1] = Type(v[1]);
			data[2] = Type(v[2]);
		}
		Type& operator [](const int index) {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		const Type& operator [](const int index)const {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		Type& at(int index) {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		const Type& at(const int index)const {
			if (index < 0 || index >= len) {
				fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			}
			return data[index];
		}
		Type* begin()
		{
			return data;
		}
		const Type* begin()const
		{
			return data;
		}
		Type* end()
		{
			if (data == nullptr)return nullptr;
			return data + len;
		}
		const Type* end()const
		{
			if (data == nullptr)return nullptr;
			return data + len;
		}
		const Vec3<Type> operator + (const Type& v)const
		{
			return Vec3<Type>(data[0] + v, data[1] + v, data[2] + v);
		}
		const Vec3<Type> operator - (const Type& v)const
		{
			return Vec3<Type>(data[0] - v, data[1] - v, data[2] - v);
		}
		const Vec3<Type> operator * (const Type& v)const
		{
			return Vec3<Type>(data[0] * v, data[1] * v, data[2] * v);
		}
		const Vec3<Type> operator / (const Type& v)const
		{
			return Vec3<Type>(data[0] / v, data[1] / v, data[2] / v);
		}
	private:
		const int len = 3;
		MatPtr<Type> data;
	};
	typedef Vec3<char> Vec3c;
	typedef Vec3<uchar> Vec3b;
	typedef Vec3<int> Vec3i;
	typedef Vec3<uint> Vec3Ui;
	typedef Vec3<float> Vec3f;
	typedef Vec3<double> Vec3d;
	typedef Vec3<mat_t> Vec3m;
	class MatCommaInitializer_;
	class  Matrix
	{
	public:
		explicit Matrix();
		/**
		生成w_的向量
		@param w		向量长度
		*/
		Matrix(int w);
		/**
		生成w_*h_的矩阵
		@param w		矩阵列数
		@param h		矩阵行数
		*/
		Matrix(int w, int h);
		/**
		生成h_*w_*depth的矩阵
		@param w		矩阵列数
		@param h		矩阵行数
		@param c		矩阵通道数
		*/
		Matrix(int w, int h, int c);
		/**
		生成size*1的矩阵
		@param size_	矩阵尺寸
		*/
		Matrix(Size size_);
		/**
		生成size的矩阵
		@param size_	矩阵尺寸
		*/
		Matrix(Size3 size_);
		/**
		拷贝函数
		@param src		拷贝对象
		*/
		Matrix(const lzh::Matrix *src);
		/**
		拷贝函数
		@param src		拷贝对象
		*/
		Matrix(const lzh::Matrix &src);
		/**
		将矩阵a和b合并(COL为按列合并|ROW为按行合并)
		@param a		输入矩阵1
		@param b		输入矩阵2
		@param merge	合并方式
		*/
		Matrix(const lzh::Matrix &a, const lzh::Matrix &b, X_Y_Z merge);
		/**
		构造函数
		深拷贝m
		@param m 矩阵
		*/
		Matrix(const MatCommaInitializer_ &m);
		/**
		生成n*n*1的矩阵,元素为matrix
		@param matrix	矩阵元素
		@param n		矩阵大小
		*/
		Matrix(int *matrix, int n);
		/**
		生成n*n*1的矩阵,元素为matrix
		@param matrix	矩阵元素
		@param n		矩阵大小
		*/
		Matrix(mat_t *matrix, int n);
		/**
		生成h_*w_*1的矩阵,元素为matrix
		@param matrix	矩阵元素
		@param w		矩阵列数
		@param h		矩阵行数
		*/
		Matrix(int *matrix, int w, int h, int c = 1);
		/**
		生成h_*w_*1的矩阵,元素为matrix
		@param matrix	矩阵元素
		@param w		矩阵列数
		@param h		矩阵行数
		*/
		Matrix(mat_t *matrix, int w, int h, int c = 1);
		/**
		生成1*w*1的向量,元素为data
		@param w		列数
		@param data		矩阵
		*/
		Matrix(int w, mat_t *data);
		/**
		生成h*w*1的矩阵,元素为data
		@param w		矩阵列数
		@param h		矩阵行数
		@param data		矩阵元素
		*/
		Matrix(int w, int h, mat_t *data);
		/**
		生成h*w*c的矩阵,元素为data
		@param w		矩阵列数
		@param h		矩阵行数
		@param c		矩阵通道数
		@param data		矩阵元素
		*/
		Matrix(int w, int h, int c, mat_t *data);
		/**
		生成h*w*c的矩阵,元素为data
		@param w		矩阵列数
		@param h		矩阵行数
		@param c		矩阵通道数
		@param step		步长
		@param data		矩阵元素
		*/
		Matrix(int w, int h, int c, int step, mat_t *data);
		/**
		向量转矩阵
		@param vec		向量
		@param dirc		矩阵方向
		*/
		template<class Type>
		Matrix(const std::vector<Type> &vec, X_Y_Z dirc)
		{
			init();
			Matrix mat;
			switch (dirc)
			{
			case lzh::ROW:mat = zeros((int)vec.size(), 1, 1);
				break;
			case lzh::COL:mat = zeros(1, (int)vec.size(), 1);
				break;
			case lzh::CHANNEL:mat = zeros(1, 1, (int)vec.size());
				break;
			default:
				return;
			}
			int idx = 0;
			for (const Type&v : vec)
				mat(idx++) = (mat_t)v;
			*this = mat;
		}
		template<class Type>
		Matrix(const Vec2<Type> &vec)
		{
			init();
			create(2);
			for (int i = 0; i < 2; i++)
				matrix[i] = _T(vec[i]);
		}
		template<class Type>
		Matrix(const Vec3<Type> &vec)
		{
			init();
			create(3);
			for (int i = 0; i < 3; i++)
				matrix[i] = _T(vec[i]);
		}
		template<class Type>
		Matrix(const Vec<Type> &vec)
		{
			init();
			create(vec.h_*vec.w_*vec.c_);
			memcpy(matrix, vec.data, memsize());
		}
		~Matrix();
		/**
		生成w的向量
		@param w		列数
		*/
		void create(int w);
		/**
		生成w*h的矩阵
		@param w		列数
		@param h		行数
		*/
		void create(int w, int h);
		/**
		生成w*h*c的张量
		@param w		列数
		@param h		行数
		@param c		通道数
		*/
		void create(int w, int h, int c);
		/**
		生成size的矩阵
		@param size		矩阵大小
		*/
		void create(Size size);
		/**
		生成size的张量
		@param size		矩阵大小
		*/
		void create(Size3 size);
		mat_t* data()const;
		mat_t* begin();
		const mat_t* begin()const;
		mat_t* end();
		const mat_t* end()const;
		/**
		@brief 内存长度
		*/
		uint memsize()const;
		/**
		@brief 检查维度
		*/
		void DimCheck()const;
		/**
		@brief 返回矩阵尺寸(h_,w_,c_)
		*/
		Size3 size3()const;
		/**
		@brief 返回矩阵偏移
		*/
		int total()const;
		/**
		@brief 返回维度
		*/
		int dims()const;
		/**
		@brief 返回行数
		*/
		int rows()const;
		/**
		@brief 返回列数
		*/
		int cols()const;
		/**
		@brief 返回通道数
		*/
		int channels()const;
		/**
		@brief 返回行秩
		*/
		int rank()const;
		/**
		@brief 返回矩阵大小(h_*w_*c_)
		*/
		uint size()const;
		/**
		@brief 返回矩阵大小Size(h_,w_)
		*/
		Size mSize()const;
		/**
		保存矩阵
		@param file		保存文件名
		@param binary	选择文本还是二进制
		binary = false	选择文本
		binary = true	选择二进制
		*/
		void save(std::string file, bool binary = true)const;
		/**
		二进制保存矩阵
		@param file		保存文件指针
		*/
		void save(FILE *file)const;
		/**
		读取矩阵
		@param file		读取文件名
		只支持二进制读取
		*/
		void load(std::string file);
		/**
		读取矩阵
		@param file		读取文件指针
		只支持二进制读取
		*/
		void load(FILE *file);
		/**
		@brief 返回矩阵大小(h_*w_*c_)
		*/
		int len()const;
		/**
		@brief 返回矩阵状态
		0为矩阵
		-1为空矩阵
		-2为非方阵
		*/
		int enable()const;
		/**
		@brief 返回矩阵是否为空
		*/
		bool empty()const;
		/**
		@brief 返回矩阵是否为矩阵
		*/
		bool Square()const;
		/**
		拷贝数据到mat
		@param mat		输入
		*/
		void copyTo(Matrix& mat)const;
		/**
		拷贝数据到mat
		@param mat		输入
		*/
		void copyTo(Matrix&& mat)const;
		/**
		在矩阵最左边或最右边添加一列1
		@param dire		选择添加方式
		*/
		void setAddOnes(Dire dire = RIGHT);
		/**
		在矩阵最左边或最右边添加一列0
		@param dire		选择添加方式
		*/
		void setAddZeros(Dire dire = RIGHT);
		/**
		释放内存
		*/
		void release();
		/**
		@brief 按索引返回矩阵元素
		@param w		索引
		*/
		mat_t& at(int w)const;
		/**
		@brief 按索引返回矩阵元素
		@param w		列索引
		@param h		行索引
		*/
		mat_t& at(int w, int h)const;
		/**
		@brief 按索引返回矩阵元素
		@param w		列索引
		@param h		行索引
		@param c		通道索引
		*/
		mat_t& at(int w, int h, int c)const;
		/**
		@brief 将索引转换为对应矩阵列索引
		@param index	索引
		*/
		int toX(int index)const;
		/**
		@brief 将索引转换为对应矩阵行索引
		@param index	索引
		*/
		int toY(int index)const;
		/**
		@brief 将索引转换为对应矩阵通道索引
		@param index	索引
		*/
		int toZ(int index)const;

		/**
		@brief 矩阵第一个元素
		*/
		mat_t frist()const;
		/**
		@brief 返回矩阵与value相等的第一个元素索引
		@param value	元素
		*/
		int find(mat_t value)const;
		/**
		@brief 返回矩阵元素最大值的索引
		*/
		int maxAt()const;
		/**
		@brief 返回矩阵元素最小值的索引
		*/
		int minAt()const;
		/**
		@brief 返回矩阵是否包含value
		@param value	元素
		*/
		bool contains(mat_t value)const;
		/**
		@brief 返回矩阵与value相等的第一个元素
		@param value	元素
		*/
		mat_t& findAt(mat_t value)const;
		mat_t max(bool is_abs = false)const;
		mat_t min(bool is_abs = false)const;
		/**
		@brief 返回矩阵元素最大值
		*/
		mat_t& findmax()const;
		/**
		@brief 返回矩阵元素最小值
		*/
		mat_t& findmin()const;
		/**
		@brief 将矩阵按索引区域拷贝元素到src矩阵中
		@param src			被拷贝矩阵
		@param Row_Start	截取行初始索引值
		@param Row_End		截取行结束索引值
		@param Col_Start	截取列初始索引值
		@param Col_End		截取列结束索引值
		*/
		void copy(Matrix &src, int Row_Start, int Row_End, int Col_Start, int Col_End)const;
		/**
		@brief 将矩阵拷贝到src
		@param src 被拷贝矩阵
		*/
		void swap(Matrix &src)const;
		/**
		@brief mChannel 将src覆盖到第c_通道
		@param src		矩阵
		@param c_	通道数
		*/
		void mChannel(const lzh::Matrix &src, int c);
		/**
		@brief mChannel 将src覆盖到第c_通道
		@param src		矩阵
		@param c_	通道数
		*/
		void mChannel(const lzh::Matrix &src, int w, int h);
		/**
		@brief 设置矩阵维度
		不允许改变矩阵长度
		@param size		矩阵大小
		*/
		void reshape(Size3 size = Size3(0, 0, 0));
		/**
		@brief 设置矩阵维度
		不允许改变矩阵长度
		@param h_		矩阵行数
		@param w_		矩阵列数
		@param c_	矩阵通道
		*/
		void reshape(int w, int h = 0, int c = 0);
		/**
		@brief 设置矩阵大小
		如果矩阵原大小不等于h_*w_*1则元素全部重置为0
		@param h_	矩阵行数
		@param w_	矩阵列数
		*/
		bool setSize(int w, int h, int c);
		/**
		@brief 拷贝矩阵src
		@param src	拷贝矩阵
		*/
		void setvalue(const lzh::Matrix &src);
		/**
		@brief 修改矩阵对应索引元素
		@param number	元素
		@param index	索引
		*/
		void setNum(mat_t number, int index);
		/**
		@brief 修改矩阵对应索引元素
		@param number	元素
		@param index_y	行索引
		@param index_x	列索引
		*/
		void setNum(mat_t number, int index_y, int index_x);
		/**
		@brief 重置矩阵
		@param mat	矩阵元素
		@param h_	矩阵行数
		@param w_	矩阵列数
		*/
		void setMat(mat_t *mat, int w, int h);
		/**
		@brief 设置逆矩阵
		*/
		void setInv();
		/**
		@brief 设置矩阵的num次幂
		@param num	次幂
		*/
		void setPow(mat_t num);
		/**
		@brief 设置取反
		*/
		void setOpp();
		/**
		@brief 设置单位矩阵
		*/
		void setIden();
		/**
		@brief 设置伴随矩阵
		*/
		void setAdj();
		/**
		@brief 设置转置矩阵
		*/
		void setTran();

		/**
		@brief 命令行输出矩阵
		*/
		void show()const;
		/**
		@brief 输出矩阵
		*/
		std::ostream & show(std::ostream & out)const;

		/**
		@brief 返回h行矩阵
		@param 索引
		*/
		Matrix Row(int h);
		/**
		@brief 返回h行矩阵
		@param 索引
		*/
		const Matrix Row(int h)const;
		/**
		@brief 返回w列矩阵
		@param 索引
		*/
		Matrix Col(int w);
		/**
		@brief 返回w列矩阵
		@param 索引
		*/
		const Matrix Col(int w)const;
		/**
		@brief 返回c通道矩阵
		@param 通道索引
		*/
		Matrix Channel(int c);
		/**
		@brief 返回c通道矩阵
		@param 通道索引
		*/
		const Matrix Channel(int c)const;
		/**
		@brief 在矩阵最左边或最右边添加一列1
		@param dire		选择添加方式
		*/
		const Matrix addones(Dire dire = RIGHT)const;
		/**
		@brief 在矩阵最左边或最右边添加一列0
		@param dire		选择添加方式
		*/
		const Matrix addzeros(Dire dire = RIGHT)const;
		/**
		@brief 返回start到end矩阵
		@param start	开始索引
		@param end		结束索引
		*/
		const Matrix Range(int start, int end);
		/**
		@brief 返回h_start到h_end行矩阵
		@param h_start	行开始索引
		@param h_end	行结束索引
		*/
		const Matrix rowRange(int h_start, int h_end);
		/**
		@brief 返回w_start到w_end列矩阵
		@param w_start	列开始索引
		@param w_end	列结束索引
		*/
		const Matrix colRange(int w_start, int w_end);
		/**
		@brief 返回c_start到c_end通道矩阵
		@param c_start	通道开始索引
		@param c_end	通道结束索引
		*/
		const Matrix channelRange(int c_start, int c_end);
		/**
		@brief 返回深拷贝矩阵
		*/
		const lzh::Matrix clone()const;
		/**
		@brief 返回取反矩阵
		*/
		const lzh::Matrix opp()const;
		/**
		@brief 返回绝对值矩阵
		*/
		const lzh::Matrix abs()const;
		/**
		@brief 返回按num次幂矩阵
		@param num 次幂
		*/
		const lzh::Matrix mPow(int num)const;
		/**
		@brief 返回按num次幂矩阵
		@param num 次幂
		*/
		const lzh::Matrix pow(mat_t num)const;
		/**
		@brief 返回按元素取指数矩阵
		*/
		const lzh::Matrix exp()const;
		/**
		@brief 返回按元素取对数矩阵
		*/
		const lzh::Matrix log()const;
		/**
		@brief 返回按元素取开方矩阵
		*/
		const lzh::Matrix sqrt()const;
		/**
		@brief 返回伴随矩阵
		*/
		const lzh::Matrix adj()const;
		/**
		@brief 返回转置矩阵
		*/
		const lzh::Matrix t()const;
		/**
		@brief 返回逆矩阵
		*/
		const lzh::Matrix inv()const;
		/**
		@brief 返回逆矩阵
		*/
		const lzh::Matrix diag(int k = 0)const;
		/**
		@brief 返回逆序矩阵
		矩阵必须是向量
		*/
		const lzh::Matrix reverse()const;
		const lzh::Matrix EigenvectorsMax(mat_t offset = 1e-8)const;
		/**
		@brief sigmoid函数
		*/
		const lzh::Matrix sigmoid()const;
		/**
		@brief tanh函数
		*/
		const lzh::Matrix tanh()const;
		/**
		@brief relu函数
		*/
		const lzh::Matrix relu()const;
		/**
		@brief elu函数
		*/
		const lzh::Matrix elu()const;
		/**
		@brief selu函数
		*/
		const lzh::Matrix selu()const;
		/**
		@brief leaky_relu函数
		*/
		const lzh::Matrix leaky_relu()const;
		/**
		@brief softmax函数
		*/
		const lzh::Matrix softmax()const;
		/**
		@brief 返回行列式
		*/
		mat_t Det();
		/**
		@brief 返回num范数
		@param num 几范数
		*/
		mat_t Norm(int num = 1)const;
		/**
		@brief 返回对应索引的余子式
		@param x 列索引
		@param y 行索引
		*/
		mat_t Cof(int x, int y);
		mat_t EigenvalueMax(mat_t offset = 1e-8)const;
		/**
		@brief 返回随机抽取的矩阵元素
		*/
		mat_t RandSample();
		/**
		@brief 返回矩阵元素和
		@param num	设置次幂
		@param _abs 是否取绝对值
		*/
		mat_t sum(int num = 1, bool _abs = false)const;
		/**
		@brief 返回平均值
		*/
		mat_t mean()const;
		/**
		@brief 返回标准差
		*/
		mat_t S()const;
		/**
		@brief 返回标准差
		*/
		mat_t D()const;
		/**
		@brief 重载运算符+
		对应元素相加
		*/
		const lzh::Matrix operator + (const mat_t val)const;
		/**
		@brief 重载运算符+
		对应元素相加
		*/
		const lzh::Matrix operator + (const lzh::Matrix &a)const;
		/**
		@brief 重载运算符+=
		按元素相加
		*/
		void operator += (const mat_t val);
		/**
		@brief 重载运算符+=
		按元素相加
		*/
		void operator += (const lzh::Matrix &a);
		/**
		@brief 友元重载运算符+
		按元素相加
		*/
		friend const lzh::Matrix operator + (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief 重载运算符-
		按元素取相反数
		*/
		const lzh::Matrix operator - (void)const;
		/**
		@brief 重载运算符-
		按元素相减
		*/
		const lzh::Matrix operator - (const mat_t val)const;
		/**
		@brief 重载运算符-
		对应元素相减
		*/
		const lzh::Matrix operator - (const lzh::Matrix &a)const;
		/**
		@brief 重载运算符-=
		按元素相减
		*/
		void operator -= (const mat_t val);
		/**
		@brief 重载运算符-=
		对应元素相减
		*/
		void operator -= (const lzh::Matrix &a);
		/**
		@brief 友元重载运算符-
		按元素相减
		*/
		friend const lzh::Matrix operator - (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief 重载运算符*
		按元素相乘
		*/
		const lzh::Matrix operator * (const mat_t val)const;
		/**
		@brief 重载运算符*
		对应元素相乘
		*/
		const lzh::Matrix operator * (const lzh::Matrix &a)const;
		/**
		@brief 重载运算符*=
		按元素相乘
		*/
		void operator *= (const mat_t val);
		/**
		@brief 重载运算符*=
		对应元素相乘
		*/
		void operator *= (const lzh::Matrix &a);
		/**
		@brief 友元重载运算符*
		按元素相乘
		*/
		friend const lzh::Matrix operator * (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief 重载运算符/
		按元素相除
		*/
		const lzh::Matrix operator / (const mat_t val)const;
		/**
		@brief 重载运算符/
		矩阵乘法
		*/
		const lzh::Matrix operator / (const lzh::Matrix &a)const;
		/**
		@brief 重载运算符/=
		按元素相除
		*/
		void operator /= (const mat_t val);
		/**
		@brief 重载运算符/=
		对应元素相除
		*/
		void operator /= (const lzh::Matrix &a);
		/**
		@brief 友元重载运算符/
		按元素相乘
		*/
		friend const lzh::Matrix operator / (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief 重载运算符=
		深拷贝
		*/
		const Matrix & operator = (const lzh::Matrix &temp);
		/**
		@brief 重载运算符==
		判断矩阵是否相等
		*/
		bool operator == (const lzh::Matrix &a)const;
		/**
		@brief 重载运算符!=
		判断矩阵是否不相等
		*/
		bool operator != (const lzh::Matrix &a)const;
		/**
		@brief 返回对应索引元素
		@param w 索引
		*/
		mat_t& operator () (const int w)const;
		/**
		@brief 返回对应索引元素
		@param h 行索引
		@param w 列索引
		*/
		mat_t& operator () (const int h, const int w)const;
		/**
		@brief 返回对应索引元素
		@param h		行索引
		@param w		列索引
		@param c		通道索引
		*/
		mat_t& operator () (const int h, const int w, const int c)const;
		/**
		@brief 返回对应索引元素
		@param p 索引
		*/
		mat_t& operator () (Point p)const;
		/**
		@brief 返回对应索引元素
		@param p 索引
		*/
		mat_t& operator () (Point3i p)const;
		/**
		@brief 返回矩阵对应索引的列或行
		@param index	索引
		@param rc		索引方式
		*/
		const Matrix operator () (const int index, X_Y_Z rc)const;
		/**
		@brief 返回矩阵对应索引的列或行
		@param index	索引
		@param rc		索引方式
		*/
		const Matrix operator () (const int v1, const int v2, X_Y_Z rc)const;
		operator mat_t *() {
			return matrix;
		}
		operator const mat_t *() const {
			return matrix;
		}

		const Matrix operator [] (const int idx) const;

		friend std::ostream & operator << (std::ostream &out, const lzh::Matrix &ma);

		/**
		设置打印宽度
		@param w 宽度
		*/
		static void setPrintW(lzh::uint w);
		/**
		设置打印有效数字
		@param n 位数
		*/
		static void setPrintSignificantDigits(lzh::uint n);
		/**
		设置打印方法
		@param t 打印方法
		*/
		static void setPrintType(lzh::PrintType t);	
	protected:
		static lzh::uint print_width;
		static lzh::uint print_precision;
		static lzh::PrintType print_type;

		void init();
		void checkSquare();
#ifdef MAT_DEBUG
		void checkindex(int index);
		void checkindex(int index_x, int index_y);
#endif // MAT_DEBUG
		void setsize(int w, int h = 1, int c = 1);
	private:
		int h_;
		int w_;
		int c_;
		int dim;
		bool square;
		int step;
		MatPtr<mat_t> matrix;

	};

	typedef Matrix Mat;
	typedef mat_t(*NF)(mat_t);
	typedef mat_t(*Fun)(const Mat & x);
	typedef const Mat(*Activate)(const Mat &x);
	typedef const Mat(*F)(const Mat & a, const Mat & x);
	typedef const Mat(*DF)(const Mat & a, const Mat & x, const Mat & dy);
	typedef const Mat(*LF)(const Mat & y, const Mat & y0);

	/**
	@brief Mat_ 工具类
	继承Mat类，用于实现
	Mat mat = (Mat_(3, 3) <<
		-1, -1, -1,
		-1,  9, -1,
		-1, -1, -1);
	*/
	class Mat_ : public Mat
	{
	public:
		explicit Mat_() {}
		/**
		@brief 生成w的向量
		@param w		长度
		*/
		Mat_(int w) : Mat(w) {}
		/**
		@brief 生成h*w的矩阵
		@param w		列数
		@param h		行数
		*/
		Mat_(int w, int h) : Mat(w, h) {}
		/**
		@brief 生成h*w*c的张量
		@param w		列数
		@param h		行数
		@param depth	通道数
		*/
		Mat_(int w, int h, int c) : Mat(w, h, c) {}
		/**
		@brief 生成size的矩阵
		@param size		尺寸
		*/
		Mat_(Size size_) : Mat(size_) {}
		/**
		@brief 生成size的矩阵
		@param size		尺寸
		*/
		Mat_(Size3 size_) : Mat(size_) {}
		operator mat_t *() {
			return *this;
		}
		operator const mat_t *() const {
			return *this;
		}
	};
	/**
	@brief MatCommaInitializer_ 工具类
	作为迭代器，用于实现
	Mat mat = (Mat_(3, 3) <<
		-1, -1, -1,
		-1,  9, -1,
		-1, -1, -1);
	*/
	class MatCommaInitializer_
	{
	public:
		explicit MatCommaInitializer_() {}
		MatCommaInitializer_(const Mat_& m) {
			head = m.data();
			it = head;
			h = m.rows();
			w = m.cols();
			c = m.channels();
		}
		template<typename Tp_>
		MatCommaInitializer_ operator , (Tp_ v);
		int rows()const { return h; }
		int cols()const { return w; }
		int channels()const { return c; }
		mat_t * matrix()const { return head; }
	private:
		int h;
		int w;
		int c;
		mat_t *it;
		mat_t *head;
	};
	template<typename Tp_>
	inline MatCommaInitializer_ MatCommaInitializer_::operator , (Tp_ v)
	{
#ifdef MAT_DEBUG
		if (this->it == this->head + h * w*c) {
			fprintf(stderr, "%s\n", errinfo[ERR_INFO_MEMOUT]);
			throw std::exception(errinfo[ERR_INFO_MEMOUT]);
		}
#endif
		*this->it = _T(v);
		++this->it;
		return *this;
	}

	template<typename Tp_>
	static inline MatCommaInitializer_ operator << (const Mat_ &m, Tp_ val)
	{
		MatCommaInitializer_ commaInitializer(m);
		return (commaInitializer, val);
	}

	/**
	设置随机数种子
	*/
	extern void Srandom();
	/**
	@brief 返回元素为0的w向量
	@param w		向量长度
	*/
	extern const Mat zeros(int w);
	/**
	@brief 返回元素为0的h*w矩阵
	@param w		矩阵列数
	@param h		矩阵行数
	*/
	extern const Mat zeros(int w, int h);
	/**
	@brief 返回元素为0的h*w*c矩阵
	@param w		矩阵列数
	@param h		矩阵行数
	@param c		矩阵通道数
	*/
	extern const Mat zeros(int w, int h, int c);
	/**
	@brief 返回元素为0的size矩阵
	@param size 矩阵大小
	*/
	extern const Mat zeros(Size size);
	/**
	@brief 返回元素为0的size矩阵
	@param size 矩阵大小
	*/
	extern const Mat zeros(Size3 size);
	/**
	@brief 返回元素为v的w向量
	@param v		填充元素
	@param w		向量长度
	*/
	extern const Mat value(mat_t v, int w);
	/**
	@brief 返回元素为v的h*w矩阵
	@param v		填充元素
	@param w		矩阵列数
	@param h		矩阵行数
	*/
	extern const Mat value(mat_t v, int w, int h);
	/**
	@brief 返回元素为v的h_*w_*c_矩阵
	@param v		填充元素
	@param w		矩阵列数
	@param h		矩阵行数
	@param c		矩阵通道数
	*/
	extern const Mat value(mat_t v, int w, int h, int c);
	/**
	@brief 返回元素为1的w向量
	@param w		向量长度
	*/
	extern const Mat ones(int w);
	/**
	@brief 返回元素为1的h*w矩阵
	@param w	矩阵列数
	@param h	矩阵行数
	*/
	extern const Mat ones(int w, int h);
	/**
	@brief 返回元素为1的h*w*c矩阵
	@param w		矩阵列数
	@param h		矩阵行数
	@param c		矩阵通道数
	*/
	extern const Mat ones(int w, int h, int c);
	/**
	@brief 返回元素为1的size矩阵
	@param size 矩阵大小
	*/
	extern const Mat ones(Size size);
	/**
	@brief 返回元素为1的size矩阵
	@param size 矩阵大小
	*/
	extern const Mat ones(Size3 size);
	/**
	@brief 返回从low到top等分成的1*len的矩阵
	@param low 下界
	@param top 上界
	@param gap 间距
	*/
	extern const Mat range(int low, int top, mat_t gap = _T(1));
	/**
	@brief 返回从low到top等分成的1*len的矩阵
	@param low 下界
	@param top 上界
	@param gap 间距
	*/
	extern const Mat range(mat_t low, mat_t top, mat_t gap = _T(1));
	/**
	@brief 返回从low到top等分成的1*len的矩阵
	@param low 下界
	@param top 上界
	@param len 等分个数
	*/
	extern const Mat linspace(int low, int top, int len);
	/**
	@brief 返回从low到top等分成的1*len的矩阵
	@param low 下界
	@param top 上界
	@param len 等分个数
	*/
	extern const Mat linspace(mat_t low, mat_t top, int len);
	/**
	@brief 返回随机矩阵
	符合高斯分布
	大小为size
	@param size		矩阵大小
	@param n1		先前通道
	@param n2		当前通道
	*/
	extern const Mat Xavier(Size3 size, int n1, int n2);
	/**
	@brief 返回随机矩阵
	符合高斯分布
	大小为h_*w_*c_
	@param h_		矩阵行数
	@param w_		矩阵列数
	@param c_	矩阵通道数
	@param n1		先前通道
	@param n2		当前通道
	*/
	extern const Mat Xavier(int w, int h, int c, int n1, int n2);
	/**
	@brief 返回随机矩阵
	符合高斯分布
	大小为size, 元素范围[0, 1]
	@param size		矩阵大小
	*/
	extern const Mat Random(Size3 size);
	/**
	@brief 返回随机矩阵
	符合高斯分布
	大小为h_*w_*c_, 元素范围[0, 1]
	@param h_		矩阵行数
	@param w_		矩阵列数
	@param c_	矩阵通道数
	*/
	extern const Mat Random(int w, int h = 1, int c = 1);
	/**
	@brief 返回矩阵的最大值
	@param src		矩阵
	@param isAbs	是否取绝对值
	*/
	extern mat_t Max(const Mat &src, bool isAbs = false);
	/**
	@brief 返回矩阵的最小值
	@param src		矩阵
	@param isAbs	是否取绝对值
	*/
	extern mat_t Min(const Mat &src, bool isAbs = false);
	/**
	@brief 返回矩阵的行列式
	@param src	矩阵
	*/
	extern mat_t det(const Mat &src);
	/**
	@brief 返回矩阵的迹
	@param src	矩阵
	*/
	extern mat_t trace(const Mat &src);
	/**
	@brief 返回矩阵对应索引的余子式
	@param src	矩阵
	@param x	列索引
	@param y	行索引
	*/
	extern mat_t cof(const Mat &src, int x, int y);
	/**
	@brief 返回随机数
	@param min		最小值
	@param max		最大值
	@param isdouble 是否随机非整数
	*/
	extern mat_t getRandData(int min, int max, bool isdouble = false);
	/**
	@brief 返回矩阵范数
	@param src 矩阵
	@param num 几范数
	*/
	extern mat_t mNorm(const Mat& src, int num = 1);
	/**
	@brief 返回矩阵的距离
	@param a	矩阵
	@param b	矩阵
	@param num	几范数
	*/
	extern mat_t mDistance(const Mat& a, const Mat& b, int num = 2);
	/**
	@brief 返回随机的矩阵元素
	@param src 矩阵
	*/
	extern mat_t mRandSample(const Mat &src);
	/**
	@brief 返回生成的n*n*1单位矩阵
	@param n 矩阵大小
	*/
	extern const Mat eye(int n);
	/**
	@brief 返回矩阵的第c_的通道
	@param src		矩阵
	@param c_	通道索引
	*/
	extern const Mat mSplit(const Mat &src, int c);
	/**
	@brief 返回按矩阵的通道数复制
	@param src 输入矩阵
	@param dst 输入矩阵通道数的矩阵数组
	*/
	extern void mSplit(const Mat &src, Mat *dst);
	/**
	@brief 返回按通道合并的矩阵
	@param src		矩阵数组
	@param channels 通道数
	*/
	extern const Mat mMerge(const Mat *src, int channels);
	/**
	@brief 返回按索引区域切割的矩阵
	@param src			矩阵
	@param Row_Start	截取行初始索引值
	@param Row_End		截取行结束索引值
	@param Col_Start	截取列初始索引值
	@param Col_End		截取列结束索引值
	*/
	extern const Mat Block(const Mat &src, int Row_Start, int Row_End, int Col_Start, int Col_End, int Chennel_Start = 0, int Chennel_End = 0);
	/**
	@brief 返回随机生成元素n*n*1矩阵
	@param n		矩阵大小
	@param low		下界
	@param top		上界
	@param isdouble 是否生成浮点数
	*/
	extern const Mat mRand(int low, int top, int n, bool isdouble = false);
	/**
	@brief 返回随机生成元素size.x*size.y*size.z矩阵
	@param low		下界
	@param top		上界
	@param size		矩阵大小
	@param isdouble 是否生成浮点数
	*/
	extern const Mat mRand(int low, int top, Size3 size, bool isdouble = false);
	/**
	@brief 返回随机生成元素h_*w_*c_矩阵
	@param h_		矩阵行数
	@param w_		矩阵列数
	@param low		下界
	@param top		上界
	@param isdouble 是否生成浮点数
	*/
	extern const Mat mRand(int low, int top, int w, int h, int c = 1, bool isdouble = false);
	/**
	@brief 返回大小为w向量
	@param w		向量长度
	*/
	extern const Mat mcreate(int w);
	/**
	@brief 返回大小为h*w矩阵
	@param w		矩阵列数
	@param h		矩阵行数
	*/
	extern const Mat mcreate(int w, int h);
	/**
	@brief 返回大小为h*w*c矩阵
	@param w		矩阵列数
	@param h		矩阵行数
	@param c		矩阵通道数
	*/
	extern const Mat mcreate(int w, int h, int c);
	/**
	@brief 返回大小为size矩阵
	@param size 矩阵大小
	*/
	extern const Mat mcreate(Size size);
	/**
	@brief 返回大小为size矩阵
	@param size 矩阵大小
	*/
	extern const Mat mcreate(Size3 size);
	/**
	@brief 返回逆序矩阵，矩阵需为一维向量
	@param src	矩阵
	*/
	extern const Mat reverse(const Mat &src);
	/**
	@brief 返回生成h*w*c的矩阵，随机抽取矩阵src元素填充
	@param src	矩阵
	@param h	矩阵行数
	@param w	矩阵列数
	@param c	矩阵通道数
	*/
	extern const Mat mRandSample(const Mat &src, int w, int h, int c = 1);
	/**
	@brief 返回随机抽取num次矩阵src的行或列组成的矩阵
	@param src	矩阵
	@param rc	抽取方式
	@param num	抽取次数
	*/
	extern const Mat mRandSample(const Mat &src, X_Y_Z rc, int num = 1);
	/**
	@brief 返回矩阵的伴随矩阵
	@param src 矩阵
	*/
	extern const Mat adj(const Mat &src);
	/**
	@brief 返回矩阵的逆矩阵
	@param src 矩阵
	*/
	extern const Mat inv(const Mat &src);
	/**
	@brief 返回对角线矩阵
	@param src	向量
	@param k	第k条对角线
	*/
	extern const Mat diag(const Mat &src, int k = 0);
	/**
	@brief 返回矩阵的伪逆矩阵
	@param src	矩阵
	@param dire 伪逆矩阵的计算方式
	*/
	extern const Mat pinv(const Mat &src, Dire dire = LEFT);
	/**
	@brief 返回矩阵的转置矩阵
	@param src 矩阵
	*/
	extern const Mat tran(const Mat &src);
	/**
	@brief 返回矩阵的绝对值矩阵
	@param src 矩阵
	*/
	extern const Mat mAbs(const Mat &src);
	/**
	@brief 返回angle度2*2的旋转矩阵
	@param angle 角度
	*/
	extern const Mat Rotate(mat_t angle);
	/**
	@brief 返回矩阵num次幂
	@param src 矩阵
	@param num 次幂
	*/
	extern const Mat POW(const Mat &src, int num);
	/**
	@brief 返回矩阵取反
	@param src 矩阵
	*/
	extern const Mat mOpp(const Mat &src);
	/**
	@brief 返回矩阵按行或列之和
	@param src	矩阵
	@param rc	求和的方向
	*/
	extern const Mat mSum(const Mat &src, X_Y_Z rc);
	/**
	@brief 返回矩阵按元素取指数
	@param src 矩阵
	*/
	extern const Mat mExp(const Mat &src);
	/**
	@brief 返回矩阵按元素取对数
	@param src 矩阵
	*/
	extern const Mat mLog(const Mat &src);
	/**
	@brief 返回矩阵按元素取开方
	@param src 矩阵
	*/
	extern const Mat mSqrt(const Mat &src);
	/**
	@brief 返回矩阵按元素取num次幂
	@param src 矩阵
	@param num 次幂
	*/
	extern const Mat mPow(const Mat &src, mat_t num);
	/**
	@brief 返回矩阵val/src按元素除
	@param src 矩阵
	@param val 除数
	*/
	extern const Mat Divi(const Mat &src, mat_t val, Dire dire = RIGHT);
	/**
	@brief 返回矩阵除法
	@param a	被除矩阵
	@param b	除矩阵
	@param dire 除法方式
	*/
	extern const Mat Divi(const Mat &a, const Mat &b, Dire dire = RIGHT);
	/**
	@brief 返回哈达玛积
	@param a 矩阵
	@param b 矩阵
	*/
	extern const Mat Mult(const Mat &a, const Mat &b);
	/**
	@brief 返回矩阵乘法
	@param a 矩阵
	@param b 矩阵
	*/
	extern const Mat Dot(const Mat &a, const Mat &b);
	/**
	@brief 返回矩阵按元素取a和b之间的最大值
	@param a 比较值
	@param b 比较矩阵
	*/
	extern const Mat mMax(mat_t a, const Mat &b);
	/**
	@brief 返回矩阵按元素取a和b之间的最大值
	@param a 比较矩阵
	@param b 比较矩阵
	*/
	extern const Mat mMax(const Mat &a, const Mat &b);
	/**
	@brief 返回矩阵按元素取a和b之间的最小值
	@param a 比较值
	@param b 比较矩阵
	*/
	extern const Mat mMin(mat_t a, const Mat &b);
	/**
	@brief 返回矩阵按元素取a和b之间的最小值
	@param a 比较矩阵
	@param b 比较矩阵
	*/
	extern const Mat mMin(const Mat &a, const Mat &b);
	/**
	@brief mCalSize 计算卷积所需扩张的边界
	返回矩阵大小
	@param src		被卷积矩阵
	@param kern		卷积核
	@param anchor	像素对应卷积核坐标
	anchor默认为Point(-1,-1), 像素对应卷积核中心
	@param strides	滑动步长
	@param top		向上扩充几行
	@param bottom	向下扩充几行
	@param left		向左扩充几列
	@param right	向右扩充几列
	*/
	extern Size3 mCalSize(Size3 src, Size3 kern, Point anchor, Size strides, int &top, int &bottom, int &left, int &right);
	/**
	@brief mCalSize 计算卷积所需扩张的边界
	返回矩阵大小
	@param src		被卷积矩阵尺寸
	@param kern		卷积核尺寸
	@param anchor	像素对应卷积核坐标
	*/
	extern Size3 mCalSize(Size3 src, Size3 kern, Point &anchor, Size strides);
	/**
	@brief mCalSize 计算卷积所需扩张的边界
	返回矩阵大小
	@param src		被卷积矩阵尺寸
	@param kern		卷积核尺寸
	@param anchor	像素对应卷积核坐标
	*/
	extern Size3 mCalSize(Size3 src, Size kern, Point &anchor, Size strides);
	/**
	@brief 返回按boundary分界填充的矩阵
	返回矩阵大小等于输入矩阵大小
	@param src				输入矩阵
	@param boundary			分界值
	@param lower			小于boundary用lower填充
	@param upper			大于boundary用upper填充
	@param boundary2upper	当矩阵元素等于boundary时
	为1归upper,				为-1归lower, 为0不处理
	*/
	extern const Mat Threshold(const Mat &src, mat_t boundary, mat_t lower, mat_t upper, int boundary2upper = 1);
	/**
	@brief 返回边界扩充的矩阵
	@param src			输入矩阵
	@param top			向上扩充几行
	@param bottom		向下扩充几行
	@param left			向左扩充几列
	@param right		向右扩充几列
	@param borderType	边界像素外推的插值方法
	@param value		常量插值的数值
	**/
	extern const Mat copyMakeBorder(const Mat &src, int top, int bottom, int left, int right, BorderTypes borderType = BORDER_CONSTANT, mat_t value = 0.0);
	/**
	@brief 返回矩阵2维卷积结果(只支持二维输入)
	返回矩阵大小为(input.h_/strides_x, input.w_/strides_y, 1)
	@param input			输入矩阵
	@param kern				卷积核
	@param anchor			矩阵元素对应卷积核的位置
	以卷积核的左上角为(0,0)点, 默认(-1,-1)为中心
	@param strides			滑动步长
	Size.hei为x轴,Size.wid为y轴
	@param is_copy_border	是否要扩展边界
	*/
	extern const Mat Filter2D(const Mat & input, const Mat & kern, Point anchor = Point(-1, -1), const Size & strides = Size(1, 1), bool is_copy_border = true);
	/**
	@brief 返回矩阵2维卷积结果(支持三维输入)
	返回矩阵大小为(input.h_/strides_x, input.w_/strides_y, 1)
	@param input			输入矩阵
	@param kern				卷积核
	@param anchor			矩阵元素对应卷积核的位置
	以卷积核的左上角为(0,0)点, 默认(-1,-1)为中心
	@param strides			滑动步长
	Size.hei为x轴,Size.wid为y轴
	@param is_copy_border	是否要扩展边界
	*/
	extern void Filter2D(const Mat & in, Mat & out, const Mat & kern, Point anchor = Point(-1, -1), const Size & strides = Size(1, 1), bool is_copy_border = true);
	/**
	@brief 返回拟合结果
	最小二乘法
	@param x 自变量
	@param y 因变量
	*/
	extern const Mat LeastSquare(const Mat& x, const Mat &y);
	/**
	@brief 返回a
	线性拟合y=a(0)*x+a(1)
	y和x必须为行向量
	保证y和x的行数相同
	@param x 自变量
	@param y 因变量
	*/
	extern const Mat regress(const Mat& x, const Mat &y);
	/**
	@brief 返回P
	多项式拟合y=P(1)*x^n + P(2)*x^(n-1) +...+ P(n)*x + P(n+1)
	y和x必须为行向量
	保证y和x的行数相同
	@param x 自变量
	@param y 因变量
	*/
	extern const Mat polyfit(const Mat& x, const Mat &y, uint n);
	/**
	@brief 返回拟合结果
	拟合圆
	p必须为行向量
	@param p 点集
	*/
	extern const Mat circlefit(const Mat& p);
	/**
	@brief 返回拟合结果
	非线性最小二乘法
	y和x必须为行向量
	保证y和x的行数相同
	@param x		自变量
	@param y		因变量
	@param a0		初始参数
	@param f		函数指针 f(a, x) = y
	@param step		更新步长
	@param error	误差(小于误差结束更新)
	*/
	extern const Mat NonLinearLeastSqures(const Mat & x, const Mat & y, const Mat & a0, F f, mat_t step = _T(1e-2), mat_t error = _T(1e-6));
	/**
	@brief 返回梯度
	@param y	因变量
	@param x	自变量(缺省差值为1)
	*/
	extern const Mat gradient(const Mat & y, const Mat & x = Mat());
	/**
	@brief 返回数值梯度
	@param f		函数指针 f(x) = y
	@param x		自变量
	@param epsilon	差值
	*/
	extern mat_t NumericGradient(NF f, mat_t x, mat_t epsilon = _T(1e-3));
	/**
	@brief 返回数值梯度
	@param f		函数指针 f(x) = y
	@param x		自变量
	@param epsilon	差值
	*/
	extern const Mat NumericGradient(Fun f, const Mat & x, mat_t epsilon = _T(1e-3));
	/**
	@brief 返回数值梯度
	@param f		函数指针 f(a, x) = y
	@param a		自变量
	@param x		输入
	@param epsilon	差值
	*/
	extern const Mat NumericGradient(F f, const Mat & a, const Mat & x, mat_t epsilon = _T(1e-3));
	/**
	@brief 返回积分值
	使用欧拉法求数值积分
	@param f		函数指针 f(x) = y
	@param low		积分下限
	@param high		积分上限
	@param epsilon	采样间隔
	*/
	extern mat_t EulerInt(NF f, mat_t low, mat_t high, mat_t epsilon = _T(1e-3));
	/**
	@brief 返回积分值
	使用梯形法求数值积分
	@param f		函数指针 f(x) = y
	@param low		积分下限
	@param high		积分上限
	@param epsilon	采样间隔
	*/
	extern mat_t TrapezoidInt(NF f, mat_t low, mat_t high, mat_t epsilon = _T(1e-3));
	/**
	@brief 返回积分值
	使用四阶龙格-库塔法求数值积分
	@param f		函数指针 f(x) = y
	@param low		积分下限
	@param high		积分上限
	@param epsilon	采样间隔
	*/
	extern mat_t RungeKuttaInt(NF f, mat_t low, mat_t high, mat_t epsilon = _T(1e-3));
	extern const Mat LinearIntersection(const Mat & line_1, const Mat & line_2);
	/**
	@brief 返回y
	多项式拟合y=P(1)*x^n + P(2)*x^(n-1) +...+ P(n)*x + P(n+1)
	x必须为行向量
	@param a 参数
	@param x 自变量
	*/
	extern const Mat polynomial(const Mat& a, const Mat &x);
	/**
	@brief 返回修改维度的矩阵
	不允许改变矩阵长度
	@param src	输入
	@param size 矩阵大小
	*/
	extern const Mat Reshape(const Mat& src, Size3 size);
	/**
	@brief 返回矩阵
	按通道求和
	@param src	输入
	*/
	extern const Mat SumChannel(const Mat& src);
	/**
	@brief 返回旋转矩阵
	@param src	输入
	@param dice 旋转角度
	*/
	extern const Mat rotate(const Mat& src, RotateAngle dice);
	/**
	按比例缩放矩阵
	@param src		输入
	@param dst		输出
	@param xRatio	x缩放比例
	@param yRatio	y缩放比例
	@param mothed	缩放方法
	*/
	extern void resize(const Mat & src, Mat & dst, mat_t xRatio, mat_t yRatio, ReductionMothed mothed);
	/**
	按大小缩放矩阵
	@param src		输入
	@param dst		输出
	@param newSize	新的矩阵大小
	@param mothed	缩放方法
	*/
	extern void resize(const Mat& src, Mat& dst, Size newSize, ReductionMothed mothed);
	/**
	@brief 返回随机数
	符合高斯分布的随机数
	@param mu		期望值
	@param sigma	标准差
	*/
	extern mat_t generateGaussianNoise(mat_t mu = 0, mat_t sigma = 1);
	/**
	@brief 返回聚类结果
	K均值聚类
	@param P			输入
	@param k			对点集的分类结果
	@param K			聚类个数
	@param iteration	迭代次数
	@param error		误差(小于误差结束迭代)
	*/
	extern const Mat kmeans(const Mat &P, Mat &k, const uint K, const uint iteration, const mat_t error = _T(1e-7));
	/**
	@brief 返回行秩
	计算矩阵的行最简矩阵
	@param src		输入
	@param dst		行最简矩阵
	*/
	extern int RowSimplest(const Mat & src, Mat & dst);	
	/**
	@brief 返回行最简矩阵
	计算矩阵的行最简矩阵
	@param src		输入
	*/
	extern const Mat RowSimplest(const Mat & src);
	/**
	@brief 返回列秩
	计算矩阵的列最简矩阵
	@param src		输入
	@param dst		列最简矩阵
	*/
	extern int ColSimplest(const Mat & src, Mat & dst);
	/**
	@brief 返回列最简矩阵
	计算矩阵的列最简矩阵
	@param src		输入
	@param dst		列最简矩阵
	*/
	extern const Mat ColSimplest(const Mat & src);
	/**
	@brief 返回求解状态
	求解非齐次线性方程组
	齐次线性方程组使用addzeros手动添加0
	@param src		输入线性方程组
	@param dst		基础解系
	@param simplest	最简行矩阵
	@param mark		解的状态(自由解为0,特解为1,无解为nullptr)
	*/
	extern EQUATION SolveLinearEquation(const Mat & src, Mat & dst, Mat *simplest = nullptr, Mat *mark = nullptr);
	/**
	排序算法类
	提供了
	冒泡排序	=> bubble;
	插入排序	=> insert;
	选择排序	=> select;
	梳排序		=> comb;
	地精排序	=> gnome;
	堆排序		=> heap;
	希尔排序	=> shell;
	快速排序	=> quick;
	归并排序	=> merge;
	*/
	class sort
	{
	public:
		/**
		冒泡排序
		@param m		排序数据
		@param order	顺序
		*/
		static void bubble(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		冒泡排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat bubble(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		插入排序
		@param m		排序数据
		@param order	顺序
		*/
		static void insert(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		插入排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat insert(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		选择排序
		@param m		排序数据
		@param order	顺序
		*/
		static void select(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		选择排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat select(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		梳排序
		@param m		排序数据
		@param order	顺序
		*/
		static void comb(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		梳排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat comb(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		地精排序
		@param m		排序数据
		@param order	顺序
		*/
		static void gnome(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		地精排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat gnome(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		堆排序
		@param m		排序数据
		@param order	顺序
		*/
		static void heap(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		堆排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat heap(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		希尔排序
		@param m		排序数据
		@param order	顺序
		*/
		static void shell(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		希尔排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat shell(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		快速排序
		@param m		排序数据
		@param order	顺序
		*/
		static void quick(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		快速排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat quick(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		归并排序
		@param m		排序数据
		@param order	顺序
		*/
		static void merge(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief 返回排序结果
		归并排序
		@param m		排序数据
		@param order	顺序
		*/
		static const Mat merge(const Mat & m, ORDER order = MIN_TO_MAX);
		template<typename Tp_>
		static inline void swap(Tp_ &a, Tp_ &b)
		{
			Tp_ temp = a;
			a = b;
			b = temp;
		}
	protected:
		static void heapdown(mat_t *m, int i, int n, ORDER order);
		static void heapup(mat_t *m, int i, ORDER order);
		static void makeheap(mat_t *m, int length, ORDER order);
		static void heaparray(mat_t *m, int i, int n, ORDER order);
		static void mergearray(mat_t*a, mat_t *b, int start, int mid, int end, ORDER order);
		static void _merge(mat_t*a, mat_t *b, int start, int end, ORDER order);
		static void _quick(mat_t *m, int start, int end, ORDER order);
	};
	/**
	B-样条曲线类
	*/
	class BSpline
	{
	public:
		/**
		UNIFORM			均匀B样条曲线
		QUASI_UNIFORM	准均匀B样条曲线
		*/
		enum BSplineType
		{
			UNIFORM	= 0,
			QUASI_UNIFORM
		};

		explicit BSpline();
		BSpline(BSplineType type, int k = 1, const Mat &p = Mat());
		/**
		@brief 返回B样条曲线在t的点
		t的定义域在[0,1]
		@param t	参数方程自变量
		*/
		void setCtrlPoint(const Mat &p);
		/**
		注册B样条曲线节点向量
		@param node	节点向量
		缺省时按曲线类型自动生成
		*/
		void NodeVector(const Mat &node = Mat());
		/**
		@brief 返回B样条曲线控制点
		*/
		const Mat CtrlPoint()const;
		/**
		@brief 返回B样条曲线节点向量
		*/
		const Mat Node()const;
		/**
		@brief 返回B样条曲线点集
		@param T	参数方程自变量
		T为定义域在[0,1]的递增序列
		*/
		const Mat BPiont(const Mat & T)const;
		/**
		@brief 返回B样条曲线点集
		@param number	[0,1]等间隔number等分
		*/
		const Mat BPoint(int number)const;
		/**
		@brief 返回B样条曲线在t的点
		t的定义域在[0,1]
		@param t	参数方程自变量
		*/
		const Mat operator ()(mat_t t)const;
		/**
		@brief 返回B样条曲线
		点集拟合B样条曲线
		@param P	点集
		@param n	控制点数量
		@param k	k阶B样条曲线
		*/
		static BSpline fitBSpline(const Mat &P, int n, int k);
	protected:
		mat_t BF(int i, int k, mat_t t)const;
		const Mat BaseFunction(mat_t t)const;
	private:
		int k;
		Mat P;
		Mat nodevector;
		BSplineType type;
	};

	/**
	@brief softmax函数
	Si = exp(y - max(y)) / sum(exp(y - max(y)))
	*/
	extern const Mat Softmax(const Mat &y);
	/**
	@brief L2函数
	E = (y - y0)^2
	*/
	extern const Mat L2(const Mat & y, const Mat & y0);
	/**
	@brief quadratic函数
	E = 1/2 * (y - y0)^2
	*/
	extern const Mat Quadratic(const Mat & y, const Mat & y0);
	/**
	@brief crossentropy函数
	E = -(y * log(y0))
	*/
	extern const Mat CrossEntropy(const Mat & y, const Mat & y0);
	/**
	@brief softmax + crossentropy函数
	E = -(y * log(softmax(y0)))
	*/
	extern const Mat SoftmaxCrossEntropy(const Mat & y, const Mat & y0);
	/**
	@brief sigmoid函数
	y = 1/(1 + exp(-x))
	*/
	extern const Mat Sigmoid(const Mat & x);
	/**
	@brief tanh函数
	y = (exp(x) - exp(-x)) / (exp(x) + exp(-x))
	*/
	extern const Mat Tanh(const Mat & x);
	/**
	@brief relu函数
	y = {x if x > 0; 0 if x < 0}
	*/
	extern const Mat ReLU(const Mat & x);
	/**
	@brief elu函数
	y = {x if x > 0; a*(exp(x) - 1) if x < 0}
	*/
	extern const Mat ELU(const Mat & x);
	/**
	@brief selu函数
	y = scale*Elu(x)
	*/
	extern const Mat SELU(const Mat & x);
	/**
	@brief leaky relu函数
	y = {x if x > 0; a*x if x < 0}
	*/
	extern const Mat LReLU(const Mat & x);

	class tools
	{
	public:
		/**
		暂停程序
		*/
		static void pause();
		/**
		命令行输出
		*/
		template<typename Tp_>
		static void print(const Tp_ &t) {
			std::cout << t << std::endl;
		}
		static void Wait(uint ms);
		static void Frequency();
		static void StartCounter();
		static mat_t EndCounter();
		static std::vector<std::string> strsplit(std::string &str, char ch);
		static void getFiles(std::string path, std::vector<std::string>& files);
		static std::vector<lzh::mat_t> str2mat_t(const std::vector<std::string> &str);
		static const Mat readcsv(std::string filename);
		/**
		@brief 返回将向量转换成1*point.size()*1的矩阵
		@param point 向量
		*/
		static const Mat Vec2Mat(std::vector<mat_t> &point);
		/**
		@brief 返回将向量转换成point.size()*point[0].size()*1的矩阵
		@param points 向量
		*/
		static const Mat Vec2Mat(std::vector<std::vector<mat_t>> &points);
		template<class Type>
		static inline std::vector<Point2<Type>> Mat2Point(const Mat & m)
		{
			std::vector<Point2<Type>> vec;
			int l1, l2, flag = -1;
			if (m.rows() == 2) {
				l1 = m.cols();
				l2 = m.channels();
				flag = 0;
			}
			else if (m.cols() == 2) {
				l1 = m.rows();
				l2 = m.channels();
				flag = 1;
			}
			else if (m.channels() == 2) {
				l1 = m.rows();
				l2 = m.cols();
				flag = 2;
			}
			else
				return vec;
			vec.resize(l1*l2);
			for (int j = 0; j < l1; j++)
				for (int k = 0; k < l2; k++)
					switch (flag)
					{
					case 0:vec[j*l2 + k] = Point2<Type>(m(0, j, k), m(1, j, k)); break;
					case 1:vec[j*l2 + k] = Point2<Type>(m(j, 0, k), m(j, 1, k)); break;
					case 2:vec[j*l2 + k] = Point2<Type>(m(j, k, 0), m(j, k, 1)); break;
					default:return vec;
					}
			return vec;
		}
		template<class Type>
		static inline std::vector<Point3<Type>> Mat2Point3(const Mat & m)
		{
			std::vector<Point3<Type>> vec;
			int l1, l2, flag = -1;
			if (m.rows() == 3) {
				l1 = m.cols();
				l2 = m.channels();
				flag = 0;
			}
			else if (m.cols() == 3) {
				l1 = m.rows();
				l2 = m.channels();
				flag = 1;
			}
			else if (m.channels() == 3) {
				l1 = m.rows();
				l2 = m.cols();
				flag = 2;
			}
			else
				return vec;
			vec.resize(l1*l2);
			for (int j = 0; j < l1; j++)
				for (int k = 0; k < l2; k++)
					switch (flag)
					{
					case 0:vec[j*l2 + k] = Point3<Type>(m(0, j, k), m(1, j, k), m(2, j, k)); break;
					case 1:vec[j*l2 + k] = Point3<Type>(m(j, 0, k), m(j, 1, k), m(j, 2, k)); break;
					case 2:vec[j*l2 + k] = Point3<Type>(m(j, k, 0), m(j, k, 1), m(j, k, 2)); break;
					default:return vec;
					}
			return vec;
		}
		template<class Type>
		static inline const Mat Point2Mat(const std::vector<Point3<Type>> & ps)
		{
			if (ps.empty())return Mat();
			Mat m(ps.size(), 3);
			int idx = 0;
			for (const Point3<Type> &p : ps)
			{
				Mat(Mat_(3) << p.x, p.y, p.z).copyTo(m.Col(idx++));
			}
		}
		template<class Type>
		static inline const Mat Point2Mat(const std::vector<Point2<Type>> & ps)
		{
			if (ps.empty())return Mat();
			Mat m(ps.size(), 2);
			int idx = 0;
			for (const Point2<Type> &p : ps)
			{
				Mat(Mat_(2) << p.x, p.y).copyTo(m.Col(idx++));
			}
		}
		/**
		@brief 返回将矩阵转换成一维向量
		@param src 矩阵
		*/
		static std::vector<mat_t> Mat2Vec(const Mat &src);
		/**
		@brief 返回将矩阵转换成向量
		@param src 矩阵
		*/
		static std::vector<std::vector<mat_t>> Mat2Vecs(const Mat &src);
		static void check(int w, int h = 1, int depth = 1);
		static std::string createfile(std::string filename);
		static std::string createtype(std::string filename);
		static std::string show_path();
		static std::string Enum2String(PrintType type);
		static std::string Enum2String(BorderTypes type);
		static std::string Enum2String(MatErrorInfo type);
		static std::string Enum2String(EQUATION type);
		static std::string Enum2String(ORDER type);
		static std::string Enum2String(X_Y_Z type);
		static std::string Enum2String(Dire type);
		static std::string Enum2String(ReductionMothed type);
		static std::string Enum2String(RotateAngle type);
		static std::string Enum2String(BSpline::BSplineType type);
	};
}
#endif //  __MAT_H__
