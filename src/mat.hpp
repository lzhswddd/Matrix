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
	//LeakyReLU�ĳ���
	static const mat_t LReLU_alpha = _T(0.2);
	//ELU�ĳ���
	static const mat_t ELU_alpha = _T(1.6732632423543772848170429916717);
	//SELU�ĳ���
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
		"error 0: ����Ϊ��!\0",
		"error 1: �����Ƿ���!\0",
		"error 2: �����Ƿ��󣬲������ð������!\0",
		"error 3: �����Ƿ��󣬲������������!\0",
		"error 4: �����Ƿ��󣬲��ܽ��д�������!\0",
		"error 5: �����Ƿ��󣬲�������Ϊ��λ����!\0",
		"error 6: ��������!\0",
		"error 7: ����û��ʵ������ֵ!\0",
		"error 8: ����ά��Ϊ0!\0",
		"error 9: ������������!\0",
		"error 10: ����������Ч!\0",
		"error 11: ��������ά�Ȳ�һ��!\0",
		"error 12: ��������ά�Ȳ�����˷�����!\0",
		"error 13: ����ά�Ȳ�Ϊ1����������!\0",
		"error 14: ����Υ��!\0",
		"error 15: ���������ʧ��!\0",
		"error 16: ����ʽΪ0!\0",
		"error 17: ��֧����ά����!\0",
		"error 18: ָ��Ϊ��!\0",
		"error 19: �������ά�ȱ���Ϊ2D!\0",
		"error 20: û���ҵ��ļ�!\0"
	};
	/**
	SPECIAL_SOLUTION	�������ؽ�
	GENERAL_SOLUTION	������ͨ��
	NO_SOLUTION			�����޽�
	*/
	enum EQUATION
	{
		SPECIAL_SOLUTION = 0,	//�������ؽ�
		GENERAL_SOLUTION,		//������ͨ��
		NO_SOLUTION				//�����޽�
	};
	/**
	MIN_TO_MAX	��С����
	MAX_TO_MIN	�Ӵ�С
	*/
	enum ORDER
	{
		MIN_TO_MAX = 0,
		MAX_TO_MIN
	};
	/**
	ROW		��
	COL		��	
	CHANNEL	ͨ��
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
	EqualIntervalSampling �ȼ������
	LocalMean �ֲ���ֵ
	*/
	enum ReductionMothed
	{
		EqualIntervalSampling = 0,
		LocalMean
	};
	/**
	��ת����˳ʱ��
	ROTATE_90_ANGLE 90��
	ROTATE_180_ANGLE 180��
	ROTATE_270_ANGLE 270��
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
	�ڴ������
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
		����w_������
		@param w		��������
		*/
		Matrix(int w);
		/**
		����w_*h_�ľ���
		@param w		��������
		@param h		��������
		*/
		Matrix(int w, int h);
		/**
		����h_*w_*depth�ľ���
		@param w		��������
		@param h		��������
		@param c		����ͨ����
		*/
		Matrix(int w, int h, int c);
		/**
		����size*1�ľ���
		@param size_	����ߴ�
		*/
		Matrix(Size size_);
		/**
		����size�ľ���
		@param size_	����ߴ�
		*/
		Matrix(Size3 size_);
		/**
		��������
		@param src		��������
		*/
		Matrix(const lzh::Matrix *src);
		/**
		��������
		@param src		��������
		*/
		Matrix(const lzh::Matrix &src);
		/**
		������a��b�ϲ�(COLΪ���кϲ�|ROWΪ���кϲ�)
		@param a		�������1
		@param b		�������2
		@param merge	�ϲ���ʽ
		*/
		Matrix(const lzh::Matrix &a, const lzh::Matrix &b, X_Y_Z merge);
		/**
		���캯��
		���m
		@param m ����
		*/
		Matrix(const MatCommaInitializer_ &m);
		/**
		����n*n*1�ľ���,Ԫ��Ϊmatrix
		@param matrix	����Ԫ��
		@param n		�����С
		*/
		Matrix(int *matrix, int n);
		/**
		����n*n*1�ľ���,Ԫ��Ϊmatrix
		@param matrix	����Ԫ��
		@param n		�����С
		*/
		Matrix(mat_t *matrix, int n);
		/**
		����h_*w_*1�ľ���,Ԫ��Ϊmatrix
		@param matrix	����Ԫ��
		@param w		��������
		@param h		��������
		*/
		Matrix(int *matrix, int w, int h, int c = 1);
		/**
		����h_*w_*1�ľ���,Ԫ��Ϊmatrix
		@param matrix	����Ԫ��
		@param w		��������
		@param h		��������
		*/
		Matrix(mat_t *matrix, int w, int h, int c = 1);
		/**
		����1*w*1������,Ԫ��Ϊdata
		@param w		����
		@param data		����
		*/
		Matrix(int w, mat_t *data);
		/**
		����h*w*1�ľ���,Ԫ��Ϊdata
		@param w		��������
		@param h		��������
		@param data		����Ԫ��
		*/
		Matrix(int w, int h, mat_t *data);
		/**
		����h*w*c�ľ���,Ԫ��Ϊdata
		@param w		��������
		@param h		��������
		@param c		����ͨ����
		@param data		����Ԫ��
		*/
		Matrix(int w, int h, int c, mat_t *data);
		/**
		����h*w*c�ľ���,Ԫ��Ϊdata
		@param w		��������
		@param h		��������
		@param c		����ͨ����
		@param step		����
		@param data		����Ԫ��
		*/
		Matrix(int w, int h, int c, int step, mat_t *data);
		/**
		����ת����
		@param vec		����
		@param dirc		������
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
		����w������
		@param w		����
		*/
		void create(int w);
		/**
		����w*h�ľ���
		@param w		����
		@param h		����
		*/
		void create(int w, int h);
		/**
		����w*h*c������
		@param w		����
		@param h		����
		@param c		ͨ����
		*/
		void create(int w, int h, int c);
		/**
		����size�ľ���
		@param size		�����С
		*/
		void create(Size size);
		/**
		����size������
		@param size		�����С
		*/
		void create(Size3 size);
		mat_t* data()const;
		mat_t* begin();
		const mat_t* begin()const;
		mat_t* end();
		const mat_t* end()const;
		/**
		@brief �ڴ泤��
		*/
		uint memsize()const;
		/**
		@brief ���ά��
		*/
		void DimCheck()const;
		/**
		@brief ���ؾ���ߴ�(h_,w_,c_)
		*/
		Size3 size3()const;
		/**
		@brief ���ؾ���ƫ��
		*/
		int total()const;
		/**
		@brief ����ά��
		*/
		int dims()const;
		/**
		@brief ��������
		*/
		int rows()const;
		/**
		@brief ��������
		*/
		int cols()const;
		/**
		@brief ����ͨ����
		*/
		int channels()const;
		/**
		@brief ��������
		*/
		int rank()const;
		/**
		@brief ���ؾ����С(h_*w_*c_)
		*/
		uint size()const;
		/**
		@brief ���ؾ����СSize(h_,w_)
		*/
		Size mSize()const;
		/**
		�������
		@param file		�����ļ���
		@param binary	ѡ���ı����Ƕ�����
		binary = false	ѡ���ı�
		binary = true	ѡ�������
		*/
		void save(std::string file, bool binary = true)const;
		/**
		�����Ʊ������
		@param file		�����ļ�ָ��
		*/
		void save(FILE *file)const;
		/**
		��ȡ����
		@param file		��ȡ�ļ���
		ֻ֧�ֶ����ƶ�ȡ
		*/
		void load(std::string file);
		/**
		��ȡ����
		@param file		��ȡ�ļ�ָ��
		ֻ֧�ֶ����ƶ�ȡ
		*/
		void load(FILE *file);
		/**
		@brief ���ؾ����С(h_*w_*c_)
		*/
		int len()const;
		/**
		@brief ���ؾ���״̬
		0Ϊ����
		-1Ϊ�վ���
		-2Ϊ�Ƿ���
		*/
		int enable()const;
		/**
		@brief ���ؾ����Ƿ�Ϊ��
		*/
		bool empty()const;
		/**
		@brief ���ؾ����Ƿ�Ϊ����
		*/
		bool Square()const;
		/**
		�������ݵ�mat
		@param mat		����
		*/
		void copyTo(Matrix& mat)const;
		/**
		�������ݵ�mat
		@param mat		����
		*/
		void copyTo(Matrix&& mat)const;
		/**
		�ھ�������߻����ұ����һ��1
		@param dire		ѡ����ӷ�ʽ
		*/
		void setAddOnes(Dire dire = RIGHT);
		/**
		�ھ�������߻����ұ����һ��0
		@param dire		ѡ����ӷ�ʽ
		*/
		void setAddZeros(Dire dire = RIGHT);
		/**
		�ͷ��ڴ�
		*/
		void release();
		/**
		@brief ���������ؾ���Ԫ��
		@param w		����
		*/
		mat_t& at(int w)const;
		/**
		@brief ���������ؾ���Ԫ��
		@param w		������
		@param h		������
		*/
		mat_t& at(int w, int h)const;
		/**
		@brief ���������ؾ���Ԫ��
		@param w		������
		@param h		������
		@param c		ͨ������
		*/
		mat_t& at(int w, int h, int c)const;
		/**
		@brief ������ת��Ϊ��Ӧ����������
		@param index	����
		*/
		int toX(int index)const;
		/**
		@brief ������ת��Ϊ��Ӧ����������
		@param index	����
		*/
		int toY(int index)const;
		/**
		@brief ������ת��Ϊ��Ӧ����ͨ������
		@param index	����
		*/
		int toZ(int index)const;

		/**
		@brief �����һ��Ԫ��
		*/
		mat_t frist()const;
		/**
		@brief ���ؾ�����value��ȵĵ�һ��Ԫ������
		@param value	Ԫ��
		*/
		int find(mat_t value)const;
		/**
		@brief ���ؾ���Ԫ�����ֵ������
		*/
		int maxAt()const;
		/**
		@brief ���ؾ���Ԫ����Сֵ������
		*/
		int minAt()const;
		/**
		@brief ���ؾ����Ƿ����value
		@param value	Ԫ��
		*/
		bool contains(mat_t value)const;
		/**
		@brief ���ؾ�����value��ȵĵ�һ��Ԫ��
		@param value	Ԫ��
		*/
		mat_t& findAt(mat_t value)const;
		mat_t max(bool is_abs = false)const;
		mat_t min(bool is_abs = false)const;
		/**
		@brief ���ؾ���Ԫ�����ֵ
		*/
		mat_t& findmax()const;
		/**
		@brief ���ؾ���Ԫ����Сֵ
		*/
		mat_t& findmin()const;
		/**
		@brief �������������򿽱�Ԫ�ص�src������
		@param src			����������
		@param Row_Start	��ȡ�г�ʼ����ֵ
		@param Row_End		��ȡ�н�������ֵ
		@param Col_Start	��ȡ�г�ʼ����ֵ
		@param Col_End		��ȡ�н�������ֵ
		*/
		void copy(Matrix &src, int Row_Start, int Row_End, int Col_Start, int Col_End)const;
		/**
		@brief �����󿽱���src
		@param src ����������
		*/
		void swap(Matrix &src)const;
		/**
		@brief mChannel ��src���ǵ���c_ͨ��
		@param src		����
		@param c_	ͨ����
		*/
		void mChannel(const lzh::Matrix &src, int c);
		/**
		@brief mChannel ��src���ǵ���c_ͨ��
		@param src		����
		@param c_	ͨ����
		*/
		void mChannel(const lzh::Matrix &src, int w, int h);
		/**
		@brief ���þ���ά��
		������ı���󳤶�
		@param size		�����С
		*/
		void reshape(Size3 size = Size3(0, 0, 0));
		/**
		@brief ���þ���ά��
		������ı���󳤶�
		@param h_		��������
		@param w_		��������
		@param c_	����ͨ��
		*/
		void reshape(int w, int h = 0, int c = 0);
		/**
		@brief ���þ����С
		�������ԭ��С������h_*w_*1��Ԫ��ȫ������Ϊ0
		@param h_	��������
		@param w_	��������
		*/
		bool setSize(int w, int h, int c);
		/**
		@brief ��������src
		@param src	��������
		*/
		void setvalue(const lzh::Matrix &src);
		/**
		@brief �޸ľ����Ӧ����Ԫ��
		@param number	Ԫ��
		@param index	����
		*/
		void setNum(mat_t number, int index);
		/**
		@brief �޸ľ����Ӧ����Ԫ��
		@param number	Ԫ��
		@param index_y	������
		@param index_x	������
		*/
		void setNum(mat_t number, int index_y, int index_x);
		/**
		@brief ���þ���
		@param mat	����Ԫ��
		@param h_	��������
		@param w_	��������
		*/
		void setMat(mat_t *mat, int w, int h);
		/**
		@brief ���������
		*/
		void setInv();
		/**
		@brief ���þ����num����
		@param num	����
		*/
		void setPow(mat_t num);
		/**
		@brief ����ȡ��
		*/
		void setOpp();
		/**
		@brief ���õ�λ����
		*/
		void setIden();
		/**
		@brief ���ð������
		*/
		void setAdj();
		/**
		@brief ����ת�þ���
		*/
		void setTran();

		/**
		@brief �������������
		*/
		void show()const;
		/**
		@brief �������
		*/
		std::ostream & show(std::ostream & out)const;

		/**
		@brief ����h�о���
		@param ����
		*/
		Matrix Row(int h);
		/**
		@brief ����h�о���
		@param ����
		*/
		const Matrix Row(int h)const;
		/**
		@brief ����w�о���
		@param ����
		*/
		Matrix Col(int w);
		/**
		@brief ����w�о���
		@param ����
		*/
		const Matrix Col(int w)const;
		/**
		@brief ����cͨ������
		@param ͨ������
		*/
		Matrix Channel(int c);
		/**
		@brief ����cͨ������
		@param ͨ������
		*/
		const Matrix Channel(int c)const;
		/**
		@brief �ھ�������߻����ұ����һ��1
		@param dire		ѡ����ӷ�ʽ
		*/
		const Matrix addones(Dire dire = RIGHT)const;
		/**
		@brief �ھ�������߻����ұ����һ��0
		@param dire		ѡ����ӷ�ʽ
		*/
		const Matrix addzeros(Dire dire = RIGHT)const;
		/**
		@brief ����start��end����
		@param start	��ʼ����
		@param end		��������
		*/
		const Matrix Range(int start, int end);
		/**
		@brief ����h_start��h_end�о���
		@param h_start	�п�ʼ����
		@param h_end	�н�������
		*/
		const Matrix rowRange(int h_start, int h_end);
		/**
		@brief ����w_start��w_end�о���
		@param w_start	�п�ʼ����
		@param w_end	�н�������
		*/
		const Matrix colRange(int w_start, int w_end);
		/**
		@brief ����c_start��c_endͨ������
		@param c_start	ͨ����ʼ����
		@param c_end	ͨ����������
		*/
		const Matrix channelRange(int c_start, int c_end);
		/**
		@brief �����������
		*/
		const lzh::Matrix clone()const;
		/**
		@brief ����ȡ������
		*/
		const lzh::Matrix opp()const;
		/**
		@brief ���ؾ���ֵ����
		*/
		const lzh::Matrix abs()const;
		/**
		@brief ���ذ�num���ݾ���
		@param num ����
		*/
		const lzh::Matrix mPow(int num)const;
		/**
		@brief ���ذ�num���ݾ���
		@param num ����
		*/
		const lzh::Matrix pow(mat_t num)const;
		/**
		@brief ���ذ�Ԫ��ȡָ������
		*/
		const lzh::Matrix exp()const;
		/**
		@brief ���ذ�Ԫ��ȡ��������
		*/
		const lzh::Matrix log()const;
		/**
		@brief ���ذ�Ԫ��ȡ��������
		*/
		const lzh::Matrix sqrt()const;
		/**
		@brief ���ذ������
		*/
		const lzh::Matrix adj()const;
		/**
		@brief ����ת�þ���
		*/
		const lzh::Matrix t()const;
		/**
		@brief ���������
		*/
		const lzh::Matrix inv()const;
		/**
		@brief ���������
		*/
		const lzh::Matrix diag(int k = 0)const;
		/**
		@brief �����������
		�������������
		*/
		const lzh::Matrix reverse()const;
		const lzh::Matrix EigenvectorsMax(mat_t offset = 1e-8)const;
		/**
		@brief sigmoid����
		*/
		const lzh::Matrix sigmoid()const;
		/**
		@brief tanh����
		*/
		const lzh::Matrix tanh()const;
		/**
		@brief relu����
		*/
		const lzh::Matrix relu()const;
		/**
		@brief elu����
		*/
		const lzh::Matrix elu()const;
		/**
		@brief selu����
		*/
		const lzh::Matrix selu()const;
		/**
		@brief leaky_relu����
		*/
		const lzh::Matrix leaky_relu()const;
		/**
		@brief softmax����
		*/
		const lzh::Matrix softmax()const;
		/**
		@brief ��������ʽ
		*/
		mat_t Det();
		/**
		@brief ����num����
		@param num ������
		*/
		mat_t Norm(int num = 1)const;
		/**
		@brief ���ض�Ӧ����������ʽ
		@param x ������
		@param y ������
		*/
		mat_t Cof(int x, int y);
		mat_t EigenvalueMax(mat_t offset = 1e-8)const;
		/**
		@brief ���������ȡ�ľ���Ԫ��
		*/
		mat_t RandSample();
		/**
		@brief ���ؾ���Ԫ�غ�
		@param num	���ô���
		@param _abs �Ƿ�ȡ����ֵ
		*/
		mat_t sum(int num = 1, bool _abs = false)const;
		/**
		@brief ����ƽ��ֵ
		*/
		mat_t mean()const;
		/**
		@brief ���ر�׼��
		*/
		mat_t S()const;
		/**
		@brief ���ر�׼��
		*/
		mat_t D()const;
		/**
		@brief ���������+
		��ӦԪ�����
		*/
		const lzh::Matrix operator + (const mat_t val)const;
		/**
		@brief ���������+
		��ӦԪ�����
		*/
		const lzh::Matrix operator + (const lzh::Matrix &a)const;
		/**
		@brief ���������+=
		��Ԫ�����
		*/
		void operator += (const mat_t val);
		/**
		@brief ���������+=
		��Ԫ�����
		*/
		void operator += (const lzh::Matrix &a);
		/**
		@brief ��Ԫ���������+
		��Ԫ�����
		*/
		friend const lzh::Matrix operator + (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief ���������-
		��Ԫ��ȡ�෴��
		*/
		const lzh::Matrix operator - (void)const;
		/**
		@brief ���������-
		��Ԫ�����
		*/
		const lzh::Matrix operator - (const mat_t val)const;
		/**
		@brief ���������-
		��ӦԪ�����
		*/
		const lzh::Matrix operator - (const lzh::Matrix &a)const;
		/**
		@brief ���������-=
		��Ԫ�����
		*/
		void operator -= (const mat_t val);
		/**
		@brief ���������-=
		��ӦԪ�����
		*/
		void operator -= (const lzh::Matrix &a);
		/**
		@brief ��Ԫ���������-
		��Ԫ�����
		*/
		friend const lzh::Matrix operator - (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief ���������*
		��Ԫ�����
		*/
		const lzh::Matrix operator * (const mat_t val)const;
		/**
		@brief ���������*
		��ӦԪ�����
		*/
		const lzh::Matrix operator * (const lzh::Matrix &a)const;
		/**
		@brief ���������*=
		��Ԫ�����
		*/
		void operator *= (const mat_t val);
		/**
		@brief ���������*=
		��ӦԪ�����
		*/
		void operator *= (const lzh::Matrix &a);
		/**
		@brief ��Ԫ���������*
		��Ԫ�����
		*/
		friend const lzh::Matrix operator * (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief ���������/
		��Ԫ�����
		*/
		const lzh::Matrix operator / (const mat_t val)const;
		/**
		@brief ���������/
		����˷�
		*/
		const lzh::Matrix operator / (const lzh::Matrix &a)const;
		/**
		@brief ���������/=
		��Ԫ�����
		*/
		void operator /= (const mat_t val);
		/**
		@brief ���������/=
		��ӦԪ�����
		*/
		void operator /= (const lzh::Matrix &a);
		/**
		@brief ��Ԫ���������/
		��Ԫ�����
		*/
		friend const lzh::Matrix operator / (const mat_t value, const lzh::Matrix &mat);
		/**
		@brief ���������=
		���
		*/
		const Matrix & operator = (const lzh::Matrix &temp);
		/**
		@brief ���������==
		�жϾ����Ƿ����
		*/
		bool operator == (const lzh::Matrix &a)const;
		/**
		@brief ���������!=
		�жϾ����Ƿ����
		*/
		bool operator != (const lzh::Matrix &a)const;
		/**
		@brief ���ض�Ӧ����Ԫ��
		@param w ����
		*/
		mat_t& operator () (const int w)const;
		/**
		@brief ���ض�Ӧ����Ԫ��
		@param h ������
		@param w ������
		*/
		mat_t& operator () (const int h, const int w)const;
		/**
		@brief ���ض�Ӧ����Ԫ��
		@param h		������
		@param w		������
		@param c		ͨ������
		*/
		mat_t& operator () (const int h, const int w, const int c)const;
		/**
		@brief ���ض�Ӧ����Ԫ��
		@param p ����
		*/
		mat_t& operator () (Point p)const;
		/**
		@brief ���ض�Ӧ����Ԫ��
		@param p ����
		*/
		mat_t& operator () (Point3i p)const;
		/**
		@brief ���ؾ����Ӧ�������л���
		@param index	����
		@param rc		������ʽ
		*/
		const Matrix operator () (const int index, X_Y_Z rc)const;
		/**
		@brief ���ؾ����Ӧ�������л���
		@param index	����
		@param rc		������ʽ
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
		���ô�ӡ���
		@param w ���
		*/
		static void setPrintW(lzh::uint w);
		/**
		���ô�ӡ��Ч����
		@param n λ��
		*/
		static void setPrintSignificantDigits(lzh::uint n);
		/**
		���ô�ӡ����
		@param t ��ӡ����
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
	@brief Mat_ ������
	�̳�Mat�࣬����ʵ��
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
		@brief ����w������
		@param w		����
		*/
		Mat_(int w) : Mat(w) {}
		/**
		@brief ����h*w�ľ���
		@param w		����
		@param h		����
		*/
		Mat_(int w, int h) : Mat(w, h) {}
		/**
		@brief ����h*w*c������
		@param w		����
		@param h		����
		@param depth	ͨ����
		*/
		Mat_(int w, int h, int c) : Mat(w, h, c) {}
		/**
		@brief ����size�ľ���
		@param size		�ߴ�
		*/
		Mat_(Size size_) : Mat(size_) {}
		/**
		@brief ����size�ľ���
		@param size		�ߴ�
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
	@brief MatCommaInitializer_ ������
	��Ϊ������������ʵ��
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
	�������������
	*/
	extern void Srandom();
	/**
	@brief ����Ԫ��Ϊ0��w����
	@param w		��������
	*/
	extern const Mat zeros(int w);
	/**
	@brief ����Ԫ��Ϊ0��h*w����
	@param w		��������
	@param h		��������
	*/
	extern const Mat zeros(int w, int h);
	/**
	@brief ����Ԫ��Ϊ0��h*w*c����
	@param w		��������
	@param h		��������
	@param c		����ͨ����
	*/
	extern const Mat zeros(int w, int h, int c);
	/**
	@brief ����Ԫ��Ϊ0��size����
	@param size �����С
	*/
	extern const Mat zeros(Size size);
	/**
	@brief ����Ԫ��Ϊ0��size����
	@param size �����С
	*/
	extern const Mat zeros(Size3 size);
	/**
	@brief ����Ԫ��Ϊv��w����
	@param v		���Ԫ��
	@param w		��������
	*/
	extern const Mat value(mat_t v, int w);
	/**
	@brief ����Ԫ��Ϊv��h*w����
	@param v		���Ԫ��
	@param w		��������
	@param h		��������
	*/
	extern const Mat value(mat_t v, int w, int h);
	/**
	@brief ����Ԫ��Ϊv��h_*w_*c_����
	@param v		���Ԫ��
	@param w		��������
	@param h		��������
	@param c		����ͨ����
	*/
	extern const Mat value(mat_t v, int w, int h, int c);
	/**
	@brief ����Ԫ��Ϊ1��w����
	@param w		��������
	*/
	extern const Mat ones(int w);
	/**
	@brief ����Ԫ��Ϊ1��h*w����
	@param w	��������
	@param h	��������
	*/
	extern const Mat ones(int w, int h);
	/**
	@brief ����Ԫ��Ϊ1��h*w*c����
	@param w		��������
	@param h		��������
	@param c		����ͨ����
	*/
	extern const Mat ones(int w, int h, int c);
	/**
	@brief ����Ԫ��Ϊ1��size����
	@param size �����С
	*/
	extern const Mat ones(Size size);
	/**
	@brief ����Ԫ��Ϊ1��size����
	@param size �����С
	*/
	extern const Mat ones(Size3 size);
	/**
	@brief ���ش�low��top�ȷֳɵ�1*len�ľ���
	@param low �½�
	@param top �Ͻ�
	@param gap ���
	*/
	extern const Mat range(int low, int top, mat_t gap = _T(1));
	/**
	@brief ���ش�low��top�ȷֳɵ�1*len�ľ���
	@param low �½�
	@param top �Ͻ�
	@param gap ���
	*/
	extern const Mat range(mat_t low, mat_t top, mat_t gap = _T(1));
	/**
	@brief ���ش�low��top�ȷֳɵ�1*len�ľ���
	@param low �½�
	@param top �Ͻ�
	@param len �ȷָ���
	*/
	extern const Mat linspace(int low, int top, int len);
	/**
	@brief ���ش�low��top�ȷֳɵ�1*len�ľ���
	@param low �½�
	@param top �Ͻ�
	@param len �ȷָ���
	*/
	extern const Mat linspace(mat_t low, mat_t top, int len);
	/**
	@brief �����������
	���ϸ�˹�ֲ�
	��СΪsize
	@param size		�����С
	@param n1		��ǰͨ��
	@param n2		��ǰͨ��
	*/
	extern const Mat Xavier(Size3 size, int n1, int n2);
	/**
	@brief �����������
	���ϸ�˹�ֲ�
	��СΪh_*w_*c_
	@param h_		��������
	@param w_		��������
	@param c_	����ͨ����
	@param n1		��ǰͨ��
	@param n2		��ǰͨ��
	*/
	extern const Mat Xavier(int w, int h, int c, int n1, int n2);
	/**
	@brief �����������
	���ϸ�˹�ֲ�
	��СΪsize, Ԫ�ط�Χ[0, 1]
	@param size		�����С
	*/
	extern const Mat Random(Size3 size);
	/**
	@brief �����������
	���ϸ�˹�ֲ�
	��СΪh_*w_*c_, Ԫ�ط�Χ[0, 1]
	@param h_		��������
	@param w_		��������
	@param c_	����ͨ����
	*/
	extern const Mat Random(int w, int h = 1, int c = 1);
	/**
	@brief ���ؾ�������ֵ
	@param src		����
	@param isAbs	�Ƿ�ȡ����ֵ
	*/
	extern mat_t Max(const Mat &src, bool isAbs = false);
	/**
	@brief ���ؾ������Сֵ
	@param src		����
	@param isAbs	�Ƿ�ȡ����ֵ
	*/
	extern mat_t Min(const Mat &src, bool isAbs = false);
	/**
	@brief ���ؾ��������ʽ
	@param src	����
	*/
	extern mat_t det(const Mat &src);
	/**
	@brief ���ؾ���ļ�
	@param src	����
	*/
	extern mat_t trace(const Mat &src);
	/**
	@brief ���ؾ����Ӧ����������ʽ
	@param src	����
	@param x	������
	@param y	������
	*/
	extern mat_t cof(const Mat &src, int x, int y);
	/**
	@brief ���������
	@param min		��Сֵ
	@param max		���ֵ
	@param isdouble �Ƿ����������
	*/
	extern mat_t getRandData(int min, int max, bool isdouble = false);
	/**
	@brief ���ؾ�����
	@param src ����
	@param num ������
	*/
	extern mat_t mNorm(const Mat& src, int num = 1);
	/**
	@brief ���ؾ���ľ���
	@param a	����
	@param b	����
	@param num	������
	*/
	extern mat_t mDistance(const Mat& a, const Mat& b, int num = 2);
	/**
	@brief ��������ľ���Ԫ��
	@param src ����
	*/
	extern mat_t mRandSample(const Mat &src);
	/**
	@brief �������ɵ�n*n*1��λ����
	@param n �����С
	*/
	extern const Mat eye(int n);
	/**
	@brief ���ؾ���ĵ�c_��ͨ��
	@param src		����
	@param c_	ͨ������
	*/
	extern const Mat mSplit(const Mat &src, int c);
	/**
	@brief ���ذ������ͨ��������
	@param src �������
	@param dst �������ͨ�����ľ�������
	*/
	extern void mSplit(const Mat &src, Mat *dst);
	/**
	@brief ���ذ�ͨ���ϲ��ľ���
	@param src		��������
	@param channels ͨ����
	*/
	extern const Mat mMerge(const Mat *src, int channels);
	/**
	@brief ���ذ����������и�ľ���
	@param src			����
	@param Row_Start	��ȡ�г�ʼ����ֵ
	@param Row_End		��ȡ�н�������ֵ
	@param Col_Start	��ȡ�г�ʼ����ֵ
	@param Col_End		��ȡ�н�������ֵ
	*/
	extern const Mat Block(const Mat &src, int Row_Start, int Row_End, int Col_Start, int Col_End, int Chennel_Start = 0, int Chennel_End = 0);
	/**
	@brief �����������Ԫ��n*n*1����
	@param n		�����С
	@param low		�½�
	@param top		�Ͻ�
	@param isdouble �Ƿ����ɸ�����
	*/
	extern const Mat mRand(int low, int top, int n, bool isdouble = false);
	/**
	@brief �����������Ԫ��size.x*size.y*size.z����
	@param low		�½�
	@param top		�Ͻ�
	@param size		�����С
	@param isdouble �Ƿ����ɸ�����
	*/
	extern const Mat mRand(int low, int top, Size3 size, bool isdouble = false);
	/**
	@brief �����������Ԫ��h_*w_*c_����
	@param h_		��������
	@param w_		��������
	@param low		�½�
	@param top		�Ͻ�
	@param isdouble �Ƿ����ɸ�����
	*/
	extern const Mat mRand(int low, int top, int w, int h, int c = 1, bool isdouble = false);
	/**
	@brief ���ش�СΪw����
	@param w		��������
	*/
	extern const Mat mcreate(int w);
	/**
	@brief ���ش�СΪh*w����
	@param w		��������
	@param h		��������
	*/
	extern const Mat mcreate(int w, int h);
	/**
	@brief ���ش�СΪh*w*c����
	@param w		��������
	@param h		��������
	@param c		����ͨ����
	*/
	extern const Mat mcreate(int w, int h, int c);
	/**
	@brief ���ش�СΪsize����
	@param size �����С
	*/
	extern const Mat mcreate(Size size);
	/**
	@brief ���ش�СΪsize����
	@param size �����С
	*/
	extern const Mat mcreate(Size3 size);
	/**
	@brief ����������󣬾�����Ϊһά����
	@param src	����
	*/
	extern const Mat reverse(const Mat &src);
	/**
	@brief ��������h*w*c�ľ��������ȡ����srcԪ�����
	@param src	����
	@param h	��������
	@param w	��������
	@param c	����ͨ����
	*/
	extern const Mat mRandSample(const Mat &src, int w, int h, int c = 1);
	/**
	@brief ���������ȡnum�ξ���src���л�����ɵľ���
	@param src	����
	@param rc	��ȡ��ʽ
	@param num	��ȡ����
	*/
	extern const Mat mRandSample(const Mat &src, X_Y_Z rc, int num = 1);
	/**
	@brief ���ؾ���İ������
	@param src ����
	*/
	extern const Mat adj(const Mat &src);
	/**
	@brief ���ؾ���������
	@param src ����
	*/
	extern const Mat inv(const Mat &src);
	/**
	@brief ���ضԽ��߾���
	@param src	����
	@param k	��k���Խ���
	*/
	extern const Mat diag(const Mat &src, int k = 0);
	/**
	@brief ���ؾ����α�����
	@param src	����
	@param dire α�����ļ��㷽ʽ
	*/
	extern const Mat pinv(const Mat &src, Dire dire = LEFT);
	/**
	@brief ���ؾ����ת�þ���
	@param src ����
	*/
	extern const Mat tran(const Mat &src);
	/**
	@brief ���ؾ���ľ���ֵ����
	@param src ����
	*/
	extern const Mat mAbs(const Mat &src);
	/**
	@brief ����angle��2*2����ת����
	@param angle �Ƕ�
	*/
	extern const Mat Rotate(mat_t angle);
	/**
	@brief ���ؾ���num����
	@param src ����
	@param num ����
	*/
	extern const Mat POW(const Mat &src, int num);
	/**
	@brief ���ؾ���ȡ��
	@param src ����
	*/
	extern const Mat mOpp(const Mat &src);
	/**
	@brief ���ؾ����л���֮��
	@param src	����
	@param rc	��͵ķ���
	*/
	extern const Mat mSum(const Mat &src, X_Y_Z rc);
	/**
	@brief ���ؾ���Ԫ��ȡָ��
	@param src ����
	*/
	extern const Mat mExp(const Mat &src);
	/**
	@brief ���ؾ���Ԫ��ȡ����
	@param src ����
	*/
	extern const Mat mLog(const Mat &src);
	/**
	@brief ���ؾ���Ԫ��ȡ����
	@param src ����
	*/
	extern const Mat mSqrt(const Mat &src);
	/**
	@brief ���ؾ���Ԫ��ȡnum����
	@param src ����
	@param num ����
	*/
	extern const Mat mPow(const Mat &src, mat_t num);
	/**
	@brief ���ؾ���val/src��Ԫ�س�
	@param src ����
	@param val ����
	*/
	extern const Mat Divi(const Mat &src, mat_t val, Dire dire = RIGHT);
	/**
	@brief ���ؾ������
	@param a	��������
	@param b	������
	@param dire ������ʽ
	*/
	extern const Mat Divi(const Mat &a, const Mat &b, Dire dire = RIGHT);
	/**
	@brief ���ع������
	@param a ����
	@param b ����
	*/
	extern const Mat Mult(const Mat &a, const Mat &b);
	/**
	@brief ���ؾ���˷�
	@param a ����
	@param b ����
	*/
	extern const Mat Dot(const Mat &a, const Mat &b);
	/**
	@brief ���ؾ���Ԫ��ȡa��b֮������ֵ
	@param a �Ƚ�ֵ
	@param b �ȽϾ���
	*/
	extern const Mat mMax(mat_t a, const Mat &b);
	/**
	@brief ���ؾ���Ԫ��ȡa��b֮������ֵ
	@param a �ȽϾ���
	@param b �ȽϾ���
	*/
	extern const Mat mMax(const Mat &a, const Mat &b);
	/**
	@brief ���ؾ���Ԫ��ȡa��b֮�����Сֵ
	@param a �Ƚ�ֵ
	@param b �ȽϾ���
	*/
	extern const Mat mMin(mat_t a, const Mat &b);
	/**
	@brief ���ؾ���Ԫ��ȡa��b֮�����Сֵ
	@param a �ȽϾ���
	@param b �ȽϾ���
	*/
	extern const Mat mMin(const Mat &a, const Mat &b);
	/**
	@brief mCalSize �������������ŵı߽�
	���ؾ����С
	@param src		���������
	@param kern		�����
	@param anchor	���ض�Ӧ���������
	anchorĬ��ΪPoint(-1,-1), ���ض�Ӧ���������
	@param strides	��������
	@param top		�������伸��
	@param bottom	�������伸��
	@param left		�������伸��
	@param right	�������伸��
	*/
	extern Size3 mCalSize(Size3 src, Size3 kern, Point anchor, Size strides, int &top, int &bottom, int &left, int &right);
	/**
	@brief mCalSize �������������ŵı߽�
	���ؾ����С
	@param src		���������ߴ�
	@param kern		����˳ߴ�
	@param anchor	���ض�Ӧ���������
	*/
	extern Size3 mCalSize(Size3 src, Size3 kern, Point &anchor, Size strides);
	/**
	@brief mCalSize �������������ŵı߽�
	���ؾ����С
	@param src		���������ߴ�
	@param kern		����˳ߴ�
	@param anchor	���ض�Ӧ���������
	*/
	extern Size3 mCalSize(Size3 src, Size kern, Point &anchor, Size strides);
	/**
	@brief ���ذ�boundary�ֽ����ľ���
	���ؾ����С������������С
	@param src				�������
	@param boundary			�ֽ�ֵ
	@param lower			С��boundary��lower���
	@param upper			����boundary��upper���
	@param boundary2upper	������Ԫ�ص���boundaryʱ
	Ϊ1��upper,				Ϊ-1��lower, Ϊ0������
	*/
	extern const Mat Threshold(const Mat &src, mat_t boundary, mat_t lower, mat_t upper, int boundary2upper = 1);
	/**
	@brief ���ر߽�����ľ���
	@param src			�������
	@param top			�������伸��
	@param bottom		�������伸��
	@param left			�������伸��
	@param right		�������伸��
	@param borderType	�߽��������ƵĲ�ֵ����
	@param value		������ֵ����ֵ
	**/
	extern const Mat copyMakeBorder(const Mat &src, int top, int bottom, int left, int right, BorderTypes borderType = BORDER_CONSTANT, mat_t value = 0.0);
	/**
	@brief ���ؾ���2ά������(ֻ֧�ֶ�ά����)
	���ؾ����СΪ(input.h_/strides_x, input.w_/strides_y, 1)
	@param input			�������
	@param kern				�����
	@param anchor			����Ԫ�ض�Ӧ����˵�λ��
	�Ծ���˵����Ͻ�Ϊ(0,0)��, Ĭ��(-1,-1)Ϊ����
	@param strides			��������
	Size.heiΪx��,Size.widΪy��
	@param is_copy_border	�Ƿ�Ҫ��չ�߽�
	*/
	extern const Mat Filter2D(const Mat & input, const Mat & kern, Point anchor = Point(-1, -1), const Size & strides = Size(1, 1), bool is_copy_border = true);
	/**
	@brief ���ؾ���2ά������(֧����ά����)
	���ؾ����СΪ(input.h_/strides_x, input.w_/strides_y, 1)
	@param input			�������
	@param kern				�����
	@param anchor			����Ԫ�ض�Ӧ����˵�λ��
	�Ծ���˵����Ͻ�Ϊ(0,0)��, Ĭ��(-1,-1)Ϊ����
	@param strides			��������
	Size.heiΪx��,Size.widΪy��
	@param is_copy_border	�Ƿ�Ҫ��չ�߽�
	*/
	extern void Filter2D(const Mat & in, Mat & out, const Mat & kern, Point anchor = Point(-1, -1), const Size & strides = Size(1, 1), bool is_copy_border = true);
	/**
	@brief ������Ͻ��
	��С���˷�
	@param x �Ա���
	@param y �����
	*/
	extern const Mat LeastSquare(const Mat& x, const Mat &y);
	/**
	@brief ����a
	�������y=a(0)*x+a(1)
	y��x����Ϊ������
	��֤y��x��������ͬ
	@param x �Ա���
	@param y �����
	*/
	extern const Mat regress(const Mat& x, const Mat &y);
	/**
	@brief ����P
	����ʽ���y=P(1)*x^n + P(2)*x^(n-1) +...+ P(n)*x + P(n+1)
	y��x����Ϊ������
	��֤y��x��������ͬ
	@param x �Ա���
	@param y �����
	*/
	extern const Mat polyfit(const Mat& x, const Mat &y, uint n);
	/**
	@brief ������Ͻ��
	���Բ
	p����Ϊ������
	@param p �㼯
	*/
	extern const Mat circlefit(const Mat& p);
	/**
	@brief ������Ͻ��
	��������С���˷�
	y��x����Ϊ������
	��֤y��x��������ͬ
	@param x		�Ա���
	@param y		�����
	@param a0		��ʼ����
	@param f		����ָ�� f(a, x) = y
	@param step		���²���
	@param error	���(С������������)
	*/
	extern const Mat NonLinearLeastSqures(const Mat & x, const Mat & y, const Mat & a0, F f, mat_t step = _T(1e-2), mat_t error = _T(1e-6));
	/**
	@brief �����ݶ�
	@param y	�����
	@param x	�Ա���(ȱʡ��ֵΪ1)
	*/
	extern const Mat gradient(const Mat & y, const Mat & x = Mat());
	/**
	@brief ������ֵ�ݶ�
	@param f		����ָ�� f(x) = y
	@param x		�Ա���
	@param epsilon	��ֵ
	*/
	extern mat_t NumericGradient(NF f, mat_t x, mat_t epsilon = _T(1e-3));
	/**
	@brief ������ֵ�ݶ�
	@param f		����ָ�� f(x) = y
	@param x		�Ա���
	@param epsilon	��ֵ
	*/
	extern const Mat NumericGradient(Fun f, const Mat & x, mat_t epsilon = _T(1e-3));
	/**
	@brief ������ֵ�ݶ�
	@param f		����ָ�� f(a, x) = y
	@param a		�Ա���
	@param x		����
	@param epsilon	��ֵ
	*/
	extern const Mat NumericGradient(F f, const Mat & a, const Mat & x, mat_t epsilon = _T(1e-3));
	/**
	@brief ���ػ���ֵ
	ʹ��ŷ��������ֵ����
	@param f		����ָ�� f(x) = y
	@param low		��������
	@param high		��������
	@param epsilon	�������
	*/
	extern mat_t EulerInt(NF f, mat_t low, mat_t high, mat_t epsilon = _T(1e-3));
	/**
	@brief ���ػ���ֵ
	ʹ�����η�����ֵ����
	@param f		����ָ�� f(x) = y
	@param low		��������
	@param high		��������
	@param epsilon	�������
	*/
	extern mat_t TrapezoidInt(NF f, mat_t low, mat_t high, mat_t epsilon = _T(1e-3));
	/**
	@brief ���ػ���ֵ
	ʹ���Ľ�����-����������ֵ����
	@param f		����ָ�� f(x) = y
	@param low		��������
	@param high		��������
	@param epsilon	�������
	*/
	extern mat_t RungeKuttaInt(NF f, mat_t low, mat_t high, mat_t epsilon = _T(1e-3));
	extern const Mat LinearIntersection(const Mat & line_1, const Mat & line_2);
	/**
	@brief ����y
	����ʽ���y=P(1)*x^n + P(2)*x^(n-1) +...+ P(n)*x + P(n+1)
	x����Ϊ������
	@param a ����
	@param x �Ա���
	*/
	extern const Mat polynomial(const Mat& a, const Mat &x);
	/**
	@brief �����޸�ά�ȵľ���
	������ı���󳤶�
	@param src	����
	@param size �����С
	*/
	extern const Mat Reshape(const Mat& src, Size3 size);
	/**
	@brief ���ؾ���
	��ͨ�����
	@param src	����
	*/
	extern const Mat SumChannel(const Mat& src);
	/**
	@brief ������ת����
	@param src	����
	@param dice ��ת�Ƕ�
	*/
	extern const Mat rotate(const Mat& src, RotateAngle dice);
	/**
	���������ž���
	@param src		����
	@param dst		���
	@param xRatio	x���ű���
	@param yRatio	y���ű���
	@param mothed	���ŷ���
	*/
	extern void resize(const Mat & src, Mat & dst, mat_t xRatio, mat_t yRatio, ReductionMothed mothed);
	/**
	����С���ž���
	@param src		����
	@param dst		���
	@param newSize	�µľ����С
	@param mothed	���ŷ���
	*/
	extern void resize(const Mat& src, Mat& dst, Size newSize, ReductionMothed mothed);
	/**
	@brief ���������
	���ϸ�˹�ֲ��������
	@param mu		����ֵ
	@param sigma	��׼��
	*/
	extern mat_t generateGaussianNoise(mat_t mu = 0, mat_t sigma = 1);
	/**
	@brief ���ؾ�����
	K��ֵ����
	@param P			����
	@param k			�Ե㼯�ķ�����
	@param K			�������
	@param iteration	��������
	@param error		���(С������������)
	*/
	extern const Mat kmeans(const Mat &P, Mat &k, const uint K, const uint iteration, const mat_t error = _T(1e-7));
	/**
	@brief ��������
	����������������
	@param src		����
	@param dst		��������
	*/
	extern int RowSimplest(const Mat & src, Mat & dst);	
	/**
	@brief ������������
	����������������
	@param src		����
	*/
	extern const Mat RowSimplest(const Mat & src);
	/**
	@brief ��������
	����������������
	@param src		����
	@param dst		��������
	*/
	extern int ColSimplest(const Mat & src, Mat & dst);
	/**
	@brief ������������
	����������������
	@param src		����
	@param dst		��������
	*/
	extern const Mat ColSimplest(const Mat & src);
	/**
	@brief �������״̬
	����������Է�����
	������Է�����ʹ��addzeros�ֶ����0
	@param src		�������Է�����
	@param dst		������ϵ
	@param simplest	����о���
	@param mark		���״̬(���ɽ�Ϊ0,�ؽ�Ϊ1,�޽�Ϊnullptr)
	*/
	extern EQUATION SolveLinearEquation(const Mat & src, Mat & dst, Mat *simplest = nullptr, Mat *mark = nullptr);
	/**
	�����㷨��
	�ṩ��
	ð������	=> bubble;
	��������	=> insert;
	ѡ������	=> select;
	������		=> comb;
	�ؾ�����	=> gnome;
	������		=> heap;
	ϣ������	=> shell;
	��������	=> quick;
	�鲢����	=> merge;
	*/
	class sort
	{
	public:
		/**
		ð������
		@param m		��������
		@param order	˳��
		*/
		static void bubble(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		ð������
		@param m		��������
		@param order	˳��
		*/
		static const Mat bubble(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		��������
		@param m		��������
		@param order	˳��
		*/
		static void insert(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		��������
		@param m		��������
		@param order	˳��
		*/
		static const Mat insert(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		ѡ������
		@param m		��������
		@param order	˳��
		*/
		static void select(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		ѡ������
		@param m		��������
		@param order	˳��
		*/
		static const Mat select(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		������
		@param m		��������
		@param order	˳��
		*/
		static void comb(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		������
		@param m		��������
		@param order	˳��
		*/
		static const Mat comb(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		�ؾ�����
		@param m		��������
		@param order	˳��
		*/
		static void gnome(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		�ؾ�����
		@param m		��������
		@param order	˳��
		*/
		static const Mat gnome(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		������
		@param m		��������
		@param order	˳��
		*/
		static void heap(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		������
		@param m		��������
		@param order	˳��
		*/
		static const Mat heap(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		ϣ������
		@param m		��������
		@param order	˳��
		*/
		static void shell(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		ϣ������
		@param m		��������
		@param order	˳��
		*/
		static const Mat shell(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		��������
		@param m		��������
		@param order	˳��
		*/
		static void quick(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		��������
		@param m		��������
		@param order	˳��
		*/
		static const Mat quick(const Mat & m, ORDER order = MIN_TO_MAX);
		/**
		�鲢����
		@param m		��������
		@param order	˳��
		*/
		static void merge(mat_t *begin, mat_t *end, ORDER order = MIN_TO_MAX);
		/**
		@brief ����������
		�鲢����
		@param m		��������
		@param order	˳��
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
	B-����������
	*/
	class BSpline
	{
	public:
		/**
		UNIFORM			����B��������
		QUASI_UNIFORM	׼����B��������
		*/
		enum BSplineType
		{
			UNIFORM	= 0,
			QUASI_UNIFORM
		};

		explicit BSpline();
		BSpline(BSplineType type, int k = 1, const Mat &p = Mat());
		/**
		@brief ����B����������t�ĵ�
		t�Ķ�������[0,1]
		@param t	���������Ա���
		*/
		void setCtrlPoint(const Mat &p);
		/**
		ע��B�������߽ڵ�����
		@param node	�ڵ�����
		ȱʡʱ�����������Զ�����
		*/
		void NodeVector(const Mat &node = Mat());
		/**
		@brief ����B�������߿��Ƶ�
		*/
		const Mat CtrlPoint()const;
		/**
		@brief ����B�������߽ڵ�����
		*/
		const Mat Node()const;
		/**
		@brief ����B�������ߵ㼯
		@param T	���������Ա���
		TΪ��������[0,1]�ĵ�������
		*/
		const Mat BPiont(const Mat & T)const;
		/**
		@brief ����B�������ߵ㼯
		@param number	[0,1]�ȼ��number�ȷ�
		*/
		const Mat BPoint(int number)const;
		/**
		@brief ����B����������t�ĵ�
		t�Ķ�������[0,1]
		@param t	���������Ա���
		*/
		const Mat operator ()(mat_t t)const;
		/**
		@brief ����B��������
		�㼯���B��������
		@param P	�㼯
		@param n	���Ƶ�����
		@param k	k��B��������
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
	@brief softmax����
	Si = exp(y - max(y)) / sum(exp(y - max(y)))
	*/
	extern const Mat Softmax(const Mat &y);
	/**
	@brief L2����
	E = (y - y0)^2
	*/
	extern const Mat L2(const Mat & y, const Mat & y0);
	/**
	@brief quadratic����
	E = 1/2 * (y - y0)^2
	*/
	extern const Mat Quadratic(const Mat & y, const Mat & y0);
	/**
	@brief crossentropy����
	E = -(y * log(y0))
	*/
	extern const Mat CrossEntropy(const Mat & y, const Mat & y0);
	/**
	@brief softmax + crossentropy����
	E = -(y * log(softmax(y0)))
	*/
	extern const Mat SoftmaxCrossEntropy(const Mat & y, const Mat & y0);
	/**
	@brief sigmoid����
	y = 1/(1 + exp(-x))
	*/
	extern const Mat Sigmoid(const Mat & x);
	/**
	@brief tanh����
	y = (exp(x) - exp(-x)) / (exp(x) + exp(-x))
	*/
	extern const Mat Tanh(const Mat & x);
	/**
	@brief relu����
	y = {x if x > 0; 0 if x < 0}
	*/
	extern const Mat ReLU(const Mat & x);
	/**
	@brief elu����
	y = {x if x > 0; a*(exp(x) - 1) if x < 0}
	*/
	extern const Mat ELU(const Mat & x);
	/**
	@brief selu����
	y = scale*Elu(x)
	*/
	extern const Mat SELU(const Mat & x);
	/**
	@brief leaky relu����
	y = {x if x > 0; a*x if x < 0}
	*/
	extern const Mat LReLU(const Mat & x);

	class tools
	{
	public:
		/**
		��ͣ����
		*/
		static void pause();
		/**
		���������
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
		@brief ���ؽ�����ת����1*point.size()*1�ľ���
		@param point ����
		*/
		static const Mat Vec2Mat(std::vector<mat_t> &point);
		/**
		@brief ���ؽ�����ת����point.size()*point[0].size()*1�ľ���
		@param points ����
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
		@brief ���ؽ�����ת����һά����
		@param src ����
		*/
		static std::vector<mat_t> Mat2Vec(const Mat &src);
		/**
		@brief ���ؽ�����ת��������
		@param src ����
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
