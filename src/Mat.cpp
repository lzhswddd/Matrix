#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <iomanip>
//#define MAT_EXPORTS
#include "mat.hpp"

#ifdef MAT_DEBUG
#define CHECK_MATRIX(matrix) if((matrix) == nullptr) {std::cerr << lzh::errinfo[lzh::ERR_INFO_EMPTY] << std::endl;throw std::exception(lzh::errinfo[lzh::ERR_INFO_EMPTY]);}
#define CHECK_MATRIX_(mat) if((mat).empty()) {std::cerr << lzh::errinfo[lzh::ERR_INFO_EMPTY] << std::endl;throw std::exception(lzh::errinfo[lzh::ERR_INFO_EMPTY]);}
#define CHECK_PTR_ORDER(start, end) if((end) - (start) < 0) {std::cerr << lzh::errinfo[lzh::ERR_INFO_UNLESS] << std::endl;throw std::exception(lzh::errinfo[lzh::ERR_INFO_UNLESS]);}
#endif // MAT_DEBUG
#define CHECK_PTR(p) if((p) == nullptr) {std::cerr << lzh::errinfo[lzh::ERR_INFO_PTR] << std::endl;throw std::exception(lzh::errinfo[lzh::ERR_INFO_PTR]);}


lzh::uint lzh::Matrix::print_width = 10;
lzh::uint lzh::Matrix::print_precision = 4;
lzh::PrintType lzh::Matrix::print_type = lzh::SCIENTIFIC;

void lzh::THROW_INFO(lzh::MatErrorInfo info)
{
	std::cerr << lzh::errinfo[info] << std::endl; 
	throw std::exception(lzh::errinfo[(info)]);
}
/****************************************************************************
矩阵类
*****************************************************************************/
lzh::Matrix::Matrix()
{
	init();
	checkSquare();
}
lzh::Matrix::Matrix(int w)
{
	init();
	create(w);
}
lzh::Matrix::Matrix(int w, int h)
{
	init();
	create(w, h);
}
lzh::Matrix::Matrix(int w, int h, int depth)
{
	init();
	create(w, h, depth);
}
lzh::Matrix::Matrix(Size size_)
{
	init();
	create(size_.w, size_.h);
}
lzh::Matrix::Matrix(Size3 size_)
{
	init();
	create(size_.w, size_.h, size_.c);
}
lzh::Matrix::Matrix(lzh::mat_t *matrix, int n)
{
	init();
	*this = Matrix(matrix, n, n);
}
lzh::Matrix::Matrix(int *matrix, int n)
{
	init();
	*this = Matrix(matrix, n, n);
}
lzh::Matrix::Matrix(int *matrix, int w, int h, int c)
{
	init();
	if (matrix != nullptr) {
		lzh::tools::tools::check(w, h, c);
		setsize(w, h, c);
		this->matrix.create(size());
#ifdef MAT_DEBUG
		CHECK_PTR(this->matrix.p);
#endif //MAT_DEBUG
		for (int index = 0; index < len(); index++)
			(*this)(index) = (lzh::mat_t)matrix[index];
	}
	checkSquare();
}
lzh::Matrix::Matrix(lzh::mat_t *matrix, int w, int h, int c)
{
	init();
	if (matrix != nullptr) {
		lzh::tools::tools::check(w, h, c);
		setsize(w, h, c);
		this->matrix.create(size());
#ifdef MAT_DEBUG
		CHECK_PTR(this->matrix.p);
#endif //MAT_DEBUG
		memcpy(this->matrix, matrix, h_*w_*c_ * sizeof(lzh::mat_t));
	}
	checkSquare();
}
lzh::Matrix::Matrix(int w, mat_t * data)
{
	init();
	setsize(w);
	matrix = data;
}
lzh::Matrix::Matrix(int w, int h, mat_t * data)
{
	init();
	setsize(w, h);
	matrix = data;
}
lzh::Matrix::Matrix(int w, int h, int c, mat_t * data)
{
	init();
	setsize(w, h, c);
	matrix = data;
}
lzh::Matrix::Matrix(int w, int h, int c, int c_offset, mat_t * data)
{
	init();
	setsize(w, h, c);
	step = c_offset;
	matrix = data;
}
lzh::Matrix::Matrix(const lzh::Matrix &src)
{
	init();
	*this = src;
	/*init();
	setvalue(src);
	checkSquare();*/
}
lzh::Matrix::Matrix(const lzh::Matrix *src)
{
	init();
	*this = *src;
	/*init();
	if (src != nullptr)
		setvalue(*src);
	checkSquare();*/
}
lzh::Matrix::Matrix(const lzh::Matrix &a, const lzh::Matrix &b, X_Y_Z merge)
{
	init();
	if (merge == ROW) {
		if (a.w_ == b.w_) {
			create(a.w_, a.h_ + b.h_, a.c_);
			memcpy(matrix, a.matrix, a.memsize());
			memcpy(matrix + a.len(), b.matrix, b.memsize());
		}
	}
	else if (merge == COL) {
		if (a.h_ == b.h_) {
			create(a.w_ + b.w_, a.h_, a.c_);
			for (int i = 0; i < h_; i++) {
				memcpy(matrix + i * w_*c_,
					a.matrix + i * a.w_*c_,
					a.w_*c_ * sizeof(lzh::mat_t));
				memcpy(matrix + i * w_*c_ + a.w_*c_,
					b.matrix + i * b.w_*c_,
					b.w_*c_ * sizeof(lzh::mat_t));
			}
		}
	}
	checkSquare();
}
lzh::Matrix::Matrix(const MatCommaInitializer_ & m)
{
	init();
	*this = Matrix(m.matrix(), m.cols(), m.rows(), m.channels());
}
lzh::Matrix::~Matrix()
{
	release();
	setsize(0, 0, 0);
	/*if (matrix != nullptr) {
		fastFree(matrix);
		matrix = nullptr;
	}*/
}

void lzh::Matrix::create(int w)
{
	release();
	lzh::tools::tools::check(w);
	setsize(w);
	matrix.create(size());
#ifdef MAT_DEBUG
	CHECK_PTR(matrix.p);
#endif //MAT_DEBUG
	checkSquare();
}
void lzh::Matrix::create(int w, int h)
{
	release();
	lzh::tools::tools::check(w, h);
	setsize(w, h);
	matrix.create(size());
#ifdef MAT_DEBUG
	CHECK_PTR(matrix.p);
#endif //MAT_DEBUG
	checkSquare();
}
void lzh::Matrix::create(int w, int h, int c)
{
	release();
	lzh::tools::tools::check(w, h, c);
	setsize(w, h, c);
	matrix.create(size());
#ifdef MAT_DEBUG
	CHECK_PTR(matrix.p);
#endif //MAT_DEBUG
	checkSquare();
}
void lzh::Matrix::create(Size size)
{
	create(size.w, size.h);
}
void lzh::Matrix::create(Size3 size)
{
	create(size.w, size.h, size.c);
}

lzh::mat_t * lzh::Matrix::data() const
{
	return matrix.p;
}

lzh::mat_t * lzh::Matrix::begin()
{
	return matrix;
}

const lzh::mat_t * lzh::Matrix::begin() const
{
	return matrix;
}

lzh::mat_t * lzh::Matrix::end()
{
	if (matrix == nullptr)return nullptr;
	return matrix + h_ * w_*c_;
}

const lzh::mat_t * lzh::Matrix::end() const
{
	if (matrix == nullptr)return nullptr;
	return matrix + h_ * w_*c_;
}

lzh::uint lzh::Matrix::memsize() const
{
	return sizeof(lzh::mat_t)*h_*w_*c_;
}

void lzh::Matrix::DimCheck()const
{
	if (c_ != 1) {
		THROW_INFO(ERR_INFO_DIM);
	}
}

void lzh::Matrix::copy(Matrix &src, int Row_Start, int Row_End, int Col_Start, int Col_End)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	int hei = Row_End - Row_Start + 1;
	int wid = Col_End - Col_Start + 1;
	lzh::tools::tools::check(hei, wid);
	if (src.matrix == nullptr) {
		src = lzh::zeros(wid, hei, c_);
	}
	for (int y = Row_Start, j = 0; y <= Row_End; y++, j++)
		for (int x = Col_Start, i = 0; x <= Col_End; x++, i++)
			for (int z = 0; z < c_; z++)
				src(y, x, z) = (*this)(j, i, z);
}
void lzh::Matrix::swap(Matrix &src)const
{
	src.setvalue(*this);
}

const lzh::Matrix lzh::Matrix::addones(Dire dire) const
{
	Matrix temp(w_ + 1, h_, c_);
	for (int i = 0; i < h_; i++) {
		for (int j = 0; j < w_ + 1; j++) {
			for (int z = 0; z < c_; z++) {
				if (dire == LEFT) {
					if (j == 0)
						temp(i, 0, z) = 1;
					else
						temp(i, j, z) = (*this)(i, j - 1, z);
				}
				else if (dire == RIGHT) {
					if (j == w_)
						temp(i, w_, z) = 1;
					else
						temp(i, j, z) = (*this)(i, j, z);
				}
			}
		}
	}
	return temp;
}
const lzh::Matrix lzh::Matrix::addzeros(Dire dire) const
{
	Matrix temp(w_ + 1, h_, c_);
	for (int i = 0; i < h_; i++) {
		for (int j = 0; j < w_ + 1; j++) {
			for (int z = 0; z < c_; z++) {
				if (dire == LEFT) {
					if (j == 0)
						temp(i, 0, z) = 0;
					else
						temp(i, j, z) = (*this)(i, j - 1, z);
				}
				else if (dire == RIGHT) {
					if (j == w_)
						temp(i, w_, z) = 0;
					else
						temp(i, j, z) = (*this)(i, j, z);
				}
			}
		}
	}
	return temp;
}
void lzh::Matrix::mChannel(const lzh::Matrix & src, int channels)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	CHECK_MATRIX(src.matrix);
	if (h_ != src.h_ || w_ != src.w_ || channels >= c_)
		THROW_INFO(ERR_INFO_SIZE);	
#endif // MAT_DEBUG
	for (int i = 0; i < h_; ++i) {
		for (int j = 0; j < w_; ++j) {
			(*this)(i, j, channels) = src(i, j);
		}
	}
}
void lzh::Matrix::mChannel(const lzh::Matrix & src, int w, int h)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	CHECK_MATRIX(src.matrix);
	if (this->h_ <= h_ || this->w_ <= w_ || src.channels() != c_)
		THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	for (int i = 0; i < c_; ++i) {
		(*this)(h_, w_, i) = src(i);
	}
}
void lzh::Matrix::reshape(Size3 size)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	reshape(size.w, size.h, size.c);
}
void lzh::Matrix::reshape(int w, int h, int c)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	if (h == 0) {
		h = h_;
	}
	if (w == 0) {
		w = w_;
	}
	if (c == 0) {
		c = c_;
	}
	if (h == -1 && w == -1 && c == -1) {
		THROW_INFO(ERR_INFO_UNLESS);
	}
	bool check_ = true;
	if (h == -1) {
		h = len() / (w * c);
		if (check_)check_ = false;
		else {
			THROW_INFO(ERR_INFO_UNLESS);
		}
	}
	else if (w == -1) {
		w = len() / (h * c);
		if (check_)check_ = false;
		else {
			THROW_INFO(ERR_INFO_UNLESS);
		}
	}
	else if (c == -1) {
		c = len() / (h * w);
		if (check_)check_ = false;
		else {
			THROW_INFO(ERR_INFO_UNLESS);
		}
	}
	if (len() != h * w*c) {
		THROW_INFO(ERR_INFO_UNLESS);
	}
	else {
		setsize(w, h, c);
	}
}
bool lzh::Matrix::setSize(int w, int h, int c)
{
	if (h*w*c > 0) {
		*this = Matrix(w, h, c);
		return true;
	}
	if (len() == h * w * c) {
		setsize(w, h, c);
		return true;
	}
	return false;
}
void lzh::Matrix::setNum(lzh::mat_t number, int index)
{
#ifdef MAT_DEBUG
	checkindex(index);
#endif // MAT_DEBUG
	(*this)(index) = number;
}
void lzh::Matrix::setNum(lzh::mat_t number, int index_y, int index_x)
{
#ifdef MAT_DEBUG
	checkindex(index_x, index_y);
#endif // MAT_DEBUG
	(*this)(index_y, index_x) = number;
}
void lzh::Matrix::setMat(lzh::mat_t *mat, int hei, int wid)
{
	if ((h_ <= 0 || w_ <= 0)) return;
	*this = Matrix(mat, wid, hei);
	checkSquare();
}
void lzh::Matrix::setvalue(const lzh::Matrix &src)
{
	setsize(src.w_, src.h_, src.c_);
	square = src.square;
	matrix.release();
	matrix.create(src.len());
	memcpy(matrix, src.matrix, src.memsize());
}
void lzh::Matrix::setOpp()
{
	*this = lzh::mOpp(*this);
}
void lzh::Matrix::setAdj()
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (square)
		*this = adj();
	else THROW_INFO(ERR_INFO_ADJ);
}
void lzh::Matrix::setTran()
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	*this = t();
}
void lzh::Matrix::setInv()
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (enable() == 0)
		*this = inv();
	else THROW_INFO(ERR_INFO_INV);
}
void lzh::Matrix::setPow(lzh::mat_t num)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (enable() == 0)
		*this = pow(num);
	else THROW_INFO(ERR_INFO_POW);
}
void lzh::Matrix::setIden()
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (enable() == 0)
		*this = lzh::eye(w_);
	else THROW_INFO(ERR_INFO_IND);
}

lzh::Size3 lzh::Matrix::size3() const
{
	return Size3(w_, h_, c_);
}
int lzh::Matrix::total() const
{
	return step;
}
int lzh::Matrix::dims() const
{
	return dim;
}
int lzh::Matrix::rows()const
{
	return h_;
}
int lzh::Matrix::cols()const
{
	return w_;
}
int lzh::Matrix::channels() const
{
	return c_;
}
int lzh::Matrix::rank() const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (dim != 2)
		THROW_INFO(ERR_INFO_NOT2D);
#endif //MAT_DEBUG
	int rank = h_;
	for (int i = 0; i < h_; i++) {
		int count = 0;
		for (int j = 0; j < w_; j++) {
			if (matrix[j + i * w_] == 0) {
				count++;
			}
		}
		if (count == h_)
			rank--;
	}
	return rank;
}
lzh::uint lzh::Matrix::size()const
{
	return (uint)len();
}
lzh::Size lzh::Matrix::mSize()const
{
	return Size(w_, h_);
}

int lzh::Matrix::len()const
{
	return h_ * w_ * c_;
}
int lzh::Matrix::enable()const
{
	if (matrix == nullptr)
		return -1;
	if (!square)
		return -2;
	return 0;
}
bool lzh::Matrix::empty()const
{
	if (matrix == nullptr)return true;
	else return false;
}
bool lzh::Matrix::Square()const
{
	return square;
}

void lzh::Matrix::save(std::string file, bool binary) const
{
	if (binary) {
		FILE *out = fopen((file + ".bin").c_str(), "wb");
		if (out) {
			save(out);
			fclose(out);
		}
		else {
			THROW_INFO(ERR_INFO_FILE);
		}
	}
	else {
		std::ofstream out(file);
		if (out.is_open()) {
			show(out);
			out.close();
		}
		else {
			THROW_INFO(ERR_INFO_FILE);
		}
	}
}
void lzh::Matrix::save(FILE * out) const
{
	if (out) {
		int param[3] = { w_, h_, c_ };
		fwrite(param, sizeof(int) * 3, 1, out);
		fwrite(matrix, memsize(), 1, out);
	}
	else {
		THROW_INFO(ERR_INFO_FILE);
	}
}
void lzh::Matrix::load(std::string file)
{
	FILE *in = fopen((file + ".bin").c_str(), "rb");
	if (in) {
		load(in);
		fclose(in);
	}
	else {
		THROW_INFO(ERR_INFO_FILE);
	}
}
void lzh::Matrix::load(FILE * in)
{
	if (in) {
		int param[3] = { 0 };
		fread(param, sizeof(int) * 3, 1, in);
		create(param[0], param[1], param[2]);
		fread(matrix, memsize(), 1, in);
	}
	else {
		THROW_INFO(ERR_INFO_FILE);
	}
}
void lzh::Matrix::copyTo(Matrix & mat) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	mat.create(size3());
	for (int i = 0; i < len(); i++)
		mat(i) = (*this)(i);
}
void lzh::Matrix::copyTo(Matrix && mat) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if(mat.empty())
		mat.create(size3());
	for (int i = 0; i < len(); i++)
		mat(i) = (*this)(i);
}
void lzh::Matrix::setAddOnes(Dire dire)
{
	*this = addones(dire);
}
void lzh::Matrix::setAddZeros(Dire dire)
{
	*this = addzeros(dire);
}
void lzh::Matrix::release()
{
	matrix.release();
}

lzh::mat_t& lzh::Matrix::at(int w)const
{
	return (*this)(w);
}
lzh::mat_t& lzh::Matrix::at(int w, int h)const
{
	return (*this)(h, w);
}
lzh::mat_t & lzh::Matrix::at(int w, int h, int c) const
{
	return (*this)(h, w, c);
}
int lzh::Matrix::toX(int index)const
{
	return (index / c_) % w_;
}
int lzh::Matrix::toY(int index)const
{
	return (index / c_) / w_;
}
int lzh::Matrix::toZ(int index) const
{
	return index % c_;
}
lzh::mat_t lzh::Matrix::frist()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	return (*this)(0);
}
lzh::mat_t& lzh::Matrix::findAt(lzh::mat_t value)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	static lzh::mat_t err = NAN;
	for (int ind = 0; ind < len(); ind++)
		if ((*this)(ind) == value)
			return (*this)(ind);
	return err;
}
lzh::mat_t lzh::Matrix::max(bool is_abs) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	lzh::mat_t maxdata = is_abs ? fabs((*this)(0)) : (*this)(0);
	if (is_abs) {
		for (int ind = 1; ind < len(); ind++) {
			lzh::mat_t v = fabs((*this)(ind));
			if (maxdata < v)
				maxdata = v;
		}
	}
	else {
		for (int ind = 1; ind < len(); ind++) {
			lzh::mat_t v = (*this)(ind);
			if (maxdata < v)
				maxdata = v;
		}
	}
	return maxdata;
}
lzh::mat_t lzh::Matrix::min(bool is_abs) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	lzh::mat_t mindata = is_abs ? fabs((*this)(0)) : (*this)(0);
	if (is_abs) {
		for (int ind = 1; ind < len(); ind++) {
			lzh::mat_t v = fabs((*this)(ind));
			if (mindata > v)
				mindata = v;
		}
	}
	else {
		for (int ind = 1; ind < len(); ind++) {
			lzh::mat_t v = (*this)(ind);
			if (mindata > v)
				mindata = v;
		}
	}
	return mindata;
}
lzh::mat_t& lzh::Matrix::findmax()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	int max_adr = 0;
	for (int ind = 1; ind < len(); ind++)
		if ((*this)(max_adr) < (*this)(ind))
			max_adr = ind;
	return (*this)(max_adr);
}
lzh::mat_t& lzh::Matrix::findmin()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	int min_adr = 0;
	for (int ind = 1; ind < len(); ind++)
		if ((*this)(min_adr) < (*this)(ind))
			min_adr = ind;
	return (*this)(min_adr);
}
int lzh::Matrix::find(lzh::mat_t value)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	for (int ind = 0; ind < len(); ind++)
		if ((*this)(ind) == value)
			return ind;
	return -1;
}
int lzh::Matrix::maxAt()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	int max_adr = 0;
	for (int ind = 1; ind < len(); ind++)
		if ((*this)(max_adr) < (*this)(ind))
			max_adr = ind;
	return max_adr;
}
int lzh::Matrix::minAt()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	int min_adr = 0;
	for (int ind = 1; ind < len(); ind++)
		if ((*this)(min_adr) < (*this)(ind))
			min_adr = ind;
	return min_adr;
}
bool lzh::Matrix::contains(lzh::mat_t value)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	for (int ind = 0; ind < len(); ind++)
		if ((*this)(ind) == value)
			return true;
	return false;
}

void lzh::Matrix::show()const
{
	show(std::cout);
}
std::ostream & lzh::Matrix::show(std::ostream & out) const
{
	if(matrix == nullptr)return out;
	switch (print_type)
	{
	case lzh::FIXED:
		out.unsetf(std::ios::fixed);
		break;
	case lzh::SCIENTIFIC:
		out.setf(std::ios::scientific);
		break;
	default:
		break;
	}
	out.setf(std::ios::showpos);
	out.setf(std::ios::left);
	for (int i = 0; i < h_; i++) {
		out << "[";
		for (int j = 0; j < w_; j++) {
			out << "[";
			for (int k = 0; k < c_; k++) {
				if (j == w_ - 1 && k == c_ - 1)
				{
					out << std::setw(print_width) 
						<< std::setprecision(print_precision) 
						<< std::setfill(' ') << (*this)(i, j, k) << "]]";
				}
				else if (k == c_ - 1)
				{
					out << std::setw(print_width) 
						<< std::setprecision(print_precision) 
						<< std::setfill(' ') << (*this)(i, j, k) << "]";
				}
				else {
					out << std::setw(print_width) 
						<< std::setprecision(print_precision) 
						<< std::setfill(' ') << (*this)(i, j, k) << ", ";
				}
			}
		}
		/*if (i != h_ - 1)*/
		out << std::endl;
	}
	switch (print_type)
	{
	case lzh::FIXED:
		out.unsetf(std::ios::fixed);
		break;
	case lzh::SCIENTIFIC:
		out.unsetf(std::ios::scientific);
		break;
	default:
		break;
	}
	out.unsetf(std::ios::showpos);
	out.unsetf(std::ios::left);
	out << std::defaultfloat << std::setprecision(6);
	return out;
}

lzh::Matrix lzh::Matrix::Row(int h)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (h >= h_) THROW_INFO(ERR_INFO_MEMOUT);
#endif // DEBUG_MAT
	return Matrix(w_, 1, c_, c_, matrix + h * w_ * c_);
}
const lzh::Matrix lzh::Matrix::Row(int h) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (h >= h_) THROW_INFO(ERR_INFO_MEMOUT);
#endif // DEBUG_MAT
	return Matrix(w_, 1, c_, c_, matrix.p + h * w_ * c_);
}
lzh::Matrix lzh::Matrix::Col(int w)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (w < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (w >= w_) THROW_INFO(ERR_INFO_MEMOUT);
#endif // DEBUG_MAT
	return Matrix(1, h_, c_, c_*w_, matrix.p + w * c_);
}
const lzh::Matrix lzh::Matrix::Col(int w) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (w < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (w >= w_) THROW_INFO(ERR_INFO_MEMOUT);
#endif // DEBUG_MAT
	return Matrix(1, h_, c_, c_*w_, matrix.p + w * c_);
}
lzh::Matrix lzh::Matrix::Channel(int c)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (c < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (c >= c_) THROW_INFO(ERR_INFO_MEMOUT);
#endif // DEBUG_MAT
	return Matrix(w_, h_, 1, c_, matrix.p + c);
}
const lzh::Matrix lzh::Matrix::Channel(int c) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (c < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (c >= c_) THROW_INFO(ERR_INFO_MEMOUT);
#endif // DEBUG_MAT
	return Matrix(w_, h_, 1, c_, matrix.p + c);
}
const lzh::Matrix lzh::Matrix::Range(int start, int end)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (start < 0 || end < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (start >= h_ || end >= h_) THROW_INFO(ERR_INFO_MEMOUT);
	if (start >= end) THROW_INFO(ERR_INFO_UNLESS);
#endif // DEBUG_MAT
	return Matrix(end - start, 1, 1, matrix);
}
const lzh::Matrix lzh::Matrix::rowRange(int h_start, int h_end)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h_start < 0 || h_end < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (h_start >= h_ || h_end >= h_) THROW_INFO(ERR_INFO_MEMOUT);
	if (h_start >= h_end) THROW_INFO(ERR_INFO_UNLESS);
#endif // DEBUG_MAT
	/*return Matrix(h_end - h_start, w_, c_, c_, matrix + h_start * w_*c_);*/
	return lzh::Block(*this, h_start, h_end, 0, w_ - 1, 0, c_ - 1);
}
const lzh::Matrix lzh::Matrix::colRange(int w_start, int w_end)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (w_start < 0 || w_end < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (w_start >= h_ || w_end >= h_) THROW_INFO(ERR_INFO_MEMOUT);
	if (w_start >= w_end) THROW_INFO(ERR_INFO_UNLESS);
#endif // DEBUG_MAT
	//int newcol = w_end - w_start;
	//return Matrix(h_, newcol, c_, c_*w_, matrix + w_start * c_);
	return lzh::Block(*this, 0, h_ - 1, w_start, w_end, 0, c_ - 1);
}
const lzh::Matrix lzh::Matrix::channelRange(int c_start, int c_end)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (c_start < 0 || c_end < 0) THROW_INFO(ERR_INFO_UNLESS);
	if (c_start >= h_ || c_end >= h_) THROW_INFO(ERR_INFO_MEMOUT);
	if (c_start >= c_end) THROW_INFO(ERR_INFO_UNLESS);
#endif // DEBUG_MAT
	//int newcol = w_end - w_start;
	//return Matrix(h_, newcol, c_, c_*w_, matrix + w_start * c_);
	return lzh::Block(*this, 0, h_ - 1, 0, w_ - 1, c_start, c_end);
}
const lzh::Matrix lzh::Matrix::clone() const
{
	Matrix dst;
	dst.setvalue(*this);
	return dst;
}
const lzh::Matrix lzh::Matrix::opp()const
{
	return lzh::mOpp(*this);
}
const lzh::Matrix lzh::Matrix::adj()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (square) {
		return lzh::adj(*this);
	}
	else THROW_INFO(ERR_INFO_ADJ);
	return Matrix();
}
const lzh::Matrix lzh::Matrix::t()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	return lzh::tran(*this);
}
const lzh::Matrix lzh::Matrix::inv()const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (!square) 
		THROW_INFO(ERR_INFO_PINV);
	return lzh::inv(*this);
}
const lzh::Matrix lzh::Matrix::diag(int k) const
{
	return lzh::diag(*this, k);
}
const lzh::Matrix lzh::Matrix::abs()const
{
	return lzh::mAbs(*this);
}
const lzh::Matrix lzh::Matrix::mPow(int num) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	if (square) {
		return lzh::POW(*this, num);
	}
	else THROW_INFO(ERR_INFO_POW);
	return Matrix();
}
const lzh::Matrix lzh::Matrix::pow(lzh::mat_t num)const
{
	return lzh::mPow(*this, num);
}
const lzh::Matrix lzh::Matrix::reverse()const
{
	return lzh::reverse(*this);
}
lzh::mat_t lzh::Matrix::sum(int num, bool _abs)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	lzh::mat_t sum = 0;
	if (num == 0) {
		return (lzh::mat_t)len();
	}
	else if (num == 1) {
		for (int ind = 0; ind < len(); ind++)
			if (_abs)
				sum += std::abs((*this)(ind));
			else
				sum += (*this)(ind);
	}
	else
		for (int ind = 0; ind < len(); ind++)
			if (_abs)
				sum += std::pow(std::abs((*this)(ind)), num);
			else
				sum += std::pow((*this)(ind), num);
	return sum;
}
lzh::mat_t lzh::Matrix::mean() const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	lzh::mat_t sum = 0;
	for (int ind = 0; ind < len(); ind++)
		sum += (*this)(ind);
	return sum / (lzh::mat_t)len();
}
lzh::mat_t lzh::Matrix::S() const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	return std::sqrt(D());
}
lzh::mat_t lzh::Matrix::D() const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif //MAT_DEBUG
	lzh::mat_t m = mean();
	lzh::mat_t d = 0;
	lzh::mat_t n = (lzh::mat_t)len();
	for (lzh::mat_t v : *this)
	{
		d += std::pow(v - m, 2.0f);
	}
	return d / n;
}
lzh::mat_t lzh::Matrix::Det()
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	if (square)
		return lzh::det(*this);
	else
		return NAN;
}
lzh::mat_t lzh::Matrix::Norm(int num)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (w_ != 1 && h_ != 1) THROW_INFO(ERR_INFO_NORM);
	if (num < 0) THROW_INFO(ERR_INFO_VALUE);
#endif // MAT_DEBUG
	if (num == 0)
		return sum();
	else if (num == 1)
		return sum(1, true);
	else if (num == 2)
		return std::sqrt(sum(2, true));
	//else if (isinf(num) == 1)
	//	return abs(matrix[find(findmax())]);
	//else if (isinf(num) == -1)
	//	return abs(matrix[find(findmin())]);
	else
		return std::pow(sum(num, true), 1 / _T(num));
}
lzh::mat_t lzh::Matrix::Cof(int x, int y)
{
	return lzh::cof(*this, x, y);
}
lzh::mat_t lzh::Matrix::EigenvalueMax(lzh::mat_t offset)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	if (square) {
		int count = 0;
		lzh::mat_t err = 100 * offset;
		Matrix v;
		Matrix u0 = lzh::ones(w_, 1);
		while (err > offset) {
			v = *this*u0;
			Matrix u1 = v * (1 / v.findmax());
			err = (u1 - u0).abs().findmax();
			u0 = u1; count += 1;
			if (count >= 1e+3) THROW_INFO(ERR_INFO_EIGEN);
		}
		return v.findmax();
	}
	else THROW_INFO(ERR_INFO_SQUARE);
	return _T(NAN);
}
lzh::mat_t lzh::Matrix::RandSample()
{
	return lzh::mRandSample(*this);
}
const lzh::Matrix lzh::Matrix::EigenvectorsMax(lzh::mat_t offset)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	if (square) {
		int count = 0;
		lzh::mat_t err = 100 * offset;
		Matrix v;
		Matrix u0 = lzh::ones(h_, 1);
		while (err > offset) {
			v = *this*u0;
			Matrix u1 = v * (1 / v.findmax());
			err = (u1 - u0).abs().findmax();
			u0 = u1; count += 1;
			if (count >= 1e+3) THROW_INFO(ERR_INFO_EIGEN);
		}
		return u0;
	}
	else THROW_INFO(ERR_INFO_SQUARE);
	return Matrix();
}

const lzh::Matrix lzh::Matrix::sigmoid() const
{
	return lzh::Sigmoid(*this);
}
const lzh::Matrix lzh::Matrix::tanh() const
{
	return lzh::Tanh(*this);
}
const lzh::Matrix lzh::Matrix::relu() const
{
	return lzh::ReLU(*this);
}
const lzh::Matrix lzh::Matrix::elu() const
{
	return lzh::ELU(*this);
}
const lzh::Matrix lzh::Matrix::selu() const
{
	return lzh::SELU(*this);
}
const lzh::Matrix lzh::Matrix::leaky_relu() const
{
	return lzh::LReLU(*this);
}
const lzh::Matrix lzh::Matrix::softmax() const
{
	return lzh::Softmax(*this);
}
const lzh::Matrix lzh::Matrix::exp()const
{
	return lzh::mExp(*this);
}
const lzh::Matrix lzh::Matrix::log()const
{
	return lzh::mLog(*this);
}
const lzh::Matrix lzh::Matrix::sqrt()const
{
	return lzh::mSqrt(*this);
}

void lzh::Matrix::init()
{
	dim = 0;
	h_ = 0;
	w_ = 0;
	c_ = 0;
	matrix.release();
	step = 1;
}

void lzh::Matrix::checkSquare()
{
	if (h_ == w_)
		square = true;
	else
		square = false;
}
#ifdef MAT_DEBUG
void lzh::Matrix::checkindex(int index)
{
	if (h_ == 0 || w_ == 0) THROW_INFO(ERR_INFO_LEN);
	if (index > len() - 1) THROW_INFO(ERR_INFO_MEMOUT);
	if (index < 0) THROW_INFO(ERR_INFO_VALUE);
	if (!matrix) THROW_INFO(ERR_INFO_EMPTY);
}
void lzh::Matrix::checkindex(int index_x, int index_y)
{
	if (h_ == 0 || w_ == 0)THROW_INFO(ERR_INFO_LEN);
	if (index_x < 0 || index_y < 0) THROW_INFO(ERR_INFO_MEMOUT);
	if (index_x*w_ + index_y > h_*w_ - 1 || index_x >= h_ || index_y >= w_) THROW_INFO(ERR_INFO_VALUE);
	if (!matrix) THROW_INFO(ERR_INFO_EMPTY);
}
#endif // MAT_DEBUG
void lzh::Matrix::setsize(int w, int h, int c)
{
	w_ = w;
	h_ = h;
	c_ = c;
	step = c;
	if (w == 1 && h == 1 && c == 1)dim = 0;
	else if (w != 1 && h == 1 && c == 1)dim = 1;
	else if (w == 1 && h != 1 && c == 1)dim = 1;
	else if (w != 1 && h != 1 && c == 1)dim = 2;
	else if (w != 1 && h == 1 && c != 1)dim = 2;
	else if (w == 1 && h != 1 && c != 1)dim = 2;
	else dim = 3; 
	checkSquare();
}

const lzh::Matrix lzh::Matrix::operator + (const lzh::mat_t val)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	Matrix mark(size3());
	for (int i = 0; i < len(); i++)
		mark(i) = (*this)(i) + val;
	return mark;
}
const lzh::Matrix lzh::Matrix::operator + (const lzh::Matrix &a)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	if (h_ == 1 && w_ == 1 && c_ == 1) {
		return (*this)(0) + a;
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == 1) {
		return *this + a(0);
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == c_) {
		Matrix mat(size3());
		for (int i = 0; i < h_; i++)
			for (int j = 0; j < w_; j++)
				for (int z = 0; z < c_; z++)
					mat(i, j, z) = (*this)(i, j, z) + a(z);
		return mat;
	}
#ifdef MAT_DEBUG
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_) 
		THROW_INFO(ERR_INFO_SIZE);		
#endif // MAT_DEBUG
	Matrix mark(size3());
	for (int i = 0; i < len(); i++)
		mark(i) = (*this)(i) + a(i);
	return mark;
}
void lzh::Matrix::operator += (const lzh::mat_t val)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) += val;
}
void lzh::Matrix::operator += (const lzh::Matrix & a)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_)
		THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) += a(i);
}
const lzh::Matrix lzh::Matrix::operator-(void) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	Matrix mark(size3());
	for (int i = 0; i < len(); i++)
		mark(i) = -(*this)(i);
	return mark;
}
const lzh::Matrix lzh::Matrix::operator - (const lzh::mat_t val)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	Matrix mark(size3());
	for (int i = 0; i < len(); i++)
		mark(i) = (*this)(i) - val;
	return mark;
}
const lzh::Matrix lzh::Matrix::operator - (const lzh::Matrix &a)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	if (h_ == 1 && w_ == 1 && c_ == 1) {
		return (*this)(0) - a;
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == 1) {
		return *this - a(0);
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == c_) {
		Matrix mat(size3());
		for (int i = 0; i < h_; i++)
			for (int j = 0; j < w_; j++)
				for (int z = 0; z < c_; z++)
					mat(i, j, z) = (*this)(i, j, z) - a(z);
		return mat;
	}
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_)
		THROW_INFO(ERR_INFO_SIZE);
	Matrix mark(size3());
	for (int i = 0; i < len(); i++)
		mark(i) = (*this)(i) - a(i);
	return mark;
}
void lzh::Matrix::operator-=(const lzh::mat_t val)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) -= val;
}
void lzh::Matrix::operator-=(const lzh::Matrix & a)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_)
		THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) -= a(i);
}
const lzh::Matrix lzh::Matrix::operator * (const lzh::mat_t val)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	Matrix mark(size3());
	for (int i = 0; i < len(); i++)
		mark(i) = (*this)(i) * val;
	return mark;
}
const lzh::Matrix lzh::Matrix::operator * (const lzh::Matrix &a)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	CHECK_MATRIX(a.matrix);
#endif // MAT_DEBUG
	if (h_ == 1 && w_ == 1 && c_ == 1) {
		return (*this)(0) * a;
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == 1) {
		return *this * a(0);
	}
#ifdef MAT_DEBUG
	if (w_ != a.h_) THROW_INFO(ERR_INFO_MULT);
	if (c_ != a.c_) THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	Matrix mark(a.w_, h_, c_);
	for (int z = 0; z < c_; z++)
		for (int i = 0; i < h_; i++)
			for (int j = 0; j < a.w_; j++) {
				lzh::mat_t temp = 0;
				for (int d = 0; d < w_; d++)
					temp = temp + (*this)(i, d, z) * a(d, j, z);
				mark(i, j, z) = temp;
			}
	return mark;
}
void lzh::Matrix::operator*=(const lzh::mat_t val)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) *= val;
}
void lzh::Matrix::operator*=(const lzh::Matrix & a)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_) THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) *= a(i);
}
const lzh::Matrix lzh::Matrix::operator / (const lzh::mat_t val)const
{
	return (*this) * (1.0f / val);
}
const lzh::Matrix lzh::Matrix::operator / (const lzh::Matrix &a)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	CHECK_MATRIX(a.matrix);
#endif // MAT_DEBUG
	if (h_ == 1 && w_ == 1 && c_ == 1) {
		return (*this)(0) / a;
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == 1) {
		return *this / a(0);
	}
	else if (a.h_ == 1 && a.w_ == 1 && a.c_ == c_) {
		Matrix mat(size3());
		for (int i = 0; i < h_; i++)
			for (int j = 0; j < w_; j++)
				for (int z = 0; z < c_; z++)
					mat(i, j, z) = (*this)(i, j, z) / a(z);
		return mat;
	}
#ifdef MAT_DEBUG
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_)
		THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	Matrix mark(size3());
	for (int i = 0; i < h_; i++)
		for (int j = 0; j < w_; j++)
			for (int z = 0; z < c_; z++)
				mark(i, j, z) = (*this)(i, j, z) / a(i, j, z);
	return mark;
}
void lzh::Matrix::operator/=(const lzh::mat_t val)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) /= val;
}
void lzh::Matrix::operator/=(const lzh::Matrix & a)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h_ != a.h_ || w_ != a.w_ || c_ != a.c_) 
		THROW_INFO(ERR_INFO_SIZE);
#endif // MAT_DEBUG
	for (int i = 0; i < len(); i++)
		(*this)(i) = (*this)(i) / a(i);
}
const lzh::Matrix& lzh::Matrix::operator = (const lzh::Matrix &temp)
{
	if (this == &temp)
		return *this;
	matrix = temp.matrix;
	h_ = temp.h_;
	w_ = temp.w_;
	c_ = temp.c_;
	step = temp.step;
	dim = temp.dim;
	square = temp.square;
	return *this;
}
bool lzh::Matrix::operator == (const lzh::Matrix &a)const
{
	if (w_ != a.w_) {
		return false;
	}
	if (h_ != a.h_) {
		return false;
	}
	if (c_ != a.c_) {
		return false;
	}
	if (memcmp(matrix, a.matrix, memsize()) == 0)
		return true;
	return false;
}
bool lzh::Matrix::operator != (const lzh::Matrix & a)const
{
	return !(*this == a);
}
lzh::mat_t & lzh::Matrix::operator()(const int w) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (w > len() - 1) THROW_INFO(ERR_INFO_MEMOUT);
	if (w < 0) THROW_INFO(ERR_INFO_VALUE);
#endif // MAT_DEBUG
	return matrix[w * step / c_];
}
lzh::mat_t& lzh::Matrix::operator()(const int h, const int w) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h > h_ - 1 || w > w_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
	if (h < 0 || w < 0) THROW_INFO(ERR_INFO_VALUE);
#endif // MAT_DEBUG
	return matrix[(h*w_ + w)*step];
}
lzh::mat_t & lzh::Matrix::operator()(const int h, const int w, const int c) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (h > h_ - 1 || w > w_ - 1 || c > c_ - 1) 
		THROW_INFO(ERR_INFO_MEMOUT);
	if (h < 0 || w < 0 || c < 0) THROW_INFO(ERR_INFO_VALUE);
#endif // MAT_DEBUG
	return matrix[(h*w_ + w)*step + c];
}
lzh::mat_t & lzh::Matrix::operator()(Point p) const
{
	return this->at(p.x, p.y);
}
lzh::mat_t & lzh::Matrix::operator()(Point3i p) const
{
	return this->at(p.x, p.y, p.z);
}
const lzh::Matrix lzh::Matrix::operator()(const int index, X_Y_Z rc) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (index < 0) THROW_INFO(ERR_INFO_VALUE);
#endif // MAT_DEBUG
	switch (rc) {
	case ROW:
		if (index > this->h_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
		return lzh::Block(*this, index, index, 0, w_ - 1);
	case COL:
		if (index > this->w_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
		return lzh::Block(*this, 0, h_ - 1, index, index);
	case CHANNEL:
		if (index > c_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
		return lzh::mSplit(*this, index);
	default:return Matrix();
	}
}
const lzh::Matrix lzh::Matrix::operator()(const int v1, const int v2, X_Y_Z rc) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX(matrix);
	if (v1 < 0 || v2 < 0) THROW_INFO(ERR_INFO_VALUE);
#endif // MAT_DEBUG
	Matrix m;
	switch (rc) {
	case ROW:
#ifdef MAT_DEBUG
		if (v1 > w_ - 1 || v2 > c_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
#endif // MAT_DEBUG
		m.create(h_, 1, 1);
		for (int i = 0; i < h_; i++)
			m(i) = (*this)(i, v1, v2);
		break;
	case COL:
#ifdef MAT_DEBUG
		if (v1 > h_ - 1 || v2 > c_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
#endif // MAT_DEBUG
		m.create(1, w_, 1);
		for (int i = 0; i < w_; i++)
			m(i) = (*this)(v1, i, v2);
		break;
	case CHANNEL:
#ifdef MAT_DEBUG
		if (v1 > h_ - 1 || v2 > w_ - 1) THROW_INFO(ERR_INFO_MEMOUT);
#endif // MAT_DEBUG
		return Matrix(1, 1, c_, matrix.p + (v1*w_ + v2)*step);
		break;
	default:break;
	}
	return m;
}
const lzh::Matrix lzh::Matrix::operator[](const int idx) const
{
	if (dim == 3) {
		return Channel(idx);
	}
	else if (dim == 2) {
		return Row(idx);
	}
	else if (dim == 1) {
#ifdef MAT_DEBUG
		if (idx < 0) THROW_INFO(ERR_INFO_VALUE);
		if (idx >= len()) THROW_INFO(ERR_INFO_MEMOUT);
#endif // MAT_DEBUG
		return Matrix(1, 1, 1, matrix.p + idx);
	}
	else {
		THROW_INFO(ERR_INFO_UNLESS);
	}
	return Matrix();
}

void lzh::Matrix::setPrintW(lzh::uint w)
{
	print_width = w;
}
void lzh::Matrix::setPrintSignificantDigits(lzh::uint n)
{
	print_precision = n;
	//print_width = 2 + n + 3;
}
void lzh::Matrix::setPrintType(lzh::PrintType t)
{
	print_type = t;
}

const lzh::Matrix lzh::operator + (const lzh::mat_t value, const lzh::Matrix &mat)
{
	return mat + value;
}
const lzh::Matrix lzh::operator - (const lzh::mat_t value, const lzh::Matrix &mat)
{
	return value + (-mat);
}
const lzh::Matrix lzh::operator * (const lzh::mat_t value, const lzh::Matrix &mat)
{
	return mat * value;
}
const lzh::Matrix lzh::operator / (const lzh::mat_t value, const lzh::Matrix &mat)
{
	return lzh::Divi(mat, value, LEFT);
}
std::ostream & lzh::operator << (std::ostream &out, const lzh::Matrix &ma)
{
	return ma.show(out);
}

/****************************************************************************
矩阵操作
*****************************************************************************/
const lzh::Mat lzh::mSplit(const lzh::Mat & src, int c)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
	tools::check(src.cols(), src.rows(), src.channels());
	if (c > src.channels() - 1) 
		THROW_INFO(ERR_INFO_MEMOUT);
	if (c < 0)
		THROW_INFO(ERR_INFO_MEMOUT);
#endif //MAT_DEBUG
	Mat mat(src.size3());
	for (int i = 0; i < src.rows(); i++)
		for (int j = 0; j < src.cols(); j++)
			mat(i, j) = src(i, j, c);
	return mat;
}
void lzh::mSplit(const lzh::Mat & src, Mat * dst)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
#endif //MAT_DEBUG
	for (int c = 0; c < src.channels(); ++c)
		dst[c] = src.Channel(c);
}
const lzh::Mat lzh::Reshape(const lzh::Mat & src, Size3 size)
{
	Mat dst;
	src.swap(dst);
	dst.reshape(size);
	return dst;
}
const lzh::Mat lzh::mMerge(const lzh::Mat * src, int channels)
{
#ifdef MAT_DEBUG
	if (channels < 0)
		THROW_INFO(ERR_INFO_MEMOUT);
	if (src == nullptr)
		THROW_INFO(ERR_INFO_PTR);
	CHECK_MATRIX_(src[channels - 1]);
#endif //MAT_DEBUG
	Mat mat(src[0].cols(), src[0].rows(), channels);
	for (int z = 0; z < channels; z++) {
		for (int i = 0; i < src[z].rows(); i++)
			for (int j = 0; j < src[z].cols(); j++) {
				mat(i, j, z) = src[z](i, j);
			}
	}
	return mat;
}
const lzh::Mat lzh::reverse(const lzh::Mat &m)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(m);
	if (!(m.cols() == 1 || m.rows() == 1))
		THROW_INFO(ERR_INFO_MEMOUT);
#endif //MAT_DEBUG
	Mat temp = m.clone();
	for (uint ind = 0; ind < m.size() / 2; ind++) {
		lzh::mat_t val = temp(ind);
		temp(ind) = temp((int)m.size() - 1 - ind);
		temp((int)m.size() - 1 - ind) = val;
	}
	return temp;
}
const lzh::Mat lzh::copyMakeBorder(const lzh::Mat & src, int top, int bottom, int left, int right, BorderTypes borderType, lzh::mat_t value)
{
	Size3 size = src.size3();
	size.h += (top + bottom);
	size.w += (left + right);
	Mat mat(size);
	switch (borderType)
	{
	case BORDER_CONSTANT:
		for (int i = 0; i < top; i++) {
			for (int j = 0; j < mat.cols(); j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = value;
				}
			}
		}
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < left; j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = value;
				}
			}
		}
		for (int i = top + src.rows(); i < mat.rows(); i++) {
			for (int j = 0; j < mat.cols(); j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = value;
				}
			}
		}
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = left + src.cols(); j < mat.cols(); j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = value;
				}
			}
		}
		break;
	case BORDER_REPLICATE:
		for (int i = 0; i < top; i++) {
			for (int j = left; j < mat.cols() - right; j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(0, j - left, z);
				}
			}
		}
		for (int i = top; i < mat.cols() - bottom; i++) {
			for (int j = 0; j < left; j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(i - top, 0, z);
				}
			}
		}
		for (int i = top + src.rows(); i < mat.rows(); i++) {
			for (int j = left; j < mat.cols() - right; j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(src.rows() - 1, j - left, z);
				}
			}
		}
		for (int i = top; i < mat.cols() - bottom; i++) {
			for (int j = left + src.cols(); j < mat.cols(); j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(i - top, src.cols() - 1, z);
				}
			}
		}
		for (int i = 0; i < top; i++) {
			for (int j = 0; j < left; j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(0, 0, z);
				}
			}
		}
		for (int i = 0; i < top; i++) {
			for (int j = mat.cols() - right; j < mat.cols(); j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(0, src.cols() - 1, z);
				}
			}
		}
		for (int i = mat.rows() - bottom; i < mat.rows(); i++) {
			for (int j = 0; j < left; j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(src.rows() - 1, 0, z);
				}
			}
		}
		for (int i = mat.rows() - bottom; i < mat.rows(); i++) {
			for (int j = mat.cols() - right; j < mat.cols(); j++) {
				for (int z = 0; z < mat.channels(); z++) {
					mat(i, j, z) = src(src.rows() - 1, src.cols() - 1, z);
				}
			}
		}
		break;
	case BORDER_REFLECT:
		break;
	case BORDER_WRAP:
		break;
	case BORDER_REFLECT_101:
		break;
	case BORDER_TRANSPARENT:
		break;
	case BORDER_ISOLATED:
		break;
	default:
		break;
	}
	for (int i = 0; i < src.rows(); i++) {
		for (int j = 0; j < src.cols(); j++) {
			for (int z = 0; z < src.channels(); z++) {
				mat(i + top, j + left, z) = src(i, j, z);
			}
		}
	}
	return mat;
}
const lzh::Mat lzh::Block(const lzh::Mat&a, int Row_Start, int Row_End, int Col_Start, int Col_End, int Chennel_Start, int Chennel_End)
{
	int h = Row_End - Row_Start + 1;
	int w = Col_End - Col_Start + 1;
	int c = Chennel_End - Chennel_Start + 1;
	tools::check(w, h, c);
	Mat mark(w, h, c);
	int i = 0;
	for (int y = Row_Start, j = 0; y <= Row_End; y++, j++)
		for (int x = Col_Start, i = 0; x <= Col_End; x++, i++)
			for (int z = Chennel_Start, k = 0; z <= Chennel_End; z++, k++)
				mark(j, i, k) = a(y, x, z);
	return mark;
}

/****************************************************************************
生成矩阵
*****************************************************************************/
const lzh::Mat lzh::mcreate(int w)
{
	return zeros(w);
}
const lzh::Mat lzh::mcreate(int w, int h)
{
	return zeros(w, h);
}
const lzh::Mat lzh::mcreate(int w, int h, int c)
{
	return zeros(w, h, c);
}
const lzh::Mat lzh::mcreate(Size size)
{
	return zeros(size);
}
const lzh::Mat lzh::mcreate(Size3 size)
{
	return zeros(size);
}
const lzh::Matrix lzh::zeros(int w)
{
	tools::check(w);
	Mat mat(w);
	memset(mat, 0, sizeof(lzh::mat_t)*mat.len());
	return mat;
}
const lzh::Matrix lzh::zeros(int w, int h)
{
	tools::check(w, h);
	Mat mat(w, h);
	memset(mat, 0, sizeof(lzh::mat_t)*mat.len());
	return mat;
}
const lzh::Matrix lzh::zeros(int w, int h, int c)
{
	tools::check(w, h, c);
	Mat mat(w, h, c);
	memset(mat, 0, sizeof(lzh::mat_t)*mat.len());
	return mat;
}
const lzh::Matrix lzh::zeros(Size size)
{
	tools::check(size.h, size.w);
	Mat mat(size.h, size.w);
	memset(mat, 0, sizeof(lzh::mat_t)*mat.len());
	return mat;
}
const lzh::Matrix lzh::zeros(Size3 size)
{
	Mat mat(size);
	memset(mat, 0, sizeof(lzh::mat_t)*mat.len());
	return mat;
}
const lzh::Matrix lzh::value(lzh::mat_t v, int w)
{
	tools::check(w);
	Mat mark(w);
	for (int ind = 0; ind < mark.len(); ind++)
		mark(ind) = v;
	return mark;
}
const lzh::Matrix lzh::value(lzh::mat_t v, int w, int h)
{
	tools::check(w, h);
	Mat mark(w, h);
	for (int ind = 0; ind < mark.len(); ind++)
		mark(ind) = v;
	return mark;
}
const lzh::Matrix lzh::value(lzh::mat_t v, int w, int h, int c)
{
	tools::check(w, h, c);
	Mat mark(w, h, c);
	for (int ind = 0; ind < mark.len(); ind++)
		mark(ind) = v;
	return mark;
}
const lzh::Matrix lzh::ones(int w)
{
	return value(1, w);
}
const lzh::Matrix lzh::ones(int w, int h)
{
	return value(1, w, h);
}
const lzh::Matrix lzh::ones(int w, int h, int c)
{
	return value(1, w, h, c);
}
const lzh::Matrix lzh::ones(Size3 size)
{
	return value(1, size.w, size.h, size.c);
}
const lzh::Matrix lzh::ones(Size size)
{
	tools::check(size.w, size.h);
	return value(1, size.w, size.h);
}
const lzh::Matrix lzh::range(int low, int top, mat_t gap)
{
	return range(_T(low), _T(top), gap);
}
const lzh::Matrix lzh::range(lzh::mat_t low, lzh::mat_t top, lzh::mat_t gap)
{
#ifdef MAT_DEBUG
	if (low >= top)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	int len = (int)((top - low) / gap + 1);
	Mat mark(len);
	mark = low + linspace(0, len - 1, len)*gap;
	if (mark.enable() != -1) {
		mark(0) = low;
		mark(len - 1) = top;
	}
	return mark;
}
const lzh::Matrix lzh::linspace(int low, int top, int len)
{
#ifdef MAT_DEBUG
	if (low >= top)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	tools::check(len, len);
	Mat mark(len);
	mark(0) = (lzh::mat_t)low;
	lzh::mat_t gap = _T(std::abs(low) + std::abs(top)) / (len - 1);
	for (int ind = 1; ind < len; ind++)
		mark(ind) = mark(ind - 1) + gap;
	return mark;
}
const lzh::Matrix lzh::linspace(lzh::mat_t low, lzh::mat_t top, int len)
{
#ifdef MAT_DEBUG
	if (low >= top)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	tools::check(len, len);
	Mat mark(len);
	lzh::mat_t gap = (top - low) / (len - 1);
	mark = low + linspace(0, len - 1, len)*gap;
	if (mark.enable() != -1) {
		mark(0) = low;
		mark(len - 1) = top;
	}
	return mark;
}
const lzh::Mat lzh::eye(int n)
{
	tools::check(n, n);
	Mat mark = zeros(n, n);
	for (int ind = 0; ind < n; ind++)
		mark(ind + ind * n) = 1;
	return mark;
}

/****************************************************************************
矩阵数学工具
*****************************************************************************/
lzh::mat_t lzh::Max(const lzh::Mat &temp, bool isAbs)
{
	if (isAbs)
		return mAbs(temp).findmax();
	else
		return temp.findmax();
}
lzh::mat_t lzh::Min(const lzh::Mat &temp, bool isAbs)
{
	if (isAbs)
		return mAbs(temp).findmin();
	else
		return temp.findmin();
}
lzh::mat_t lzh::trace(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	if (temp.enable() == -2) THROW_INFO(ERR_INFO_SQUARE);
#endif //MAT_DEBUG
	lzh::mat_t sum = 0;
	for (int index = 0; index < temp.rows(); index++) {
		sum += temp((index + index * temp.cols())*temp.channels());
	}
	return sum;
}
lzh::mat_t lzh::cof(const lzh::Mat &temp, int x, int y)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	if (x >= temp.cols() || y >= temp.rows())
		THROW_INFO(ERR_INFO_MEMOUT);
#endif //MAT_DEBUG
	temp.DimCheck();
	Mat a(temp.cols() - 1, temp.rows() - 1);
	int n = temp.rows();
	for (int i = 0, k = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if ((i != x) && (j != y)) {
				a(k) = temp(i*n + j);
				k++;
			}
	return det(a);
}
lzh::mat_t lzh::det(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	temp.DimCheck();
	if (temp.enable() == -2)
		THROW_INFO(ERR_INFO_SQUARE);
#endif //MAT_DEBUG
	int n = temp.rows();
	if (n == 1)
		return temp(0);
	else {
		Mat a = temp.clone();
		for (int j = 0; j < n; j++)
			for (int i = 0; i < n; i++) {
				if (a(j + j * n) == 0) {
					lzh::mat_t m;
					for (int d = j + 1; d < n; d++)
						if (a(j + d * n) != 0) {
							for (int f = j; f < n; f++)
								a(f + j * n) += a(f + d * n);
							m = -a(j + d * n) / a(j + j * n);
							for (int f = j; f < n; f++)
								a(f + d * n) += a(f + j * n) * m;
						}
				}
				else if (i != j) {
					lzh::mat_t w = -a(j + i * n) / a(j + j * n);
					for (int f = j; f < n; f++)
						a(f + i * n) += a(f + j * n) * w;
				}
			}
		lzh::mat_t answer = 1;
		for (int i = 0; i < n; i++)
			answer *= a(i + i * n);
		return answer;
	}
}
lzh::mat_t lzh::mNorm(const lzh::Mat &temp, int num)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	if (temp.cols() != 1 && temp.rows() != 1)
		THROW_INFO(ERR_INFO_NORM);
	temp.DimCheck();
	if (num < 0)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	if (num == 1)
		return temp.sum(1, true);
	else if (num == 2)
		return sqrt(temp.sum(2, true));
	//else if (isinf(num) == 1)
	//	return abs(matrix[find(findmax())]);
	//else if (isinf(num) == -1)
	//	return abs(matrix[find(findmin())]);
	else
		return pow(temp.sum(num, true), 1 / _T(num));
}
lzh::mat_t lzh::mDistance(const lzh::Mat &a, const lzh::Mat &b, int num)
{
	return (a - b).Norm(num);
}
const lzh::Mat lzh::adj(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	if (temp.enable() == -2)
		THROW_INFO(ERR_INFO_SQUARE);
#endif //MAT_DEBUG
	int n = temp.rows();
	int depth = temp.channels();
	Mat a(n, n, depth);
	for (int z = 0; z < depth; z++) {
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				lzh::mat_t m = cof(temp, i, j);
				a((i*n + j)*depth + z) = (lzh::mat_t)pow(-1, i + j + 2)*m;
			}
	}
	return tran(a);
}
const lzh::Mat lzh::inv(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	if (temp.channels() != 1)
		THROW_INFO(ERR_INFO_DIM);
#endif //MAT_DEBUG
	Mat m;
	temp.swap(m);
	lzh::mat_t answer = det(m);
	if (answer != 0 && answer == answer) {
		m = adj(m);
		int n = m.rows();
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				m(i, j) = (1 / answer)*m(i, j);
	}
	else
		THROW_INFO(ERR_INFO_DET);
	return m;
}
const lzh::Mat lzh::diag(const lzh::Mat & src, int k)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
	if (src.channels() != 1) THROW_INFO(ERR_INFO_DIM);
	if (src.dims() != 1) THROW_INFO(ERR_INFO_NORM);
#endif //MAT_DEBUG
	int pos = k < 0 ? -k : k;
	Mat mat = zeros(src.len() + pos, src.len() + pos);
	const mat_t *data = src;
	mat_t *p = mat;
	for (int h = 0; h < mat.rows(); h++)
	{
		for (int w = 0; w < mat.cols(); w++)
		{
			if (h + k == w)
				*p = *data;
			p += 1;
		}
		data += 1;
	}
	return mat;
}
const lzh::Mat lzh::pinv(const lzh::Mat &temp, Dire direc)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	switch (direc)
	{
	case LEFT:return (temp.t()*temp).inv()*temp.t();
	case RIGHT: {
		Mat m = temp.t();
		return lzh::pinv(m, LEFT).t();
	}
	default:
		THROW_INFO(ERR_INFO_VALUE);
	}
	return Mat();
}
const lzh::Mat lzh::tran(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat a(temp.rows(), temp.cols(), temp.channels());
	int n = temp.rows(),
		m = temp.cols();
	for (int z = 0; z < temp.channels(); z++)
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				a(j, i, z) = temp(i, j, z);
	return a;
}
const lzh::Mat lzh::mAbs(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat m(temp.size3());
	for (int ind = 0; ind < temp.len(); ind++)
		m(ind) = fabs(temp(ind));
	return m;
}
const lzh::Mat lzh::POW(const lzh::Mat &temp, int num)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
	if (temp.enable() == -2) 
		THROW_INFO(ERR_INFO_SQUARE);
#endif //MAT_DEBUG
	Mat m;
	temp.swap(m);
	if (num > 0) {
		for (int i = 1; i < num; i++)
			m = m * temp;
		return m;
	}
	else if (num < 0) {
		Mat a(temp);
		temp.swap(a);
		m.setInv();
		a.setInv();
		for (int i = -1; i > num; i--)
			a = a * m;
		return a;
	}
	else
		return eye(temp.rows());

}
const lzh::Mat lzh::mPow(const lzh::Mat &temp, lzh::mat_t num)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat m(temp.size3());
	for (int ind = 0; ind < temp.len(); ind++)
		m(ind) = pow(temp(ind), num);
	return m;
}
const lzh::Mat lzh::mSum(const lzh::Mat &temp, X_Y_Z r_c)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	if (r_c == COL) {
		Mat m = zeros(temp.cols(), 1, 1);
		for (int i = 0; i < temp.cols(); i++)
			for (int z = 0; z < temp.channels(); z++)
				for (int j = 0; j < temp.rows(); j++)
					m(i) += temp(j, i, z);
		return m;
	}
	else if (r_c == ROW) {
		Mat m = zeros(1, temp.rows(), 1);
		for (int i = 0; i < temp.rows(); i++)
			for (int z = 0; z < temp.channels(); z++)
				for (int j = 0; j < temp.cols(); j++)
					m(i) += temp(i, j, z);
		return m;
	}
	else if (r_c == CHANNEL) {
		Mat m = zeros(1, 1, temp.channels());
		for (int z = 0; z < temp.channels(); z++)
			for (int i = 0; i < temp.rows(); i++)
				for (int j = 0; j < temp.cols(); j++)
					m(z) += temp(i, j, z);
		return m;
	}
	else
		THROW_INFO(ERR_INFO_VALUE);
	return Mat();
}
const lzh::Mat lzh::mExp(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat m(temp.size3());
	for (int ind = 0; ind < temp.len(); ind++) {
		m(ind) = exp(temp(ind));
		if (m(ind) == 0)
			m(ind) = (std::numeric_limits<lzh::mat_t>::min)();
	}
	return m;
}
const lzh::Mat lzh::mLog(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat m(temp.size3());
	for (int ind = 0; ind < temp.len(); ind++)
		if (temp(ind) == 0)
			m(ind) = (std::numeric_limits<lzh::mat_t>::min)();
		else
			m(ind) = log(temp(ind));
	return m;
}
const lzh::Mat lzh::mSqrt(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat m(temp.size3());
	for (int ind = 0; ind < temp.len(); ind++)
		m(ind) = sqrt(temp(ind));
	return m;
}
const lzh::Mat lzh::mOpp(const lzh::Mat &temp)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(temp);
#endif //MAT_DEBUG
	Mat m(temp.size3());
	for (int ind = 0; ind < temp.len(); ind++)
		m(ind) = -temp(ind);
	return m;
}
const lzh::Mat lzh::Divi(const lzh::Mat &a, lzh::mat_t val, Dire dire)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(a);
#endif //MAT_DEBUG
	Mat mark(a.size3());
	for (int ind = 0; ind < mark.len(); ind++)
		if (dire == LEFT)
			mark(ind) = val / a(ind);
		else if (dire == RIGHT)
			mark(ind) = a(ind) / val;
	return mark;
}
const lzh::Mat lzh::Divi(const lzh::Mat &a, const lzh::Mat &b, Dire dire)
{
	switch (dire)
	{
	case LEFT:return a.inv()*b;
	case RIGHT:return a / b;
	default:THROW_INFO(ERR_INFO_VALUE);
	}
	return Mat();
}
const lzh::Mat lzh::Mult(const lzh::Mat &a, const lzh::Mat &b)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(a);
	CHECK_MATRIX_(b);
#endif //MAT_DEBUG
	if (b.rows() == 1 && b.cols() == 1 && a.channels() == b.channels()) {
		Mat mat(a.size3());
		for (int z = 0; z < mat.channels(); z++) {
			mat_t v = b(z);
			for (int i = 0; i < mat.rows(); i++)
				for (int j = 0; j < mat.cols(); j++)
					mat(i, j, z) = a(i, j, z) * v;
		}
		return mat;
	}
	else if (a.rows() == 1 && a.cols() == 1 && a.channels() == b.channels())
	{
		Mat mat(b.size3());
		for (int z = 0; z < mat.channels(); z++) {
			mat_t v = a(z);
			for (int i = 0; i < mat.rows(); i++)
				for (int j = 0; j < mat.cols(); j++)
					mat(i, j, z) = b(i, j, z) * v;
		}
		return mat;
	}
#ifdef MAT_DEBUG
	else if (a.rows() != b.rows() || a.cols() != b.cols() || a.channels() != b.channels())
		THROW_INFO(ERR_INFO_SIZE);
#endif //MAT_DEBUG
	Mat temp(a.size3());
	for (int ind = 0; ind < a.len(); ind++)
		temp(ind) = a(ind) * b(ind);
	return temp;
}
const lzh::Mat lzh::Dot(const lzh::Mat & a, const lzh::Mat & b)
{
	return a * b;
}
const lzh::Mat lzh::mMax(lzh::mat_t a, const lzh::Mat &b)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(b);
#endif //MAT_DEBUG
	Mat mark(b.size3());
	for (int ind = 0; ind < b.len(); ind++)
		mark(ind) = a > b(ind) ? a : b(ind);
	return mark;
}
const lzh::Mat lzh::mMax(const lzh::Mat &a, const lzh::Mat &b)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(a);
	CHECK_MATRIX_(b);
	if (a.rows() != b.rows() || a.cols() != b.cols() || a.channels() != b.channels())
		THROW_INFO(ERR_INFO_SIZE);
#endif //MAT_DEBUG
	Mat mark(b.size3());
	for (int ind = 0; ind < b.len(); ind++)
		mark(ind) = a(ind) > b(ind) ? a(ind) : b(ind);
	return mark;
}
const lzh::Mat lzh::mMin(lzh::mat_t a, const lzh::Mat &b)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(b);
#endif //MAT_DEBUG
	Mat mark(b.size3());
	for (int ind = 0; ind < b.len(); ind++)
		mark(ind) = a < b(ind) ? a : b(ind);
	return mark;
}
const lzh::Mat lzh::mMin(const lzh::Mat &a, const lzh::Mat &b)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(a);
	CHECK_MATRIX_(b);
	if (a.rows() != b.rows() || a.cols() != b.cols() || a.channels() != b.channels())
		THROW_INFO(ERR_INFO_SIZE);
#endif //MAT_DEBUG
	Mat mark(b.size3());
	for (int ind = 0; ind < b.len(); ind++)
		mark(ind) = a(ind) < b(ind) ? a(ind) : b(ind);
	return mark;
}

/****************************************************************************
数值运算
*****************************************************************************/
const lzh::Mat lzh::LeastSquare(const lzh::Mat & x, const lzh::Mat & y)
{
	return ((x.t()*x).inv()*x.t()*y).t();
}
const lzh::Mat lzh::regress(const lzh::Mat & x, const lzh::Mat & y)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(x);
	CHECK_MATRIX_(y);
	if (x.rows() != y.rows())
		THROW_INFO(ERR_INFO_SIZE);
#endif //MAT_DEBUG
	return LeastSquare(x.addones(RIGHT), y);
}
const lzh::Mat lzh::polyfit(const lzh::Mat & x, const lzh::Mat & y, uint n)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(x);
	CHECK_MATRIX_(y);
	if (x.rows() != y.rows())
		THROW_INFO(ERR_INFO_SIZE);
#endif //MAT_DEBUG
	Mat param(n + 1, x.rows());
	for (uint idx = n, i = 0; idx > 0; idx--, i++)
	{
		x.pow((lzh::mat_t)idx).copyTo(param.Col(i));
	}
	ones(1, x.rows()).copyTo(param.Col(n));
	return LeastSquare(param, y);
}
const lzh::Mat lzh::circlefit(const Mat & p)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(p);
#endif //MAT_DEBUG
	Mat y = zeros(1, p.rows());
	Mat x = zeros(p.cols() + 1, p.rows());
	for (int w = 0; w < p.cols(); ++w)
	{
		y += p.Col(w).pow(2);
		p.Col(w).copyTo(x.Col(w));
	}
	ones(1, p.rows()).copyTo(x.Col(p.cols()));
	Mat param = LeastSquare(y, x);
	return param;
}
const lzh::Mat lzh::polynomial(const lzh::Mat & a, const lzh::Mat & x)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(a);
	CHECK_MATRIX_(x);
#endif //MAT_DEBUG
	Mat y = zeros(x.size3());
	for (int idx = 0, power = a.len() - 1; idx < a.len(); ++idx, --power)
	{
		y += a(idx)*x.pow(_T(power));
	}
	return y;
}
const lzh::Mat lzh::NonLinearLeastSqures(const Mat & x, const Mat & y, const Mat & a0, F f, mat_t step, mat_t error)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
	CHECK_MATRIX_(x);
	CHECK_MATRIX_(y);
	CHECK_MATRIX_(a0);
#endif //MAT_DEBUG
	Mat a = a0;
	Mat da = a;
	while (da.Norm() > error) {
		try {
			Mat p = NumericGradient(f, a, x, _T(1e-3));
			Mat q = f(a, x) - y;
			da = step * LeastSquare(p, q);
			a = a - da;
		}
		catch (std::exception e)
		{
			std::cout << "拟合失败!" << std::endl;
			break;
		}
	}
	return a;
}
const lzh::Mat lzh::gradient(const Mat & y, const Mat & x)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(y);
#endif //MAT_DEBUG
	int dim = y.dims();
	Mat df(y.cols(), y.rows(), dim*y.channels());
	for (int c = 0; c < y.channels(); ++c) {
		if (x.empty()) {
			if (y.dims() >= 2) {
				for (int row = 0; row < y.rows(); ++row)
					for (int col = 0; col < y.cols(); ++col)
						if (col == 0) {
							df(row, col, dim*c) = y(row, col + 1, c) - y(row, col, c);
						}
						else if (col == y.cols() - 1) {
							df(row, col, dim*c) = y(row, col, c) - y(row, col - 1, c);
						}
						else {
							df(row, col, dim*c) = ((y(row, col + 1, c) - y(row, col, c)) + (y(row, col, c) - y(row, col - 1, c))) / _T(2.0);
						}
				for (int row = 0; row < y.rows(); ++row)
					for (int col = 0; col < y.cols(); ++col)
						if (row == 0) {
							df(row, col, dim*c + 1) = y(row + 1, col, c) - y(row, col, c);
						}
						else if (row == y.rows() - 1) {
							df(row, col, dim*c + 1) = y(row, col, c) - y(row - 1, col, c);
						}
						else {
							df(row, col, dim*c + 1) = ((y(row + 1, col, c) - y(row, col, c)) + (y(row, col, c) - y(row - 1, col, c))) / _T(2.0);
						}
			}
			else {
				for (int idx = 0; idx < y.len(); ++idx)
					if (idx == 0) {
						df(idx) = y(idx + 1) - y(idx);
					}
					else if (idx == y.len() - 1) {
						df(idx) = y(idx) - y(idx - 1);
					}
					else {
						df(idx) = ((y(idx + 1) - y(idx)) + (y(idx) - y(idx - 1))) / _T(2.0);
					}
			}
		}
		else {
			if (y.dims() >= 2) {
				for (int row = 0; row < y.rows(); ++row)
					for (int col = 0; col < y.cols(); ++col)
						if (col == 0) {
							df(row, col, dim*c) = (y(row, col + 1, c) - y(row, col, c)) / (x(row, col + 1, c) - x(row, col, c));
						}
						else if (col == y.cols() - 1) {
							df(row, col, dim*c) = (y(row, col, c) - y(row, col - 1, c)) / (x(row, col, c) - x(row, col - 1, c));
						}
						else {
							df(row, col, dim*c) = ((y(row, col + 1, c) - y(row, col, c)) + (y(row, col, c) - y(row, col - 1, c))) / ((x(row, col + 1, c) - x(row, col, c)) + (x(row, col, c) - x(row, col - 1, c)));
						}
				for (int row = 0; row < y.rows(); ++row)
					for (int col = 0; col < y.cols(); ++col)
						if (row == 0) {
							df(row, col, dim*c + 1) = (y(row + 1, col, c) - y(row, col, c)) / (x(row + 1, col, c) - x(row, col, c));
						}
						else if (row == y.rows() - 1) {
							df(row, col, dim*c + 1) = (y(row, col, c) - y(row - 1, col, c)) / (x(row, col, c) - x(row - 1, col, c));
						}
						else {
							df(row, col, dim*c + 1) = ((y(row + 1, col, c) - y(row, col, c)) + (y(row, col, c) - y(row - 1, col, c))) / ((x(row, col + 1, c) - x(row, col, c)) + (x(row, col, c) - x(row, col - 1, c)));
						}
			}
			else {
				for (int idx = 0; idx < y.len(); ++idx)
					if (idx == 0) {
						df(idx) = (y(idx + 1) - y(idx)) / (x(idx + 1) - x(idx));
					}
					else if (idx == y.len() - 1) {
						df(idx) = (y(idx) - y(idx - 1)) / (x(idx) - x(idx - 1));
					}
					else {
						df(idx) = ((y(idx + 1) - y(idx)) + (y(idx) - y(idx - 1))) / ((x(idx + 1) - x(idx)) + (x(idx) - x(idx - 1)));
					}
			}
		}
	}
	return df;
}
lzh::mat_t lzh::NumericGradient(NF f, mat_t x, mat_t epsilon)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
#endif //MAT_DEBUG
	mat_t tmp1 = x;
	mat_t tmp2 = x;
	tmp1 = tmp1 + epsilon;
	tmp2 = tmp2 - epsilon;
	return (f(tmp1) - f(tmp2)) / (_T(2) * epsilon);
}
const lzh::Mat lzh::NumericGradient(Fun f, const Mat & x, mat_t epsilon)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
#endif //MAT_DEBUG
	int n = x.len();
	Mat numgrad = zeros(n, 1);
	for (int i = 0; i < n; ++i) {
		Mat tmp1 = x.clone();
		Mat tmp2 = x.clone();
		tmp1(i) = tmp1(i) + epsilon;
		tmp2(i) = tmp2(i) - epsilon;
		numgrad(i) = (f(tmp1) - f(tmp2)) / (_T(2) * epsilon);
	}
	return numgrad;
}
const lzh::Mat lzh::NumericGradient(F f, const Mat & a, const Mat & x, mat_t epsilon)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
	CHECK_MATRIX_(a);
	CHECK_MATRIX_(x);
#endif //MAT_DEBUG
	int n = a.len();
	Mat numgrad = zeros(n, x.rows());
	for (int i = 0; i < n; ++i) {
		Mat tmp1 = a.clone();
		Mat tmp2 = a.clone();
		tmp1(i) = tmp1(i) + epsilon;
		tmp2(i) = tmp2(i) - epsilon;
		((f(tmp1, x) - f(tmp2, x)) / (_T(2) * epsilon)).copyTo(numgrad.Col(i));
	}
	return numgrad;
}
lzh::mat_t lzh::EulerInt(NF f, mat_t low, mat_t high, mat_t epsilon)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
#endif //MAT_DEBUG
	mat_t y = 0;
	for (mat_t x = low; x <= high; x += epsilon)
	{
		y += f(x) * epsilon;
	}
	return y;
}
lzh::mat_t lzh::TrapezoidInt(NF f, mat_t low, mat_t high, mat_t epsilon)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
#endif //MAT_DEBUG
	mat_t y = 0;
	for (mat_t x = low; x <= high; x += epsilon)
	{
		y += _T(0.5) * epsilon * (f(x) + f(x + epsilon));
	}
	return y;
}
lzh::mat_t lzh::RungeKuttaInt(NF f, mat_t low, mat_t high, mat_t epsilon)
{
#ifdef MAT_DEBUG
	CHECK_PTR(f);
#endif //MAT_DEBUG	
	mat_t y = 0;
	for (mat_t x = low; x <= high; x += epsilon)
	{
		mat_t k1 = f(x);
		mat_t k2 = f(x + _T(0.5)*epsilon);
		mat_t k3 = f(x + _T(0.5)*epsilon);
		mat_t k4 = f(x + epsilon);
		y += epsilon * (k1 + 2 * k2 + 2 * k3 + k4) / _T(6);
	}
	return y;
}
const lzh::Mat lzh::LinearIntersection(const lzh::Mat & line_1, const lzh::Mat & line_2)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(line_1);
	CHECK_MATRIX_(line_2);
	if (line_1.size3() != line_2.size3())
		THROW_INFO(ERR_INFO_SIZE);
#endif //MAT_DEBUG
	return Mat();
}
/****************************************************************************
图像运算
*****************************************************************************/
lzh::Size3 lzh::mCalSize(Size3 src, Size3 kern, Point anchor, Size strides, int & top, int & bottom, int & left, int & right)
{
	int kern_row = kern.h;
	int kern_col = kern.w;
	top = anchor.x;
	bottom = kern_row - anchor.x - 1;
	left = anchor.y;
	right = kern_col - anchor.y - 1;
	return Size3(int(_T(src.w - kern_col + 1) / _T(strides.w) + 0.5f), int(_T(src.h - kern_row + 1) / _T(strides.h) + 0.5f), kern.c / src.c);
}
lzh::Size3 lzh::mCalSize(Size3 src, Size3 kern, Point &anchor, Size strides)
{
	int kern_row = kern.h;
	int kern_col = kern.w;
	if (anchor == Point(-1, -1)) {
		anchor.x = kern_row / 2;
		anchor.y = kern_col / 2;
	}
	return Size3(int(_T(src.w - kern_col + 1) / _T(strides.w) + 0.5f), int(_T(src.h - kern_row + 1) / _T(strides.h) + 0.5f), kern.c / src.c);
}
lzh::Size3 lzh::mCalSize(Size3 src, Size kern, Point & anchor, Size strides)
{
	int kern_row = kern.h;
	int kern_col = kern.w;
	if (anchor == Point(-1, -1)) {
		anchor.x = kern_row / 2;
		anchor.y = kern_col / 2;
	}
	return Size3(int(_T(src.w - kern_col + 1) / _T(strides.w) + 0.5f), int(_T(src.h - kern_row + 1) / _T(strides.h) + 0.5f), src.c);
}
const lzh::Mat lzh::Threshold(const lzh::Mat & src, lzh::mat_t boundary, lzh::mat_t lower, lzh::mat_t upper, int boundary2upper)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
#endif //MAT_DEBUG
	Mat mark;
	src.swap(mark);
	switch (boundary2upper)
	{
	case -1:
		for (int ind = 0; ind < mark.len(); ind++)
			mark(ind) = mark(ind) <= boundary ? lower : upper;
		break;
	case 0:
		for (int ind = 0; ind < mark.len(); ind++)
			mark(ind) = mark(ind) >= boundary ? upper : lower;
		break;
	case 1:
		for (int ind = 0; ind < mark.len(); ind++)
			mark(ind) = mark(ind) < boundary ? lower : (mark(ind) == boundary ? boundary : upper);
		break;
	default:
		THROW_INFO(ERR_INFO_UNLESS);
	}
	return mark;
}
const lzh::Mat lzh::Filter2D(const lzh::Mat & input, const lzh::Mat & kern, Point anchor, const Size & strides, bool is_copy_border)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(input);
	CHECK_MATRIX_(kern);
	if (input.dims() != 2)
		THROW_INFO(ERR_INFO_NOT2D);
	if (kern.dims() != 2)
		THROW_INFO(ERR_INFO_NOT2D);
#endif //MAT_DEBUG
	int kern_row = kern.rows();
	int kern_col = kern.cols();
	int left, right, top, bottom;
	Size3 size = mCalSize(input.size3(), kern.size3(), anchor, strides, left, right, top, bottom);
	Mat dst, src;
	if (is_copy_border) {
		src = copyMakeBorder(input, top, bottom, left, right);
		dst = zeros(input.cols() / strides.h, input.rows() / strides.w);
	}
	else {
		src = input;
		dst = zeros(size.w, size.h);
	}
	for (int h = top, x = 0; h < src.rows() - bottom; h += strides.h, x++)
		for (int w = left, y = 0; w < src.cols() - right; w += strides.w, y++) {
			lzh::mat_t value = 0;
			for (int i = 0; i < kern_row; ++i) {
				for (int j = 0; j < kern_col; ++j) {
					value += src(h + i - anchor.x, w + j - anchor.y)*kern(i, j);
				}
			}
			dst(x, y) = value;
		}
	return dst;
}
const lzh::Mat lzh::SumChannel(const lzh::Mat & src)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
#endif //MAT_DEBUG
	Mat m = zeros(src.size3());
	for (int i = 0; i < src.rows(); i++)
		for (int j = 0; j < src.cols(); j++) {
			lzh::mat_t sum = 0;
			for (int z = 0; z < src.channels(); z++)
				sum += src(i, j, z);
			m(i, j) = sum;
		}
	return m;
}
const lzh::Mat lzh::rotate(const lzh::Mat & src, RotateAngle dice)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
#endif //MAT_DEBUG
	Mat dst;
	switch (dice)
	{
	case lzh::ROTATE_90_ANGLE:
	{
		Mat mat(src.rows(), src.cols(), src.channels());
		for (int h = 0; h < src.rows(); ++h)
			for (int w = 0; w < src.cols(); ++w)
				for (int depth = 0; depth < src.channels(); ++depth)
					mat(w, src.rows() - 1 - h, depth) = src(h, w, depth);
		dst = mat;
	}
	break;
	case lzh::ROTATE_180_ANGLE:
	{
		Mat mat(src.cols(), src.rows(), src.channels());
		for (int h = 0, y = src.rows() - 1; h < src.rows() && y >= 0; ++h, --y)
			for (int w = 0, x = src.cols() - 1; w < src.cols() && x >= 0; ++w, --x)
				for (int depth = 0; depth < src.channels(); ++depth)
					mat(h, w, depth) = src(y, x, depth);
		dst = mat;
	}
	break;
	case lzh::ROTATE_270_ANGLE:
	{
		Mat mat(src.rows(), src.cols(), src.channels());
		for (int h = 0; h < src.rows(); ++h)
			for (int w = 0; w < src.cols(); ++w)
				for (int depth = 0; depth < src.channels(); ++depth)
					mat(src.cols() - 1 - w, h, depth) = src(h, w, depth);
		dst = mat;
	}
	break;
	default:
		break;
	}
	return dst;
}
const lzh::Mat lzh::Rotate(lzh::mat_t angle)
{
	Mat rotate_mat(2, 2);
	angle = angle * PI / 180.0f;
	rotate_mat(0) = cos(angle);
	rotate_mat(1) = -sin(angle);
	rotate_mat(2) = sin(angle);
	rotate_mat(3) = cos(angle);
	return rotate_mat;
}
void lzh::Filter2D(const lzh::Mat & input, Mat & out, const lzh::Mat & kern, Point anchor, const Size & strides, bool is_copy_border)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(input);
	CHECK_MATRIX_(kern);
#endif //MAT_DEBUG
	Mat in;
	int kern_row = kern.rows();
	int kern_col = kern.cols();
	int left, right, top, bottom;
	Size3 size = mCalSize(input.size3(), kern.size3(), anchor, strides, left, right, top, bottom);
	if (is_copy_border) {
		in = copyMakeBorder(input, top, bottom, left, right);
		out = zeros(input.cols() / strides.w, input.rows() / strides.h, size.c);
	}
	else {
		in = input;
		out = zeros(size);
	}
	for (int kern_c = 0; kern_c < size.c; kern_c++) {
		for (int in_c = 0; in_c < input.channels(); in_c++) {
			int c = kern_c * input.channels() + in_c;
			for (int h = top, x = 0; h < in.rows() - bottom; h += strides.h, x++)
				for (int w = left, y = 0; w < in.cols() - right; w += strides.w, y++) {
					int px = w - anchor.y;
					int py = h - anchor.x;
					lzh::mat_t value = 0;
					for (int i = 0; i < kern_row; ++i) {
						for (int j = 0; j < kern_col; ++j) {
							value += in(py + i, px + j, in_c)*kern(i, j, c);
						}
					}
					out(x, y, kern_c) += value;
				}
		}
	}
}
void lzh::resize(const Mat & src, Mat & dst, mat_t xRatio, mat_t yRatio, ReductionMothed mothed)
{
	if (src.empty())return;
	int rows = static_cast<int>(src.rows() * yRatio);
	int cols = static_cast<int>(src.cols() * xRatio);
	Mat img(cols, rows, src.channels());
	switch (mothed)
	{
	case lzh::EqualIntervalSampling:
		for (int i = 0; i < rows; i++) {
			int h = static_cast<int>((i + 1) / yRatio + _T(0.5)) - 1;
			for (int j = 0; j < cols; j++) {
				int w = static_cast<int>((j + 1) / xRatio + _T(0.5)) - 1;
				for (int z = 0; z < src.channels(); z++) {
					img(i, j, z) = src(h, w, z); //取得采样像素
				}
			}
		}
		break;
	case lzh::LocalMean:
	{
		int lastRow = 0;
		int lastCol = 0;

		for (int i = 0; i < rows; i++) {
			int h = static_cast<int>((i + 1) / yRatio + _T(0.5)) - 1;
			for (int j = 0; j < cols; j++) {
				int w = static_cast<int>((j + 1) / xRatio + _T(0.5)) - 1;
				Vec3m temp;
				for (int idx = lastCol; idx <= w; idx++) {
					for (int jdx = lastRow; jdx <= h; jdx++) {
						temp[0] += src(jdx, idx, 0);
						temp[1] += src(jdx, idx, 1);
						temp[2] += src(jdx, idx, 2);
					}
				}

				int count = (w - lastCol + 1) * (h - lastRow + 1);
				img(i, j, 0) = temp[0] / count;
				img(i, j, 1) = temp[1] / count;
				img(i, j, 2) = temp[2] / count;

				lastCol = w + 1; //下一个子块左上角的列坐标，行坐标不变
			}
			lastCol = 0; //子块的左上角列坐标，从0开始
			lastRow = h + 1; //子块的左上角行坐标
		}
	}
	break;
	default:
		break;
	}
	dst = img;
}
void lzh::resize(const Mat & src, Mat & dst, Size newSize, ReductionMothed mothed)
{
	resize(src, dst, newSize.w / _T(src.cols()), newSize.h / _T(src.rows()), mothed);
}

/****************************************************************************
随机数
*****************************************************************************/
void lzh::Srandom()
{
	srand(uint(time(NULL)));
}
lzh::mat_t lzh::generateGaussianNoise(mat_t mu, mat_t sigma)
{
	const mat_t epsilon = std::numeric_limits<mat_t>::min();
	const mat_t two_pi = _T(2.0)*PI;


	static mat_t z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	mat_t u1, u2;
	do
	{
		u1 = rand() * (_T(1.0) / RAND_MAX);
		u2 = rand() * (_T(1.0) / RAND_MAX);
	} while (u1 <= epsilon);

	z0 = sqrt(_T(-2.0) * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(_T(-2.0) * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}
const lzh::Matrix lzh::Xavier(Size3 size, int n1, int n2)
{
	Mat m(size);
	float *p = m;
	for (int i = 0; i < m.len(); ++i)
	{
		*p = generateGaussianNoise() * _T(1.0) / std::sqrt(_T(n1 + n2));
		//*p = generateGaussianNoise(0, 1) * sqrt(6.0f / (w_ + h_));
		p++;
	}
	return m;
}
const lzh::Matrix lzh::Xavier(int w, int h, int c, int n1, int n2)
{
	Mat m(w, h, c);
	float *p = m;
	for (int i = 0; i < m.len(); ++i)
	{
		*p = generateGaussianNoise() * _T(1.0) / std::sqrt(_T(n1 + n2));
		//*p = generateGaussianNoise(0, 1) * sqrt(6.0f / (w_ + h_));
		p++;
	}
	return m;
}
const lzh::Matrix lzh::Random(Size3 size)
{
	Mat m(size);
	float *p = m;
	for (int i = 0; i < m.len(); ++i)
	{
		*p = generateGaussianNoise();
		p++;
	}
	return m;
}
const lzh::Matrix lzh::Random(int w, int h, int c)
{
	Mat m(w, h, c);
	float *p = m;
	for (int i = 0; i < m.len(); ++i)
	{
		*p = generateGaussianNoise();
		p++;
	}
	return m;
}
const lzh::Mat lzh::mRandSample(const lzh::Mat &src, int w, int h, int c)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
#endif //MAT_DEBUG
	tools::check(w, h, c);
	Mat dst(w, h, c);
	for (uint ind = 0; ind < src.size(); ind++)
		dst(ind) = mRandSample(src);
	return dst;
}
const lzh::Mat lzh::mRandSample(const lzh::Mat& m, X_Y_Z rc, int num)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(m);
#endif //MAT_DEBUG
	Mat dst = m(rand() % m.rows(), rc);
	for (int i = 1; i < num; i++) {
		dst = Mat(dst, m(rand() % m.rows(), rc), rc);
	}
	return dst;
}
const lzh::Mat lzh::mRand(int low, int top, int n, bool isdouble)
{
	return mRand(low, top, n, n, 1, isdouble);
}
const lzh::Mat lzh::mRand(int low, int top, Size3 size, bool isdouble)
{
#ifdef MAT_DEBUG
	tools::check(size.w, size.h, size.c);
	if (low >= top)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	Mat m(size);
	for (int index = 0; index < m.len(); index++)
		m(index) = getRandData(low, top, isdouble);
	return m;
}
const lzh::Mat lzh::mRand(int low, int top, int w, int h, int c, bool isdouble)
{
#ifdef MAT_DEBUG
	tools::check(w, h, c);
	if (low >= top)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	Mat m(w, h, c);
	for (int index = 0; index < m.len(); index++)
		m(index) = getRandData(low, top, isdouble);
	return m;
}
lzh::mat_t lzh::getRandData(int min, int max, bool isdouble)
{
#ifdef MAT_DEBUG
	if (min > max)
		THROW_INFO(ERR_INFO_VALUE);
#endif //MAT_DEBUG
	if (isdouble) {
		lzh::mat_t m1 = (lzh::mat_t)(rand() % 101) / 101;
		min++;
		lzh::mat_t m2 = (lzh::mat_t)((rand() % (max - min + 1)) + min);
		m2 = m2 - 1;
		return m1 + m2;
	}
	else {
		int m = rand() % (max - min) + 1 + min;
		return (lzh::mat_t)m;
	}
}
lzh::mat_t lzh::mRandSample(const lzh::Mat &m)
{
	int h = rand() % m.rows();
	int w = rand() % m.cols();
	int depth = rand() % m.channels();
	return m(h, w, depth);
}

/****************************************************************************
k均值聚类
*****************************************************************************/
const lzh::Mat lzh::kmeans(const Mat & P, Mat & k, const uint K, const uint iteration, const mat_t error)
{
	uint iterate = 0;
	Mat center(P.cols(), K);
	Mat group(P);
	Mat dist(P.rows(), K);
	Mat kclass(1, P.rows());

	//随机生成聚点
	for (uint i = 0; i < K; ++i) {
		Mat p = P.Row(rand() % P.rows());
		p.copy(center, i, i, 0, P.cols() - 1);
		for (uint j = i; j != 0; j--)
			if (center.Row(i) == center.Row(i - 1)) {
				i--;
				break;
			}
	}
	//center.show();

	mat_t error_ = _T(1);
	Mat point = center;
	Mat a0 = zeros(1, K);

	while (iterate++ != iteration && error_ > error) {
		center = point;

		//计算聚点与所有点的距离
		for (uint i = 0; i < K; ++i) {
			for (uint j = 0; j < uint(P.rows()); ++j) {
				dist(i, j) = (center(i, ROW) - group(j, ROW)).Norm(2);
			}
		}
		//数据点按最小距离分组
		for (uint j = 0; j < uint(P.rows()); ++j) {
			int min = 0;
			for (uint i = 1; i < K; ++i) {
				if (dist(i, j) < dist(min, j)) {
					min = i;
				}
			}
			kclass(j) = _T(min);
		}
		//计算聚点与所属数据点的距离
		Mat a = zeros(1, K);
		for (uint i = 0; i < K; ++i) {
			for (uint j = 0; j < uint(P.rows()); ++j) {
				if (kclass(j) == i)
					a(i) += dist(i, j);
			}
		}
		//a.show();
		//计算与上一次距离的二范数
		error_ = (a0 - a).Norm(2);
		a0 = a;

		point = zeros(P.cols(), K);
		Mat sum = zeros(K);
		//计算聚点所属数据点坐标平均值作为新的聚点
		for (uint i = 0; i < K; ++i) {
			Mat p = zeros(P.cols());
			for (uint j = 0; j < uint(P.rows()); ++j) {
				if (kclass(j) == i) {
					sum(i) += 1;
					p += group(j, ROW);
				}
			}
			p /= sum(i);
			p.copy(point, i, i, 0, P.cols() - 1);
		}
	}
	kclass.swap(k);
	return center;
}

/****************************************************************************
求解方程
*****************************************************************************/
int lzh::RowSimplest(const Mat & src, Mat & dst)
{
	dst = RowSimplest(src);
	return dst.rank();
}
const lzh::Mat lzh::RowSimplest(const Mat & src)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(src);
	if (src.dims() != 2)
		THROW_INFO(ERR_INFO_NOT2D);
#endif //MAT_DEBUG 
	Mat dst = src.clone();
	Mat t(dst.size3());
	int height = dst.rows();
	int width = dst.cols();
	Mat zeronum(2 * height);
	mat_t *matrix = dst;
	int *ZeroNum = new int[2 * height];
	for (int count = 0; count < height; count++) {
		for (int i = 0; i < height; i++) {
			int count = 0;
			for (int j = 0; j < width; j++) {
				if (matrix[i*width + j] == 0)
					count++;
				else break;
			}
			ZeroNum[i * 2] = i;
			ZeroNum[i * 2 + 1] = count;
		}
		for (int i = 0; i < height; i++) {
			for (int j = i; j < height; j++) {
				if (ZeroNum[i * 2 + 1] > ZeroNum[j * 2 + 1]) {
					sort::swap(ZeroNum[i * 2], ZeroNum[j * 2]);
					sort::swap(ZeroNum[i * 2 + 1], ZeroNum[j * 2 + 1]);
				}
			}
		}
		mat_t *temp = t;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				temp[j + i * width] = matrix[j + ZeroNum[i * 2] * width];
			}
		}
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				matrix[j + i * width] = temp[j + i * width];
			}
		}
		if ((width - ZeroNum[count * 2 + 1]) != (width - count))continue;
		if (matrix[count + count * width] != 1) {
			mat_t m = _T(1) / matrix[count + count * width];
			for (int i = 0; i < width; i++) {
				matrix[i + count * width] *= m;
			}
		}
		for (int i = count + 1; i < height; i++) {
			if ((width - ZeroNum[i * 2 + 1]) != (width - count))continue;
			mat_t m = matrix[count + i * width] / matrix[count + count * width];
			for (int j = ZeroNum[i * 2 + 1]; j < width; j++) {
				matrix[j + i * width] -= (m*matrix[j + count * width]);
			}
		}
	}
	FREE_PTR(ZeroNum);
	for (int i = 1; i < height; i++) {
		bool flag = false; 
		int j;
		for (j = 0; j < width; j++) {
			if (matrix[i*width + j] != 0) {
				flag = true; 
				break;
			}
		}
		if (!flag)continue;
		for (int k = i - 1; k >= 0; k--) {
			mat_t m = matrix[j + k * width] / matrix[j + i * width];
			for (int l = j; l < width; l++) {
				matrix[l + k * width] = matrix[l + k * width] - m * matrix[l + i * width];
			}
		}
	}
	return dst;
}
int lzh::ColSimplest(const Mat & src, Mat & dst)
{
	return RowSimplest(src.t(), dst);
}
const lzh::Mat lzh::ColSimplest(const Mat & src)
{
	return RowSimplest(src.t());
}
lzh::EQUATION lzh::SolveLinearEquation(const Mat & src, Mat & dst, Mat *simplest, Mat *mark)
{
	Mat answer = RowSimplest(src);
	mat_t * matrix = answer;
	int height = answer.rows();
	int width = answer.cols();
	int rank = height;
	int augmentedRank = rank;
	for (int i = 0; i < height; i++) {
		int augmented_count = 0;
		int rank_count = 0;
		for (int j = 0; j < width; j++) {
			if (j < width - 1) {
				if (matrix[j + i * width] == 0) {
					rank_count++;
					augmented_count++;
				}
			}
			else {
				if (matrix[j + i * width] == 0) {
					augmented_count++;
				}
			}
		}
		if (rank_count == width - 1)
			rank--;
		if (augmented_count == width)
			augmentedRank--;
	}
	EQUATION state;
	if (rank >= augmentedRank) {
		if (rank == height){
			answer.Col(width - 1).copyTo(dst);
			if(mark!=nullptr)
				*mark = ones(height);
			state = SPECIAL_SOLUTION;
		}
		else {
			int freeparam = height - rank; 
			int *concer = new int[rank];
			int *markfree = new int[freeparam];
			for (int h = 0; h < rank; h++) {
				for (int w = h; w < width; w++) {
					if (answer(h, w) == 1)
					{
						concer[h] = w; break;
					}
				}
			}
			int count = 0;
			for (int i = 0; i < height; i++) {
				bool flag = false;
				for (int j = 0; j < rank; j++)
					if (i == concer[j]) {
						flag = true;
						break;
					}
				if (flag)continue;
				markfree[count++] = i;
			}
			dst = zeros(freeparam, height);
			Mat param = eye(freeparam);
			for (int h = 0; h < param.rows(); h++) {
				for (int i = 0; i < rank; i++) {
					mat_t v = answer(i, width - 1);
					for (int j = 0; j < freeparam; j++) {
						v -= answer(i, markfree[j])*param(h, j);
					}
					dst(concer[i], h) = v;
				}
				for (int i = 0; i < freeparam; i++) {
					dst(markfree[i], h) = param(h, i);
				}
			}
			if (mark != nullptr) {
				*mark = ones(height);
				for (int i = 0; i < freeparam; i++) {
					mark->at(markfree[i]) = _T(0);
				}
			}
			FREE_PTR(concer);
			FREE_PTR(markfree);
			state = GENERAL_SOLUTION;
		}
	}
	else {
		state = NO_SOLUTION;
	}
	if (simplest != nullptr)
		*simplest = answer;
	return state;
}

/****************************************************************************
排序
*****************************************************************************/
void lzh::sort::heapdown(mat_t *m, int i, int n, ORDER order)
{
	for (int j = 2 * i + 1; j < n; i = j, j = 2 * j + 1) {
		if (j + 1 < n && (order == MAX_TO_MIN ? m[j + 1] < m[j] : m[j + 1] > m[j]))
			j++;
		if (order == MAX_TO_MIN ? m[j] >= m[i] : m[j] <= m[i])
			break;
		swap(m[j], m[i]);
	}
}
void lzh::sort::heapup(mat_t *m, int i, ORDER order)
{
	for (int j = (i - 1) / 2; j >= 0 && order == MAX_TO_MIN ? m[i] < m[j] : m[i] > m[j]; i = j, j = (i - 1) / 2)
		swap(m[i], m[j]);
}
void lzh::sort::makeheap(mat_t *m, int length, ORDER order)
{
	for (int i = 0; i < length; i++)
		heapup(m, i, order);
}
void lzh::sort::heaparray(mat_t *m, int i, int n, ORDER order)
{
	swap(m[i], m[n]);
	heapdown(m, i, n, order);
}
void lzh::sort::mergearray(mat_t * a, mat_t * b, int start, int mid, int end, ORDER order)
{
	int i = start,
		j = mid + 1,
		k = 0;
	while (i <= mid && j <= end) {
		if (order == MAX_TO_MIN ? a[i] >= a[j] : a[i] <= a[j])
			b[k++] = a[i++];
		else
			b[k++] = a[j++];
	}
	while (i <= mid)
		b[k++] = a[i++];
	while (j <= end)
		b[k++] = a[j++];
	for (i = 0; i < k; i++)
		a[start + i] = b[i];
}
void lzh::sort::_merge(mat_t * a, mat_t * b, int start, int end, ORDER order)
{
	int mid;
	if (start < end) {
		mid = (start + end) / 2;
		_merge(a, b, start, mid, order);
		_merge(a, b, mid + 1, end, order);
		mergearray(a, b, start, mid, end, order);
	}
}
void lzh::sort::_quick(mat_t *m, int low, int high, ORDER order)
{
	int start, end;
	mat_t mark;
	if (low >= high)return;
	start = low;
	end = high;
	mark = m[start];
	while (start < end) {
		while (start < end && (order == MAX_TO_MIN ? m[end] <= mark : m[end] >= mark)) {
			end -= 1;
		}
		m[start] = m[end];
		while (start < end && (order == MAX_TO_MIN ? m[start] >= mark : m[start] <= mark)) {
			start += 1;
		}
		m[end] = m[start];
	}
	m[start] = mark;
	_quick(m, low, start - 1, order);
	_quick(m, start + 1, high, order);
}
void lzh::sort::bubble(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
CHECK_PTR(begin);
CHECK_PTR(end);
CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int i, j, 
		length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0)
		for (i = 0; i < length - 1; i++)
			for (j = i + 1; j < length; j++) {
				if (order == MAX_TO_MIN ? m[i] < m[j] : m[i] > m[j])
					swap(m[i], m[j]);
			}
}
const lzh::Mat lzh::sort::bubble(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	bubble(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::insert(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int i, j, k;
	int length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0)
		for (i = 0; i < length; i++)
			for (j = 0; j < i; j++) {
				if (order == MAX_TO_MIN ? m[i] > m[j] : m[i] < m[j])
					for (k = i; k > j; k--)
						swap(m[k], m[k - 1]);
			}
}
const lzh::Mat lzh::sort::insert(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	insert(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::select(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int i, j,
		index = 0,
		length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0)
		for (i = 0; i < length - 1; index = ++i) {
			for (j = i + 1; j < length; j++) {
				if (order == lzh::MAX_TO_MIN ? m[index] < m[j] : m[index] > m[j])
					index = j;
			}
			if (i != index)
				swap(m[i], m[index]);
		}
}
const lzh::Mat lzh::sort::select(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	select(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::comb(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int length = ARRAY_LEN(begin, end);
	int j, i, k;
	mat_t * m = begin;
	if (length > 0)
		for (k = 0, i = (int)(length / shrink_factor); i > 1 || k; i = (i > 1) ? (int)(i / shrink_factor) : i) {
			k = 0;
			for (j = 0; j < length - i; j++) {
				if ((order == MAX_TO_MIN ? m[j] < m[j + i] : m[j] > m[j + i])) {
					swap(m[j], m[j + i]);
					k = 1;
				}
			}
		}
}
const lzh::Mat lzh::sort::comb(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	comb(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::gnome(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int i, flag,
		length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0)
		for (i = 0; i < length&&i >= 0; flag ? i++ : i--, flag = 0) {
			if (i == 0 || (order == MAX_TO_MIN ? m[i - 1] >= m[i] : m[i - 1] <= m[i]))
				flag = 1;
			else
				swap(m[i - 1], m[i]);
		}
}
const lzh::Mat lzh::sort::gnome(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	gnome(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::heap(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int length = ARRAY_LEN(begin, end);
	if (length > 0) {
		makeheap(begin, length, order);
		for (int i = length - 1; i >= 1; i--)
			heaparray(begin, 0, i, order);
	}
}
const lzh::Mat lzh::sort::heap(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	heap(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::shell(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int j, i, k,
		length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0)
		for (i = length >> 1; i > 0; i >>= 1)
			for (j = i; j < length; j++) {
				for (k = j - i; k >= 0 && (order == lzh::MAX_TO_MIN ? m[k] < m[k + i] : m[k] > m[k + i]); k -= i)
					swap(m[k], m[k + i]);
			}
}
const lzh::Mat lzh::sort::shell(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	shell(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::merge(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	mat_t *p;
	int length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0) {
		p = (mat_t *)malloc(length * sizeof(mat_t));
		CHECK_PTR(p);
		_merge(m, p, 0, length - 1, order);
		free(p);
	}
}
const lzh::Mat lzh::sort::merge(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	merge(dst.begin(), dst.end(), order);
	return dst;
}
void lzh::sort::quick(mat_t * begin, mat_t * end, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_PTR(begin);
	CHECK_PTR(end);
	CHECK_PTR_ORDER(begin, end);
#endif //MAT_DEBUG
	int length = ARRAY_LEN(begin, end);
	mat_t * m = begin;
	if (length > 0)
		_quick(m, 0, length - 1, order);
}
const lzh::Mat lzh::sort::quick(const Mat & mat, ORDER order)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(mat);
#endif //MAT_DEBUG
	Mat dst = mat.clone();
	quick(dst.begin(), dst.end(), order);
	return dst;
}

/****************************************************************************
B-样条曲线类
*****************************************************************************/
lzh::BSpline::BSpline() :
	type(UNIFORM), k(0),
	P(), nodevector()
{

}
lzh::BSpline::BSpline(BSplineType type, int k, const Mat &p) :
	type(type), k(k),
	P(p), nodevector()
{

}
void lzh::BSpline::setCtrlPoint(const Mat &p)
{
	P = p;
}
const lzh::Mat lzh::BSpline::CtrlPoint()const
{
	return P;
}
const lzh::Mat lzh::BSpline::Node()const
{
	return nodevector;
}
lzh::mat_t lzh::BSpline::BF(int i, int k, mat_t t)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(nodevector);
#endif //MAT_DEBUG
	if (k == 0)
		if (t >= nodevector(i) && t < nodevector(i + 1))
			return _T(1);
		else
			return _T(0);
	else {
		mat_t len1 = nodevector(i + k) - nodevector(i);
		mat_t len2 = nodevector(i + k + 1) - nodevector(i + 1);
		if (len1 == _T(0))
			len1 = _T(1);
		if (len2 == _T(0))
			len2 = _T(1);
		mat_t s1 = (t - nodevector(i)) / len1 * BF(i, k - 1, t);
		mat_t s2 = (nodevector(i + k + 1) - t) / len2 * BF(i + 1, k - 1, t);
		return s1 + s2;
	}
}
const lzh::Mat lzh::BSpline::BaseFunction(mat_t t)const
{
	Mat mark(1, P.cols());
	for (int i = 0; i < P.cols(); i++)
		mark(i) = BF(i, k, t);
	return mark;
}
void lzh::BSpline::NodeVector(const Mat &node)
{
	if (!node.empty())
		nodevector = node;
	else {
		switch (type)
		{
		case UNIFORM:
			nodevector = linspace(0, 1, P.cols() + k + 1);
			break;
		case QUASI_UNIFORM: {
			nodevector = zeros(P.cols() + k + 1);
			int linepage = P.cols() - k;
			if (linepage == 1)
				for (int ind = P.cols() + 1; ind < P.cols() + k + 1; ind++)
					nodevector(ind) = 1;
			else {
				int judge = 1;
				while (judge != linepage) {
					nodevector(k + judge) = nodevector(k + judge - 1) + 1 / linepage;
					judge++;
				}
				for (int ind = P.cols() + 1; ind < P.cols() + k + 1; ind++)
					nodevector(ind) = 1;
			}
			break;
		}
		default:
			THROW_INFO(ERR_INFO_UNLESS);
		}
	}
}
const lzh::Mat lzh::BSpline::BPiont(const Mat & T)const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(P);
	CHECK_MATRIX_(nodevector);
	CHECK_MATRIX_(T);
#endif //MAT_DEBUG
	Mat m(2, T.len());
	for (int i = 0; i < T.len(); ++i) {
		Mat p_u = P * BaseFunction(T(i));
		//p_u.show();
		m(2 * i) = p_u(0);
		m(2 * i + 1) = p_u(1);
	}
	return m;
}
const lzh::Mat lzh::BSpline::BPoint(int number)const
{
	int n = P.rows() - 1;
	mat_t start = _T(k) / (n + k + 1);
	mat_t end = _T(n + 1) / (n + k + 1);
	Mat t = linspace(start, end, number);
	return BPiont(t);
}
const lzh::Mat lzh::BSpline::operator()(mat_t t) const
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(P);
	CHECK_MATRIX_(nodevector);
#endif //MAT_DEBUG
	return P * BaseFunction(t);
}
lzh::BSpline lzh::BSpline::fitBSpline(const Mat &P, int n, int k)
{
	int r = P.rows() - 1;
	n = n + 2;
	mat_t start = _T(k) / (n + k + 1),
		end = _T(n + 1) / (n + k + 1);
	Mat t = linspace(start, end, P.rows());
	Mat Node = linspace(0, 1, n + k + 1);
	//t.show();
	BSpline bs(BSpline::UNIFORM, k);
	bs.NodeVector(Node);
	Mat N(n, r + 1);
	//n = n - 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < P.rows(); j++) {
			N(j*n + i) = bs.BF(i, k, t(j));
		}
	}
	Mat m = LeastSquare(N, P).t();
	Mat d = Block(m, 1, m.rows() - 2, 0, 1);
	Mat p(d, Block(d, 0, k, 0, 1), ROW);
	bs.setCtrlPoint(p.t());
	bs.NodeVector();
	return bs;
}

/****************************************************************************
激活函数/损失函数
*****************************************************************************/
const lzh::Mat lzh::Softmax(const Mat &y)
{
	Mat out;
	if (y.dims() == 3) {
		out = y.clone();
		for (int i = 0; i < y.rows(); ++i)
		{
			for (int j = 0; j < y.cols(); ++j)
			{
				Mat y_ = out(i, j, CHANNEL);
				y_ -= Max(y_);
				Mat y_exp = mExp(y_);
				lzh::mat_t y_sum = y_exp.sum();
				(y_exp / y_sum).copyTo(y_);
			}
		}
	}
	else
	{
		Mat y_ = y.clone();
		y_ -= Max(y_);
		Mat y_exp = mExp(y_);
		lzh::mat_t y_sum = y_exp.sum();
		out = y_exp / y_sum;
	}
	return out;
}
const lzh::Mat lzh::L2(const Mat & y, const Mat & y0)
{
	return (y - y0).pow(2);
}
const lzh::Mat lzh::Quadratic(const Mat &y, const Mat &y0)
{
	return 0.5 * mPow(y - y0, 2);
}
const lzh::Mat lzh::CrossEntropy(const Mat &y, const Mat &y0)
{
	return -Mult(y, mLog(y0));
}
const lzh::Mat lzh::SoftmaxCrossEntropy(const Mat & y, const Mat & y0)
{
	return CrossEntropy(y, Softmax(y0));
}
const lzh::Mat lzh::Sigmoid(const Mat &x)
{
	return 1.0 / (1.0 + mExp(-x));
}
const lzh::Mat lzh::Tanh(const Mat &x)
{
	return 2.0 * Sigmoid(2 * x) - 1.0;
}
const lzh::Mat lzh::ReLU(const Mat &x)
{
	return mMax(0, x);
}
const lzh::Mat lzh::ELU(const Mat & x)
{
	Mat x1(x.size3());
	lzh::mat_t *p = x1;
	const lzh::mat_t *mat = x;
	for (int i = 0; i < x.len(); ++i) {
		if (*mat <= 0)
			*p = ELU_alpha * (exp(*mat) - 1);
		else
			*p = *mat;
		p++;
		mat++;
	}
	return x1;
}
const lzh::Mat lzh::SELU(const Mat & x)
{
	return SELU_scale * ELU(x);
}
const lzh::Mat lzh::LReLU(const Mat & x)
{
	Mat x1(x.size3());
	lzh::mat_t *p = x1;
	const lzh::mat_t *mat = x;
	for (int i = 0; i < x.len(); ++i) {
		if (*mat <= 0)
			*p = *mat*LReLU_alpha;
		else
			*p = *mat;
		p++;
		mat++;
	}
	return x1;
}

/****************************************************************************
其他工具
*****************************************************************************/
std::string lzh::tools::createfile(std::string filename)
{
	return filename.substr(filename.rfind('\\') + 1);
}
std::string lzh::tools::createtype(std::string filename)
{
	return filename.substr(filename.rfind('.') + 1);
}
void lzh::tools::pause()
{
	fprintf(stderr, "waitting press enter key...\n");
	while (getchar() != '\n');
}
void lzh::tools::check(int w, int h, int c)
{
	if (w <= 0 || h <= 0 || c <= 0)
		THROW_INFO(ERR_INFO_VALUE);
}
const lzh::Mat lzh::tools::Vec2Mat(std::vector<lzh::mat_t> &p)
{
	if (p.empty())return Mat();
	Mat m(int(p.size()));
	for (int iter = 0; iter != int(p.size()); ++iter) {
		m(iter) = p[iter];
	}
	return m;
}
const lzh::Mat lzh::tools::Vec2Mat(std::vector<std::vector<lzh::mat_t>> &ps)
{
	if (ps.empty())return Mat();
	int size = 0;
	for (int i = 0; i < int(ps.size() - 1); ++i) {
		for (int j = i + 1; j < int(ps.size()); ++j) {
			if (ps[i].size() != ps[j].size())
				THROW_INFO(ERR_INFO_SIZE);
		}
	}
	int hei = int(ps.size());
	int wid = int(ps[0].size());
	Mat m(wid, hei);
	for (int i = 0; i < hei; ++i) {
		for (int j = 0; j < wid; ++j) {
			m(i, j) = ps[i][j];
		}
	}
	return m;
}
std::vector<lzh::mat_t> lzh::tools::Mat2Vec(const lzh::Mat & m)
{
	if (m.empty())return std::vector<lzh::mat_t>();
	std::vector<lzh::mat_t> p(m.size());
	for (int iter = 0; iter != m.size(); ++iter) {
		p[iter] = m(iter);
	}
	return p;
}
std::vector<std::vector<lzh::mat_t>> lzh::tools::Mat2Vecs(const lzh::Mat & m)
{
	if (m.empty())return std::vector<std::vector<lzh::mat_t>>();
	std::vector<std::vector<lzh::mat_t>> ps;
	for (int h = 0; h != m.rows(); ++h) {
		std::vector<lzh::mat_t> p;
		for (int w = 0; w != m.cols(); ++w) {
			p.push_back(m(h, w));
		}
		ps.push_back(p);
	}
	return ps;
}
std::vector<std::string> lzh::tools::strsplit(std::string &str, char ch)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(str);
#endif //MAT_DEBUG
	size_t idx = 0;
	size_t offset = 0;
	std::vector<std::string> spl;
	while (true) {
		offset = str.find(ch, idx);
		if (offset == std::string::npos)break;
		spl.push_back(str.substr(idx, offset - idx));
		idx = offset + 1;
	}
	str = str.substr(idx);
	if (str != "")
		spl.push_back(str);
	return spl;
}
std::vector<lzh::mat_t> lzh::tools::str2mat_t(const std::vector<std::string> &str)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(str);
#endif //MAT_DEBUG
	std::vector<mat_t> v;
	for (const std::string &s : str)
		if(sizeof(mat_t) == sizeof(float))
			v.push_back(_T(stof(s)));
		else if (sizeof(mat_t) == sizeof(double))
			v.push_back(_T(stod(s)));
		else if (sizeof(mat_t) == sizeof(long double))
			v.push_back(_T(stold(s)));
	return v;
}
const lzh::Mat lzh::tools::readcsv(std::string filename)
{
#ifdef MAT_DEBUG
	CHECK_MATRIX_(filename);
#endif //MAT_DEBUG
	std::ifstream in(filename);
	if (!in.is_open())
		THROW_INFO(ERR_INFO_FILE);
	std::vector<std::vector<lzh::mat_t>> vec;
	int len = -1;
	std::string str;
	while (std::getline(in, str))
	{
		std::vector<lzh::mat_t> data = str2mat_t(strsplit(str, ','));
		if (len == -1)
			len = (int)data.size();
		else
			if (len != data.size())
				THROW_INFO(ERR_INFO_SIZE);
		vec.push_back(data);
	}
	in.close();
	if (vec.empty())return Mat();
	return Vec2Mat(vec);
}
std::string lzh::tools::Enum2String(PrintType type)
{
	std::string str;
	switch (type)
	{
	case lzh::FIXED:
		str = "FIXED INFO:以小数形式格式化输出";
		break;
	case lzh::SCIENTIFIC:
		str = "SCIENTIFIC INFO:以科学记数法形式格式化输出";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(BorderTypes type)
{
	std::string str;
	switch (type)
	{
	case lzh::BORDER_CONSTANT:		
		str = "BORDER_CONSTANT //!< `iiiiii|abcdefgh|iiiiiii`  with some specified `i`";
		break;
	case lzh::BORDER_REPLICATE:
		str = "BORDER_CONSTANT INFO://!< `iiiiii|abcdefgh|iiiiiii`  with some specified `i`";
		break;
	case lzh::BORDER_REFLECT:
		str = "BORDER_REPLICATE INFO://!< `aaaaaa|abcdefgh|hhhhhhh`";
		break;
	case lzh::BORDER_WRAP:
		str = "BORDER_REFLECT INFO://!< `fedcba|abcdefgh|hgfedcb`";
		break;
	case lzh::BORDER_REFLECT_101:
		str = "BORDER_REFLECT_101 INFO://!< `gfedcb|abcdefgh|gfedcba`";
		break;
	case lzh::BORDER_TRANSPARENT:
		str = "BORDER_TRANSPARENT INFO://!< `uvwxyz|abcdefgh|ijklmno`";
		break;
	case lzh::BORDER_ISOLATED:
		str = "BORDER_ISOLATED INFO://!< do not look outside of ROI";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(MatErrorInfo type)
{
	return errinfo[type];
}
std::string lzh::tools::Enum2String(EQUATION type)
{
	std::string str;
	switch (type)
	{
	case lzh::SPECIAL_SOLUTION:
		str = "SPECIAL_SOLUTION INFO:方程有特解";
		break;
	case lzh::GENERAL_SOLUTION:
		str = "GENERAL_SOLUTION INFO:方程有通解";
		break;
	case lzh::NO_SOLUTION:
		str = "NO_SOLUTION INFO:方程无解";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(ORDER type)
{
	std::string str;
	switch (type)
	{
	case lzh::MIN_TO_MAX:	
		str = "MIN_TO_MAX INFO:从小到大";
		break;
	case lzh::MAX_TO_MIN:
		str = "MAX_TO_MIN INFO:从大到小";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(X_Y_Z type)
{
	std::string str;
	switch (type)
	{
	case lzh::ROW:
		str = "ROW INFO:行方向";
		break;
	case lzh::COL:
		str = "COL INFO:列方向";
		break;
	case lzh::CHANNEL:
		str = "CHANNEL INFO:通道方向";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(Dire type)
{
	std::string str;
	switch (type)
	{
	case lzh::LEFT:
		str = "LEFT INFO:左";
		break;
	case lzh::RIGHT:
		str = "RIGHT INFO:右";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(ReductionMothed type)
{
	std::string str;
	switch (type)
	{
	case lzh::EqualIntervalSampling:
		str = "EqualIntervalSampling INFO:等间隔采样";
		break;
	case lzh::LocalMean:
		str = "LocalMean INFO:局部均值";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(RotateAngle type)
{
	std::string str;
	switch (type)
	{
	case lzh::ROTATE_90_ANGLE:
		str = "ROTATE_90_ANGLE INFO:顺时针旋转90度";
		break;
	case lzh::ROTATE_180_ANGLE:
		str = "ROTATE_180_ANGLE INFO:顺时针旋转180度";
		break;
	case lzh::ROTATE_270_ANGLE:
		str = "ROTATE_270_ANGLE INFO:顺时针旋转270度";
		break;
	default:
		break;
	}
	return str;
}
std::string lzh::tools::Enum2String(BSpline::BSplineType type)
{
	std::string str;
	switch (type)
	{
	case BSpline::UNIFORM:
		str = "UNIFORM INFO:均匀B样条曲线90度"; break;
	case BSpline::QUASI_UNIFORM:
		str = "QUASI_UNIFORM INFO:准均匀B样条曲线90度"; break;
	default:
		break;
	}
	return std::string();
}
/****************************************************************************
计时器
*****************************************************************************/
#if defined(__linux__)
#include <sys/time.h>
#include <unistd.h>
static struct timeval t1, t2;
void lzh::tools::StartCounter()
{
	gettimeofday(&t1, NULL);
}
lzh::mat_t lzh::tools::EndCounter()
{
	gettimeofday(&t2, NULL);
	return _T((t2.tv_sec - t1.tv_sec) * 1000.0 + t2.tv_usec - t1.tv_usec)
}
void lzh::tools::Wait(uint ms)
{
	sleep(ms);
}
#elif defined(_WIN32)
#include <windows.h>  
#include <io.h>
#include <direct.h>  
static LARGE_INTEGER cpuFreq;
static LARGE_INTEGER startTime;
static LARGE_INTEGER endTime;
void lzh::tools::Frequency()
{
	QueryPerformanceFrequency(&cpuFreq);
}
void lzh::tools::StartCounter()
{
	QueryPerformanceCounter(&startTime);
}
lzh::mat_t lzh::tools::EndCounter()
{
	QueryPerformanceCounter(&endTime);
	return _T(((endTime.QuadPart - startTime.QuadPart) * 1000.0) / cpuFreq.QuadPart);
}
/**
@brief getFiles 得到路径下所有文件的路径
@param path 文件夹路径
@param files 保存path下的所有文件路径
*/
void lzh::tools::getFiles(std::string path, std::vector<std::string>& files)
{
	//文件句柄  
	intptr_t hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	std::string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之  
			//如果不是,加入列表  
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}
std::string lzh::tools::show_path()
{
	char buffer[260];
	_getcwd(buffer, 260);
	return std::string(buffer);
}
void lzh::tools::Wait(uint ms)
{
	Sleep(ms);
}
#endif
