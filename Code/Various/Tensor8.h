/*
 * Tensor8.h
 *
 *  Created on: Apr 14, 2016
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_TENSOR8_H_
#define CODE_VARIOUS_TENSOR8_H_

using namespace std;
using namespace dealii;

template<int dim>
void print(std::ostream &o, const Tensor<4, dim, double> &t) {

	int w = 14;
	o << "--------------------------------------------" << endl;
	o << setw(w) << t[0][0][0][0] << setw(w) << t[0][0][1][1] << setw(w)
			<< t[0][0][2][2] << setw(w) << t[0][0][1][2] << setw(w)
			<< t[0][0][0][2] << setw(w) << t[0][0][0][1] << endl;

	o << setw(w) << t[1][1][0][0] << setw(w) << t[1][1][1][1] << setw(w)
			<< t[1][1][2][2] << setw(w) << t[1][1][1][2] << setw(w)
			<< t[1][1][0][2] << setw(w) << t[1][1][0][1] << endl;

	o << setw(w) << t[2][2][0][0] << setw(w) << t[2][2][1][1] << setw(w)
			<< t[2][2][2][2] << setw(w) << t[2][2][1][2] << setw(w)
			<< t[2][2][0][2] << setw(w) << t[2][2][0][1] << endl;

	o << setw(w) << t[1][2][0][0] << setw(w) << t[1][2][1][1] << setw(w)
			<< t[1][2][2][2] << setw(w) << t[1][2][1][2] << setw(w)
			<< t[1][2][0][2] << setw(w) << t[1][2][0][1] << endl;

	o << setw(w) << t[0][2][0][0] << setw(w) << t[0][2][1][1] << setw(w)
			<< t[0][2][2][2] << setw(w) << t[0][2][1][2] << setw(w)
			<< t[0][2][0][2] << setw(w) << t[2][0][0][1] << endl;

	o << setw(w) << t[0][1][0][0] << setw(w) << t[0][1][1][1] << setw(w)
			<< t[0][1][2][2] << setw(w) << t[0][1][1][2] << setw(w)
			<< t[0][1][0][2] << setw(w) << t[0][1][0][1] << endl;
	o << "--------------------------------------------" << endl;
}

template<int dim>
void print(LogStream &o, const Tensor<4, dim, double> &t) {
	print(o.get_console(), t);
}

template<int dim>
void print(std::ostream &o, const SymmetricTensor<4, dim, double> &t) {
	int w = 14;
	o << "--------------------------------------------" << endl;
	o << setw(w) << t[0][0][0][0] << setw(w) << t[0][0][1][1] << setw(w)
			<< t[0][0][2][2] << setw(w) << t[0][0][1][2] << setw(w)
			<< t[0][0][0][2] << setw(w) << t[0][0][0][1] << endl;

	o << setw(w) << t[1][1][0][0] << setw(w) << t[1][1][1][1] << setw(w)
			<< t[1][1][2][2] << setw(w) << t[1][1][1][2] << setw(w)
			<< t[1][1][0][2] << setw(w) << t[1][1][0][1] << endl;

	o << setw(w) << t[2][2][0][0] << setw(w) << t[2][2][1][1] << setw(w)
			<< t[2][2][2][2] << setw(w) << t[2][2][1][2] << setw(w)
			<< t[2][2][0][2] << setw(w) << t[2][2][0][1] << endl;

	o << setw(w) << t[1][2][0][0] << setw(w) << t[1][2][1][1] << setw(w)
			<< t[1][2][2][2] << setw(w) << t[1][2][1][2] << setw(w)
			<< t[1][2][0][2] << setw(w) << t[1][2][0][1] << endl;

	o << setw(w) << t[0][2][0][0] << setw(w) << t[0][2][1][1] << setw(w)
			<< t[0][2][2][2] << setw(w) << t[0][2][1][2] << setw(w)
			<< t[0][2][0][2] << setw(w) << t[2][0][0][1] << endl;

	o << setw(w) << t[0][1][0][0] << setw(w) << t[0][1][1][1] << setw(w)
			<< t[0][1][2][2] << setw(w) << t[0][1][1][2] << setw(w)
			<< t[0][1][0][2] << setw(w) << t[0][1][0][1] << endl;
	o << "--------------------------------------------" << endl;
}

template<int dim>
void print(LogStream &o, const SymmetricTensor<4, dim, double> &t) {
	print(o.get_console(), t);
}

template<int dim>
inline Tensor<4, dim, double> contract4(const Tensor<8, dim, double> &t1,
		const Tensor<4, dim, double> &t2) {
	Tensor<4, dim, double> tmp;
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = 0; l < dim; ++l)
					for (unsigned int m = 0; m < dim; ++m)
						for (unsigned int n = 0; n < dim; ++n)
							for (unsigned int o = 0; o < dim; ++o)
								for (unsigned int p = 0; p < dim; ++p)
									tmp[i][j][k][l] +=
											t1[i][j][k][l][m][n][o][p]
													* t2[m][n][o][p];
	return tmp;
}

template<int dim>
inline Tensor<4, dim, double> contract2(const Tensor<4, dim, double> &t1,
		const Tensor<4, dim, double> &t2) {
	Tensor<4, dim, double> tmp;
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = 0; l < dim; ++l)
					for (unsigned int o = 0; o < dim; ++o)
						for (unsigned int p = 0; p < dim; ++p)
							tmp[i][j][o][p] += t1[i][j][k][l] * t2[k][l][o][p];
	return tmp;
}

template<int dim>
inline Tensor<2, dim, double> unit_identity2() {
	Tensor<2, dim, double> tmp;

	for (unsigned int d = 0; d < dim; ++d)
		tmp[d][d] = 1;

	return tmp;
}

template<int dim>
inline Tensor<4, dim, double> unit_identity_sym() {
	Tensor<4, dim, double> tmp;

	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = 0; l < dim; ++l)
					tmp[i][j][k][l] = 0.5
							* ((i == k ? 1 : 0) * (j == l ? 1 : 0)
									+ (i == l ? 1 : 0) * (j == k ? 1 : 0));

	return tmp;
}

template<int dim, typename Number>
inline Tensor<2, dim, Number> outer_product(const Tensor<1, dim, Number> &t1,
		const Tensor<1, dim, Number> &t2) {
	Tensor<2, dim, Number> tmp;

	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			tmp[i][j] = t1[i] * t2[j];

	return tmp;
}

template<int dim, typename Number>
inline Tensor<4, dim, Number> outer_product(const Tensor<2, dim, Number> &t1,
		const Tensor<2, dim, Number> &t2) {
	Tensor<4, dim, Number> tmp;

	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = 0; l < dim; ++l)
					tmp[i][j][k][l] = t1[i][j] * t2[k][l];

	return tmp;
}

template<int dim, typename Number>
inline Tensor<8, dim, Number> outer_product(const Tensor<4, dim, Number> &t1,
		const Tensor<4, dim, Number> &t2) {
	Tensor<8, dim, Number> tmp;

	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = 0; l < dim; ++l)
					for (unsigned int m = 0; m < dim; ++m)
						for (unsigned int n = 0; n < dim; ++n)
							for (unsigned int o = 0; o < dim; ++o)
								for (unsigned int p = 0; p < dim; ++p)
									tmp[i][j][k][l][m][n][o][p] = t1[i][j][k][l]
											* t2[m][n][o][p];

	return tmp;
}

template<int dim>
Tensor<8, dim> calculate_transversal_isotropic_projector() {
	Tensor<2, dim> identity2 = unit_identity2<dim>();
//	Tensor<4, dim> identity_sym4 = unit_identity_sym<dim>();
//	Tensor<4, dim> J = outer_product(identity2, identity2) / 3;
	Tensor<1, dim> a( { 1, 0, 0 });
	Tensor<1, dim> b( { 0, 1, 0 });
	Tensor<1, dim> c( { 0, 0, 1 });
	Tensor<2, dim> P = outer_product(c, c);
	Tensor<2, dim> Q = (identity2 - P) / std::sqrt(2);
	Tensor<2, dim> U = (outer_product(a, c) + outer_product(c, a))
			/ std::sqrt(2);
	Tensor<2, dim> V = (outer_product(b, c) + outer_product(c, b))
			/ std::sqrt(2);
	Tensor<2, dim> W = (outer_product(a, b) + outer_product(b, a))
			/ std::sqrt(2);
	Tensor<2, dim> Z = (outer_product(a, a) - outer_product(b, b))
			/ std::sqrt(2);

	Tensor<4, dim> E1 = outer_product(P, P);
	Tensor<4, dim> E2 = outer_product(Q, Q);
	Tensor<4, dim> E3 = outer_product(P, Q);
	Tensor<4, dim> E4 = outer_product(Q, P);
	Tensor<4, dim> F = outer_product(W, W) + outer_product(Z, Z);
	Tensor<4, dim> G = outer_product(U, U) + outer_product(V, V);

	Tensor<8, dim> tmp = outer_product(E1, E1) + outer_product(E2, E2)
			+ 0.5 * outer_product(E3 + E4, E4 + E3) + 0.5 * outer_product(F, F)
			+ 0.5 * outer_product(G, G);

//	cout << "-------E1-------";
//	print(cout, E1);
//	cout << "-------E2-------";
//	print(cout, E2);
//	cout << "-------E3-------";
//	print(cout, E3);
//	cout << "-------E4-------";
//	print(cout, E4);
//	cout << "-------F--------";
//	print(cout, F);
//	cout << "-------G--------";
//	print(cout, G);
//
//	exit(0);
	return tmp;
}

#endif /* CODE_VARIOUS_TENSOR8_H_ */
