/**
 * -------------------------------------- Metody Numeryczne ----------------------------------------
 * ------------------------------------------ Zadanie 13 -------------------------------------------
 * ------------------------------- autor: Grzegorz Poreda, 21.01.2021 ------------------------------
 * 
 * Stosujac metode Laguerre'a wraz ze strategia obnizania stopnia wielomianu i wygladzania, znajdz
 * wszystkie rozwiazania rownan.
 * 
 * -------------------------------------------------------------------------------------------------
 */

/*

Liczby zespolone wypisywane sa w formacie ( Re , Im )

a) Miejsca zerowe:
z[0] == ( 0.666671 , 2.35823e-006 )
z[1] == ( 0.666664 , -5.81474e-008 )
z[2] == ( 0.666668 , -2.83893e-006 )
z[3] == ( -4.61717e-017 , 1.41421 )
z[4] == ( 0.333333 , 9.86076e-032 )
z[5] == ( -0.333333 , 0 )
z[6] == ( -4.3832e-017 , -1.41421 )

b) Miejsca zerowe:
z[0] == ( 6.6268e-008 , 1 )
z[1] == ( -1.61974e-017 , 1.41421 )
z[2] == ( -1.97319e-008 , 1 )
z[3] == ( 1.41421 , -3.67342e-040 )
z[4] == ( -0.5 , 0.866025 )
z[5] == ( 4.24132e-008 , -1 )
z[6] == ( 2.57067e-008 , -1 )
z[7] == ( 4.354e-016 , -1.41421 )
z[8] == ( -0.5 , -0.866025 )
z[9] == ( -1.41421 , -3.76158e-037 )

c) Miejsca zerowe:
z[0] == ( 0.951057 , 0.309017 )
z[1] == ( 0.587785 , -0.809017 )
z[2] == ( -0.951057 , 0.309017 )
z[3] == ( -0.587785 , -0.809017 )

*/


#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

#define PRECISION 0.000000000001


class Complex
{
  private:
    double real;
    double imaginary;

  public:
    explicit Complex( double = 0.0, double = 0.0 );

	void setReal(double Re);
	double getReal() const;

	void setImaginary(double Im);
	double getImaginary() const;

	double getAbsoluteValue();
	Complex getSquareRoot();

    Complex operator+( const Complex & ) const;
    Complex operator-( const Complex & ) const;

    Complex operator*( const Complex & ) const;
	Complex operator*( const double & ) const;

	Complex operator/( const Complex & ) const;

    bool operator==( const Complex & ) const;
    bool operator!=( const Complex & ) const;

    friend std::ostream &operator<<( std::ostream &, const Complex & );
    friend std::istream &operator>>( std::istream &, Complex & );
};


double countError(Complex newI, Complex oldI)
{
    return (newI - oldI).getAbsoluteValue();
}

std::vector<Complex> getDerivativeCoeff(const std::vector<Complex>& coeff)
{
	size_t n = coeff.size() - 1;

	if(n == 0)
	{
		std::vector<Complex> derivCoeff(1);
		derivCoeff[0] = Complex(0.0,0.0);
		return derivCoeff;
	}
		
	std::vector<Complex> derivCoeff(n);
	for(size_t i = 0; i < n; i++)
		derivCoeff[i] = coeff[i+1] * (i+1);
	return derivCoeff;
}


Complex getValue(const std::vector<Complex>& coeff, Complex x)
{
	size_t n = coeff.size();
	Complex value = coeff[n-1];
	for(size_t i = n-2; i!=-1 ; i--)
		value = value*x + coeff[i];
	return value;
}


Complex laguerre(const std::vector<Complex>& coeff, Complex guess)
{
	double n = coeff.size() - 1;
	std::vector<Complex> coeffDeriv1 = getDerivativeCoeff(coeff);
	std::vector<Complex> coeffDeriv2 = getDerivativeCoeff(coeffDeriv1);
	Complex oldRoot;
	Complex newRoot = guess;

	do
	{
		oldRoot = newRoot;
		Complex nPn_zi = getValue(coeff, oldRoot) * n;
		Complex derPn_zi = getValue(coeffDeriv1, oldRoot);
		Complex squareRoot = ( ( ((derPn_zi*derPn_zi)*(n-1))-(nPn_zi*getValue(coeffDeriv2,oldRoot)) ) * (n-1) ).getSquareRoot();
		Complex denominator1 = derPn_zi + squareRoot;
		Complex denominator2 = derPn_zi - squareRoot;
		if(denominator1.getAbsoluteValue() >= denominator2.getAbsoluteValue())
			newRoot = oldRoot - (nPn_zi / denominator1);
		else
			newRoot = oldRoot - (nPn_zi / denominator2);
		
	} while( countError( getValue(coeff,newRoot),getValue(coeff,oldRoot) ) > PRECISION );

	return newRoot;
}


std::vector<Complex> deflate(const std::vector<Complex>& A, Complex z0)
{
	size_t n = A.size()-1;
	std::vector<std::vector<Complex>> Z(n-1);
	std::vector<Complex> B(n);

	B[n-1] = A[n];
	for(size_t i = n-2; i!=-1 ; i--)
	{
		B[i] = A[i+1] + B[i+1]*z0;
	}
	return B;
}


std::vector<Complex> getRoots(const std::vector<Complex>& A)
{
	size_t n = A.size() - 1;
	std::vector<Complex> roots(n);
	
	roots[0] = laguerre(A,Complex(1.0,1.0));
	std::vector<Complex> B = deflate(A,roots[0]);

	for(size_t i = 1; i < n; i++)
	{
		roots[i] = laguerre(B,Complex(1.0,1.0));
		roots[i] = laguerre(A,roots[i]);
		B = deflate(B,roots[i]);
	}

	return roots;
}


void printRoots(const std::vector<Complex>& x)
{
    for(int i=0; i<x.size(); i++)
        std::cout << "z[" << i << "] == " << x[i] << std::endl;
    std::cout << std::endl;
}



int main()
{
	std::cout << "Liczby zespolone wypisywane sa w formacie ( Re , Im )" << std::endl << std::endl;

	std::vector<Complex> A(8);
	A[0] = Complex(  16.0  , 0 );
	A[1] = Complex( -72.0  , 0 );
	A[2] = Complex( -28.0  , 0 );
	A[3] = Complex(  558.0 , 0 );
	A[4] = Complex( -990.0 , 0 );
	A[5] = Complex(  783.0 , 0 );
	A[6] = Complex( -486.0 , 0 );
	A[7] = Complex(  243.0 , 0 );
	std::vector<Complex> aRoots = getRoots(A);
	std::cout << "a) Miejsca zerowe:" << std::endl;
	printRoots(aRoots);

	std::vector<Complex> B(11);
	B[0]  = Complex( -4.0  , 0 );
	B[1]  = Complex( -4.0  , 0 );
	B[2]  = Complex( -12.0 , 0 );
	B[3]  = Complex( -8.0  , 0 );
	B[4]  = Complex( -11.0 , 0 );
	B[5]  = Complex( -3.0  , 0 );
	B[6]  = Complex( -1.0  , 0 );
	B[7]  = Complex(  2.0  , 0 );
	B[8]  = Complex(  3.0  , 0 );
	B[9]  = Complex(  1.0  , 0 );
	B[10] = Complex(  1.0  , 0 );
	std::vector<Complex> bRoots = getRoots(B);
	std::cout << "b) Miejsca zerowe:" << std::endl;
	printRoots(bRoots);

	std::vector<Complex> C(5);
	C[0] = Complex(  1.0 ,  0   );
	C[1] = Complex(  0.0 , -1.0 );
	C[2] = Complex( -1.0 ,  0   );
	C[3] = Complex(  0.0 ,  1.0 );
	C[4] = Complex(  1.0 ,  0   );
	std::vector<Complex> cRoots = getRoots(C);
	std::cout << "c) Miejsca zerowe:" << std::endl;
	printRoots(cRoots);

	return 0;
}



std::ostream &operator<<( std::ostream &output, const Complex &com )
{
    output << "( " << com.real << " , " << com.imaginary << " )";
    return output;
}

std::istream &operator>>( std::istream &input, Complex &com )
{
    input.ignore(2);
    input >> com.real;
    input.ignore(3);
    input >> com.imaginary;
    input.ignore(2);
    return input;
}

Complex::Complex( double realPart, double imaginaryPart )
    : real( realPart ),
      imaginary( imaginaryPart )
{

}

void Complex::setReal(double Re)
{
	real = Re;
}

double Complex::getReal() const
{
	return real;
}

void Complex::setImaginary(double Im)
{
    imaginary = Im;
}

double Complex::getImaginary() const
{
    return imaginary;
}

double Complex::getAbsoluteValue()
{
	return sqrt( (real*real) + (imaginary*imaginary) );
}

Complex Complex::getSquareRoot()
{
	double r = sqrt(getAbsoluteValue());
    double theta = atan2(imaginary, real) / 2.0;
    return Complex(r*cos(theta), r*sin(theta));
}

Complex Complex::operator+( const Complex &operand2 ) const
{
    return Complex( real + operand2.real,
                   imaginary + operand2.imaginary );
}

Complex Complex::operator-( const Complex &operand2 ) const
{
    return Complex( real - operand2.real,
                   imaginary - operand2.imaginary );
}

Complex Complex::operator*( const Complex &operand2 ) const
{
    return Complex( real * operand2.real - imaginary * operand2.imaginary,
                   real * operand2.imaginary + imaginary * operand2.real );
}

Complex Complex::operator*( const double &num ) const
{
    return Complex( real * num, imaginary * num );
}

Complex Complex::operator/( const Complex &com ) const
{
	double Re = (getReal()*com.getReal() + getImaginary()*com.getImaginary())
    	/ (com.getReal()*com.getReal() + com.getImaginary()*com.getImaginary());
    double Im = (getImaginary()*com.getReal() - getReal()*com.getImaginary())
    	/ (com.getReal()*com.getReal() + com.getImaginary()*com.getImaginary());
	return Complex(Re,Im);
}

bool Complex::operator==( const Complex &operand2 ) const
{
    if( real == operand2.real && imaginary == operand2.imaginary ) return true;
    return false;
}

bool Complex::operator!=( const Complex &operand2 ) const
{
    if( real != operand2.real || imaginary != operand2.imaginary ) return true;
    return false;
}
