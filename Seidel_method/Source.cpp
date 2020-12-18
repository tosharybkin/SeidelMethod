#include "NumberSolution.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

using namespace std;


//U"xx + U"yy = -f(x, y)
//a <= x <= b, c <= y <= d
//U(a, y) = M1(y), U(b, y) = M2(y)
//U(x, c) = M3(x), U(x, d) = M3(x)

int main(void)
{
	setlocale(LC_ALL, "rus");
	system("title ������ �����: ������ ������� ��� ��������� ��������");

	int		Nmax = 10000;		//������������ ����� ��������
	int		S = 0;				//������� ��������
	double	eps = 0.0000001;	//�������� ��������� ��������
	double	epsMax = 0;			//����������� ��������
	int		n = 4, m = 4;		//����������� �����
	double** V = NULL;			//������� ������ 
	double** F = NULL;			//f(x, y) �� ������������������ ��������� � ����� �����
	double	a, b, c, d;			//������� ������� ���������� ���������
	int Exit = 1, Show = 0;

	//������� �� x � y
	a = 0;
	b = 2;
	c = 0;
	d = 1;

	Nmax = 0;

	system("cls");
	cout << "����������� ����� 4 �� 4";

	while (Nmax <= 0)
	{
		cout << endl << "������� ������������ ���������� ��������: ";
		cin >> Nmax;
	}
	system("cls");

	V = MemoryAllocator(n + 1, m + 1);
	FillStartSolution(V, n, m, a, b, c, d);
	ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
	cout << "�������� ��������� ������� ����� ��� " << endl;
	ShowSolution(V, n, m);

	//�������
	cout << endl << "---------------------------------------------" << endl;
	cout << "����������� �����: (" << n << ", " << m << ")" << endl;
	cout << "������� ����������� X: [" << a << "; " << b << "]" << endl;
	cout << "������� ����������� Y: [" << c << "; " << d << "]" << endl;
	cout << "��� ����� �� ��� Ox: h = " << (b - a) / n << endl;
	cout << "��� ����� �� ��� Oy: k = " << (d - c) / m << endl;
	cout << "�������� ��������: " << eps << endl;
	cout << "����������� ��������: " << epsMax << endl;
	cout << "��������� ��������: " << S << endl;
	cout << "������� �������: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;
	cout << "����������� �������: " << CheckNumSolution(V, n, m, a, b, c, d) << endl;
	cout << "---------------------------------------------" << endl << endl;

	MemoryCleaner(V, n);
	
	cout << endl;
	return 0;
}