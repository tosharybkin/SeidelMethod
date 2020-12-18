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
	system("title Рыбкин Антон: задача Дирихле для уравнения Пуассона");

	int		Nmax = 10000;		//Максимальное число итераций
	int		S = 0;				//Счетчик итераций
	double	eps = 0.0000001;	//Параметр требуемой точности
	double	epsMax = 0;			//Достигнутая точность
	int		n = 4, m = 4;		//Размерность сетки
	double** V = NULL;			//Искомый вектор 
	double** F = NULL;			//f(x, y) из дифференециального уравнения в узлах сетки
	double	a, b, c, d;			//Границы области определния уравнения
	int Exit = 1, Show = 0;

	//Границы по x и y
	a = 0;
	b = 2;
	c = 0;
	d = 1;

	Nmax = 0;

	system("cls");
	cout << "Размерность сетки 4 на 4";

	while (Nmax <= 0)
	{
		cout << endl << "Введите максимальное количество итераций: ";
		cin >> Nmax;
	}
	system("cls");

	V = MemoryAllocator(n + 1, m + 1);
	FillStartSolution(V, n, m, a, b, c, d);
	ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
	cout << "Итоговое численное решение имеет вид " << endl;
	ShowSolution(V, n, m);

	//Справка
	cout << endl << "---------------------------------------------" << endl;
	cout << "Размерность сетки: (" << n << ", " << m << ")" << endl;
	cout << "Область определения X: [" << a << "; " << b << "]" << endl;
	cout << "Область определения Y: [" << c << "; " << d << "]" << endl;
	cout << "Шаг сетки по оси Ox: h = " << (b - a) / n << endl;
	cout << "Шаг сетки по оси Oy: k = " << (d - c) / m << endl;
	cout << "Заданная точность: " << eps << endl;
	cout << "Достигнутая точность: " << epsMax << endl;
	cout << "Проведено итераций: " << S << endl;
	cout << "Невязка решения: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;
	cout << "Погрешность решения: " << CheckNumSolution(V, n, m, a, b, c, d) << endl;
	cout << "---------------------------------------------" << endl << endl;

	MemoryCleaner(V, n);
	
	cout << endl;
	return 0;
}