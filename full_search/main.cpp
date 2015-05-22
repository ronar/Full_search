#include <windows.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include "math.h"

using namespace std;

// Отбор признаков с использованием метода полного перебора

// Шаг на спуске
void gradient_step (double** &points, double &theta0, double &theta1, float a, int m) {
	double grad0 = 0, grad1 = 0;
	double x, y;

	for (int i = 0; i < m; i++)
	{
		x = points[i][0];
		y = points[i][1];

		grad0 += -(2.0 / m) * (y - ((theta1 * x) + theta0));
		grad1 += -(2.0 / m) * x * (y - ((theta1 * x) + theta0));
	}

	theta0 = theta0 - (a * grad0);
	theta1 = theta1 - (a * grad1);
}

// m - количество точек
double full_search(double** &dataset, int m, int num_iterations, double &theta0, double &theta1, double learning_rate) {
	double err = 0;
			
	// Градиентный метод
	for(int i = 0; i < num_iterations; i++) {
		gradient_step(dataset, theta0, theta1, learning_rate, m);
	}

	float sum = 0;

	// Находим ошибку на данном множестве точек
	for (int s = 0; s < m; s++)
	{
		sum += pow(((theta0 + theta1 * dataset[s][0]) - dataset[s][1]), 2);
	}

	err = sum / m; // Ошибка

	return err;
}

int** get_first_last_test(int num_data_items, int num_folds)
{
	// return[fold][firstIndex][lastIndex] для тестовых данных k-fold cross validation
	int interval = num_data_items / num_folds; 
	int** result = new int *[num_folds]; // Парные индексы для каждого разбиения
	for (int i = 0; i < num_folds; ++i)
		result[i] = new int[2];

	for (int k = 0; k < num_folds; ++k)
	{
		int first = k * interval;
		int last = (k+1) * interval - 1;
		result[k][0] = first;
		result[k][1] = last;
	}

	result[num_folds-1][1] = result[num_folds-1][1] + num_data_items % num_folds;

	return result;
}

double** get_train_data(double** dataset, int length, int num_folds, int fold)
{
	int** first_and_last_test = get_first_last_test(length, num_folds); // Первый и последний индекс строк помеченных тестовыми
	int num_train = length - (first_and_last_test[fold][1] - first_and_last_test[fold][0] + 1); // общее количество строк - число тестовых строк
	double** result = new double * [num_train];
	
	for (int i = 0; i < num_train; i++) result[i] = new double[];

	int i = 0; // Индекс результирующих/тестовых данных
	int ia = 0; // Индекс во всем наборе

	while (i < num_train)
	{
		if (ia < first_and_last_test[fold][0] || ia > first_and_last_test[fold][1]) // Это строка тестового набора
		{
			result[i] = dataset[ia];
			++i;
		}
		++ia;
	}

	return result;
}

double** get_test_data(double** dataset, int length, int num_folds, int fold)
{
	// Возвращает указатель на тестовые данные
	int** first_and_last_test = get_first_last_test(length, num_folds); // Первый и последний индекс строк помеченных тестовыми
	int num_test = first_and_last_test[fold][1] - first_and_last_test[fold][0] + 1;
	double** result = new double *[num_test];
	int ia = first_and_last_test[fold][0]; // Индеск во всем наборе
	
	for (int i = 0; i < num_test; ++i)
	{
		result[i] = dataset[ia]; // Индексы тестовых данных смежные
		++ia;
	}

	return result;
}

double cross_validate(double** &dataset, int length, int num_folds, double learn_rate, int features_num) {
	double theta0 = 0.0, theta1 = 0.0;
	double err;
	int* cum_err = new int[2];
	int num_iterations = 1000; // Число итераций

	cum_err[0] = 0;
	cum_err[1] = 0;
	
	for (int k = 0; k < num_folds; ++k)
	{
		double** train_data = get_train_data(dataset, length, num_folds, k); // Получаем обучающую выборку для разбиения
		double** test_data = get_test_data(dataset, length, num_folds, k); // Получаем тестовую выборку для разбиения
		
		err = full_search(train_data, 36, num_iterations, theta0, theta1, 0.0001);

		cout << "Обучающая выборка ****************" << endl;
		cout << "theta0 = " << theta0 << endl;
		cout << "theta1 = " << theta1 << endl;

		cum_err[0] += err;

		err = full_search(test_data, 12, num_iterations, theta0, theta1, 0.0001);

		cout << "Тестовая выборка *****************" << endl;
		cout << "theta0 = " << theta0 << endl;
		cout << "theta1 = " << theta1 << endl;

		cum_err[1] += err;
	}
	
	return (cum_err[0] * 1.0) / (cum_err[0] + cum_err[1]); // mean classification error
}

int main()
{
	setlocale(0,"Russian");

	ifstream input_set; // Поток с данными
	string line; // Строка в файле
	double **starting_set; // Массив в данными
	int features_num = 10; // Число признаков
	double j_best = 100000;
	int l = 0;

	// Создаем двумерный массив
	starting_set = new double *[20000];
	for (int i = 0; i < 20000; i++) starting_set[i] = new double [10];

	input_set.open("magic04.data");

	// Читаем файл и загружаем данные в массив
	while(!input_set.eof())
	{
		getline(input_set, line);

		for (int i = 0; i < line.length(); i++)
		{
			if (line[i] == ',')
				line[i] = ' ';
		}
		
		stringstream ss(line);
		float temp;
		vector<float> array;

		while (ss >> temp)
			array.push_back(temp); // splitted array
		
		for (int i = 0; i < 10; i++) {
			starting_set[l][i] = array[i];
		}

		l++;
	}

	float a = 0.0001; // Learning rate

	int m = 50; // Количество точек

	//int num_nodes[3] = { 4, 7, 3 };
	int num_folds = 4;
	
	// Создаем вектор подмножества признаков
	double** subset = new double *[m];
	for (int i = 0; i < m; i++) subset[i] = new double [2];

	for(int n = 0; n < features_num - 1; n++)
		for (int k = n + 1; k < features_num; k++) {
			// Выделяем признаки
			for(int i = 0; i < m; i++) {
				subset[i][0] = starting_set[i][n];
				subset[i][1] = starting_set[i][k];
			}
			double mce = cross_validate(subset, m, num_folds, a, features_num); // mean classification error across folds
			double mca = 1.0 - mce;

			cout << "CV error = " << mce << endl;
		}

	input_set.close();

	system("PAUSE");

	return(0);
}