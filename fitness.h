#pragma once
#include<iostream>
#include<vector>
#include<random>
#include<ctime>
using namespace std;
#define pi 3.1415926535
#define E  2.71828182845904523536

double do1t(vector<double>a, vector<double>b)
{
	double result = 0;
	for (int i = 0; i < a.size(); i++)
	{
		result += a[i] * b[i];
	}
	return result;
}

vector<double> divided(double b, vector<double>a)
{
	vector<double> temp = a;
	for (int i = 0; i < a.size(); i++)
	{
		temp[i] = b / a[i];
	}
	return temp;
}

vector<double> divided(vector<double>a, double b)
{
	vector<double> temp = a;
	for (int i = 0; i < a.size(); i++)
	{
		temp[i] = a[i] / b;
	}
	return temp;
}

vector<double> multi(double b, vector<double>a)
{
	vector<double> temp = a;
	for (int i = 0; i < a.size(); i++)
	{
		temp[i] = b * a[i];
	}
	return temp;
}

vector<double> sub(vector<double>a, vector<double> b)
{
	vector<double> temp = a;
	for (int i = 0; i < a.size(); i++)
	{
		temp[i] = a[i] - b[i];
	}
	return temp;
}

void print(vector<double>p)
{
	for (auto a : p)
		cout << a << " ";
	cout << endl;
}

class settings
{
public:
	int num;
	double min;
	double max;
	int D;
	double global_opt;
	vector<vector<double>>orthogonal;


	settings(int num)
	{
		this->num = num;
		cout << num << endl;
		setParameters();
	}

	void preProcessing()
	{
		static default_random_engine engine(time(nullptr));
		uniform_real_distribution<double> dis(1, 10);
		vector<double>a;


		a.resize(D);
		for (int i = 0; i < D; i++)
			orthogonal.push_back(a);

		for (int i = 0; i < D; i++)
			for (int j = 0; j < D; j++)
			{
				orthogonal[i][j] = dis(engine);
			}

		vector<vector<double>>o_copy = orthogonal;
		for (int i = 0; i < D; i++)
		{
			o_copy[i] = orthogonal[i];
			for (int j = 0; j < i; j++)
			{
				o_copy[i] = sub(o_copy[i], divided(do1t(orthogonal[i], o_copy[j]),
					multi(do1t(o_copy[j], o_copy[j]), o_copy[j])));
			}
		}
		for (int i = 0; i < D; i++)
		{
			print(o_copy[i]);
		}

		for (int i = 0; i < D; i++)
		{
			o_copy[i] = divided(o_copy[i], sqrt(do1t(o_copy[i], o_copy[i])));
		}
		for (int i = 0; i < D; i++)
		{
			print(o_copy[i]);
		}


	}

	void setParameters()
	{
		switch (num)
		{
		case 1:
			min = -100;
			max = 100;
			D = 30;
			global_opt = 0;
			break;
		case 2:
			min = -10;
			max = 10;
			D = 30;
			global_opt = 0;
			break;
		case 3:
			min = -10;
			max = 10;
			D = 30;
			global_opt = 0;
			break;
		case 4:
			min = -1.28;
			max = 1.28;
			D = 30;
			global_opt = 0;
			break;
		case 5:
			min = -500;
			max = 500;
			D = 30;
			global_opt = 0;
			break;
		case 6:
			min = -5.12;
			max = 5.12;
			D = 30;
			global_opt = 0;
			break;
		case 7:
			min = -32;
			max = 32;
			D = 30;
			global_opt = 0;
			break;
		case 8:
			min = -600;
			max = 600;
			D = 30;
			global_opt = 0;
			break;
		case 9:
			min = -50;
			max = 50;
			D = 30;
			global_opt = 0;
			break;
		case 10:
			min = -50;
			max = 50;
			D = 30;
			global_opt = 0;
			break;
		case 11:
			min = -500;
			max = 500;
			D = 30;
			global_opt = 0;
			break;
		case 12:

			min = -5.12;
			max = 5.12;
			D = 30;
			global_opt = 0;
			break;


		case 13:
			
			min = -32;
			max = 32;
			D = 30;
			global_opt = 0;
			break;
		case 14:
			min = -600;
			max = 600;
			D = 30;
			global_opt = 0;
			break;
		case 15:
			min = -100;
			max = 100;
			D = 30;
			global_opt = 390;
			break;
		case 16:
			min = -5;
			max = 5;
			D = 30;
			global_opt = -330;
			break;

		default:
			cout << 13 << endl;
			break;
		}
	}

	double fitnessFunction(vector<double>pos)
	{
		double result = 0.0;
		double part1 = 0.0;
		double part2 = 0.0;
		switch (num)
		{
		case 1:
			for (int i = 0; i < pos.size(); i++)
			{
				result += pow(pos[i], 2);
			}
			break;
		case 2:
			break;
		case 3:
			for (int i = 0; i < pos.size() - 1; i++)
			{
				result += 100 * pow(pos[i + 1] - pow(pos[i], 2), 2) + pow(pos[i] - 1, 2);
			}
			break;
		case 4:
			break;
		case 5:
			result += 418.9829*pos.size();
			for (int i = 0; i < pos.size(); i++)
			{
				result -= pos[i] * sin(pow(abs(pos[i]), 0.5));
			}
			break;
		case 6:
			for (int i = 0; i < pos.size(); i++)
			{
				result += pow(pos[i], 2) - 10 * cos(2 * pi*pos[i]) + 10;
			}
			break;
		case 7:

			break;
		case 8:
			break;
		case 9:
			break;
		case 10:
			break;
		case 11:
			break;
		case 12:
			break;
		case 13:

			for (int i = 0; i < pos.size(); i++)
			{
				part1 += pow(pos[i], 2);
				part2 += cos(2 * pi*pos[i]);
			}

			result = -20 * exp(-0.2*sqrt(part1 / pos.size())) - exp(part2 / pos.size()) + 20 + E;
			break;
		case 14:
			part2 = 1;
			for (int i = 0; i < pos.size(); i++)
			{
				part1 += pow(pos[i], 2);
				part2 *= cos(pos[i] / sqrt(i + 1));
			}

			result = part1 / 4000 - part2 + 1;
			break;
		case 15:
			break;
		case 16:
			break;

		}
		return result;
	}
};

