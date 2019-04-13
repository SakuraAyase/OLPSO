#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#include"particles.h"
#include"fitness.h"


using namespace std;


settings set = settings(14);

double fitnessFunction(vector<double>pos)
{
	return set.fitnessFunction(pos);
}


bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class Particle
{
public:
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;
	vector<double> P0;
	int stagnated = 0;

	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;
	}
};



double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}


class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c,double dt, double percent,int G)
	{
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c = c;
		this->dt = dt;
		this->percent = percent;
		this->G = G;
		particles.resize(m);
	}
	

	vector<vector<int>> getOA()
	{
		int M = pow(2, ceil(log2(dim + 1)));
		int N = M - 1;
		int u = ceil(log2(dim + 1));
		vector<vector<int>>OA;
		for (int i = 0; i < M; i++)
		{
			vector<int>temp;
			temp.resize(N);
			OA.push_back(temp);
		}

		int k = 1;
		for (int b = 1; b <= N; b++)
		{
			for (int a = 1; a <= M; a++)
			{
				if (b == pow(2, k - 1))
					OA[a - 1][b - 1] = (int)(floor((double)(a - 1) / pow(2, u - k))) % 2;
				else
					OA[a - 1][b - 1] = (OA[a - 1][pow(2, k - 2) - 1] + OA[a - 1][b - pow(2, k - 2) - 1]) % 2;
			}
			if (b == pow(2, k - 1))
				k++;
		}

		return OA;
	}

	void orthogonal(int n)
	{
		//cout << gBestIndex;
		//print(particles[gBestIndex].pBest);
		//cout << 'v' << " ";
		static default_random_engine engine(time(nullptr));
		uniform_int_distribution<int> dis(0, m-1);
		vector<double>Pb;
		vector<double>Po;
		Pb.resize(dim);
		Po.resize(dim);
		int g_index = gBestIndex;
		while(n == g_index)
			g_index = dis(engine);


		vector<vector<int>>OA = getOA();
		vector<vector<double>>Fs;
		int index = 0;

		//cout << fitnessFunction(particles[n].pBest) << " " << fitnessFunction(particles[g_index].pBest) << endl;
		for (int j = 0; j < OA.size(); j++)
		{
			vector<double>p;
			p.resize(dim);
			for (int k = 0; k < dim; k++)
			{
				if (OA[j][k] == 0)
				{ 
					p[k] = particles[n].pBest[k];
				}
				else
				{
					p[k] = particles[g_index].pBest[k];
				}
			}
			//cout <<"fur"<< fitnessFunction(p) << endl;
			Fs.push_back(p);
			if (better(fitnessFunction(Fs[j]),fitnessFunction(Fs[index])))
				index = j;
		}

		Pb = Fs[index];
		//cout << fitnessFunction(Pb) << endl;
		//cout<<Fs.size(
		
		for (int j = 0; j < dim; j++)
		{
			double f1 = 0.0;
			int num1 = 0;
			double f2 = 0.0;
			int num2 = 0;

			for (int k = 0; k < OA.size(); k++)
			{
				if (OA[k][j] == 0)
				{
					num1++;
					f1 += fitnessFunction(Fs[k]);
					//cout <<"1 :"<< fitnessFunction(Fs[k]) << endl;
				}
				else
				{
					num2++;
					f2 += fitnessFunction(Fs[k]);
					//cout << "2 :" << fitnessFunction(Fs[k]) << endl;
				}
			}
			//cout << (num1) << " " << (num2) << endl;
			if (better((f1 / num1),(f2 / num2)))
				Po[j] = particles[n].pBest[j];
			else
				Po[j] = particles[g_index].pBest[j];
		}

		if (better(fitnessFunction(Po), fitnessFunction(Pb)))
			particles[n].P0 = Po;
		else
			particles[n].P0 = Pb;
		//cout << "v: " << fitnessFunction(Po)<<" "<< fitnessFunction(Pb) << endl;
		//cout << "s: " << fitnessFunction(particles[n].P0)  << endl;
	}

	void initialParticles(int i)
	{
		particles[i].stagnated = 0;
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		particles[i].P0.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(this->min, this->max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].pBest[j] = particles[i].position[j];
		}
		
	}
	//
	

	void initialAllParticles()
	{
		for (int i = 0; i < m; i++)
			initialParticles(i);
		for (int i = 0; i < m; i++)
		{
			if (better(fitnessFunction(particles[i].pBest), fitnessFunction(particles[gBestIndex].pBest)))
				gBestIndex = i;
		}
		
		
		for (int i = 0; i < m; i++)
			orthogonal(i);
	}

	void inertiaWeight()
	{
		double t = T / ((double)Tmax);
		w = 0.9 - 0.5*t;
	}

	void updateParticle(int i)
	{
		bool legal = true;

		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);

			particles[i].velocity[j] = w * particles[i].velocity[j] +
				c * randDouble(0, 1) * (particles[i].P0[j] - particles[i].position[j]);

			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = range;

			if (particles[i].velocity[j] < -range)
				particles[i].velocity[j] = -range;
			particles[i].position[j] += dt * particles[i].velocity[j];
			if (particles[i].position[j] > max || particles[i].position[j] < min)
				legal = false;
		}

		if (legal)
		{
			if (better(fitnessFunction(particles[i].position), fitnessFunction(particles[i].pBest)))
			{
				particles[i].pBest = particles[i].position;
				particles[i].stagnated = 0;

				if (better(fitnessFunction(particles[i].pBest), fitnessFunction(particles[gBestIndex].pBest)))
					gBestIndex = i;
			}
			else
			{
				particles[i].stagnated++;
				if (particles[i].stagnated > G)
				{
					orthogonal(i);
					particles[i].stagnated = 0;
				}
			}
		}

	}


	void updateAllParticles()
	{
		inertiaWeight();
		for (int i = 0; i < m; i++)
			updateParticle(i);
		T++;
	}

	double getFitness()
	{
		int index = 0;
		for (int i = 0; i < m; i++)
		{
			if (better(fitnessFunction(particles[i].pBest),fitnessFunction(particles[index].pBest)))
				index = i;
		}
		return fitnessFunction(particles[index].pBest);
	}
private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	int G;

	double w;
	double max;
	double min;
	double c;

	double dt;//时间步长
	double percent;

	int gBestIndex;

	vector<Particle> particles;


};

void run(vector<double>& result1)
{
	int dim = set.D;
	int m = 40;
	int Tmax = 2000;
	double max = set.max;
	double min = set.min;
	double c = 2;
	double dt = 1.0;
	double percent = 0.2;
	int G = 5;

	PSO pso = PSO(dim, m, Tmax, max, min, c, dt, percent,G);
	pso.getOA();
	pso.initialAllParticles();

	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		pso.updateAllParticles();
		cout << ":";
		//fitness.push_back(pso.getFitness());
		fitness.push_back(pso.getFitness());
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;
	}

	result1 = fitness;
}

int main()
{

	//set.preProcessing();
	int times = 1;
	int interval = 10;
	vector<double> result1;

	run(result1);

	for (int i = 1; i < times; i++)
	{
		vector<double> result1_temp;
		run(result1_temp);
		for (int j = 0; j < result1_temp.size(); j++)
		{
			result1[j] += result1_temp[j];
		}
	}
	for (int j = 0; j < result1.size(); j++)
	{
		result1[j] /= times;
	}

	for (int j = 0; j < result1.size(); j++)
	{
		if (j%interval == 0)
			cout << result1[j] << " ";
	}

	system("pause");
}