#include <iostream>
#include <omp.h>
#include <fstream>
#include <iomanip>

using namespace std;

const int N = 500;
const int T = 500;
double *u = new double[N];
double *uNext = new double[N];

void writeSnapshot(ofstream &out, int t, const double u[], int N)
{
    out << t;
    for (int i = 0; i < N; ++i)
    {
        out << ',' << u[i];
    }
    out << "\n";
}

void initialize()
{
    for (int i = 0; i < N; i++)
    {
        u[i] = 0.0;
    }
    u[N / 2] = 9999.0;
}

void heat_seq()
{
    // ofstream outFileSeq("heat_output_seq.csv");
    // outFileSeq << fixed << setprecision(6);

    // writeSnapshot(outFileSeq, 0, u, N);

    for (int i = 1; i <= T; i++)
    {
        for (int j = 1; j < N - 1; j++)
        {
            uNext[j] = 0.5 * (u[j - 1] + u[j + 1]);
        }

        double *temp = u;
        u = uNext;
        uNext = temp;

        // if (i % 50 == 0)
        // {
        //     writeSnapshot(outFileSeq, i, u, N);
        // }
    }

    // outFileSeq.close();
}

void heat_par()
{
    ofstream outFilePar("heat_output.csv");

    writeSnapshot(outFilePar, 0, u, N);

    for (int i = 1; i <= T; i++)
    {
#pragma omp parallel for
        for (int j = 1; j < N - 1; j++)
        {
            uNext[j] = 0.5 * (u[j - 1] + u[j + 1]);
        }

        double *temp = u;
        u = uNext;
        uNext = temp;

        if (i % 50 == 0)
        {
            writeSnapshot(outFilePar, i, u, N);
        }
    }

    outFilePar.close();
}




void heat_par(int x)
{
	omp_set_num_threads(x);
    
    for (int i = 1; i <= T; i++)
    {
#pragma omp parallel for
        for (int j = 1; j < N - 1; j++)
        {
            uNext[j] = 0.5 * (u[j - 1] + u[j + 1]);
        }

        double* temp = u;
        u = uNext;
        uNext = temp;

        
    }

}

void heat_par_conv()
{
	double ep = 1e-4 ;
	double diff;

    do {
        diff = 0;
#pragma omp parallel for reduction (max:diff)
        for (int j = 1; j < N - 1; j++)
        {

            uNext[j] = 0.5 * (u[j - 1] + u[j + 1]);
			
            curr_diff = abs(uNext[j] - u[j]);
            if (curr_diff > diff)
				diff = curr_diff;
        }
        
        double* temp = u;
        u = uNext;
        uNext = temp;

        while (diff > ep);
    }

}



void test_num_threads(double seq_time) {

    for (i = 8; i >= 1; i / 2) {
		initialize();
		double t1 = omp_get_wtime();
		heat_par(i);
		double t2 = omp_get_wtime();
		cout << "Parallel time with " << i << " threads = " << t2 - t1 << "\n";
		cout << "Speedup " << seq_time / (t2 - t1) << "\n";
		cout << "Efficiency " << (seq_time / (t2 - t1)) / i << "\n";
    }

}



int main()
{
    initialize();
    double t1 = omp_get_wtime();
    heat_seq();
    double t2 = omp_get_wtime();
    cout << "Sequential time = " << t2 - t1 << "\n";
    initialize();
    double t3 = omp_get_wtime();
    heat_par();
    double t4 = omp_get_wtime();
    cout << "Parallel time = " << t4 - t3 << "\n";
	double speedup = (t2 - t1) / (t4 - t3);
	cout << "Speedup = " << speedup << "\n";
	cout << "Efficiency = " << speedup / omp_get_num_threads() << "\n";

	test_num_threads(t2-t1);

    cout << "advanced version " << endl;
	initialize();
	double t5 = omp_get_wtime();
    heat_par_conv();
    double t6 = omp_get_wtime();
    cout << "parallel time advanced version: " << t6 - t5 << endl;


    return 0;
}