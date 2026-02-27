#include <iostream>
#include <omp.h>
#include <fstream>
#include <iomanip>

using namespace std;

const int N = 50000;
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
    ofstream outFileSeq("heat_output_seq.csv");
    outFileSeq << fixed << setprecision(6);

    writeSnapshot(outFileSeq, 0, u, N);

    for (int i = 1; i <= T; i++)
    {
        for (int j = 1; j < N - 1; j++)
        {
            uNext[j] = 0.5 * (u[j - 1] + u[j + 1]);
        }

        double *temp = u;
        u = uNext;
        uNext = temp;

        if (i % 50 == 0)
        {
            writeSnapshot(outFileSeq, i, u, N);
        }
    }

    outFileSeq.close();
}

void heat_par()
{
    ofstream outFilePar("heat_output_par.csv");

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

int main()
{
    initialize();
    double t1 = omp_get_wtime();
    heat_seq();
    double t2 = omp_get_wtime();
    cout << "Seq time = " << t2 - t1 << "\n";
    initialize();
    double t3 = omp_get_wtime();
    heat_par();
    double t4 = omp_get_wtime();
    cout << "Par time = " << t4 - t3 << "\n";
    return 0;
}