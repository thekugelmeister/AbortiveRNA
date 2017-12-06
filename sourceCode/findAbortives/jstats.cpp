#include "jstats.h"
#include "errno.h"

//calculate the mean of an array of integers
double jMean(int a[], int len)
{
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += a[i];
    }
    return sum / len;
}

//calculate the sample standard deviation of an array of integers
double jSSD(int a[], int len)
{
    double mean = jMean(a, len);
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += ((a[i] - mean) * (a[i] - mean));
    }
    if (sum == 0) {
        if (mean != 0) {
            cerr << "STRANGE MEAN ALERT!!!!" << endl;
            cerr << "Sum: " << sum << endl
                 << "Mean: " << mean << endl;
        }
        sum += ((1 - mean) * (1 - mean));
    }
    // if (isnan(sqrt(abs(sum / (len - 1))))) {
    //     cerr << "NAN ALERT!!!!" << endl;
    // }
    errno = 0;
    //double temp = sqrt(abs(sum / (len - 1)));
    if (errno == EDOM) {
        cerr << "EDOM ALERT!!!!" << endl;
        cerr << "Sum: " << sum << endl
             << "Len - 1: " << len - 1 << endl
             << "Sum / Len - 1: " << sum / (len - 1) << endl
             << "Abs(Sum / Len - 1): " << abs(sum / (len - 1)) << endl;
        exit(1);
    }
    return (sqrt(abs(sum / (len - 1))));
}

double jSTT(int a[], int len, int val)
{
    //return ((jMean(a, len) - val) / (jSSD(a, len) / sqrt(len)));
    return ((val - jMean(a, len)) / (jSSD(a, len)));
}
