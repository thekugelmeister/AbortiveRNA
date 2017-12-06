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

double jStringProbability(string substring, string s, string alphabet) {
   // Calculate the distribution of characters in the string, and use that to calculate the likelihood of finding the given substring in that string
   double charDistribution [alphabet.length()];
   for (int i = 0; i < alphabet.length(); i++) {
      int charMatches = count(s.begin(), s.end(), alphabet[i]);
      charDistribution[i] = double(charMatches) / s.length();
   }
   // QUICK AND DIRTY: For mitigating lack of ready-made dict-like functionality, brute force determining the mapping for rapid prototyping.
   double substringProbability = 1;
   for (int substrIdx = 0; substrIdx < substring.length(); substrIdx++) {
      // for (int alphaIdx = 0; alphaIdx < alphabet.length(); alphaIdx++) {
      //    if (substring[substrIdx] == alphabet[alphaIdx]) {
      //       substringProbability *= charDistribution[alphaIdx];
      //    }
      // }
     switch (substring[substrIdx]) {
     case 'A':
       substringProbability *= charDistribution[1];
       break;
     case 'T':
       substringProbability *= (charDistribution[0] + charDistribution[3]);
       break;
     case 'C':
       substringProbability *= charDistribution[3];
       break;
     case 'G':
       substringProbability *= (charDistribution[1] + charDistribution[2]);
     }
   }
   return substringProbability;
}

