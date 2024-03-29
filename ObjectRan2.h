#include <stdio.h>

constexpr long IM1 = 2147483563;
constexpr long IM2 = 2147483399;
constexpr float AM = 1.0f / (float) IM1;
constexpr long IMM1 = IM1-1;
constexpr long IA1 = 40014;
constexpr long IA2 = 40692;
constexpr long IQ1 = 53668;
constexpr long IQ2 = 52774;
constexpr long IR1 = 12211;
constexpr long IR2 = 3791;
constexpr int NTAB = 32;
constexpr long NDIV  = 1 + IMM1/NTAB;
constexpr float EPS = 1.2e-7f;
constexpr float RNMX = 1.0f - EPS;

class ran2_gen {
    long idum = 0, idum2 = 12345, iy = 0, iv[NTAB];

public:

    ran2_gen(long seed = 0) {
        if (seed < 0) idum = -seed;
        else idum = seed + 1;
        idum2 = idum;
        int j; long k;
        for (j = NTAB + 7; j >= 0; j--) {
            k = idum / IQ1;
            idum = IA1 * (idum - k*IQ1) - k*IR1;
            if (idum < 0) idum += IM1;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
    }

    ran2_gen(const ran2_gen& gen) {
        idum = gen.idum;
        idum2 = gen.idum2;
        iy = gen.iy;
        for (int j = 0; j < NTAB; j++) {
            iv[j] = gen.iv[j];
        }
    }

    ran2_gen(FILE *in) {
        int sum = 0;
        sum += fscanf(in, "%ld\t%ld\t%ld\n", &idum, &idum2, &iy);
        for(int j = 0; j < NTAB; j++) sum += fscanf (in, "%ld\n", iv + j);
        if (sum != NTAB + 3) perror("BUILDING RAN2 FROM FILE: wrong file layout.\n");
    }

    float ran2() {
        long k = idum/IQ1;
        idum = IA1 * (idum - k*IQ1) - k*IR1;
        if (idum < 0) idum += IM1;
        k = idum2/IQ2;
        idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
        if (idum2 < 0) idum2 += IM2;
        int j= (int) (iy/NDIV);
        iy=iv[j]-idum2;
        iv[j] = idum;
        if (iy < 1) iy += IMM1;
        float temp;
        if ((temp = AM * (float)iy) > RNMX) return RNMX;
        else return temp;
    }

    long ran2int(int rmax) {
        long k = idum/IQ1;
        idum = IA1 * (idum - k*IQ1) - k*IR1;
        if (idum < 0) idum += IM1;
        k = idum2/IQ2;
        idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
        if (idum2 < 0) idum2 += IM2;
        int j = (int) (iy/NDIV);
        iy=iv[j]-idum2;
        iv[j] = idum;
        if (iy < 1) iy += IMM1;
        return iy % rmax;
    }

    void print_status(FILE *out = stdout) {
        fprintf(out, "%ld\t%ld\t%ld\n", idum, idum2, iy);
        for (int j = 0; j < NTAB; j++) fprintf(out, "%ld\n", iv[j]);
    }
};
