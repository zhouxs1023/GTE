// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/QFNumber.h>

// The conversion functions here are used to obtain arbitrary-precision
// approximations to rational numbers and to quadratic field numbers.
// The arbitrary-precision arithmetic is described in
// https://www.geometrictools.com/Documentation/ArbitraryPrecision.pdf
// The quadratic field numbers and conversions are described in
// https://www.geometrictools.com/Documentation/QuadraticFields.pdf

namespace gte
{
    template <typename Rational>
    class APConversion
    {
    public:
        using QFN1 = QFNumber<Rational, 1>;
        using QFN2 = QFNumber<Rational, 2>;

        // Construction and destruction. The default precision is 53 to match
        // 'double' and the default maximum iterations is 8.  TODO: The maximum
        // iterations should be small to avoid the geometric growth in the
        // number of bits for the rational numbers generated during Newton's
        // method. Implement support for truncating rational numbers to a
        // user-specified number of bits and use this mechanism during the
        // Newton iterations.
        APConversion(
            int32_t numBitsPrecision = std::numeric_limits<double>::digits,
            uint32_t maxIterations = 8)
            :
            mZero(0),
            mOne(1),
            mThree(3),
            mFive(5),
            mNumBitsPrecision(numBitsPrecision),
            mMaxIterations(maxIterations),
            mEpsilon(std::ldexp(mOne, -std::numeric_limits<double>::digits)),
            mThreshold(std::ldexp(mOne, -mNumBitsPrecision))
        {
            LogAssert(numBitsPrecision > 0, "Invalid precision.");
            LogAssert(maxIterations > 0, "Invalid maximum iterations.");
        }

        ~APConversion()
        {
        }

        // Member access.
        void SetPrecision(int32_t numBitsPrecision)
        {
            LogAssert(numBitsPrecision > 0, "Invalid precision.");
            mNumBitsPrecision = numBitsPrecision;
            mThreshold = std::ldexp(mOne, -mNumBitsPrecision);
        }

        void SetMaxIterations(uint32_t maxIterations)
        {
            LogAssert(maxIterations > 0, "Invalid maximum iterations.");
            mMaxIterations = maxIterations;
        }

        inline int32_t GetPrecision() const
        {
            return mNumBitsPrecision;
        }

        inline uint32_t GetMaxIterations() const
        {
            return mMaxIterations;
        }

        // Disallow copying and moving.
        APConversion(APConversion const&) = delete;
        APConversion(APConversion&&) = delete;
        APConversion& operator=(APConversion const&) = delete;
        APConversion& operator=(APConversion&&) = delete;

        // The input a^2 is rational, but a itself is usually irrational,
        // although a rational value is allowed.
        uint32_t EstimateSqrt(Rational const& aSqr, Rational& aMin, Rational& aMax)
        {
            // Factor a^2 = r^2 * 2^e , where r^2 in [1/2,1). Compute s^2 and
            // the exponent used to generate the estimate of sqrt(a^2).
            Rational sSqr;
            int exponentA;
            PreprocessSqr(aSqr, sSqr, exponentA);

            // Use the FPU to estimate s = sqrt(sSqr) to 53-bit precision and
            // then round up.  Multiply by the appropriate exponent to obtain
            // bound aMax > a.
            aMax = GetMax(sSqr, exponentA);

            // Compute a lower bound aMin < a.
            aMin = aSqr / aMax;

            // Compute Newton iterates until convergence. The estimate closest
            // to a is aMin with aMin <= a <= aMax and a - aMin <= aMax - a.
            uint32_t iterate;
            for (iterate = 1; iterate <= mMaxIterations; ++iterate)
            {
                if (aMax - aMin < mThreshold)
                {
                    break;
                }
                // aMax = (aMin + aMax) / 2
                aMax = std::ldexp(aMin + aMax, -1);
                aMin = aSqr / aMax;
            }
            return iterate;
        }

        uint32_t EstimateApB(Rational const& aSqr, Rational const& bSqr,
            Rational& tMin, Rational& tMax)
        {
            // Factor a^2 = r^2 * 2^e , where r^2 in [1/2,1). Compute u^2 and
            // the exponent used to generate the estimate of sqrt(a^2).
            Rational uSqr;
            int32_t exponentA;
            PreprocessSqr(aSqr, uSqr, exponentA);

            // Factor b^2 = s^2 * 2^e , where s^2 in [1/2,1). Compute v^2 and
            // the exponent used to generate the estimate of sqrt(b^2).
            Rational vSqr;
            int32_t exponentB;
            PreprocessSqr(bSqr, vSqr, exponentB);

            // Use the FPU to estimate u = sqrt(u^2) and v = sqrt(v^2) to
            // 53 bits of precision and then round up.  Multiply by the
            // appropriate exponents to obtain bounds aMax > a and bMax > b.
            // This ensures tMax = aMax + bMax > a + b.
            Rational aMax = GetMax(uSqr, exponentA);
            Rational bMax = GetMax(vSqr, exponentB);
            tMax = aMax + bMax;

            // Compute a lower bound tMin < a + b.
            Rational a2pb2 = aSqr + bSqr;
            Rational a2mb2 = aSqr - bSqr;
            Rational a2mb2Sqr = a2mb2 * a2mb2;
            Rational tMaxSqr = tMax * tMax;
            tMin = (a2pb2 * tMaxSqr - a2mb2Sqr) / (tMax * (tMaxSqr - a2pb2));

            // Compute Newton iterates until convergence. The estimate closest
            // to a + b is tMin with tMin < a + b < tMax and
            // (a + b) - tMin < tMax - (a + b).
            uint32_t iterate;
            for (iterate = 1; iterate <= mMaxIterations; ++iterate)
            {
                if (tMax - tMin < mThreshold)
                {
                    break;
                }
                // tMax = (3/4) * tMax + (1/4) * tMin
                tMax = std::ldexp(mThree * tMax + tMin, -2);
                tMaxSqr = tMax * tMax;
                tMin = (a2pb2 * tMaxSqr - a2mb2Sqr) / (tMax * (tMaxSqr - a2pb2));
            }
            return iterate;
        }

        uint32_t EstimateAmB(Rational const& aSqr, Rational const& bSqr,
            Rational& tMin, Rational& tMax)
        {
            // Compute various quantities that are used later in the code.
            Rational a2tb2 = aSqr * bSqr;       // a^2 * b^2
            Rational a2pb2 = aSqr + bSqr;       // a^2 + b^2
            Rational a2mb2 = aSqr - bSqr;       // a^2 - b^2
            Rational a2mb2Sqr = a2mb2 * a2mb2;  // (a^2 - b^2)^2

            // Factor a^2 = r^2 * 2^e , where r^2 in [1/2,1). Compute u^2 and
            // the exponent used to generate the estimate of sqrt(a^2).
            Rational uSqr;
            int32_t exponentA;
            PreprocessSqr(aSqr, uSqr, exponentA);

            // Factor b^2 = s^2 * 2^e , where s^2 in [1/2,1). Compute v^2 and
            // the exponent used to generate the estimate of sqrt(b^2).
            Rational vSqr;
            int32_t exponentB;
            PreprocessSqr(bSqr, vSqr, exponentB);

            // Compute the sign of f''(a-b)/8 = a^2 - 3*a*b + b^2.  To test
            // a^2 + b^2 $ 3*a*b, where $ is one of <, =, >, square the test
            // to obtain (a^2)^2 - 7 * a^2 * b^2 + (b^2)^2 $ 0. The left-hand
            // side is/ computed as (a^2 - b^2)^2 - 5 * a^2 * b^2 to reduce
            // the number of arithmetic operations by re-using (a^2 - b^2)^2.
            Rational test = a2mb2Sqr - mFive * a2tb2;

            if (test > mZero)
            {
                // Choose an initial guess tMin < a - b. Use the FPU to
                // estimate u = sqrt(u^2) and v = sqrt(v^2) to 53 bits of
                // precision. Round and then multiply by the appropriate
                // exponents to obtain tMin = aMin - bMax < a - b.
                Rational aMin = GetMin(uSqr, exponentA);
                Rational bMax = GetMax(vSqr, exponentB);
                tMin = aMin - bMax;

                // Test whether tMin is in the positive f"(t) basin containing
                // a - b. If it is not, compute a tMin that is in the basis.
                Rational tMinSqr = tMin * tMin;
                Rational tInflectionSqr = a2pb2 / mThree;
                uint32_t iterate;
                if (tMinSqr < tInflectionSqr)
                {
                    Rational twoa2pb2 = std::ldexp(a2pb2, 1);
                    tMax = tInflectionSqr / tMin;
                    for (iterate = 1; iterate <= mMaxIterations; ++iterate)
                    {
                        Rational tMaxSqr = tMax * tMax;
                        Rational f = tMaxSqr * (tMaxSqr - twoa2pb2) + a2mb2Sqr;
                        if (f >= mZero && tMaxSqr <= a2pb2)
                        {
                            tMin = tMax;
                            tMinSqr = tMaxSqr;
                            break;
                        }
                        tMax = std::ldexp(tMin + tMax, -1);
                        tMin = tInflectionSqr / tMax;
                    }
                    if (iterate > mMaxIterations)
                    {
                        LogWarning("Exceeded maximum iterations. No guarantee Newton's converges.");
                        tMinSqr = tMin * tMin;
                    }
                }

                // Compute an upper bound tMax > a - b.
                tMax = (a2pb2 * tMinSqr - a2mb2Sqr) / (tMin * (tMinSqr - a2pb2));

                // Compute Newton iterates until convergence. The estimate
                // closest to a - b is tMax with tMin > a - b < tMax and
                // tMax - (a - b) < (a - b) - tMin.
                for (iterate = 1; iterate <= mMaxIterations; ++iterate)
                {
                    if (tMax - tMin < mThreshold)
                    {
                        break;
                    }
                    // tMin = (3 * tMin + tMax) / 4
                    tMin = std::ldexp(mThree * tMin + tMax, -2);
                    tMinSqr = tMin * tMin;
                    tMax = (a2pb2 * tMinSqr - a2mb2Sqr) / (tMin * (tMinSqr - a2pb2));
                }
                return iterate;
            }

            if (test < mZero)
            {
                // Choose an initial guess tMax > a - b. Use the FPU to
                // estimate u = sqrt(u^2) and v = sqrt(v^2) to 53 bits of
                // precision. Round and then multiply by the appropriate
                // exponents to obtain tMax = aMax - bMin > a - b.
                Rational aMax = GetMax(uSqr, exponentA);
                Rational bMin = GetMin(vSqr, exponentB);
                tMax = aMax - bMin;

                // Test whether tMax is in the negative f"(t) basin containing
                // a - b. If it is not, compute a tMax that is in the basis.
                Rational tMaxSqr = tMax * tMax;
                Rational tInflectionSqr = a2pb2 / mThree;
                uint32_t iterate;
                if (tMaxSqr > tInflectionSqr)
                {
                    Rational twoa2pb2 = std::ldexp(a2pb2, 1);
                    tMin = tInflectionSqr / tMax;
                    for (iterate = 1; iterate <= mMaxIterations; ++iterate)
                    {
                        Rational tMinSqr = tMin * tMin;
                        Rational f = tMinSqr * (tMinSqr - twoa2pb2) + a2mb2Sqr;
                        if (f <= mZero)
                        {
                            tMax = tMin;
                            tMaxSqr = tMinSqr;
                            break;
                        }
                        // tMax = (tMin + tMax) / 2
                        tMax = std::ldexp(tMin + tMax, -1);
                        tMin = tInflectionSqr / tMax;
                    }
                    if (iterate > mMaxIterations)
                    {
                        LogWarning("Exceeded maximum iterations. No guarantee Newton's converges.");
                        tMaxSqr = tMax * tMax;
                    }
                }

                // Compute a lower bound tMin < a - b.
                tMin = (a2pb2 * tMaxSqr - a2mb2Sqr) / (tMax * (tMaxSqr - a2pb2));

                // Compute Newton iterates until convergence. The estimate
                // closest to a - b is tMin with tMin < a - b < tMax and
                // (a - b) - tMin < tMax - (a - b).
                for (iterate = 1; iterate <= mMaxIterations; ++iterate)
                {
                    if (tMax - tMin < mThreshold)
                    {
                        break;
                    }
                    // tMax = (3 * tMax + tMin) / 4
                    tMax = std::ldexp(mThree * tMax + tMin, -2);
                    tMaxSqr = tMax * tMax;
                    tMin = (a2pb2 * tMaxSqr - a2mb2Sqr) / (tMax * (tMaxSqr - a2pb2));
                }
                return iterate;
            }

            // The 'test' value is a^4 - 7*a^2*b^2 + b^4 for rational a^2 and
            // b^2 and is itself rational. It cannot be zero. Define rational
            // r = a^2/b^2 so that test == 0 implies r^2 - 7*r^2 + 1 = 0.  The
            // roots of the quadratic are r = (7 +- sqrt(45))/2, which are
            // irrational, a contradiction to r being rational.
            LogError("This case cannot happen theoretically.");
        }

        uint32_t Estimate(QFN1 const& q, Rational& qMin, Rational& qMax, bool& qMinIsClosest)
        {
            Rational const& x = q.x[0];
            Rational const& y = q.x[1];
            Rational const& d = q.d;

            uint32_t numIterates;
            if (d != mZero && y != mZero)
            {
                Rational aSqr = y * y * d;
                numIterates = EstimateSqrt(aSqr, qMin, qMax);
                if (y > mZero)
                {
                    qMin = x + qMin;
                    qMax = x + qMax;
                    qMinIsClosest = true;
                }
                else
                {
                    Rational diff = x - qMax;
                    qMax = x - qMin;
                    qMin = diff;
                    qMinIsClosest = false;
                }
            }
            else
            {
                numIterates = 0;
                qMin = x;
                qMax = x;
                qMinIsClosest = true;
            }

            return numIterates;
        }

    private:
        void PreprocessSqr(Rational const& aSqr, Rational& rSqr, int& exponentA)
        {
            // Factor a^2 = r^2 * 2^e , where r^2 in [1/2,1).
            int32_t exponentASqr;
            rSqr = std::frexp(aSqr, &exponentASqr);
            if (exponentASqr & 1)  // odd exponent
            {
                // a = sqrt(2*r^2) * 2^{(e-1)/2}
                exponentA = (exponentASqr - 1) / 2;
                rSqr = std::ldexp(rSqr, 1);  // = 2*rSqr
                // rSqr in [1,2)
            }
            else  // even exponent
            {
                // a = sqrt(r^2) * 2^{e/2}
                exponentA = exponentASqr / 2;
                // rSqr in [1/2,1)
            }
        }

        Rational GetMin(Rational const& aSqr, int const exponent)
        {
            Rational aMin = Rational(std::sqrt((double)aSqr)) - mEpsilon;
            aMin = std::ldexp(aMin, exponent);
            return aMin;
        }

        Rational GetMax(Rational const& sqr, int const exponent)
        {
            Rational aMax = Rational(std::sqrt((double)sqr)) + mEpsilon;
            aMax = std::ldexp(aMax, exponent);
            return aMax;
        }

        Rational const mZero, mOne, mThree, mFive;
        int32_t mNumBitsPrecision;
        uint32_t mMaxIterations;
        Rational mEpsilon;
        Rational mThreshold;
    };
}
