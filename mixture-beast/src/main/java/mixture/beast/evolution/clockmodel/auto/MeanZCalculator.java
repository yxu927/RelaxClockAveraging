package mixture.beast.evolution.clockmodel.auto;


public class MeanZCalculator {

    private static final double MAX_RATE_MULTIPLIER = 1.0e4;


    public static double computeMeanZ(double r0, double rt,
                                      double t, double phi,
                                      int order) {

        if (t <= 0.0) {
            throw new IllegalArgumentException("MeanZCalculator.computeMeanZ: t must be > 0, got " + t);
        }
        if (order < 1) {
            throw new IllegalArgumentException("MeanZCalculator.computeMeanZ: order must be >= 1, got " + order);
        }
        if (r0 <= 0.0 || rt <= 0.0) {
            throw new IllegalArgumentException("MeanZCalculator.computeMeanZ: rates must be positive. r0="
                    + r0 + ", rt=" + rt);
        }

        // v0, vT
        final double v0 = Math.log(r0);
        final double vT = Math.log(rt);

        final double a = (vT - v0) / t;
        final double b = phi / (2.0 * t);

        final double A = a + b * t;
        final double B = b;

        double sumAll = 0.0;

        for (int k = 0; k < order; k++) {
            double term_k = 0.0;

            for (int m = 0; m <= k; m++) {
                long c = binomial(k, m);
                double sign = ((m & 1) == 0) ? 1.0 : -1.0;

                double cAB = c * sign * Math.pow(A, k - m) * Math.pow(B, m);

                int p = k + m;
                double factor = cAB * (1.0 / (p + 1.0)) * Math.pow(t, p + 1.0);
                term_k += factor;
            }

            double factInv = 1.0 / factorial(k);
            sumAll += factInv * term_k;
        }

        double result = Math.exp(v0) * sumAll / t;


        boolean bad = false;

        if (!Double.isFinite(result) || result <= 0.0) {
            bad = true;
        } else {

            double maxEndpoint = Math.max(r0, rt);
            double maxReasonable = MAX_RATE_MULTIPLIER * maxEndpoint;
            if (result > maxReasonable) {
                bad = true;
            }
        }

        if (bad) {

            result = simpsonMeanZ(v0, vT, t, phi, 256);
        }

        if (!Double.isFinite(result) || result <= 0.0) {
            result = Math.sqrt(r0 * rt);
        }

        return result;
    }



    private static double simpsonMeanZ(double v0, double vT,
                                       double t, double phi,
                                       int steps) {

        if (steps <= 0) {
            throw new IllegalArgumentException("Simpson steps must be positive, got " + steps);
        }

        final double h = t / steps;

        double sum = integrand(v0, vT, 0.0, t, phi)
                + integrand(v0, vT, t,   t, phi);

        for (int i = 1; i < steps; i++) {
            double s = i * h;
            double val = integrand(v0, vT, s, t, phi);
            sum += ((i & 1) == 1) ? 4.0 * val : 2.0 * val;
        }

        double integral = sum * (h / 3.0);
        return integral / t;
    }


    private static double integrand(double v0, double vT,
                                    double s, double t, double phi) {
        double meanVs = v0 + (vT - v0) * (s / t);
        double varVs  = phi * (s * (t - s) / t);
        return Math.exp(meanVs + 0.5 * varVs);
    }

    private static double factorial(int n) {
        double f = 1.0;
        for (int i = 1; i <= n; i++) {
            f *= i;
        }
        return f;
    }

    private static long binomial(int n, int k) {
        if (k < 0 || k > n) return 0;
        long c = 1;
        for (int i = 0; i < k; i++) {
            c = c * (n - i) / (i + 1);
        }
        return c;
    }

    public static void main(String[] args) {
        double r0   = 1.0;
        double rt   = 1.0;
        double tLen = 10;
        double phi  = 0;
        int order   = 16;

        double valSeries = computeMeanZ(r0, rt, tLen, phi, order);
        System.out.println("E[Z] (safe series+Simpson) = " + valSeries);

        double v0 = Math.log(r0);
        double vT = Math.log(rt);
        double valSimpson = simpsonMeanZ(v0, vT, tLen, phi, 200000);
        System.out.println("Pure Simpson(200000 steps) = " + valSimpson);
        System.out.println("difference = " + Math.abs(valSeries - valSimpson));
    }
}
