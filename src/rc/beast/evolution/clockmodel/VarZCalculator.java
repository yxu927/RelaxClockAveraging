package rc.beast.evolution.clockmodel;

public class VarZCalculator {

    /**
     * Approximate Var(Z) = (1/t^2) * double integral_{0..t,0..t} Cov( e^{V_a}, e^{V_b} ) da db,
     * where Cov( e^{V_a}, e^{V_b} ) = E[e^{V_a + V_b}] - E[e^{V_a}] * E[e^{V_b}].
     *
     * Using Simpson's rule in 2D. steps must be even as well.
     *
     * @param r0    rate at branch start (> 0)
     * @param rt    rate at branch end (> 0)
     * @param t     branch time length (> 0)
     * @param phi   autocorrelation parameter (Brownian variance scale)
     * @param steps must be an even number (e.g. 20, 40, etc.)
     * @return approximate Var(Z)
     */
    public static double computeVarZ(double r0, double rt, double t, double phi, int steps) {
        if (t <= 0.0 || steps < 2 || (steps % 2) != 0) {
            throw new IllegalArgumentException("Invalid t or steps. steps must be even >= 2.");
        }
        if (r0 <= 0.0 || rt <= 0.0) {
            throw new IllegalArgumentException("Rates must be positive.");
        }

        double v0 = Math.log(r0);
        double vT = Math.log(rt);

        double h = t / steps; // step size in each dimension

        // We'll do a 2D Simpson sum: sum_{i=0..steps} sum_{j=0..steps} W(i)*W(j) * Cov(e^{V_a}, e^{V_b})
        // Where a = i*h, b = j*h

        double sum2D = 0.0;

        for (int i = 0; i <= steps; i++) {
            double a = i * h;
            // Simpson weight in 1D:
            double w_i = ( (i==0 || i==steps) ? 1.0 : ( (i%2==1) ? 4.0 : 2.0 ) );

            for (int j = 0; j <= steps; j++) {
                double b = j * h;
                double w_j = ( (j==0 || j==steps) ? 1.0 : ( (j%2==1) ? 4.0 : 2.0 ) );

                double covVal = covExpVaVb(a, b, v0, vT, t, phi);
                sum2D += w_i * w_j * covVal;
            }
        }

        // Multiply by (h/3) for i dimension and (h/3) for j dimension => (h^2 / 9) if we strictly do Simpson
        double factorSimpson = (h / 3.0) * (h / 3.0);
        double integral2D = sum2D * factorSimpson;

        // multiply by (1/t^2)
        double varZ = integral2D / (t * t);
        return varZ;
    }

    /**
     * Cov( e^{V_a}, e^{V_b} ) = E[ e^{V_a + V_b} ] - E[e^{V_a}] * E[e^{V_b}].
     */
    private static double covExpVaVb(double a, double b,
                                     double v0, double vT,
                                     double t, double phi) {
        double eAB = EexpVaVb(a, b, v0, vT, t, phi);
        double eA  = EexpVa(a, v0, vT, t, phi);
        double eB  = EexpVa(b, v0, vT, t, phi);
        return eAB - (eA * eB);
    }

    /**
     * E[ e^{V_a} ] = exp( E[V_a] + 0.5 * Var[V_a] ).
     */
    private static double EexpVa(double a,
                                 double v0, double vT,
                                 double t, double phi) {
        // E(V_a) = v0 + (vT - v0)*(a/t)
        double meanVa = v0 + (vT - v0)*(a/t);
        // Var(V_a) = phi*(a*(t-a)/t)
        double varVa  = phi * (a*(t-a)/t);
        return Math.exp(meanVa + 0.5*varVa);
    }

    /**
     * E[ e^{V_a + V_b} ] = exp( E[V_a]+E[V_b] + 0.5*(Var(V_a)+Var(V_b)+2Cov(V_a,V_b)) ).
     */
    private static double EexpVaVb(double a, double b,
                                   double v0, double vT,
                                   double t, double phi) {
        // E[V_a], Var(V_a)
        double meanVa = v0 + (vT - v0)*(a/t);
        double varVa  = phi*(a*(t-a)/t);

        // E[V_b], Var(V_b)
        double meanVb = v0 + (vT - v0)*(b/t);
        double varVb  = phi*(b*(t-b)/t);

        // Cov(V_a, V_b): see Guindon eq. (16)-(20).
        //   Cov(V_a, V_b) = phi * (min(a,b)*(t-max(a,b)))/t
        //   Or piecewise. But an easy way:
        double covVab = 0.0;
        if (a <= b) {
            covVab = phi*(a*(t-b)/t);
        } else {
            covVab = phi*(b*(t-a)/t);
        }
        // E(V_a + V_b) = meanVa + meanVb
        double meanSum = meanVa + meanVb;
        double varSum  = varVa + varVb + 2.0*covVab;

        return Math.exp(meanSum + 0.5*varSum);
    }

    // quick test
    public static void main(String[] args) {
        double r0 = 2.0;
        double rt = 3.5;
        double t = 1.0;
        double phi = 0.2;
        int steps = 2000;

        double varZ = computeVarZ(r0, rt, t, phi, steps);
        System.out.println("VarZ (2D Simpson, steps="+steps+"): " + varZ);

        // compare with bigger steps
        int biggerSteps = 40;
        double varZ2 = computeVarZ(r0, rt, t, phi, biggerSteps);
        System.out.println("VarZ (2D Simpson, steps="+biggerSteps+"): " + varZ2);
    }

}
